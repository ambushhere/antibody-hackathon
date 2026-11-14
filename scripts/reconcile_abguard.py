"""Reconcile saved abguard evidence with exact input rows; never submit designs."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

FIELDS = ['sequence_id', 'kind', 'heavy_sequence', 'light_sequence', 'design_description']


def reconcile(raw: Path, output: Path, accepted: Path):
    original = list(csv.DictReader((raw / 'candidates.csv').open()))
    annotations = json.loads((raw / 'annotated.json').read_text())
    revalidated = list(csv.DictReader((raw / 'revalidated.csv').open()))
    ids = [r['sequence_id'] for r in original]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate input IDs')
    by_id = {r['sequence_id']: r for r in annotations}
    if len(by_id) != len(annotations) or set(by_id) != set(ids):
        raise ValueError('Annotation coverage differs from input')
    for row in original:
        annotated = by_id[row['sequence_id']]
        if any(row[f] != annotated[f] for f in FIELDS):
            raise ValueError('Annotation input mismatch')
        if type(annotated['is_novel_enough']) is not bool or type(annotated['is_valid']) is not bool:
            raise ValueError('Unchecked novelty or validity cannot be accepted')
        if annotated['is_valid'] and not annotated['is_novel_enough']:
            raise ValueError('Validity and novelty disagree')
    passed = [r for r in original if by_id[r['sequence_id']]['is_valid']]
    # Strict CLI revalidation must preserve every accepted row, including description.
    if revalidated != passed:
        raise ValueError('Strict revalidation differs from accepted input rows')
    wt = json.loads((raw / 'wt-control.json').read_text())[0]
    if wt['is_valid'] or wt['is_novel_enough'] or wt['cdrh3_min_distance'] != 0:
        raise ValueError('Known WT must fail novelty at distance zero')
    generated = json.loads((raw / 'generator-all64.json').read_text())
    generator_input = list(csv.DictReader((raw / 'generator-all64.csv').open()))
    if len(generated) != 64 or len({r['sequence_id'] for r in generated}) != 64:
        raise ValueError('Incomplete generator benchmark')
    generated_by_id = {r['sequence_id']: r for r in generated}
    if len(generator_input) != 64 or {r['sequence_id'] for r in generator_input} != set(generated_by_id):
        raise ValueError('Generator input coverage mismatch')
    for row in generator_input:
        observed = generated_by_id.get(row['sequence_id'], {})
        if (any(row[f] != observed.get(f) for f in FIELDS)
                or type(observed.get('is_novel_enough')) is not bool
                or type(observed.get('is_valid')) is not bool):
            raise ValueError('Generator evidence mismatch or unchecked novelty')
    output.mkdir(parents=True, exist_ok=True)
    accepted.parent.mkdir(parents=True, exist_ok=True)
    with accepted.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator='\n')
        w.writeheader()
        w.writerows(passed)
    columns = ['sequence_id', 'novelty_class', 'is_valid', 'is_novel_enough',
               'identity_to_trastuzumab', 'cdrh3', 'cdrh3_min_distance',
               'required_cdrh3_mutations', 'cdrh3_nearest_source', 'cdrh3_nearest_id',
               'n_high_severity', 'n_liabilities', 'novelty_reason']
    with (output / 'designs.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(by_id[i] for i in ids)
    summary = {
        'status': 'installed_tool_eligibility_only_not_submitted',
        'input_rows': len(original), 'passed': len(passed), 'rejected': len(original)-len(passed),
        'accepted_classes': dict(Counter(by_id[r['sequence_id']]['novelty_class'] for r in passed)),
        'accepted_with_high_severity_motifs': sum(by_id[r['sequence_id']]['n_high_severity'] > 0 for r in passed),
        'rejected_ids': [i for i in ids if not by_id[i]['is_valid']],
        'generation_denovo_but_tool_novel': [i for i in ids if i.startswith('denovo') and by_id[i]['novelty_class']=='novel'],
        'within_20_design_final_limit': len(passed)<=20,
        'binding_epitope_thermal_status': 'not_measured',
        'wt_control': {k: wt[k] for k in ('is_valid', 'is_novel_enough', 'cdrh3_min_distance', 'n_high_severity')},
        'generator_arms': {
            arm: {'input': len(rows), 'passed': sum(r['is_valid'] for r in rows),
                  'passed_classes': dict(Counter(r['novelty_class'] for r in rows if r['is_valid']))}
            for arm in sorted({r['sequence_id'].rsplit('_', 1)[0] for r in generated})
            for rows in [[r for r in generated if r['sequence_id'].rsplit('_', 1)[0] == arm]]
        },
        'raw_files_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(raw.iterdir()) if p.is_file()},
        'accepted_csv_sha256': hashlib.sha256(accepted.read_bytes()).hexdigest(),
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, default=Path('results/abguard-validation/raw'))
    parser.add_argument('--output', type=Path, default=Path('results/abguard-validation'))
    parser.add_argument('--accepted', type=Path, default=Path('output/submission-draft/team-3-abguard-pass.csv'))
    a = parser.parse_args()
    result = reconcile(a.raw, a.output, a.accepted)
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},indent=2))
