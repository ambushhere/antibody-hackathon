"""Audit installed abguard against the entire frozen optimization proposal pool."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def audit(folder=Path('results/abguard-expanded')):
    source_path = Path('results/screen-inputs/expanded-preflight.json')
    source = json.loads(source_path.read_text())
    submitted = list(csv.DictReader((folder/'input.csv').open()))
    annotated = json.loads((folder/'annotations.json').read_text())
    if len(source) != 743 or len(submitted) != len(source) or len(annotated) != len(source):
        raise ValueError('Incomplete full-pool evidence')
    by_id = {r['sequence_id']: r for r in annotated}
    if len(by_id) != len(source):
        raise ValueError('Duplicate annotation IDs')
    records = []
    for i, (original, row) in enumerate(zip(source, submitted)):
        identity = f"pool_{i:04d}_{original['cdrh3']}"
        if row['sequence_id'] != identity or row['track'] != 'optimization':
            raise ValueError('Changed pool order or declaration')
        a = by_id[identity]
        for column, key in [('heavy_sequence', 'vh_sequence'), ('light_sequence', 'vl_sequence')]:
            if original[key] != row[column] or row[column] != a[column]:
                raise ValueError('Changed input sequence')
        if type(a['is_valid']) is not bool or type(a['is_novel_enough']) is not bool:
            raise ValueError('Unchecked tool result')
        records.append({'sequence_id': identity, 'cdrh3': original['cdrh3'],
                        'origin': original['origin'], 'local_pass': bool(original.get('strict_local_gate')),
                        'abguard_pass': a['is_valid'], 'novelty_class': a['novelty_class'],
                        'cdrh3_min_distance': a['cdrh3_min_distance'], 'reason': a['novelty_reason']})
    with (folder/'designs.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(records)
    result = {
        'endpoint': 'installed_tool_eligibility_not_biological_performance',
        'rows': len(records), 'passed': sum(r['abguard_pass'] for r in records),
        'local_tool_disagreements': [r['sequence_id'] for r in records if r['local_pass'] != r['abguard_pass']],
        'errors': dict(Counter(e for r in annotated for e in r['errors'])),
        'by_origin': {origin: {'rows': len(rows), 'passed': sum(r['abguard_pass'] for r in rows)}
                      for origin in sorted({r['origin'] for r in records})
                      for rows in [[r for r in records if r['origin'] == origin]]},
        'sources': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in (source_path, folder/'input.csv', folder/'annotations.json')},
    }
    (folder/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    print(json.dumps(audit(), indent=2))
