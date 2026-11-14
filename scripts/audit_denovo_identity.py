"""Compare exploratory Fv designs with two references, not the Adaptyv database.

IMGT-position identity uses the union of numbered positions as denominator.
This is an explicit local diagnostic, not the unpublished official alignment.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure import get_residues
from biotite.sequence import ProteinSequence
from immunum import Annotator

from src.her2_candidates import VH_PREFIX, VH_SUFFIX, WT_CDRH3, TRASTUZUMAB_VL


def numbered_identity(left, right):
    positions = set(left) | set(right)
    matches = sum(left.get(p) == right.get(p) for p in positions)
    return matches, len(positions)


def audit(campaign, reference_path):
    annotator = Annotator(chains=['H', 'K', 'L'], scheme='imgt', min_confidence=0.1)
    references = {'trastuzumab': {
        'H': annotator.number(VH_PREFIX + WT_CDRH3 + VH_SUFFIX).numbering,
        'L': annotator.number(TRASTUZUMAB_VL).numbering}}
    adalimumab = {}
    for entity in json.loads(reference_path.read_text())['entities']:
        result = annotator.number(entity['entity_poly']['pdbx_seq_one_letter_code_can'])
        assert not result.error
        adalimumab['H' if result.chain == 'H' else 'L'] = result.numbering
    assert set(adalimumab) == {'H', 'L'}
    references['adalimumab'] = adalimumab
    rows = []
    for arm in ['denovo-pilot', 'denovo-fullseq-pilot', 'denovo-fullseq-k1']:
        files = sorted((campaign / arm).glob('design/generation/inference_samples/*/*_generated.cif'))
        if not files:
            raise ValueError(f'No generated structures for {arm}')
        for path in files:
            atoms = get_structure(CIFFile.read(path), model=1)
            chains = {}
            sequences = {}
            # Organizer outputs preserve H/L; HER2 itself has Ig-like domains
            # and must never be admitted by antibody numbering alone.
            for chain in ['H', 'L']:
                _, names = get_residues(atoms[atoms.chain_id == chain])
                seq = ''.join(ProteinSequence.convert_letter_3to1(name) for name in names)
                number = annotator.number(seq)
                if not number.error and number.confidence >= 0.1:
                    kind = 'H' if number.chain == 'H' else 'L'
                    if kind != chain:
                        raise ValueError(f'Chain type mismatch {path}: {chain}')
                    if kind in chains:
                        raise ValueError(f'Duplicate antibody chain {path}: {kind}')
                    chains[kind] = number
                    sequences[kind] = seq
            if set(chains) != {'H', 'L'}:
                raise ValueError(f'Expected paired Fv: {path}')
            row = {'arm': arm, 'design': path.parent.name,
                   'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                   'vh_sequence': sequences['H'], 'vl_sequence': sequences['L'],
                   'paired_length': sum(map(len, sequences.values())),
                   'vh_confidence': chains['H'].confidence,
                   'vl_confidence': chains['L'].confidence,
                   'canonical_cysteines_only': all(
                       sequences[k].count('C') == 2 and chains[k].numbering.get('23') == 'C'
                       and chains[k].numbering.get('104') == 'C' for k in chains),
                   'official_novelty': 'unverified'}
            for name, ref in references.items():
                counts = [numbered_identity(chains[k].numbering, ref[k]) for k in ['H', 'L']]
                row[f'{name}_paired_imgt_identity'] = sum(c[0] for c in counts) / sum(c[1] for c in counts)
            row['below_70_both_local_references'] = all(row[f'{name}_paired_imgt_identity'] < .7 for name in references)
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--reference', type=Path, default=Path('configs/adalimumab_reference.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    frame = audit(args.campaign, args.reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(frame.groupby('arm').agg(n=('design', 'size'), mean_wt_identity=('trastuzumab_paired_imgt_identity','mean'), mean_scaffold_identity=('adalimumab_paired_imgt_identity','mean'), below_both=('below_70_both_local_references','sum')).to_string())


if __name__ == '__main__':
    main()
