"""Report matched target-only-template refolds, not measured binding."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from biotite.sequence import ProteinSequence
from biotite.structure import get_residues
from biotite.structure.io.pdbx import CIFFile, get_structure
from src.her2_structure_audit import contacts, remap_target_residues

QUERY_MANIFEST = Path('results/experiment-20261114/denovo-refold-queries.json')
CHAIN_ORDER = ('C', 'D', 'H', 'L')
SEED = 42
SAMPLE = 1


def load_queries(manifest):
    queries = {}
    for query in json.loads(manifest.read_text())['queries']:
        query_id = query['query_id']
        # ASVS 5.3.2: query IDs become directory and file name components.
        if not re.fullmatch(r'[A-Za-z0-9_-]+', query_id) or query_id in queries:
            raise ValueError(f'Invalid or duplicate manifest query ID: {query_id}')
        chains = {chain['chain_id']: chain['sequence'] for chain in query['chains']}
        if len(query['chains']) != len(CHAIN_ORDER) or set(chains) != set(CHAIN_ORDER):
            raise ValueError(f'Expected four target and antibody chains in manifest query {query_id}')
        queries[query_id] = chains
    if not {'control_wt', 'control_facs_low'} <= queries.keys():
        raise ValueError('Query manifest is missing a required control')
    return queries


def expected_artifacts(root, queries):
    # run_denovo_refold.sh requests seed 42 and one 1-indexed diffusion sample.
    stems = {query_id: root / query_id / f'seed_{SEED}' / f'{query_id}_seed_{SEED}_sample_{SAMPLE}'
             for query_id in queries}
    for suffix in ('_scores.json', '_model.cif'):
        expected = {stem.with_name(stem.name + suffix) for stem in stems.values()}
        actual = set(root.rglob(f'*{suffix}'))
        if actual != expected:
            missing = ', '.join(str(path.relative_to(root)) for path in sorted(expected - actual))
            unexpected = ', '.join(str(path.relative_to(root)) for path in sorted(actual - expected))
            raise ValueError(f'Refold {suffix} files differ from manifest. Missing: {missing}. Unexpected: {unexpected}.')
    return stems


def validate_model(atoms, expected_chains, path):
    chains = tuple(dict.fromkeys(atoms.chain_id))
    if chains != CHAIN_ORDER:
        raise ValueError(f'Unexpected chain order: {path}: {chains}')
    for chain_id in chains:
        _, residue_names = get_residues(atoms[atoms.chain_id == chain_id])
        sequence = ''.join(ProteinSequence.convert_letter_3to1(name) for name in residue_names)
        if sequence != expected_chains[chain_id]:
            raise ValueError(f'Model sequence differs from the manifest for chain {chain_id}: {path}')


def interface_scores(score, path):
    matrix = np.asarray(score['chain_chain_interface_predicted_tm_score'], dtype=float)
    overall = np.asarray(score['interface_predicted_tm_score'], dtype=float)
    if (matrix.shape != (len(CHAIN_ORDER), len(CHAIN_ORDER)) or overall.shape != ()
            or not np.isfinite(matrix).all() or not np.isfinite(overall)):
        raise ValueError(f'Invalid or nonfinite interface scores: {path}')
    return matrix, float(overall)


def audit(campaign, manifest=QUERY_MANIFEST):
    queries = load_queries(manifest)
    stems = expected_artifacts(campaign / 'denovo-target-only-protenix-v3', queries)
    mapping = json.loads(Path('configs/her2_target_residue_map.json').read_text())
    reference = get_structure(CIFFile.read(campaign / 'her2_fv_1n8z.cif'), model=1)
    reference_contacts = set(contacts(reference, ['A', 'B'], ['C', 'D'])['target_contact_residues'])
    rows = []
    for query_id, stem in sorted(stems.items()):
        path = stem.with_name(stem.name + '_scores.json')
        score = json.loads(path.read_text())
        structure_path = stem.with_name(stem.name + '_model.cif')
        atoms = get_structure(CIFFile.read(structure_path), model=1)
        validate_model(atoms, queries[query_id], structure_path)
        matrix, overall = interface_scores(score, path)
        geometry = contacts(remap_target_residues(atoms, mapping), ['H', 'L'], ['C', 'D'])
        observed = set(geometry['target_contact_residues'])
        row = {'query_id': query_id,
               'overall_iptm': overall,
               'vh_vl_iptm': float(matrix[2, 3]),
               'max_antibody_target_pair_iptm': float(matrix[:2, 2:].max()),
               'wt_footprint_recall': len(observed & reference_contacts) / len(reference_contacts),
               'wt_footprint_jaccard': len(observed & reference_contacts) / len(observed | reference_contacts),
               'scores_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'structure_sha256': hashlib.sha256(structure_path.read_bytes()).hexdigest(),
               'evidence': 'single_seed_predicted_complex_not_measured_binding', **geometry}
        for ti, target in enumerate(CHAIN_ORDER[:2]):
            for ai, antibody in enumerate(CHAIN_ORDER[2:], start=2):
                row[f'iptm_{antibody}_{target}'] = float(matrix[ai, ti])
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = audit(args.campaign)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'geometry.json').write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n')
    frame = pd.DataFrame(rows).drop(columns=['target_contact_residues', 'antibody_contact_residues'])
    frame.to_csv(args.output / 'scores.csv', index=False)
    print(frame[['query_id', 'vh_vl_iptm', 'max_antibody_target_pair_iptm', 'wt_footprint_recall']].to_string(index=False))


if __name__ == '__main__':
    main()
