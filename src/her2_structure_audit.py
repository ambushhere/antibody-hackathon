"""Explicit antibody-to-antigen geometry, never VH/VL interface confidence.

Input chain IDs must be supplied from the structure provenance. This reports
geometry only, not affinity, competition or thermal stability.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from scipy.spatial.distance import cdist


def contacts(atoms, antibody_chains, target_chains, cutoff=4.5):
    if set(antibody_chains) & set(target_chains):
        raise ValueError('Antibody and target chains overlap')
    available = set(atoms.chain_id)
    if not set(antibody_chains + target_chains) <= available:
        raise ValueError('Requested chain missing')
    heavy = ~np.isin(atoms.element, ['H', 'D'])
    ab = atoms[heavy & np.isin(atoms.chain_id, antibody_chains)]
    ag = atoms[heavy & np.isin(atoms.chain_id, target_chains)]
    if not len(ab) or not len(ag):
        raise ValueError('Empty antibody or target')
    distances = cdist(ab.coord, ag.coord)
    if not np.isfinite(distances).all():
        raise ValueError('Nonfinite coordinates')
    aidx, tidx = np.where(distances < cutoff)
    target_residues = sorted({f'{ag.chain_id[i]}:{ag.res_id[i]}' for i in tidx})
    ab_residues = sorted({f'{ab.chain_id[i]}:{ab.res_id[i]}' for i in aidx})
    return {'target_contact_residues': target_residues,
            'antibody_contact_residues': ab_residues,
            'heavy_atom_contact_pairs': int(len(tidx)),
            'heavy_atom_clash_pairs_below_1p5A': int((distances < 1.5).sum()),
            'minimum_target_distance_A': float(distances.min()),
            'contact_cutoff_A': cutoff}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structures', nargs='+', type=Path, required=True)
    p.add_argument('--antibody-chains', nargs='+', required=True)
    p.add_argument('--target-chains', nargs='+', required=True)
    p.add_argument('--reference', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    def read(path):
        return get_structure(CIFFile.read(path), model=1)
    ref = None
    if args.reference:
        ref = set(contacts(read(args.reference), args.antibody_chains, args.target_chains)['target_contact_residues'])
    rows = []
    for path in args.structures:
        row = contacts(read(path), args.antibody_chains, args.target_chains)
        row['source'] = str(path)
        row['evidence_type'] = 'predicted_geometry_not_experimental_binding'
        if ref is not None:
            observed = set(row['target_contact_residues'])
            row['reference_contact_recall'] = len(observed & ref) / len(ref) if ref else None
            row['reference_contact_jaccard'] = len(observed & ref) / len(observed | ref) if observed | ref else None
            row['numbering_assumption'] = 'matching target chain IDs and residue numbering'
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=2) + '\n')

if __name__ == '__main__':
    main()
