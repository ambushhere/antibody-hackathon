"""Prepare HER2 targeted Fv CDR redesign using organizer antibody scaffold.

This is an exploratory de novo arm, not proof of novelty or submission eligibility.
"""
import argparse
import json
import shutil
from pathlib import Path
from biotite.structure.io.pdbx import CIFFile, get_structure, set_structure


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--docs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.docs / 'Data/12-her2/her2_fv_1n8z.cif'
    atoms = get_structure(CIFFile.read(source), model=1)
    target = atoms[(atoms.chain_id == 'C') | (atoms.chain_id == 'D')]
    assert len(target) > 0
    for chain, resid, atom in [('D',593,'NZ'),('C',573,'CD1'),('D',591,'CG2'),('C',572,'CD')]:
        assert ((target.chain_id == chain) & (target.res_id == resid) & (target.atom_name == atom)).sum() == 1
    cif = CIFFile()
    set_structure(cif, target)
    cif.write(args.output / 'her2_target.cif')
    shutil.copy2(args.docs / 'Data/03-fv/framework_adalimumab.cif', args.output)
    original = json.loads((args.docs / 'Data/03-fv/design.json').read_text())
    spec = next(iter(original.values()))
    spec['input'] = ['her2_target.cif', 'framework_adalimumab.cif']
    spec['contig'] = spec['contig'].replace('F1-76', 'C489-580, /0, D591-607')
    spec['conditions']['C_CRD']['true'] = [{'select': "sel('C/*/*/[N,CA,C,O]')"}, {'select': "sel('D/*/*/[N,CA,C,O]')"}]
    spec['conditions']['C_SEQ']['true'] = [v for v in spec['conditions']['C_SEQ']['true'] if v['select'] != "sel('F')"] + [{'select': "sel('C')"}, {'select': "sel('D')"}]
    spec['conditions']['C_HOT']['true'] = [{'select': f"sel('{c}/*/{r}/{a}')"} for c,r,a in [('D',593,'NZ'),('C',573,'CD1'),('D',591,'CG2'),('C',572,'CD')]]
    (args.output / 'design.json').write_text(json.dumps({'her2_denovo_fv': spec}, indent=2) + '\n')

if __name__ == '__main__':
    main()
