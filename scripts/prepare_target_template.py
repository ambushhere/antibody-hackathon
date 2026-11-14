"""Remove antibody chains and entities while retaining target template metadata."""
import argparse
from pathlib import Path
import gemmi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    structure = gemmi.read_structure(str(args.input))
    for model in structure:
        for name in [chain.name for chain in model if chain.name not in ['C', 'D']]:
            del model[name]
    structure.setup_entities()
    present = {residue.subchain for chain in structure[0] for residue in chain}
    for i in range(len(structure.entities) - 1, -1, -1):
        entity = structure.entities[i]
        entity.subchains = [name for name in entity.subchains if name in present]
        if not entity.subchains:
            del structure.entities[i]
    structure.assign_label_seq_id()
    assert {chain.name for chain in structure[0]} == {'C', 'D'}
    assert {name for entity in structure.entities for name in entity.subchains} == present
    structure.make_mmcif_document().write_file(str(args.output))


if __name__ == '__main__':
    main()
