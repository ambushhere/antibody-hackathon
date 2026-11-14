import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import gemmi
from immunum import Annotator
from scripts.audit_denovo_identity import numbered_identity
from scripts.prepare_target_template import main as prepare_template


class DenovoAuditTests(unittest.TestCase):
    def test_numbered_identity_counts_gaps_and_trims_constant_domains(self):
        self.assertEqual(numbered_identity({'1': 'A', '2': 'C'}, {'1': 'A', '3': 'G'}), (1, 3))
        annotator = Annotator(chains=['H', 'K', 'L'], scheme='imgt', min_confidence=0.1)
        reference = json.loads(Path('configs/adalimumab_reference.json').read_text())
        for entity in reference['entities']:
            sequence = entity['entity_poly']['pdbx_seq_one_letter_code_can']
            full = annotator.number(sequence)
            segment = annotator.segment(sequence)
            variable = ''.join(getattr(segment, field) for field in ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4'])
            self.assertLess(len(variable), len(sequence))
            self.assertEqual(full.numbering, annotator.number(variable).numbering)

    def test_target_template_removes_antibody_entities_not_only_atoms(self):
        structure = gemmi.Structure()
        model = gemmi.Model('1')
        for chain_id in ['A', 'B', 'C', 'D']:
            chain = gemmi.Chain(chain_id)
            residue = gemmi.Residue()
            residue.name = 'ALA'
            residue.seqid = gemmi.SeqId(1, ' ')
            residue.entity_type = gemmi.EntityType.Polymer
            residue.subchain = chain_id
            atom = gemmi.Atom()
            atom.name = 'CA'
            atom.element = gemmi.Element('C')
            residue.add_atom(atom)
            chain.add_residue(residue)
            model.add_chain(chain)
            entity = gemmi.Entity(chain_id)
            entity.entity_type = gemmi.EntityType.Polymer
            entity.polymer_type = gemmi.PolymerType.PeptideL
            entity.subchains = [chain_id]
            entity.full_sequence = ['ALA']
            structure.entities.append(entity)
        structure.add_model(model)
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / 'input.cif', Path(temp) / 'target.cif'
            structure.make_mmcif_document().write_file(str(source))
            with patch('sys.argv', ['prepare', '--input', str(source), '--output', str(output)]):
                prepare_template()
            result = gemmi.read_structure(str(output))
            self.assertEqual({chain.name for chain in result[0]}, {'C', 'D'})
            self.assertEqual({chain for entity in result.entities for chain in entity.subchains}, {'C', 'D'})


if __name__ == '__main__':
    unittest.main()
