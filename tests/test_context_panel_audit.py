import json

import numpy as np
import pytest
from biotite.structure import AtomArray

from scripts import audit_context_panel as audit


def test_repeated_completed_query_fails_before_scoring():
    with pytest.raises(ValueError, match='duplicate'):
        audit.exact_coverage([*audit.PANEL, audit.PANEL[0]], audit.PANEL, 'RF3 results receipt')


def test_protenix_missing_query_in_result_receipt(tmp_path, monkeypatch):
    root = tmp_path / audit.PROTENIX / 'split'
    root.mkdir(parents=True)
    monkeypatch.setattr(audit, 'CAMPAIGN', tmp_path)
    config = {'missing_msas': 'fail', 'num_diffusion_samples': 1, 'query_failures': 'fail',
              'seeds': [42], 'template_paths': ['split-target.cif'],
              'weights': '/app/weights/protenix-v2.pt', 'write_detailed_scores': True}
    (root / 'config.json').write_text(json.dumps({'schema_version': 'foundry.cli.v1', **config}))
    (root / 'result.json').write_text(json.dumps({
        'task': 'predict', 'schema_version': 'foundry.cli.v1',
        'provenance': {'module': {'name': 'apheris-protenix', 'version': '1.3.0'},
                       'inputs': [{'path': 'split-queries.json'}], 'resolved_config': config},
        'result': {'predictions': [{'query_id': 'control_wt'}]},
    }))
    with pytest.raises(ValueError, match=r"missing \['control_spr_weak'\]"):
        audit.validate_protenix_run('split', ['control_wt', 'control_spr_weak'],
                                   {'template': 'split-target.cif', 'queries': 'split-queries.json'})


def test_chain_roles_allow_rf3_reorder_and_renaming_but_reject_mismatch():
    expected = {'T': 'CYS', 'H': 'EVQ', 'L': 'DIQ'}
    observed = {'Y': 'DIQ', 'X': 'CYS', 'Z': 'EVQ'}
    assert audit.sequence_roles(observed, expected, 'RF3') == {'Y': 'L', 'X': 'T', 'Z': 'H'}
    with pytest.raises(ValueError, match='sequence mismatch'):
        audit.sequence_roles({**observed, 'Z': 'EVA'}, expected, 'RF3')
    with pytest.raises(ValueError, match='IDs or order mismatch'):
        audit.sequence_roles({'H': 'EVQ', 'T': 'CYS', 'L': 'DIQ'}, expected, 'Protenix',
                             require_query_order=True)


def test_rf3_summary_order_uses_array_values_and_sequence_roles():
    roles = {'Y': 'L', 'X': 'T', 'Z': 'H'}
    tokens = np.array(['Y', 'X', 'Z'])
    matrix = np.array([[0, 6, 1], [8, 0, 9], [2, 4, 0]], dtype=float)
    summary = {
        'chain_pair_pae': [[None, 1.5, 6.5], [None, None, 7.0], [None, None, None]],
        'chain_pair_pae_min': [[None, 2.0, 4.0], [None, None, 6.0], [None, None, None]],
    }
    order = audit.rf3_summary_order(roles, tokens, {'pae': matrix}, summary, 'RF3')
    assert order == ['H', 'L', 'T']
    scores = audit.score_pairs(order, np.array([[1, 0.8, 0.22], [0.8, 1, 0.31], [0.22, 0.31, 1]]))
    assert scores == {'vh_vl_iptm': 0.8, 'max_antibody_target_pair_iptm': 0.31,
                      'iptm_H_T': 0.22, 'iptm_L_T': 0.31}
    with pytest.raises(ValueError, match='no unique chain order'):
        audit.rf3_summary_order(roles, tokens, {'pae': matrix + 1}, summary, 'RF3')


def test_split_gap_maps_to_canonical_positions_by_role():
    c, d = audit.ARMS['split'].values()
    atoms = AtomArray(len(c) + len(d) + 1)
    atoms.chain_id = np.array(['B'] * len(d) + ['H'] + ['A'] * len(c))
    atoms.res_id = np.array([*range(1, len(d) + 1), 1, *range(1, len(c) + 1)])
    mapped = audit.canonical_target(atoms, {'B': 'D', 'H': 'H', 'A': 'C'}, {'C': list(c), 'D': list(d)})
    assert list(mapped.res_id[mapped.chain_id == 'A']) == list(range(511, 603))
    assert list(mapped.res_id[mapped.chain_id == 'B']) == list(range(613, 630))
    assert mapped.res_id[mapped.chain_id == 'H'].item() == 1
    assert not set(range(603, 613)) & set(mapped.res_id)
    assert atoms.res_id[0] == 1


def test_target_alignment_preserves_antibody_displacement():
    fixed = np.array([[0, 0, 0], [3, 0, 0], [0, 4, 0], [1, 1, 5]], dtype=float)
    axes = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    shift = np.array([2, 8, 4])
    def transform(points):
        return points @ axes + shift

    ref_ab = np.array([[5, 1, 7], [5, 2, 7], [5, 2, 8]])
    displaced = np.array([12, 0, 0])
    rotation, translation = audit.target_transform(transform(fixed), fixed)
    assert audit.rmsd(transform(fixed) @ rotation + translation, fixed) == pytest.approx(0, abs=1e-10)
    assert audit.rmsd(transform(ref_ab + displaced) @ rotation + translation, ref_ab) == pytest.approx(12)
