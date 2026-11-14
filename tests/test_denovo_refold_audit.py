import numpy as np
import pytest
from biotite.structure import AtomArray

from scripts import audit_denovo_refolds as refold


@pytest.fixture
def score_tree(tmp_path):
    root = tmp_path / 'denovo-target-only-protenix-v3'
    for query_id in refold.load_queries(refold.QUERY_MANIFEST):
        seed = root / query_id / 'seed_42'
        seed.mkdir(parents=True)
        for suffix in ('scores.json', 'model.cif'):
            (seed / f'{query_id}_seed_42_sample_1_{suffix}').touch()
    return root


def test_missing_control_cannot_be_replaced_by_another_sample(score_tree):
    (score_tree / 'control_facs_low/seed_42/control_facs_low_seed_42_sample_1_scores.json').unlink()
    (score_tree / 'control_wt/seed_42/control_wt_seed_42_sample_2_scores.json').touch()
    with pytest.raises(ValueError, match=r'Missing: .*control_facs_low.*Unexpected: .*sample_2'):
        refold.audit(score_tree.parent)


def test_model_requires_exact_seed_and_sample(score_tree):
    model = score_tree / 'control_wt/seed_42/control_wt_seed_42_sample_1_model.cif'
    wrong_seed = model.parent.with_name('seed_7')
    wrong_seed.mkdir()
    model.rename(wrong_seed / 'control_wt_seed_7_sample_1_model.cif')
    with pytest.raises(ValueError, match=r'Missing: .*seed_42.*Unexpected: .*seed_7'):
        refold.audit(score_tree.parent)


@pytest.mark.parametrize('chain_id', refold.CHAIN_ORDER)
def test_model_rejects_wrong_sequence_with_the_right_chain_order(chain_id):
    expected = dict(zip(refold.CHAIN_ORDER, 'ACDE'))
    atoms = AtomArray(4)
    atoms.chain_id = np.array(refold.CHAIN_ORDER)
    atoms.res_name = np.array(['ALA', 'CYS', 'ASP', 'GLU'])
    refold.validate_model(atoms, expected, 'fixture.cif')
    atoms.res_name[refold.CHAIN_ORDER.index(chain_id)] = 'TYR'
    with pytest.raises(ValueError, match=f'Model sequence differs from the manifest for chain {chain_id}'):
        refold.validate_model(atoms, expected, 'fixture.cif')


@pytest.mark.parametrize('matrix, overall', [(np.eye(3), 0.5), (np.eye(4), np.nan), (np.full((4, 4), np.inf), 0.5)])
def test_interface_scores_reject_invalid_shapes_and_nonfinite_values(matrix, overall):
    with pytest.raises(ValueError, match='Invalid or nonfinite interface scores'):
        refold.interface_scores({
            'chain_chain_interface_predicted_tm_score': matrix,
            'interface_predicted_tm_score': overall,
        }, 'scores.json')
