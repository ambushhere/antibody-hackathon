import pandas as pd
import pytest

from scripts import audit_spr_distance as distance
from src import her2_benchmark as facs
from src import her2_spr_benchmark as spr
from test_spr_candidate_ranking import replay_files


def run_audit(replay_files):
    return distance.audit(replay_files['data_dir'], replay_files['oof_path'], replay_files['metrics_path'])


def test_distance_uses_canonical_facs_sequences(replay_files):
    rows, report = run_audit(replay_files)
    # The padded lowercase training sequence differs from this variant by two edits.
    variant = rows.set_index('cdrh3').loc['AAAAAAAAAA']
    assert variant['nearest_facs_trainval'] == 2
    assert variant['nearest_spr_outer_train'] >= 3
    assert report['subsets']['primary']['n'] == 8
    assert report['subsets']['distance_at_least_3_to_both_training_sources']['n'] == 7


def test_distance_rejects_a_swapped_split(replay_files):
    data_dir = replay_files['data_dir']
    train, test = data_dir / 'train.csv', data_dir / 'test.csv'
    train_bytes, test_bytes = train.read_bytes(), test.read_bytes()
    train.write_bytes(test_bytes)
    test.write_bytes(train_bytes)
    with pytest.raises(ValueError, match='train does not match the completed benchmark input hash'):
        run_audit(replay_files)


@pytest.mark.parametrize(('column', 'invalid', 'error'), [
    ('cdrh3', 'PPPPPPPPPP', 'Primary cohort IDs differ'),
    ('group', '99', 'Primary sequence groups differ'),
    ('fold', '99', 'Outer fold .* differs'),
    (spr.TARGET, '', f'OOF {spr.TARGET} must contain finite numbers'),
    (spr.TARGET, '2.0', f'OOF {spr.TARGET} differs'),
    ('position_ridge', 'NaN', 'OOF position_ridge must contain finite numbers'),
    (facs.SCORE, 'inf', f'OOF {facs.SCORE} must contain finite numbers'),
])
def test_distance_rejects_inconsistent_oof_rows(replay_files, column, invalid, error):
    path = replay_files['oof_path']
    frame = pd.read_csv(path, dtype=str)
    frame.loc[0, column] = invalid
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match=error):
        run_audit(replay_files)
