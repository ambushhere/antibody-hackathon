import json
import shutil
from pathlib import Path

import pytest

from scripts.reconcile_abguard import reconcile


def test_saved_tool_evidence_and_exact_export(tmp_path):
    result = reconcile(Path('results/abguard-validation/raw'), tmp_path/'report', tmp_path/'pass.csv')
    assert result['passed'] == 52
    assert result['accepted_classes'] == {'novel': 11, 'de_novo': 41}
    assert result['explicit_submission_tracks'] == {'optimization': 10, 'denovo': 42}
    assert result['track_class_mismatch_ids'] == ['denovo-scale-her2novel-k1_30']
    assert result['rejected_ids'] == ['denovo-scale-antibody-k1_23']
    assert result['meets_minimum_28_design_count']
    assert not result['is_final_ranked_submission']


@pytest.mark.parametrize('tamper', ['unchecked', 'sequence', 'missing', 'duplicate'])
def test_reject_incomplete_or_mismatched_evidence(tmp_path, tamper):
    raw = tmp_path/'raw'
    shutil.copytree('results/abguard-validation/raw', raw)
    path = raw/'annotated.json'
    rows = json.loads(path.read_text())
    if tamper == 'unchecked':
        rows[0]['is_novel_enough'] = None
    elif tamper == 'sequence':
        rows[0]['heavy_sequence'] += 'A'
    elif tamper == 'missing':
        rows.pop()
    else:
        rows.append(rows[0])
    path.write_text(json.dumps(rows))
    with pytest.raises(ValueError):
        reconcile(raw, tmp_path/'report', tmp_path/'pass.csv')
    assert not (tmp_path/'pass.csv').exists()
