import pandas as pd
import pytest

from scripts import build_evaluation_board as board


def test_board_preserves_missing_geometry_and_scaffold_boundary(tmp_path):
    board.build(tmp_path)
    candidates = pd.read_csv(tmp_path / "candidates.csv")
    denovo = candidates[candidates.track.eq("de_novo")]
    assert len(denovo) > 0
    assert denovo.predicted_neg_log10_kd_nm.isna().all()
    assert denovo.affinity_evidence.eq("fixed-scaffold SPR model not applicable").all()
    assert candidates.lab_status.eq("unmeasured").all()
    structures = pd.read_csv(tmp_path / "structures.csv")
    assert not structures[["query_id", "model", "target_context", "seed"]].isna().any().any()
    initial = structures[structures.target_context.eq("initial_split")]
    assert len(initial) > 0 and initial.wt_contact_count.isna().all()
    assert "not in this audit" in (tmp_path / "report.md").read_text()


def test_board_rejects_double_counted_structural_run(tmp_path, monkeypatch):
    read = board.pd.read_csv

    def duplicate(path, *args, **kwargs):
        frame = read(path, *args, **kwargs)
        if "context-panel" in str(path):
            return pd.concat([frame, frame.iloc[:1]], ignore_index=True)
        return frame

    monkeypatch.setattr(board.pd, "read_csv", duplicate)
    with pytest.raises(ValueError, match="duplicate query/model/context/seed"):
        board.build(tmp_path)


def test_board_rejects_mislabelled_replica_model(tmp_path, monkeypatch):
    if not (board.ROOT / "results/structural-replicates/scores.csv").exists():
        pytest.skip("Optional replica audit absent")
    read = board.pd.read_csv

    def mislabel(path, *args, **kwargs):
        frame = read(path, *args, **kwargs)
        if "structural-replicates" in str(path):
            frame["model"] = "a-different-model"
        return frame

    monkeypatch.setattr(board.pd, "read_csv", mislabel)
    with pytest.raises(ValueError, match="Unexpected replica model"):
        board.build(tmp_path)
