import numpy as np
import pandas as pd
import pytest

from scripts import benchmark_hybrid_spr as hybrid


def test_blends_ignore_outcomes_and_preserve_row_order():
    frame = pd.DataFrame({name: np.arange(4.) + i for i, name in enumerate((*hybrid.BASELINES, *hybrid.COMPONENTS))})
    expected = hybrid.fixed_blends(frame)
    frame["neg_log10_kd_nm"] = [999, -999, 20, 0]
    actual = hybrid.fixed_blends(frame.iloc[::-1])
    for name in hybrid.RULES:
        np.testing.assert_array_equal(actual[name], expected[name][::-1])
    frame.loc[0, hybrid.COMPONENTS[1]] = np.nan
    with pytest.raises(ValueError, match="Nonfinite"):
        hybrid.fixed_blends(frame)


def test_unfrozen_predictions_rejected_before_evaluation(tmp_path):
    source = tmp_path / "predictions.csv"
    source.write_text("cdrh3,position_ridge\nAAAAAAAAAA,0\n")
    with pytest.raises(ValueError, match="fingerprint"):
        hybrid.run(source, tmp_path / "output")
    assert not (tmp_path / "output").exists()
