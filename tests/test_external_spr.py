import numpy as np
import pandas as pd
import pytest

from scripts import benchmark_external_spr as external
from src import her2_benchmark as facs


def test_reference_offsets_and_chain_verification():
    vh = facs.VH_PREFIX + facs.WT_CDRH3 + facs.VH_SUFFIX
    reference = f"""
    Heavy Chain Sequence: `{vh}AST`
    Light Chain Sequence: `{facs.WT_VL}RTV`
    HCDR1: `{external.H1_WT}`
    HCDR2: `{external.H2_WT}`
    HCDR3: `{external.H3_WT}`
    """
    audit = external.audit_reference(reference)
    assert audit["linear_vh_bounds_1_based_inclusive"] == {
        "HCDR1": [26, 33], "HCDR2": [51, 58], "HCDR3": [97, 109], "model_window": [99, 108],
    }
    assert audit["per_row_paired_sequences_supplied"] is False
    with pytest.raises(ValueError, match="WT annotation differs"):
        external.audit_reference(reference.replace(facs.WT_VL, "A" + facs.WT_VL[1:]))


@pytest.mark.parametrize(
    ("h3", "h1", "h2", "expected"),
    [
        ("SRWGGDGFYAMDY", external.H1_WT, external.H2_WT, ("WGGDGFYAMD", "compatible_window")),
        ("SRAAAAAAAAAAY", external.H1_WT, external.H2_WT, ("AAAAAAAAAA", "compatible_window")),
        ("SRWGGDGFYAMDY", None, None, ("WGGDGFYAMD", "compatible_window")),
        ("SRWGGDGFYAMDY", "GFNIKDTA", external.H2_WT, (None, "h1_h2_mismatch")),
        ("SRWGGDGFYAMDY", external.H1_WT, "", (None, "h1_h2_mismatch")),
        ("SRWGGDGFYAMY", external.H1_WT, external.H2_WT, (None, "unsupported_h3_length")),
        ("SRWGGDGFYAMDDY", external.H1_WT, external.H2_WT, (None, "unsupported_h3_length")),
        ("ARWGGDGFYAMDY", external.H1_WT, external.H2_WT, (None, "non_wt_h3_anchors")),
        ("SRWGGDGFYAMDF", external.H1_WT, external.H2_WT, (None, "non_wt_h3_anchors")),
        ("SRWGGDGFYAMXY", external.H1_WT, external.H2_WT, (None, "noncanonical_h3")),
    ],
)
def test_window_mapping(h3, h1, h2, expected):
    assert external.map_window(h3, h1, h2) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("", (None, "missing")), (">500", (None, "censored")), ("≤0.01", (None, "censored")),
        ("NB", (None, "invalid")), ("0", (None, "invalid")), ("-3", (None, "invalid")),
        ("nan", (None, "invalid")), ("inf", (None, "invalid")), ("4e-2", (0.04, "finite_positive")),
    ],
)
def test_kd_parsing_preserves_missing_and_censored_values(raw, expected):
    assert external.parse_kd(raw) == expected


def make_source(windows, kd, calls=None, source=external.CONTROLS):
    frame = pd.DataFrame({"HCDR3": ["SR" + window + "Y" for window in windows], external.KD_COLUMN: kd})
    if source == external.CONTROLS:
        frame = frame.assign(HCDR1=external.H1_WT, HCDR2=external.H2_WT, Binder=calls)
    return external.source_rows(frame, source)


def test_exact_overlap_uses_distinct_internal_sources_and_control_membership():
    rows = pd.concat([
        make_source(
            ["AAAAAAAAAA", "CCCCCCCCCC", "DDDDDDDDDD", "WWWWWWWWWW", "VVVVVVVVVV", "TTTTTTTTTT"],
            ["1", "1", "1", "1", "", ""], ["True", "True", "True", "True", "False", "False"],
        ),
        make_source(["VVVVVVVVVV", "GGGGGGGGGG"], ["1", "1"], source=external.SELECTED),
    ], ignore_index=True)
    audited = external.audit_eligibility(rows, {
        "spr_fit": {"AAAAAAAAAA"}, "facs_train": {"CCCCCCCCCC"},
        "facs_val": {"DDDDDDDDDD"}, "facs_test": {"WWWWWWWWWW"}, "spr_panel": {"TTTTTTTTTT"},
    })
    assert audited["eligibility"].tolist() == [
        "internal_exact_overlap", "internal_exact_overlap", "internal_exact_overlap", "internal_exact_overlap",
        "eligible", "internal_exact_overlap", "cross_series_control_match", "eligible",
    ]
    assert audited.loc[6, "matching_control_csv_rows"] == "6"
    assert audited.loc[5, "exact_spr_panel"]
    assert audited.loc[2, "nearest_facs_val_edits"] == 0
    assert audited.loc[2, "nearest_facs_train_edits"] == 10


def test_incompatible_row_stays_visible_after_join():
    rows = make_source(["AAAAAAAAAA", "CCCCCCCCCCX"], ["2", ""], ["True", "False"])
    audited = external.audit_eligibility(rows, {
        name: {"TTTTTTTTTT"} for name in (*external.DISTANCE_REFERENCES, "facs_test", "spr_panel")
    })
    assert audited.loc[1, "eligibility"] == "noncanonical_h3"
    assert pd.isna(audited.loc[1, "nearest_spr_fit_edits"])


def test_duplicate_measurements_preserve_conflicts_and_censoring():
    rows = make_source(
        ["AAAAAAAAAA"] * 4 + ["CCCCCCCCCC"] * 2 + ["DDDDDDDDDD", "EEEEEEEEEE", "FFFFFFFFFF"],
        ["1", "1", "100", ">500", "1", "", "", "3", "30"],
        ["True", "True", "True", "True", "True", "False", "False", "False", "unknown"],
    ).assign(eligibility="eligible")
    windows = external.combine_measurements(rows).set_index("model_cdrh3")
    assert windows.loc["AAAAAAAAAA", "source_csv_rows"] == "2|3|4|5"
    assert windows.loc["AAAAAAAAAA", "distinct_numeric_kd_count"] == 2
    assert windows.loc["AAAAAAAAAA", "evaluated_neg_log10_kd_nm"] == -1
    assert windows.loc["CCCCCCCCCC", "endpoint_conflict"]
    assert pd.isna(windows.loc["CCCCCCCCCC", "published_binder"])
    assert windows.loc["EEEEEEEEEE", "endpoint_conflict"]
    assert not windows.loc["DDDDDDDDDD", "published_binder"]
    assert windows["evaluated_neg_log10_kd_nm"].notna().sum() == 1
    assert pd.isna(windows.loc["FFFFFFFFFF", "published_binder"])


def test_joint_group_bootstrap_for_binder_comparison():
    frame = pd.DataFrame({
        "published_binder": [True, False] * 5,
        "sequence_group": np.repeat(np.arange(5), 2),
        external.RIDGE: [1, 0] * 5,
        external.FACS: [0, 1] * 5,
        external.WT_DISTANCE: [-1, -1] * 5,
    })
    result = external.grouped_summary(frame, "binder", draws=20)
    assert result["groups"] == 5
    assert result["group_sizes_desc"] == [2] * 5
    assert result["average_precision_baseline"] == 0.5
    delta = result["metrics"]["ridge_minus_facs/roc_auc"]
    assert delta == {"point": 1, "valid_draws": 20, "interval_reason": None, "interval_95": [1, 1]}
    assert result["leave_one_group_out_paired"]["ridge_minus_facs/roc_auc"] == {"valid": 5, "range": [1, 1]}
    small = external.grouped_summary(frame.iloc[:8], "binder", draws=20)
    assert small["attempted_bootstrap_draws"] == 0
    assert small["metrics"]["ridge_minus_facs/roc_auc"]["interval_reason"] == "too_few_endpoint_groups"


def test_undefined_bootstrap_draws_are_counted():
    frame = pd.DataFrame({
        "published_binder": [True, False, False, False, False],
        "sequence_group": range(5),
        external.RIDGE: [1, 0, 0, 0, 0],
        external.FACS: [0, 1, 1, 1, 1],
        external.WT_DISTANCE: [-1] * 5,
    })
    result = external.grouped_summary(frame, "binder", draws=100)
    delta = result["metrics"]["ridge_minus_facs/roc_auc"]
    assert delta["point"] == 1
    assert 0 < delta["valid_draws"] < 90
    assert delta["interval_95"] is None
    assert delta["interval_reason"] == "too_few_valid_bootstrap_draws"


def test_kd_direction_and_empty_selected_endpoint():
    frame = pd.DataFrame({
        "evaluated_neg_log10_kd_nm": [-2.0, -1.0, np.nan, 0.0],
        "sequence_group": range(4),
        external.RIDGE: [-1.0, 0.0, 30.0, 1.0],
        external.FACS: [0.0, 1.0, 3.0, 2.0],
        external.WT_DISTANCE: [-1.0] * 4,
    })
    result = external.grouped_summary(frame, "kd", draws=20)
    assert result["n"] == 3
    metrics = result["metrics"]
    assert metrics[f"{external.RIDGE}/spearman"]["point"] == 1
    assert metrics[f"{external.RIDGE}/mean_signed_error"]["point"] == 1
    assert metrics[f"{external.RIDGE}/rmse"]["point"] == 1
    assert metrics[f"{external.WT_DISTANCE}/spearman"]["point"] is None
    empty = external.grouped_summary(frame.iloc[:0], "kd")
    assert empty["n"] == 0
    assert all(item["point"] is None for item in empty["metrics"].values())


def test_changed_external_file_fails_hash_check(tmp_path):
    bad_source = tmp_path / "README.md"
    bad_source.write_text("changed content")
    with pytest.raises(ValueError, match="does not match the frozen SHA-256"):
        external.load_external(tmp_path)


def test_followup_links_original_kd_at_displayed_precision():
    selected = make_source(["AAAAAAAAAA"], ["8.05"], source=external.SELECTED)
    followup = pd.DataFrame({
        "Candidate": ["1"],
        "HCDR3": ["SRAAAAAAAAAAY"],
        external.KD_COLUMN: ["8.1"],
        "Fab KD (nM) (hHER2)": ["7.20"],
        "mAb KD (nM) (hHER2)": ["9.75"],
    })
    audited = external.audit_followup(followup, selected).iloc[0]
    assert audited["selected_rows_matching_h3"] == "2"
    assert not audited["original_kd_matches_any_selected"]
    assert audited["original_kd_matches_selected_at_displayed_precision"]
    assert audited["reported_followup_fab_kd_nm_raw"] == "7.20"
    assert audited["use"] == "provenance_only"
