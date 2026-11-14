"""Compare the frozen HER2 models on published Absci SPR outcomes."""

import argparse
import json
import platform
import re
from decimal import ROUND_HALF_UP, Decimal
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from scipy.stats import kendalltau
from sklearn.metrics import average_precision_score, roc_auc_score

from scripts import rank_spr_candidates as rank
from src import her2_benchmark as facs
from src import her2_spr_benchmark as spr

ABSCI_REVISION = "cd2bd4b6932b5301d435555dc22084d407e4b92b"
EXTERNAL_DIR = Path("data/raw/external-spr")
OUTPUT_DIR = Path("results/external-spr")
PROTOCOL = Path("docs/research/external-spr-protocol.md")
RIDGE_MODEL = Path("results/spr-candidate-ranking/model.json")
SPR_BENCHMARK = Path("results/spr-benchmark/metrics.json")
FACS_MODEL = Path("results/experiment-20261114/selected_model.json")
FACS_BENCHMARK = Path("results/experiment-20261114/benchmark_metrics.json")
IMPLEMENTATION = Path("scripts/benchmark_external_spr.py")

CONTROLS = "spr-controls.csv"
SELECTED = "zero-shot-binders.csv"
FOLLOWUP = "functionality-developability-cross_reactivity-data.csv"
SOURCE_HASHES = {
    "README.md": "698f39c41a8e5a9e27fe0ad869cf59ada8c5daf1b0bea6e94c54ccdd6abe0b0d",
    "LICENSE": "12e4c580339a117464d1ca9417c7c06026d7f535ea5cb55ef0aa482e6c7bb769",
    CONTROLS: "33939f99b489b36e9f893641f3b4013ae4898093875067db31569801b15b2489",
    SELECTED: "c12742bd36ae8881a38cd41e50fd3922cbd885b184531d1151889f5fa610b639",
    FOLLOWUP: "f67f8284ce21f86b56b17f75493b05217d8f603565487692372083f626fa08b5",
}
FROZEN_HASHES = {
    RIDGE_MODEL: "b124a240dbcd9059da47181e0c27e038a176d2a9586d01072cf3308f76dc7697",
    SPR_BENCHMARK: "9668a01ca8dcfdcc8e84b127c4d5f8de1580232090f5efd2925f608ae3c41edf",
    FACS_MODEL: "4067c8e66e703a4db828cec37086f3d157fcad38c172ab265f1aebffe9f3a487",
    FACS_BENCHMARK: "229f4e90912aad73e4923f50f59970a53b922da8d4ce1ef0b8e413773fe36a78",
    PROTOCOL: "e14c2799fc707ae66b684e6d13c68e824b594ecd62263725c0d499e26980e329",
    Path("src/her2_benchmark.py"): "47faea9e341cbfc76a43af13fb16bd694545c9814ac8de5028b4040e20243257",
    Path("src/her2_spr_benchmark.py"): "afdabf4a88357dccba743871c382f8ec1cbc4bb61cfb923c8d987ef63d1b6d94",
}
H1_WT = "GFNIKDTY"
H2_WT = "IYPTNGYT"
H3_WT = "SRWGGDGFYAMDY"
CANONICAL = re.compile(f"[{facs.AMINO_ACIDS}]+")
KD_COLUMN = "KD (nM)"
RIDGE = rank.SPR_SCORE
FACS = facs.SCORE
WT_DISTANCE = "negative_wt_edit_distance"
DISTANCE_REFERENCES = ("spr_fit", "facs_train", "facs_val")
DISTANCE_CUTOFF = 3
BOOTSTRAP_DRAWS = 2000
MIN_GROUPS = 5
MIN_VALID_FRACTION = 0.90
SEED = 0


def verify_hash(record, expected):
    if record["sha256"] != expected:
        raise ValueError(f"{record['source']} does not match the frozen SHA-256")
    return record


def audit_reference(readme):
    """Anchor the source annotation to the exact WT model reference."""
    heavy = re.search(r"Heavy Chain Sequence: `([A-Z]+)`", readme)
    light = re.search(r"Light Chain Sequence: `([A-Z]+)`", readme)
    annotations = dict(re.findall(r"(HCDR[123]): `([A-Z]+)`", readme))
    wt_vh = facs.VH_PREFIX + facs.WT_CDRH3 + facs.VH_SUFFIX
    if (
        not heavy or not light or not heavy[1].startswith(wt_vh) or not light[1].startswith(facs.WT_VL)
        or annotations != {"HCDR1": H1_WT, "HCDR2": H2_WT, "HCDR3": H3_WT}
        or H3_WT != facs.VH_PREFIX[-2:] + facs.WT_CDRH3 + facs.VH_SUFFIX[0]
    ):
        raise ValueError("Absci WT annotation differs from the frozen model reference")
    bounds = {}
    for name, seq in annotations.items():
        if wt_vh.count(seq) != 1:
            raise ValueError(f"Expected a unique WT {name} location in VH")
        offset = wt_vh.index(seq)
        bounds[name] = [offset + 1, offset + len(seq)]
    return {
        "reference_vh_prefix_matches": True,
        "reference_vl_prefix_matches": True,
        "reference_vh_length": len(wt_vh),
        "reference_vl_length": len(facs.WT_VL),
        "wt_annotations": annotations,
        "linear_vh_bounds_1_based_inclusive": {
            **bounds, "model_window": [len(facs.VH_PREFIX) + 1, len(facs.VH_PREFIX) + facs.LOOP_LENGTH],
        },
        "per_row_paired_sequences_supplied": False,
        "per_row_assay_construct_supplied": False,
    }


def map_window(h3, h1=None, h2=None):
    """Project the central ten positions only after the source sequence passes the gate."""
    if h1 is not None and (h1, h2) != (H1_WT, H2_WT):
        return None, "h1_h2_mismatch"
    if not CANONICAL.fullmatch(h3):
        return None, "noncanonical_h3"
    if len(h3) != len(H3_WT):
        return None, "unsupported_h3_length"
    if not (h3.startswith(H3_WT[:2]) and h3.endswith(H3_WT[-1])):
        return None, "non_wt_h3_anchors"
    return h3[2:-1], "compatible_window"


def parse_kd(raw):
    """Return a number only for explicit, positive and finite KD."""
    raw = raw.strip()
    if not raw:
        return None, "missing"
    if raw.startswith(("<", ">", "≤", "≥")):
        return None, "censored"
    try:
        value = float(raw)
    except ValueError:
        return None, "invalid"
    return (value, "finite_positive") if np.isfinite(value) and value > 0 else (None, "invalid")


def source_rows(frame, source):
    """Preserve published rows and distinguish inferred scaffold from measured outcome."""
    controls = source == CONTROLS
    rows = []
    for source_row, record in enumerate(frame.to_dict("records"), start=2):
        h1, h2, h3 = record.get("HCDR1"), record.get("HCDR2"), record["HCDR3"]
        window, status = map_window(h3, h1, h2)
        kd, kd_status = parse_kd(record[KD_COLUMN])
        binder = {"True": True, "False": False}.get(record.get("Binder"))
        if binder is False and kd is not None:
            outcome_status = "false_label_with_finite_kd"
        elif controls and binder is None:
            outcome_status = "unknown_binder_call"
        else:
            outcome_status = "published_binder_call" if controls else "binder_selected_series"
        rows.append({
            "source": source, "source_csv_row": source_row,
            "published_h1": h1, "published_h2": h2, "published_imgt_h3": h3,
            "h1_h2_status": (
                "observed_wt" if (h1, h2) == (H1_WT, H2_WT) else
                "observed_mismatch" if controls else "inferred_wt_from_h3_only_methods"
            ),
            "imgt_h3_length": len(h3),
            "h3_has_wt_anchors": h3.startswith(H3_WT[:2]) and h3.endswith(H3_WT[-1:]),
            "model_cdrh3": window, "sequence_status": status,
            "scaffold_status": "conditional_study_fab" if window is not None else "not_admitted",
            "published_binder_call": binder, "published_kd_nm_raw": record[KD_COLUMN],
            "published_kd_nm": kd, "kd_status": kd_status, "outcome_status": outcome_status,
        })
    return pd.DataFrame(rows)


def load_external(external_dir):
    inputs = {
        name: verify_hash(rank.source(external_dir / name), SOURCE_HASHES[name])
        for name in ("README.md", "LICENSE")
    }
    columns = {
        CONTROLS: ("HCDR1", "HCDR2", "HCDR3", KD_COLUMN, "Binder"),
        SELECTED: ("HCDR3", KD_COLUMN),
        FOLLOWUP: ("Candidate", "HCDR3", KD_COLUMN, "Fab KD (nM) (hHER2)", "mAb KD (nM) (hHER2)"),
    }
    frames = {}
    for name, required in columns.items():
        frames[name], provenance = facs.read_csv(external_dir / name, columns=required)
        inputs[name] = verify_hash(provenance, SOURCE_HASHES[name])
    reference = audit_reference((external_dir / "README.md").read_text(encoding="utf-8"))
    return frames, inputs, reference


def load_frozen(data_dir):
    inputs = {str(path): verify_hash(rank.source(path), sha) for path, sha in FROZEN_HASHES.items()}
    ridge = json.loads(RIDGE_MODEL.read_text())
    benchmark, _ = rank.load_benchmark(SPR_BENCHMARK)
    primary, source_inputs, splits = rank.load_primary(data_dir, benchmark)
    facs_benchmark = json.loads(FACS_BENCHMARK.read_text())
    facs_model = facs.load_model(FACS_MODEL)
    if (
        ridge["model"] != "position_one_hot_ridge" or ridge["target"] != spr.TARGET
        or ridge["alpha"] != 10 or ridge["fit"]["cohort"] != rank.COHORT
        or len(primary) != 158 or ridge["fit"]["fit_ids"] != primary["cdrh3"].tolist()
        or ridge["fit"]["groups"] != primary["group"].nunique()
        or benchmark["inputs"]["facs_model"]["sha256"] != FROZEN_HASHES[FACS_MODEL]
        or facs_benchmark["selection"]["fit_split"] != "train"
        or facs_model["features"] != "pairwise" or facs_benchmark["selected_model"] != "pairwise"
    ):
        raise ValueError("The frozen model metadata does not match the original benchmark fit")
    for name, record in source_inputs.items():
        verify_hash(record, ridge["inputs"][name]["sha256"])
        verify_hash(record, facs_benchmark["inputs"][name]["sha256"])
    assay, provenance = facs.read_csv(data_dir / "spr.csv", columns=("cdrh3",))
    verify_hash(provenance, source_inputs["spr"]["sha256"])
    return ridge, facs_model, {
        "spr_fit": set(primary["cdrh3"]),
        "facs_train": splits["train"], "facs_val": splits["val"], "facs_test": splits["test"],
        "spr_panel": set(facs.normalize_sequences(assay["cdrh3"], "Internal SPR CDRH3")),
    }, {**inputs, **source_inputs}


def audit_distances(windows, references):
    """Report exact nearest training sequences, with lexical tie breaking."""
    indexes = {name: sorted(references[name]) for name in DISTANCE_REFERENCES}
    if any(not sequences for sequences in indexes.values()):
        raise ValueError("Every distance reference must contain its verified input sequences")
    audit = {}
    for window in sorted(set(windows)):
        record = {f"exact_{name}": window in ref for name, ref in references.items()}
        for name, sequences in indexes.items():
            nearest, distance, _ = process.extractOne(window, sequences, scorer=Levenshtein.distance)
            record[f"nearest_{name}_cdrh3"] = nearest
            record[f"nearest_{name}_edits"] = distance
        record["wt_edit_distance"] = Levenshtein.distance(window, facs.WT_CDRH3)
        record["distance_at_least_3_to_training"] = min(
            record[f"nearest_{name}_edits"] for name in DISTANCE_REFERENCES
        ) >= DISTANCE_CUTOFF
        audit[window] = record
    return audit


def eligibility(row, controls):
    if pd.isna(row.model_cdrh3):
        return row.sequence_status
    reasons = []
    if row.model_cdrh3 == facs.WT_CDRH3:
        reasons.append("wild_type")
    if any(getattr(row, f"exact_{name}") for name in ("spr_fit", "facs_train", "facs_val", "facs_test", "spr_panel")):
        reasons.append("internal_exact_overlap")
    if row.source == SELECTED and row.model_cdrh3 in controls:
        reasons.append("cross_series_control_match")
    return "|".join(reasons) or "eligible"


def audit_eligibility(rows, references):
    audit = audit_distances(rows["model_cdrh3"].dropna(), references)
    distances = pd.DataFrame.from_dict(audit, orient="index")
    rows = rows.join(distances, on="model_cdrh3")
    controls = rows.loc[rows["source"].eq(CONTROLS)].groupby("model_cdrh3")["source_csv_row"].agg(
        lambda row_ids: "|".join(map(str, row_ids))
    )
    rows["eligibility"] = [eligibility(row, controls) for row in rows.itertuples()]
    rows["matching_control_csv_rows"] = (
        rows["model_cdrh3"].map(controls).where(rows["source"].eq(SELECTED))
    )
    return rows


def combine_measurements(rows):
    """Give each sequence one endpoint record and retain every source row ID."""
    rows = rows.loc[rows["eligibility"].eq("eligible")]
    combined = []
    for (source, window), group in rows.groupby(["source", "model_cdrh3"], sort=True):
        record = group.iloc[0].drop(labels=[
            "source_csv_row", "published_binder_call", "published_kd_nm_raw", "published_kd_nm",
            "kd_status", "outcome_status",
        ]).to_dict()
        calls = set(group["published_binder_call"].dropna())
        conflict = len(calls) > 1 or group["outcome_status"].eq("false_label_with_finite_kd").any()
        binder = bool(next(iter(calls))) if len(calls) == 1 and not conflict else None
        use_kd = source == SELECTED or binder is True
        usable = group if source == SELECTED else group.loc[group["published_binder_call"].eq(True)]
        values = sorted(usable["published_kd_nm"].dropna().unique()) if use_kd else []
        target = float(np.median(-np.log10(values))) if values else None
        record.update({
            "source_csv_rows": "|".join(map(str, group["source_csv_row"])),
            "source_row_count": len(group), "distinct_numeric_kd_count": len(values),
            "published_binder": binder, "endpoint_conflict": bool(conflict),
            "evaluated_neg_log10_kd_nm": target,
        })
        combined.append(record)
    if not combined:
        # Preserve columns needed for honest empty endpoint reports.
        return pd.DataFrame(columns=[
            "source", "model_cdrh3", "published_binder", "endpoint_conflict",
            "evaluated_neg_log10_kd_nm", "source_row_count", "distance_at_least_3_to_training",
            "wt_edit_distance",
        ])
    return pd.DataFrame(combined)


def score_windows(windows, ridge, facs_model):
    scores = windows.copy()
    if scores.empty:
        return scores.assign(sequence_group=None, **dict.fromkeys([RIDGE, FACS, WT_DISTANCE]))
    ids = scores["model_cdrh3"]
    if not ids.is_unique:
        raise ValueError("Expected disjoint model windows after cross-series exclusion")
    scores["sequence_group"] = spr.sequence_groups(ids)
    scores[RIDGE] = rank.finite_array(rank.predict_from_parameters(ridge, ids), "Frozen ridge score")
    scores[FACS] = rank.finite_array(
        facs.predict_scores(facs_model["estimator"], facs_model["features"], ids), "Frozen FACS score",
    )
    scores[WT_DISTANCE] = -scores["wt_edit_distance"]
    return scores


def audit_followup(followup, selected):
    by_h3 = {
        h3: frame for h3, frame in selected.groupby("published_imgt_h3", sort=False)
    }
    matches = []
    for source_row, entry in enumerate(followup.to_dict("records"), start=2):
        h3 = entry["HCDR3"]
        original, original_status = parse_kd(entry[KD_COLUMN])
        linked = by_h3.get(h3)
        matching_rows = [] if linked is None else linked["source_csv_row"].tolist()
        matching_values = set() if linked is None else set(linked["published_kd_nm"].dropna())
        displayed_match = False
        if original is not None:
            original_decimal = Decimal(entry[KD_COLUMN])
            displayed_match = any(
                Decimal(str(value)).quantize(original_decimal, rounding=ROUND_HALF_UP) == original_decimal
                for value in matching_values
            )
        window, sequence_status = map_window(h3)
        matches.append({
            "source_csv_row": source_row, "candidate": entry["Candidate"], "published_imgt_h3": h3,
            "is_wt": h3 == H3_WT, "model_cdrh3": window, "sequence_status": sequence_status,
            "reported_original_kd_nm_raw": entry[KD_COLUMN], "reported_original_kd_nm": original,
            "original_kd_status": original_status,
            "reported_followup_fab_kd_nm_raw": entry["Fab KD (nM) (hHER2)"],
            "reported_followup_mab_kd_nm_raw": entry["mAb KD (nM) (hHER2)"],
            "selected_rows_matching_h3": "|".join(map(str, matching_rows)),
            "selected_kd_nm_values": "|".join(map(str, sorted(matching_values))),
            "original_kd_matches_any_selected": original is not None and original in matching_values,
            "original_kd_matches_selected_at_displayed_precision": displayed_match,
            "use": "provenance_only",
        })
    return pd.DataFrame(matches)


def score_metrics(labels, scores, endpoint):
    metrics = {}
    keys = ("roc_auc", "average_precision") if endpoint == "binder" else ("spearman", "kendall_tau_b")
    for model, values in scores.items():
        first = second = None
        if len(labels) >= 2 and np.ptp(labels) > 0:
            if endpoint == "binder":
                first = float(roc_auc_score(labels, values))
                second = float(average_precision_score(labels, values))
            elif np.ptp(values) > 0:
                first = spr.spearman(labels, values)
                second = float(kendalltau(labels, values, variant="b").statistic)
        metrics.update({f"{model}/{key}": value for key, value in zip(keys, (first, second))})
    for key in keys:
        ridge, facs_score = metrics[f"{RIDGE}/{key}"], metrics[f"{FACS}/{key}"]
        metrics[f"ridge_minus_facs/{key}"] = None if ridge is None or facs_score is None else ridge - facs_score
    if endpoint == "kd":
        errors = scores[RIDGE] - labels
        metrics.update({
            f"{RIDGE}/mae": float(np.mean(np.abs(errors))) if len(errors) else None,
            f"{RIDGE}/rmse": spr.rmse(labels, scores[RIDGE]) if len(errors) else None,
            f"{RIDGE}/mean_signed_error": float(np.mean(errors)) if len(errors) else None,
        })
    return metrics


def grouped_summary(frame, endpoint, *, draws=BOOTSTRAP_DRAWS, seed=SEED):
    label = "published_binder" if endpoint == "binder" else "evaluated_neg_log10_kd_nm"
    frame = frame.loc[frame[label].notna()]
    labels = frame[label].to_numpy(dtype=int if endpoint == "binder" else float)
    scores = {model: frame[model].to_numpy(dtype=float) for model in (RIDGE, FACS, WT_DISTANCE)}
    groups = frame["sequence_group"].to_numpy()
    members = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    point = score_metrics(labels, scores, endpoint)

    def evaluate(indices):
        return score_metrics(labels[indices], {key: values[indices] for key, values in scores.items()}, endpoint)

    samples = {key: [] for key in point}
    attempted = draws if len(members) >= MIN_GROUPS else 0
    rng = np.random.default_rng(seed)
    for _ in range(attempted):
        indices = np.concatenate([members[i] for i in rng.integers(len(members), size=len(members))])
        for key, value in evaluate(indices).items():
            if value is not None:
                samples[key].append(value)
    intervals = {}
    for key, values in samples.items():
        if point[key] is None:
            reason = "undefined_point_estimate"
        elif not attempted:
            reason = "too_few_endpoint_groups"
        elif len(values) < MIN_VALID_FRACTION * attempted:
            reason = "too_few_valid_bootstrap_draws"
        else:
            reason = None
        intervals[key] = {
            "point": point[key], "valid_draws": len(values), "interval_reason": reason,
            "interval_95": None if reason else np.quantile(values, [0.025, 0.975]).tolist(),
        }

    paired_keys = [key for key in point if key.startswith("ridge_minus_facs/")]
    leave_one_out = {key: [] for key in paired_keys}
    for group in np.unique(groups):
        metrics = evaluate(np.flatnonzero(groups != group))
        for key in paired_keys:
            if metrics[key] is not None:
                leave_one_out[key].append(metrics[key])
    return {
        "n": len(frame), "groups": len(members), "group_sizes_desc": sorted(map(len, members), reverse=True),
        "positive": int(labels.sum()) if endpoint == "binder" else None,
        "average_precision_baseline": float(labels.mean()) if endpoint == "binder" and len(labels) else None,
        "kd_nm_range": [float(10 ** -max(labels)), float(10 ** -min(labels))] if endpoint == "kd" and len(labels) else None,
        "attempted_bootstrap_draws": attempted, "metrics": intervals,
        "leave_one_group_out_paired": {
            key: {"valid": len(values), "range": [min(values), max(values)] if values else None}
            for key, values in leave_one_out.items()
        },
    }


def endpoint_summaries(scores):
    controls = scores.loc[scores["source"].eq(CONTROLS)]
    selected = scores.loc[scores["source"].eq(SELECTED)]
    sections = {}
    for name, frame, endpoint in (
        ("controls_binder", controls, "binder"),
        ("controls_kd", controls, "kd"),
        ("selected_kd", selected, "kd"),
    ):
        far = frame.loc[frame["distance_at_least_3_to_training"].eq(True)]
        sections[name] = {"all_eligible": grouped_summary(frame, endpoint), "distance_at_least_3": grouped_summary(far, endpoint)}
    return sections


def count_values(series):
    return {str(key): int(value) for key, value in series.value_counts(dropna=False).items()}


def summarize_audit(rows, scores, followup, references):
    cohorts = {}
    for name, data in rows.groupby("source", sort=True):
        windows = scores.loc[scores["source"].eq(name)]
        compatible = data.loc[data["sequence_status"].eq("compatible_window")]
        cohorts[name] = {
            "published_rows": len(data), "published_unique_h3": int(data["published_imgt_h3"].nunique()),
            "published_binder_calls": count_values(data["published_binder_call"].fillna("unavailable")),
            "published_binder_by_kd_status": {
                str(call): count_values(group["kd_status"])
                for call, group in data.groupby(data["published_binder_call"].fillna("unavailable"), sort=False)
            },
            "published_outcome_status": count_values(data["outcome_status"]),
            "sequence_status": count_values(data["sequence_status"]),
            "eligibility": count_values(data["eligibility"]),
            "compatible_wt_rows": int(compatible["model_cdrh3"].eq(facs.WT_CDRH3).sum()),
            "compatible_exact_internal_matches_by_panel": {
                key: int(compatible[f"exact_{key}"].eq(True).sum()) if len(compatible) else 0
                for key in references
            },
            "eligible_published_rows": int(data["eligibility"].eq("eligible").sum()),
            "scored_windows": len(windows),
            "collapsed_eligible_rows": int(windows["source_row_count"].sum()) - len(windows),
            "conflicting_scored_windows": int(windows["endpoint_conflict"].sum()),
            "scored_windows_by_distance": {
                key: count_values(windows[f"nearest_{key}_edits"].astype(int)) if len(windows) else {}
                for key in DISTANCE_REFERENCES
            },
        }
    return {
        "internal_reference_unique_sequences": {name: len(ids) for name, ids in references.items()},
        "cohorts": cohorts,
        "followup": {
            "rows": len(followup),
            "rows_matching_any_selected_h3": int(followup["selected_rows_matching_h3"].ne("").sum()),
            "rows_matching_any_selected_h3_and_original_kd": int(followup["original_kd_matches_any_selected"].sum()),
            "rows_matching_selected_at_displayed_kd_precision": int(
                followup["original_kd_matches_selected_at_displayed_precision"].sum()
            ),
            "scored_rows": 0,
        },
    }


def run(data_dir, external_dir, output):
    ridge, facs_model, references, internal_inputs = load_frozen(data_dir)
    external, source_inputs, reference_audit = load_external(external_dir)
    rows = pd.concat([source_rows(external[name], name) for name in (CONTROLS, SELECTED)], ignore_index=True)
    rows = audit_eligibility(rows, references)
    followup = audit_followup(external[FOLLOWUP], rows.loc[rows["source"].eq(SELECTED)])
    scores = score_windows(combine_measurements(rows), ridge, facs_model)
    result = {
        "scope": "retrospective_conditional_trastuzumab_fab_transfer",
        "external_source": {
            "attribution": "Absci Corporation (2023)", "revision": ABSCI_REVISION, "inputs": source_inputs,
        },
        "internal_inputs": internal_inputs,
        "implementation": rank.source(IMPLEMENTATION),
        "shared_ranking_code": rank.source(Path("scripts/rank_spr_candidates.py")),
        "versions": {"python": platform.python_version(), **{
            package: version(package) for package in ("numpy", "pandas", "rapidfuzz", "scipy", "scikit-learn")
        }},
        "reference_audit": reference_audit,
        "frozen_models": {
            RIDGE: {"alpha": ridge["alpha"], "n_fit": len(references["spr_fit"]), "target": spr.TARGET},
            FACS: {"features": facs_model["features"], "n_fit": len(references["facs_train"]), "target": "FACS high versus mid or low"},
        },
        "uncertainty": {
            "seed": SEED, "draws": BOOTSTRAP_DRAWS, "sequence_group_levenshtein_cutoff": spr.MAX_DISTANCE,
            "minimum_groups": MIN_GROUPS, "minimum_valid_fraction": MIN_VALID_FRACTION,
            "training_distance_sensitivity_cutoff": DISTANCE_CUTOFF,
        },
        "audit": summarize_audit(rows, scores, followup, references),
        "endpoints": endpoint_summaries(scores),
    }
    output.mkdir(parents=True, exist_ok=True)
    for name, frame in (("source_rows.csv", rows), ("window_scores.csv", scores), ("followup_audit.csv", followup)):
        rank.write_ranked_csv(output / name, frame)
    (output / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True, help="Original SPR and FACS split CSV directory")
    parser.add_argument("--external-dir", type=Path, default=EXTERNAL_DIR)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    run(args.data_dir, args.external_dir, args.output)


if __name__ == "__main__":
    main()
