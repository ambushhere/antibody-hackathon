"""Evaluate whether SPR labels improve trastuzumab CDRH3 affinity ranking."""

import argparse
import hashlib
import json
import platform
from importlib.metadata import version
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

from src import her2_benchmark as facs

TARGET = "neg_log10_kd_nm"
SPR_COLUMNS = ("cdrh3", "vh_sequence", "vl_sequence", "kd_nM", "binding_result")
ALPHAS = (0.1, 1.0, 10.0, 100.0)
MAX_FOLDS = 5
MAX_DISTANCE = 2
TOP_K = (5, 10)
MODEL_NAMES = ("training_mean", facs.SCORE, "facs_scalar", "position_ridge", "facs_residual_ridge")
MODEL_LABELS = dict(zip(MODEL_NAMES, (
    "Training mean", "Frozen FACS probability", "FACS scalar calibration",
    "Position ridge", "FACS residual ridge",
)))
COMPARISONS = {
    "position_ridge": (facs.SCORE,),
    "facs_residual_ridge": (facs.SCORE, "facs_scalar"),
}
SEED = 0
BOOTSTRAP_DRAWS = 2000


def prepare_spr(frame, split_sequences):
    """Filter rows before taking the median transformed KD for each sequence."""
    frame = frame.copy()
    for column in ("cdrh3", "vh_sequence", "vl_sequence"):
        frame[column] = facs.normalize_sequences(frame[column], f"SPR {column}")
    status = frame["binding_result"].str.strip().replace("", "(blank)")
    kd_text = frame["kd_nM"].str.strip()
    kd = pd.to_numeric(kd_text, errors="coerce")
    valid_kd = np.isfinite(kd) & kd.gt(0)
    kd_state = np.select([kd_text.eq(""), valid_kd], ["missing", "finite_positive"], "invalid")
    status_counts = pd.crosstab(status, kd_state).reindex(
        columns=["missing", "invalid", "finite_positive"], fill_value=0
    )
    quantified = status.eq("quantified_binding") & valid_kd
    correct_background = frame["vl_sequence"].eq(facs.WT_VL) & frame["vh_sequence"].eq(
        facs.VH_PREFIX + frame["cdrh3"] + facs.VH_SUFFIX
    )
    in_background = quantified & correct_background
    is_wt = frame["cdrh3"].eq(facs.WT_CDRH3)
    references = kd[in_background & is_wt].sort_values().tolist()
    is_variant = in_background & ~is_wt
    supported_length = frame["cdrh3"].str.len().eq(facs.LOOP_LENGTH)
    eligible = is_variant & supported_length
    variants = (
        frame.loc[eligible, ["cdrh3"]]
        .assign(**{TARGET: -np.log10(kd[eligible])})
        .groupby("cdrh3", as_index=False)
        .agg(**{TARGET: (TARGET, "median"), "quantified_measurements": (TARGET, "size")})
    )
    for name, sequences in split_sequences.items():
        variants[f"in_facs_{name}"] = variants["cdrh3"].isin(sequences)
    audit = {
        "input_rows": len(frame),
        "kd_state_by_binding_result": status_counts.to_dict(orient="index"),
        "quantified_rows_with_finite_positive_kd": int(quantified.sum()),
        "excluded_other_background_rows": int((quantified & ~correct_background).sum()),
        "excluded_wild_type_rows": len(references),
        "excluded_other_length_rows_by_length": (
            frame.loc[is_variant & ~supported_length, "cdrh3"].str.len().value_counts().to_dict()
        ),
        "eligible_quantified_variant_rows": int(eligible.sum()),
        "eligible_unique_cdrh3": len(variants),
        "replicate_rows_collapsed": int(eligible.sum()) - len(variants),
        "overlaps_by_facs_split": {
            name: int(variants[f"in_facs_{name}"].sum()) for name in split_sequences
        },
        "wt_reference_kd_nm": references,
    }
    return variants, audit


def sequence_groups(sequences):
    """Connect every pair within the edit distance cutoff, including transitive links."""
    sequences = list(sequences)
    left, right = [], []
    for i, j in combinations(range(len(sequences)), 2):
        if Levenshtein.distance(sequences[i], sequences[j], score_cutoff=MAX_DISTANCE) <= MAX_DISTANCE:
            left.append(i)
            right.append(j)
    graph = csr_matrix(
        (np.ones(len(left), dtype=bool), (left, right)), shape=(len(sequences), len(sequences))
    )
    return connected_components(graph, directed=False)[1]


def group_splits(groups):
    return list(GroupKFold(n_splits=min(MAX_FOLDS, len(np.unique(groups)))).split(
        groups, groups=groups
    ))


def rmse(labels, predictions):
    return float(np.sqrt(np.mean((labels - predictions) ** 2)))


def calibrate(train_scores, labels, test_scores):
    return LinearRegression().fit(train_scores, labels).predict(test_scores)


def ridge_predictions(positions, labels, facs_scores, groups, test_positions, test_scores):
    """Tune only alpha, using inputs scoped to a single outer training fold."""
    observed = []
    inner_predictions = {
        model: np.empty((len(ALPHAS), len(labels)))
        for model in ("position_ridge", "facs_residual_ridge")
    }
    for train, val in group_splits(groups):
        baseline = calibrate(facs_scores[train], labels[train], facs_scores[np.r_[train, val]])
        train_baseline, val_baseline = np.split(baseline, [len(train)])
        for index, alpha in enumerate(ALPHAS):
            ridge = Ridge(alpha=alpha, solver="cholesky")
            inner_predictions["position_ridge"][index, val] = ridge.fit(
                positions[train], labels[train]
            ).predict(positions[val])
            inner_predictions["facs_residual_ridge"][index, val] = val_baseline + ridge.fit(
                positions[train], labels[train] - train_baseline
            ).predict(positions[val])
        observed.append({
            "train_n": len(train), "val_n": len(val),
            "val_groups": np.unique(groups[val]).tolist(),
        })

    baseline = calibrate(facs_scores, labels, np.concatenate((facs_scores, test_scores)))
    train_baseline, test_baseline = np.split(baseline, [len(labels)])
    predictions = {"facs_scalar": test_baseline}
    tuning = {}
    for model, fit_labels in (
        ("position_ridge", labels),
        ("facs_residual_ridge", labels - train_baseline),
    ):
        errors = [rmse(labels, predictions) for predictions in inner_predictions[model]]
        # Prefer stronger regularization when the inner errors tie.
        best = min(range(len(ALPHAS)), key=lambda index: (errors[index], -ALPHAS[index]))
        predictions[model] = Ridge(alpha=ALPHAS[best], solver="cholesky").fit(
            positions, fit_labels
        ).predict(test_positions)
        tuning[model] = {"alpha": ALPHAS[best], "inner_pooled_rmse_by_alpha": dict(zip(ALPHAS, errors))}
    predictions["facs_residual_ridge"] += test_baseline
    return predictions, tuning, observed


def cross_validate(variants):
    """Return predictions aligned to variants and enough fold detail to audit grouping."""
    oof = variants.copy()
    groups = sequence_groups(variants["cdrh3"])
    oof["group"] = groups
    if len(np.unique(groups)) < 3:
        return oof, [], "Nested GroupKFold needs at least three nonempty sequence groups."

    positions = facs.sequence_features(variants["cdrh3"], "position").toarray().astype(float)
    facs_scores = variants[[facs.SCORE]].to_numpy()
    labels = variants[TARGET].to_numpy()
    folds = []
    fold_predictions = np.full((len(variants), 4), np.nan)
    fold_ids = np.zeros(len(variants), dtype=int)
    with threadpool_limits(limits=1):
        for fold_id, (train, test) in enumerate(group_splits(groups), start=1):
            predictions, tuning, inner = ridge_predictions(
                positions[train], labels[train], facs_scores[train], groups[train],
                positions[test], facs_scores[test],
            )
            fold_ids[test] = fold_id
            fold_predictions[test, 0] = labels[train].mean()
            for index, model in enumerate(MODEL_NAMES[2:], start=1):
                fold_predictions[test, index] = predictions[model]
            folds.append({
                "fold": fold_id, "train_n": len(train), "test_n": len(test),
                "train_groups": len(np.unique(groups[train])),
                "test_groups": np.unique(groups[test]).tolist(),
                "inner_folds": inner, "tuning": tuning,
            })
    if (fold_ids == 0).any() or not np.isfinite(fold_predictions).all():
        raise ValueError("Missing or invalid out-of-fold predictions")
    oof["fold"] = fold_ids
    oof[["training_mean", *MODEL_NAMES[2:]]] = fold_predictions
    return oof, folds, None


def spearman(labels, scores):
    if len(labels) < 2 or np.ptp(labels) == 0 or np.ptp(scores) == 0:
        return None
    return float(spearmanr(labels, scores).statistic)


def enrichment(positives, scores, k):
    k = min(k, len(scores))
    if not k:
        return {"k": k, "precision": None, "enrichment": None}
    cutoff = np.partition(scores, len(scores) - k)[-k]
    above, tied = scores > cutoff, scores == cutoff
    precision = float((positives[above].sum() + (k - above.sum()) * positives[tied].mean()) / k)
    prevalence = positives.mean()
    return {"k": k, "precision": precision, "enrichment": precision / prevalence if prevalence else None}


def summarize(oof, reference_kds, seed):
    labels = oof[TARGET].to_numpy()
    scores = {model: oof[model].to_numpy() for model in MODEL_NAMES if model in oof}
    thresholds = {}
    if reference_kds:
        thresholds = {
            "better_than_best_wt": labels > -np.log10(min(reference_kds)),
            "at_or_better_than_weakest_wt": labels >= -np.log10(max(reference_kds)),
        }
    metrics = {}
    for model, values in scores.items():
        metrics[model] = {
            "spearman": spearman(labels, values),
            "rmse": rmse(labels, values) if model != facs.SCORE and len(values) else None,
            "top_k": {
                threshold: {k: enrichment(positives, values, k) for k in TOP_K}
                for threshold, positives in thresholds.items()
            },
        }
    intervals, deltas = group_bootstrap(labels, scores, thresholds, oof["group"].to_numpy(), seed)
    for model in metrics:
        metrics[model]["intervals"] = intervals.get(model, {})
    return metrics, {name: int(flags.sum()) for name, flags in thresholds.items()}, deltas


def group_bootstrap(labels, scores, thresholds, groups, seed):
    """Resample intact groups with replacement and use the same draw for all models."""
    members = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    if len(members) < 2:
        return {}, {}
    rng = np.random.default_rng(seed)
    results = {model: {"spearman": [], "rmse": []} for model in scores}
    top_k_draws = {
        model: {name: {k: [] for k in TOP_K} for name in thresholds} for model in scores
    }
    pairs = [
        (model, baseline) for model, baselines in COMPARISONS.items() for baseline in baselines
        if model in scores and baseline in scores
    ]
    differences = {f"{model} minus {baseline}": [] for model, baseline in pairs}
    for _ in range(BOOTSTRAP_DRAWS):
        draw = rng.integers(len(members), size=len(members))
        rows = np.concatenate([members[group] for group in draw])
        drawn_positives = {name: flags[rows] for name, flags in thresholds.items()}
        rho = {}
        for model, values in scores.items():
            drawn_scores = values[rows]
            rho[model] = spearman(labels[rows], drawn_scores)
            results[model]["spearman"].append(rho[model])
            if model != facs.SCORE:
                results[model]["rmse"].append(rmse(labels[rows], drawn_scores))
            for name, positives in drawn_positives.items():
                if positives.any():
                    for k in TOP_K:
                        value = enrichment(positives, drawn_scores, k)["enrichment"]
                        top_k_draws[model][name][k].append(value)
        for model, baseline in pairs:
            if rho[model] is not None and rho[baseline] is not None:
                differences[f"{model} minus {baseline}"].append(rho[model] - rho[baseline])

    def interval(values):
        defined = [value for value in values if value is not None]
        return {
            "95_percent": np.percentile(defined, [2.5, 97.5]).tolist() if defined else None,
            "valid_draws": len(defined),
        }

    intervals = {
        model: {
            **{metric: interval(values) for metric, values in model_results.items() if values},
            "top_k_enrichment": {
                name: {k: interval(values) for k, values in top_k.items()}
                for name, top_k in top_k_draws[model].items()
            },
        }
        for model, model_results in results.items()
    }
    deltas = {}
    for model, baseline in pairs:
        key = f"{model} minus {baseline}"
        left, right = spearman(labels, scores[model]), spearman(labels, scores[baseline])
        deltas[key] = {
            "spearman_delta": left - right if left is not None and right is not None else None,
            **interval(differences[key]),
        }
    return intervals, deltas


def render_report(result):
    audit = result["spr_audit"]
    references = audit["wt_reference_kd_nm"]
    lines = [
        "# HER2 SPR affinity benchmark", "",
        "The target is −log10(KD in nM), aggregated by the median across quantified replicates "
        "of each CDRH3. Rows must have finite positive KD, the exact trastuzumab VH framework "
        "and VL and a ten-residue loop. Wild-type reference rows only define the enrichment thresholds.",
        "", "## Assay audit", "",
        "| Binding result | Missing KD | Invalid KD | Finite positive KD |",
        "| --- | ---: | ---: | ---: |",
    ]
    for status, counts in audit["kd_state_by_binding_result"].items():
        lines.append(f"| {status} | {counts['missing']} | {counts['invalid']} | {counts['finite_positive']} |")
    lengths = audit["excluded_other_length_rows_by_length"]
    length_summary = ", ".join(f"{count} of length {length}" for length, count in lengths.items()) or "zero"
    lines.extend([
        "",
        f"From {audit['quantified_rows_with_finite_positive_kd']} quantified rows, exclusions were "
        f"{audit['excluded_other_background_rows']} other frameworks, "
        f"{audit['excluded_wild_type_rows']} wild-type references and "
        f"{sum(lengths.values())} unsupported lengths ({length_summary}). "
        f"Aggregation collapsed {audit['replicate_rows_collapsed']} variant replicate rows.",
        "",
        "Training and validation FACS overlaps are excluded from the primary cohort by exact CDRH3 "
        "match. FACS was fit on train and selected on validation. The broader cohort includes those "
        "overlaps. The primary cohort can still contain near neighbors of the FACS training data.",
        "",
    ])
    if references:
        lines.extend([
            f"The {len(references)} wild-type measurements span {min(references):.8g} to "
            f"{max(references):.8g} nM. Top-k enrichment uses two thresholds, strictly below "
            "the best wild type and at or below the weakest wild type. The observed reference range "
            "is descriptive and supplies no statistical equivalence margin. All results condition "
            "on reportable KD. The missing statuses supply no quantitative labels.", "",
        ])
    else:
        lines.extend(["No quantified, matching-framework wild-type reference was available.", ""])
    lines.extend([
        "Each cohort is grouped separately by connected components at Levenshtein distance ≤2. "
        "Outer GroupKFold produces pooled out-of-fold predictions. In each outer training set, "
        "inner GroupKFold minimizes pooled RMSE over ridge alpha "
        f"{', '.join(map(str, ALPHAS))}. Both splitters use at most {MAX_FOLDS} folds.",
        "",
        "Position ridge uses a fixed amino-acid one-hot encoding. The scalar baseline fits a linear "
        "calibration of the frozen FACS probability against SPR labels. Residual ridge adds a "
        "one-hot ridge prediction of that calibration's residuals. Both calibration and residual "
        "fitting are repeated inside each training fold. The FACS probability is evaluated by rank.",
        "",
        f"Intervals are 95% percentile intervals from {BOOTSTRAP_DRAWS:,} draws of whole sequence "
        f"groups with replacement (seed {result['seed']}). They summarize uncertainty in the "
        "fixed out-of-fold predictions. Paired differences use the same draws. The mean baseline "
        "ties within each fold and its pooled rank can reflect differences between fold means.",
    ])
    for name, cohort in result["cohorts"].items():
        sizes = ", ".join(f"{size} × {count}" for size, count in cohort["group_size_counts"].items())
        overlaps = ", ".join(f"{split} {count}" for split, count in cohort["facs_overlaps"].items())
        lines.extend([
            "", f"## {name}", "",
            f"{cohort['n']} unique CDRH3s in {cohort['groups']} groups. "
            f"Group size × count is {sizes or 'empty'}. FACS overlaps are {overlaps}.", "",
        ])
        if cohort["split_issue"]:
            lines.extend([cohort["split_issue"], ""])
        else:
            lines.extend([
                "| Outer fold | Train sequences / groups | Test sequences / groups | Inner folds | "
                "Position alpha | Residual alpha |",
                "| ---: | ---: | ---: | ---: | ---: | ---: |",
            ])
            for fold in cohort["folds"]:
                tuning = fold["tuning"]
                lines.append(
                    f"| {fold['fold']} | {fold['train_n']} / {fold['train_groups']} | "
                    f"{fold['test_n']} / {len(fold['test_groups'])} | {len(fold['inner_folds'])} | "
                    f"{tuning['position_ridge']['alpha']:g} | {tuning['facs_residual_ridge']['alpha']:g} |"
                )
            lines.append("")

        def format_number(value):
            return "—" if value is None else f"{value:.3f}"

        def format_interval(interval):
            endpoints = interval.get("95_percent")
            return "—" if endpoints is None else f"[{endpoints[0]:.3f}, {endpoints[1]:.3f}]"

        lines.extend([
            "| Model | Spearman ρ | 95% interval | RMSE | 95% interval |",
            "| --- | ---: | --- | ---: | --- |",
        ])
        for model, row in cohort["metrics"].items():
            intervals = row["intervals"]
            lines.append(
                f"| {MODEL_LABELS[model]} | {format_number(row['spearman'])} | "
                f"{format_interval(intervals.get('spearman', {}))} | {format_number(row['rmse'])} | "
                f"{format_interval(intervals.get('rmse', {}))} |"
            )
        for pair, delta in cohort["paired_spearman_deltas"].items():
            model, baseline = pair.split(" minus ")
            lines.extend([
                "",
                f"- {MODEL_LABELS[model]} minus {MODEL_LABELS[baseline].lower()} "
                f"Δρ {format_number(delta['spearman_delta'])} "
                f"{format_interval(delta)} ({delta['valid_draws']:,} defined draws).",
            ])
        sensitivity = cohort["without_largest_group"]
        if sensitivity and not cohort["split_issue"]:
            rho = sensitivity["spearman"]
            lines.extend([
                "", f"Removing the largest group from evaluation leaves {sensitivity['n']} "
                f"predictions. Spearman ρ is {format_number(rho['position_ridge'])} for position "
                f"ridge, {format_number(rho['facs_residual_ridge'])} for residual ridge and "
                f"{format_number(rho[facs.SCORE])} for frozen FACS. This check uses the original "
                "out-of-fold predictions.",
            ])

        label_counts = cohort["reference_label_counts"]
        if label_counts:
            strict, permissive = label_counts
            lines.extend([
                "", f"{label_counts[strict]}/{cohort['n']} variants beat the best wild type and "
                f"{label_counts[permissive]}/{cohort['n']} match or beat the weakest.", "",
                "| Model | Top-k | Below best WT precision / enrichment | "
                "At or below weakest WT precision / enrichment |",
                "| --- | ---: | ---: | ---: |",
            ])
            for model, row in cohort["metrics"].items():
                for k in TOP_K:
                    strict_top = row["top_k"][strict][k]
                    permissive_top = row["top_k"][permissive][k]
                    lines.append(
                        f"| {MODEL_LABELS[model]} | {strict_top['k']} | "
                        f"{format_number(strict_top['precision'])} / "
                        f"{format_number(strict_top['enrichment'])} | "
                        f"{format_number(permissive_top['precision'])} / "
                        f"{format_number(permissive_top['enrichment'])} |"
                    )
    lines.extend([
        "", "Top-k ties contribute expected precision. Enrichment uses each cohort's prevalence. "
        "With very few positives these values are descriptive. The JSON includes conditional "
        "enrichment intervals and valid draw counts, since draws without positives have undefined "
        "enrichment.",
        "", "## Provenance", "",
        *[f"- {name} `{value}`" for name, value in result["versions"].items()],
        *[f"- {name} `{entry['sha256']}` (`{entry['source']}`)"
          for name, entry in result["inputs"].items()],
    ])
    return "\n".join(lines) + "\n"


def run(data_dir, model_path, output):
    splits, inputs = {}, {}
    for name in facs.SPLIT_NAMES:
        splits[name], inputs[name] = facs.load_split(data_dir / f"{name}.csv")
    facs.audit_splits(splits)
    inputs["facs"], _ = facs.audit_facs(data_dir / "facs.csv", splits)
    frame, inputs["spr"] = facs.read_csv(data_dir / "spr.csv", columns=SPR_COLUMNS)
    inputs["facs_model"] = {
        "source": str(model_path),
        "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
    }
    model = facs.load_model(model_path)
    variants, audit = prepare_spr(frame, {name: set(split["seq"]) for name, split in splits.items()})
    variants[facs.SCORE] = (
        facs.predict_scores(model["estimator"], model["features"], variants["cdrh3"])
        if len(variants) else np.empty(0)
    )
    result = {
        "seed": SEED, "bootstrap_draws": BOOTSTRAP_DRAWS, "alpha_grid": ALPHAS,
        "max_folds": MAX_FOLDS, "levenshtein_cutoff": MAX_DISTANCE,
        "versions": {"python": platform.python_version(), **{
            package: version(package) for package in (
                "numpy", "pandas", "scikit-learn", "scipy", "rapidfuzz", "threadpoolctl"
            )
        }},
        "inputs": inputs, "facs_features": model["features"], "spr_audit": audit, "cohorts": {},
    }
    novel = variants[~(variants["in_facs_train"] | variants["in_facs_val"])]
    predictions = []
    for name, cohort in (("Primary novelty cohort", novel), ("Broader cohort with FACS overlaps", variants)):
        oof, folds, split_issue = cross_validate(cohort)
        metrics, label_counts, deltas = summarize(oof, audit["wt_reference_kd_nm"], SEED)
        sizes = oof["group"].value_counts()
        sensitivity = None
        if len(sizes) > 1:
            largest = sizes[sizes.eq(sizes.max())].index.min()
            rest = oof[oof["group"].ne(largest)]
            sensitivity = {
                "excluded_group": int(largest), "n": len(rest),
                "spearman": {
                    model: spearman(rest[TARGET].to_numpy(), rest[model].to_numpy())
                    for model in metrics
                },
            }
        result["cohorts"][name] = {
            "n": len(oof), "groups": len(sizes),
            "group_size_counts": sizes.value_counts().sort_index().to_dict(),
            "facs_overlaps": {split: int(oof[f"in_facs_{split}"].sum()) for split in splits},
            "folds": folds, "split_issue": split_issue, "metrics": metrics,
            "reference_label_counts": label_counts, "paired_spearman_deltas": deltas,
            "without_largest_group": sensitivity,
        }
        predictions.append(oof.assign(cohort=name))
    output.mkdir(parents=True, exist_ok=True)
    pd.concat(predictions, ignore_index=True).to_csv(output / "oof_predictions.csv", index=False)
    (output / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    report = output / "report.md"
    report.write_text(render_report(result))
    print(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--facs-model", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.data_dir, args.facs_model, args.output)


if __name__ == "__main__":
    main()
