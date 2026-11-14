"""Apply the completed HER2 SPR benchmark to the existing strict optimisation pool."""

import argparse
import hashlib
import json
import platform
import re
import shlex
from importlib.metadata import version
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits

from src import her2_benchmark as facs
from src import her2_candidates
from src import her2_spr_benchmark as spr

COHORT = "Primary novelty cohort"
SPR_SCORE = "predicted_neg_log10_kd_nm"
SHORTLIST_SIZE = 10
MIN_DISTANCE = 2
OOF_TOLERANCE = 1e-10
# Frozen candidate probabilities were serialized with a lower precision than the model replay.
FACS_TOLERANCE = 2e-7
STRUCTURAL_FLAG = "review_flag_low_vh_vl_confidence"
FEATURE_ORDER = [
    f"{position}:{residue}"
    for position in range(1, facs.LOOP_LENGTH + 1)
    for residue in facs.AMINO_ACIDS
]


def source(path):
    return {"source": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def position_features(sequences):
    return facs.sequence_features(sequences, "position").toarray().astype(float)


def load_primary(data_dir, benchmark):
    """Rebuild the primary cohort from the benchmark's source files."""
    splits, inputs = {}, {}
    for name in facs.SPLIT_NAMES:
        split, inputs[name] = facs.load_split(data_dir / f"{name}.csv")
        splits[name] = set(split["seq"])
    assay, inputs["spr"] = facs.read_csv(data_dir / "spr.csv", columns=spr.SPR_COLUMNS)
    for name, entry in inputs.items():
        if entry["sha256"] != benchmark["inputs"][name]["sha256"]:
            raise ValueError(f"{name} does not match the completed benchmark input hash")
    variants, _ = spr.prepare_spr(assay, splits)
    primary = variants.loc[~(variants["in_facs_train"] | variants["in_facs_val"])].copy()
    primary["group"] = spr.sequence_groups(primary["cdrh3"])
    return primary.reset_index(drop=True), inputs, splits


def finite_array(values, name):
    # ASVS 2.2.1: validate numerical inputs before reductions or comparisons.
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain finite numbers") from exc
    if not array.size or not np.isfinite(array).all():
        raise ValueError(f"{name} must contain finite numbers")
    return array


def max_gap(expected, observed, name, tolerance):
    expected = finite_array(expected, f"{name} reference")
    observed = finite_array(observed, name)
    if expected.shape != observed.shape:
        raise ValueError(f"{name} shape differs from the reference")
    gap = float(np.max(np.abs(expected - observed)))
    if not np.isfinite(gap) or gap > tolerance:
        raise ValueError(f"{name} differs from the reference by {gap:.6g}, limit {tolerance:g}")
    return gap


def load_oof(path, score_columns):
    frame, provenance = facs.read_csv(path)
    required = {"cdrh3", "cohort", "group", "fold", spr.TARGET, *score_columns}
    if missing := required - set(frame.columns):
        raise ValueError(f"OOF CSV is missing columns: {', '.join(sorted(missing))}")
    primary = frame.loc[frame["cohort"].eq(COHORT)].copy()
    primary = primary.astype({"group": int, "fold": int})
    for column in (spr.TARGET, *score_columns):
        primary[column] = finite_array(primary[column], f"OOF {column}")
    return primary, provenance


def original_folds(primary, oof, cohort):
    """Reproduce the original held-out groups and recover their recorded alphas."""
    aligned = oof.loc[oof["cohort"].eq(COHORT)].set_index("cdrh3")
    if not aligned.index.is_unique or set(aligned.index) != set(primary["cdrh3"]) or len(primary) != cohort["n"]:
        raise ValueError("Primary cohort IDs differ from the completed benchmark")
    aligned = aligned.loc[primary["cdrh3"]]
    groups = primary["group"].to_numpy()
    if not np.array_equal(aligned["group"], groups) or len(np.unique(groups)) != cohort["groups"]:
        raise ValueError("Primary sequence groups differ from the completed benchmark")
    splits = spr.group_splits(groups)
    if len(splits) != len(cohort["folds"]):
        raise ValueError("Outer fold count differs from the completed benchmark")
    folds = []
    for fold, ((train, test), record) in enumerate(zip(splits, cohort["folds"]), start=1):
        alpha = record["tuning"]["position_ridge"]["alpha"]
        if (
            record["fold"] != fold or len(train) != record["train_n"]
            or len(test) != record["test_n"] or not aligned["fold"].iloc[test].eq(fold).all()
            or np.unique(groups[test]).tolist() != record["test_groups"] or alpha not in spr.ALPHAS
        ):
            raise ValueError(f"Outer fold {fold} differs from the completed benchmark")
        folds.append((fold, train, test, alpha))
    return folds, aligned


def select_alpha(positions, labels, groups):
    """Use the benchmark's fixed group split, pooled RMSE and alpha tie break."""
    splits = spr.group_splits(groups)
    errors = {}
    for alpha in spr.ALPHAS:
        predictions = np.empty(len(labels))
        for train, validation in splits:
            model = Ridge(alpha=alpha, solver="cholesky").fit(positions[train], labels[train])
            predictions[validation] = model.predict(positions[validation])
        errors[alpha] = spr.rmse(labels, predictions)
    return min(errors, key=lambda alpha: (errors[alpha], -alpha)), errors


def fit_outer(positions, labels, fold):
    _, train, _, alpha = fold
    return Ridge(alpha=alpha, solver="cholesky").fit(positions[train], labels[train])


def export_parameters(model):
    return {
        "alpha": float(model.alpha),
        "intercept": float(model.intercept_),
        "coefficients": model.coef_.tolist(),
    }


def predict_from_parameters(parameters, sequences):
    """Score with JSON numbers and the declared feature order."""
    if parameters["feature_order"] != FEATURE_ORDER:
        raise ValueError("Unexpected SPR position feature ordering")
    coefficients = np.asarray(parameters["coefficients"], dtype=float)
    intercept = float(parameters["intercept"])
    if (
        coefficients.shape != (len(FEATURE_ORDER),)
        or not np.isfinite(coefficients).all() or not np.isfinite(intercept)
    ):
        raise ValueError("Invalid SPR ridge parameters")
    return position_features(sequences) @ coefficients + intercept


def fit_models(primary, folds):
    positions = position_features(primary["cdrh3"])
    labels = primary[spr.TARGET].to_numpy()
    groups = primary["group"].to_numpy()
    with threadpool_limits(limits=1):
        alpha, errors = select_alpha(positions, labels, groups)
        final = Ridge(alpha=alpha, solver="cholesky").fit(positions, labels)
        outer = [fit_outer(positions, labels, fold) for fold in folds]
    return final, outer, errors


def pareto_fronts(scores):
    """Number successive nondominated fronts starting at 1, maximizing both columns."""
    dominates = (scores[:, None] >= scores).all(axis=2) & (scores[:, None] > scores).any(axis=2)
    fronts = np.zeros(len(scores), dtype=int)
    front = 1
    while (remaining := np.flatnonzero(fronts == 0)).size:
        nondominated = remaining[~dominates[np.ix_(remaining, remaining)].any(axis=0)]
        fronts[nondominated] = front
        front += 1
    return fronts


def exploratory_shortlist(ranked):
    """Traverse fronts in order, alternating axis preference after each addition."""
    picked = []
    for _, front in ranked.groupby("pareto_front", sort=True):
        axes = [
            list(front.sort_values([axis, "cdrh3"], ascending=[False, True])["cdrh3"])
            for axis in (SPR_SCORE, facs.SCORE)
        ]
        while len(picked) < SHORTLIST_SIZE:
            candidate = next((
                seq for seq in axes[len(picked) % 2]
                if all(Levenshtein.distance(seq, previous) >= MIN_DISTANCE for previous in picked)
            ), None)
            if candidate is None:
                break
            picked.append(candidate)
        if len(picked) == SHORTLIST_SIZE:
            break
    return ranked.set_index("cdrh3", drop=False).loc[picked].reset_index(drop=True)


def validate_candidate_scores(candidates, old, facs_model):
    replay = facs.predict_scores(facs_model["estimator"], facs_model["features"], candidates["cdrh3"])
    recorded = finite_array(candidates[facs.SCORE], "Preflight FACS scores")
    if not ((recorded >= 0) & (recorded <= 1)).all():
        raise ValueError("Preflight FACS scores must be between 0 and 1")
    model_gap = max_gap(replay, recorded, "Preflight FACS scores", FACS_TOLERANCE)
    aligned = pd.DataFrame({"recorded": recorded, "replayed": replay}, index=candidates["cdrh3"]).loc[old.index]
    old_scores = finite_array(old[facs.SCORE], "Old shortlist FACS scores")
    if not ((old_scores >= 0) & (old_scores <= 1)).all():
        raise ValueError("Old shortlist FACS scores must be between 0 and 1")
    max_gap(aligned["replayed"], old_scores, "Old shortlist FACS scores", FACS_TOLERANCE)
    old_gap = max_gap(aligned["recorded"], old_scores, "Old shortlist versus preflight FACS scores", FACS_TOLERANCE)
    candidates[facs.SCORE] = replay
    return old_gap, model_gap


def rank_candidates(preflight, old, primary, final, outer, facs_model):
    candidates = pd.DataFrame(row for row in preflight if row["strict_local_gate"] is True)
    if candidates.empty:
        raise ValueError("No strict candidates to rank")
    # Headers are provided by preflight keys and the CSV writer escapes values only.
    if not all(re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", column) for column in candidates.columns):
        raise ValueError("Malformed ranked CSV column name")
    candidates = candidates.sort_values("cdrh3", ignore_index=True)
    seqs = facs.normalize_sequences(candidates["cdrh3"], "Strict candidate CDRH3", length=facs.LOOP_LENGTH)
    if not seqs.eq(candidates["cdrh3"]).all():
        raise ValueError("Strict candidate CDRH3 must be canonical")
    if not seqs.is_unique or old.empty or not old["cdrh3"].is_unique or not set(old["cdrh3"]).issubset(seqs):
        raise ValueError("Expected distinct strict candidates covering the old provisional shortlist")
    if not candidates["vh_sequence"].eq(facs.VH_PREFIX + seqs + facs.VH_SUFFIX).all():
        raise ValueError("Strict candidate VH differs from the fitted trastuzumab background")
    if not candidates["vl_sequence"].eq(facs.WT_VL).all():
        raise ValueError("Strict candidate VL differs from the fitted trastuzumab background")
    old = old.set_index("cdrh3")
    prior = candidates.set_index("cdrh3").loc[old.index]
    old_gap, model_gap = validate_candidate_scores(candidates, old, facs_model)
    for chain in ("vh_sequence", "vl_sequence"):
        if not prior[chain].eq(old[chain]).all():
            raise ValueError(f"Strict candidate {chain} differs from its old shortlist record")
    candidates["old_top20_rank"] = seqs.map(dict(zip(old.index, range(1, len(old) + 1)))).astype("Int64")
    previous_structure = seqs.map(old["structural_status"])
    candidates["structural_status_preflight"] = candidates["structural_status"]
    candidates["structural_status"] = previous_structure.fillna(candidates["structural_status"])
    candidates["structural_review_flag"] = candidates["structural_status"].eq(STRUCTURAL_FLAG)
    candidates["thermal_validation_pending"] = candidates["thermal_status"].isin(["", "missing", "pending"])

    positions = position_features(seqs)
    with threadpool_limits(limits=1):
        candidates[SPR_SCORE] = finite_array(final.predict(positions), "Candidate SPR predictions")
        sensitivity = finite_array(
            np.column_stack([model.predict(positions) for model in outer]), "Candidate SPR sensitivity predictions"
        )
    for fold, prediction in enumerate(sensitivity.T, start=1):
        candidates[f"spr_sensitivity_fold_{fold}"] = prediction
    candidates["spr_sensitivity_min"] = sensitivity.min(axis=1)
    candidates["spr_sensitivity_max"] = sensitivity.max(axis=1)
    candidates["spr_sensitivity_range"] = np.ptp(sensitivity, axis=1)
    training_ids = primary["cdrh3"].tolist()
    nearest = [process.extractOne(seq, training_ids, scorer=Levenshtein.distance) for seq in seqs]
    candidates["nearest_spr_training_cdrh3"] = [match[0] for match in nearest]
    candidates["nearest_spr_training_edit_distance"] = [match[1] for match in nearest]
    candidates["facs_rank"] = candidates[facs.SCORE].rank(method="min", ascending=False).astype(int)
    candidates["spr_rank"] = candidates[SPR_SCORE].rank(method="min", ascending=False).astype(int)
    candidates["facs_rank_minus_spr_rank"] = candidates["facs_rank"] - candidates["spr_rank"]
    candidates["pareto_front"] = pareto_fronts(candidates[[facs.SCORE, SPR_SCORE]].to_numpy())
    ranked = candidates.sort_values(["pareto_front", "spr_rank", "facs_rank", "cdrh3"])
    return ranked.drop(columns=["numbering", "cdrh3_by_scheme"]), {
        "max_facs_source_mismatch": old_gap, "max_preflight_facs_model_mismatch": model_gap,
    }


def write_ranked_csv(path, frame):
    # ASVS 1.2.10: escape text at export while retaining numeric types.
    records = frame.astype(object).where(frame.notna(), None).to_dict("records")
    her2_candidates.write_csv(path, frame.columns, records)


def render_report(ranked, shortlist, parameters, output):
    fit = parameters["fit"]
    diagnostics = parameters["diagnostics"]
    inputs = parameters["inputs"]
    old = set(ranked.loc[ranked["old_top20_rank"].notna(), "cdrh3"])
    facs_top = set(ranked.sort_values(["facs_rank", "cdrh3"]).head(20)["cdrh3"])
    spr_top = set(ranked.sort_values(["spr_rank", "cdrh3"]).head(20)["cdrh3"])
    common = shortlist.loc[shortlist["cdrh3"].isin(old), "cdrh3"].tolist()
    min_distance = min((
        Levenshtein.distance(left, right) for left, right in combinations(shortlist["cdrh3"], 2)
    ), default=None)
    front_sizes = ", ".join(f"{n} ({count})" for n, count in ranked["pareto_front"].value_counts().sort_index().items())
    disagreement = ranked.reindex(ranked["facs_rank_minus_spr_rank"].abs().nlargest(5).index)
    spread = ranked["spr_sensitivity_range"]
    distance = ranked["nearest_spr_training_edit_distance"]
    flagged = ranked.loc[ranked["structural_review_flag"], "cdrh3"].tolist()
    rho = spr.spearman(ranked[facs.SCORE].to_numpy(), ranked[SPR_SCORE].to_numpy())
    lines = [
        "# HER2 SPR candidate ranking", "",
        f"Scored all {len(ranked)} strict local passes from {diagnostics['preflight_rows']} existing candidates. "
        "This is an exploratory ranking for review. The SPR axis predicts −log10(KD in nM) "
        "from CDRH3 positions. The frozen FACS high class score remains a separate axis. "
        "Candidate KD, superiority to trastuzumab and thermal improvement remain experimentally unverified. "
        "The benchmark and fit condition on variants with reportable KD.", "",
        "## Fit and replay", "",
        f"The completed benchmark measured primary cohort OOF Spearman ρ "
        f"{diagnostics['benchmark_spr_rho']:.3f} for position ridge and "
        f"{diagnostics['benchmark_facs_rho']:.3f} for FACS. The post-benchmark fit used "
        f"{len(fit['fit_ids'])} measured primary variants in {fit['groups']} sequence groups. "
        f"Full-cohort grouped tuning selected alpha **{parameters['alpha']:g}**. "
        "Its pooled tuning RMSE is a selection diagnostic.", "",
        "| Alpha | Pooled tuning RMSE |", "| ---: | ---: |",
        *[f"| {alpha:g} | {value:.6f} |" for alpha, value in fit["pooled_rmse_by_alpha"].items()], "",
        "The sensitivity columns give the minimum, maximum and range across the five original outer "
        "training refits at their recorded alphas. They describe sensitivity to training composition "
        "and regularization. They carry no calibrated coverage.", "",
        f"Across candidates the median range is {spread.median():.3f} and the maximum is {spread.max():.3f} "
        f"score units. Nearest primary SPR training distance is {distance.min()}–{distance.max()} edits. "
        "All screened rows retain the source VH and VL, official validation and thermal status.", "",
        f"Maximum source mismatch is {diagnostics['max_oof_prediction_source_mismatch']:.3g} score units "
        "between replayed held-out predictions and the original OOF position ridge scores. "
        f"The rebuilt assay labels differ from the exported OOF target by at most "
        f"{diagnostics['max_oof_label_source_mismatch']:.3g}. "
        f"The maximum preflight versus old top 20 FACS difference is {diagnostics['max_facs_source_mismatch']:.3g}. "
        f"The verified FACS model rescore supplies the ranking axis and differs from preflight by at most "
        f"{diagnostics['max_preflight_facs_model_mismatch']:.3g}.", "",
        "## Rank disagreements", "",
        f"Candidate SPR and FACS ranks have Spearman ρ {rho:.3f}. Their individual top 20 sets share "
        f"{len(facs_top & spr_top)} rows. The old top 20 shares {len(old & facs_top)} with FACS top 20 "
        f"and {len(old & spr_top)} with SPR top 20. Equal axis scores share their best rank. "
        "Top 20 boundaries break ties by CDRH3.", "",
        "| CDRH3 | FACS rank | SPR rank | FACS rank minus SPR rank |", "| --- | ---: | ---: | ---: |",
        *[f"| {r.cdrh3} | {r.facs_rank} | {r.spr_rank} | {r.facs_rank_minus_spr_rank:+d} |"
          for r in disagreement.itertuples()], "",
        f"Pareto front sizes, listed as front (count), are {front_sizes}. Front 1 contains all "
        "nondominated rows when maximizing both axes. Later fronts apply the same rule after prior fronts are removed.", "",
        "## Exploratory shortlist", "",
        f"Selected {len(shortlist)} sequences with a realized minimum pairwise CDRH3 Levenshtein "
        f"distance of {min_distance}. The review list takes up to {SHORTLIST_SIZE} rows by Pareto front, "
        f"alternates SPR and FACS priority within the remaining eligible front and requires distance ≥{MIN_DISTANCE} "
        "from every earlier choice. Score ties use lexical CDRH3 order.", "",
        f"Overlap with the old top 20 is {len(common)} ({', '.join(common) or 'none'}). "
        f"Front 1 overlaps with the old top 20 in {len(old & set(ranked.loc[ranked['pareto_front'].eq(1), 'cdrh3']))} rows.", "",
        "| Pick | CDRH3 | Front | FACS score | SPR score | Fold spread | Nearest SPR edits | Old rank | Structural flag |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in shortlist.to_dict("records"):
        old_rank = row["old_top20_rank"] if pd.notna(row["old_top20_rank"]) else "—"
        lines.append(
            f"| {row['shortlist_pick']} | {row['cdrh3']} | {row['pareto_front']} | {row[facs.SCORE]:.6f} | "
            f"{row[SPR_SCORE]:.3f} | {row['spr_sensitivity_range']:.3f} | "
            f"{row['nearest_spr_training_edit_distance']} | {old_rank} | "
            f"{'Review' if row['structural_review_flag'] else '—'} |"
        )
    lines.extend([
        "", f"Structural review flag from the old shortlist: {', '.join(flagged) or 'none'}. "
        "Low predicted VH/VL confidence calls for review. Experimental rejection requires experimental evidence.",
        "", "## Reproduction", "",
        "Run from the repository root with the existing environment. The fit and shortlist are deterministic. "
        f"Seed {parameters['benchmark_seed']} identifies the completed benchmark's bootstrap.", "",
        "```sh", "source .venv/bin/activate", "PYTHONDONTWRITEBYTECODE=1 python -m scripts.rank_spr_candidates \\",
        f"  --data-dir {shlex.quote(str(Path(inputs['spr']['source']).parent))} \\",
        f"  --preflight {shlex.quote(inputs['preflight']['source'])} \\",
        f"  --oof {shlex.quote(inputs['oof']['source'])} \\",
        f"  --metrics {shlex.quote(inputs['benchmark_metrics']['source'])} \\",
        f"  --facs-model {shlex.quote(inputs['facs_model']['source'])} \\",
        f"  --old-top20 {shlex.quote(inputs['old_top20']['source'])} \\",
        f"  --output {shlex.quote(str(output))}", "```", "",
        "`model.json` records input paths and SHA-256 hashes, versions, feature order, numeric parameters, "
        "all fit CDRH3 IDs and original outer test IDs. `candidates.csv` includes both axes, ranks, fronts, "
        "fold predictions and flags. `shortlist.csv` contains the selected rows in pick order.", "",
    ])
    return "\n".join(lines)


def load_benchmark(metrics_path):
    benchmark = json.loads(metrics_path.read_text())
    if (
        tuple(benchmark["alpha_grid"]) != spr.ALPHAS or benchmark["max_folds"] != spr.MAX_FOLDS
        or benchmark["levenshtein_cutoff"] != spr.MAX_DISTANCE
    ):
        raise ValueError("The recorded benchmark settings differ from the loaded benchmark code")
    cohort = benchmark["cohorts"][COHORT]
    if cohort["split_issue"] is not None:
        raise ValueError("Expected a completed primary SPR benchmark")
    return benchmark, cohort


def run(data_dir, preflight_path, oof_path, metrics_path, old_path, output, facs_model_path=None):
    benchmark, cohort = load_benchmark(metrics_path)
    primary, inputs, _ = load_primary(data_dir, benchmark)
    oof, inputs["oof"] = load_oof(oof_path, ("position_ridge",))
    folds, aligned = original_folds(primary, oof, cohort)
    final, outer, errors = fit_models(primary, folds)
    outer_parameters = []
    positions = position_features(primary["cdrh3"])
    replay = np.full(len(primary), np.nan)
    for (fold, train, test, _), model in zip(folds, outer):
        with threadpool_limits(limits=1):
            replay[test] = model.predict(positions[test])
        outer_parameters.append({
            "fold": fold, "fit_ids": primary["cdrh3"].iloc[train].tolist(),
            "held_out_ids": primary["cdrh3"].iloc[test].tolist(), **export_parameters(model),
        })
    audit = {
        "max_oof_prediction_source_mismatch": max_gap(
            replay, aligned["position_ridge"], "Benchmark replay position_ridge", OOF_TOLERANCE
        ),
        "max_oof_label_source_mismatch": max_gap(
            primary[spr.TARGET], aligned[spr.TARGET], f"Benchmark replay {spr.TARGET}", OOF_TOLERANCE
        ),
    }

    recorded_model = benchmark["inputs"]["facs_model"]
    facs_model_path = facs_model_path or Path(recorded_model["source"])
    inputs["facs_model"] = source(facs_model_path)
    if inputs["facs_model"]["sha256"] != recorded_model["sha256"]:
        raise ValueError("FACS model does not match the completed benchmark input hash")
    facs_model = facs.load_model(facs_model_path)
    preflight = json.loads(preflight_path.read_text())
    old, inputs["old_top20"] = facs.read_csv(old_path, columns=(
        "cdrh3", "vh_sequence", "vl_sequence", "structural_status", facs.SCORE,
    ))
    ranked, facs_audit = rank_candidates(preflight, old, primary, final, outer, facs_model)
    shortlist = exploratory_shortlist(ranked)
    shortlist.insert(0, "shortlist_pick", range(1, len(shortlist) + 1))
    diagnostics = {
        **audit, **facs_audit, "preflight_rows": len(preflight),
        "benchmark_spr_rho": cohort["metrics"]["position_ridge"]["spearman"],
        "benchmark_facs_rho": cohort["metrics"][facs.SCORE]["spearman"],
    }
    inputs.update({
        "preflight": source(preflight_path), "benchmark_metrics": source(metrics_path),
        "ranking_script": source(Path("scripts/rank_spr_candidates.py")),
        "benchmark_code": source(Path("src/her2_spr_benchmark.py")),
        "feature_code": source(Path("src/her2_benchmark.py")),
    })
    parameters = {
        "version": 1, "model": "position_one_hot_ridge", "solver": "cholesky",
        "target": spr.TARGET, "feature_order": FEATURE_ORDER, **export_parameters(final),
        "benchmark_seed": benchmark["seed"], "randomized_steps_in_ranking": False,
        "fit": {
            "cohort": COHORT, "id_type": "cdrh3", "fit_ids": primary["cdrh3"].tolist(),
            "groups": int(primary["group"].nunique()), "group_levenshtein_cutoff": spr.MAX_DISTANCE,
            "max_group_folds": spr.MAX_FOLDS, "alpha_grid": spr.ALPHAS,
            "pooled_rmse_by_alpha": errors,
        },
        "outer_refits": outer_parameters,
        "shortlist": {"limit": SHORTLIST_SIZE, "minimum_pairwise_levenshtein": MIN_DISTANCE},
        "inputs": inputs, "diagnostics": diagnostics,
        "versions": {"python": platform.python_version(), **{
            package: version(package) for package in benchmark["versions"] if package != "python"
        }},
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "model.json").write_text(json.dumps(parameters, indent=2, allow_nan=False) + "\n")
    write_ranked_csv(output / "candidates.csv", ranked)
    write_ranked_csv(output / "shortlist.csv", shortlist)
    (output / "report.md").write_text(render_report(ranked, shortlist, parameters, output))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True, help="Original SPR and FACS split CSV directory")
    parser.add_argument("--preflight", type=Path, default=Path("results/screen-inputs/initial-preflight.json"))
    parser.add_argument("--oof", type=Path, default=Path("results/spr-benchmark/oof_predictions.csv"))
    parser.add_argument("--metrics", type=Path, default=Path("results/spr-benchmark/metrics.json"))
    parser.add_argument("--facs-model", type=Path, help="Verified model at a different path than recorded in metrics")
    parser.add_argument("--old-top20", type=Path, default=Path("results/experiment-20261114/provisional_top20.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/spr-candidate-ranking"))
    args = parser.parse_args()
    run(args.data_dir, args.preflight, args.oof, args.metrics, args.old_top20, args.output, args.facs_model)


if __name__ == "__main__":
    main()
