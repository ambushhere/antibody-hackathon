"""Evaluate prespecified CDRH3 interaction ridge on the frozen primary SPR cohort."""

import argparse
import json
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from scripts import rank_spr_candidates as frozen
from src import her2_benchmark as facs
from src import her2_spr_benchmark as spr

POSITION = "position_ridge"
PAIRWISE = "position_pairwise_ridge"
DISTANT = "distant75"
BASELINE_RHO = 0.6127403292354852
DISTANCE_SHA256 = "359330bd47638d774c8f0e6754634cb0d876ad671d61ac7338d81ccc730ab770"


def ridge_predict(features, train, test, fit_labels, alphas):
    """Fit intercept ridge using an n_train by n_train centered Gram matrix."""
    x_train, x_test = features[train], features[test]
    if issparse(features):
        gram = (x_train @ x_train.T).toarray()
        cross = (x_test @ x_train.T).toarray()
    else:
        scaler = StandardScaler().fit(x_train)
        x_train, x_test = scaler.transform(x_train), scaler.transform(x_test)
        gram, cross = x_train @ x_train.T, x_test @ x_train.T

    # These are X_centered @ X_centered.T and the corresponding test/train product.
    train_mean = gram.mean(axis=0)
    grand_mean = train_mean.mean()
    gram = gram - train_mean[:, None] - train_mean + grand_mean
    cross = cross - cross.mean(axis=1, keepdims=True) - train_mean + grand_mean
    label_mean = fit_labels.mean()
    centered_labels = fit_labels - label_mean
    identity = np.eye(len(train))
    return np.stack([
        label_mean + cross @ solve(gram + alpha * identity, centered_labels, assume_a="pos")
        for alpha in alphas
    ])


def fold_metrics(labels, scores):
    return {
        name: {
            "spearman": spr.spearman(labels, prediction),
            "rmse": None if name == facs.SCORE else spr.rmse(labels, prediction),
        }
        for name, prediction in scores.items()
    }


def nested_oof(features, labels, groups, folds):
    """Use the supplied outer folds and tune each model on matching inner group splits."""
    scores = {name: np.full(len(labels), np.nan) for name in features}
    provenance = []
    with threadpool_limits(limits=1):
        for fold, train, test, frozen_alpha in folds:
            inner = spr.group_splits(groups[train])
            tuning = {}
            for name, matrix in features.items():
                validation_scores = np.full((len(spr.ALPHAS), len(train)), np.nan)
                for fit, val in inner:
                    fit_rows = train[fit]
                    validation_scores[:, val] = ridge_predict(
                        matrix, fit_rows, train[val], labels[fit_rows], spr.ALPHAS
                    )
                frozen.finite_array(validation_scores, f"{name} inner predictions")
                errors = [spr.rmse(labels[train], prediction) for prediction in validation_scores]
                best = min(range(len(spr.ALPHAS)), key=lambda i: (errors[i], -spr.ALPHAS[i]))
                alpha = spr.ALPHAS[best]
                if name == POSITION and alpha != frozen_alpha:
                    raise ValueError(f"Position ridge selected a different alpha in outer fold {fold}")
                scores[name][test] = ridge_predict(matrix, train, test, labels[train], (alpha,))[0]
                tuning[name] = {"alpha": alpha, "inner_pooled_rmse_by_alpha": dict(zip(spr.ALPHAS, errors))}
            provenance.append({
                "fold": fold,
                "train_n": len(train),
                "test_n": len(test),
                "train_groups": np.unique(groups[train]).tolist(),
                "test_groups": np.unique(groups[test]).tolist(),
                "inner_folds": [
                    {"train_n": len(fit), "val_n": len(val), "val_groups": np.unique(groups[train[val]]).tolist()}
                    for fit, val in inner
                ],
                "tuning": tuning,
            })
    for name, predictions in scores.items():
        frozen.finite_array(predictions, f"{name} outer predictions")
    return scores, provenance


def load_distances(path, primary, aligned):
    """Join the already-viewed diagnostic without rerunning or changing its filter."""
    distances, source = frozen.load_oof(path, (POSITION, facs.SCORE))
    if source["sha256"] != DISTANCE_SHA256:
        raise ValueError("Distance rows differ from the frozen distance diagnostic")
    distances = distances.set_index("cdrh3")
    if not distances.index.is_unique or set(distances.index) != set(primary["cdrh3"]):
        raise ValueError("Distance diagnostic IDs differ from the frozen primary cohort")
    distances = distances.loc[primary["cdrh3"]]
    for column in ("fold", "group", spr.TARGET, POSITION, facs.SCORE):
        frozen.max_gap(aligned[column], distances[column], f"Distance diagnostic {column}", frozen.OOF_TOLERANCE)
    diagnostics = distances[["nearest_facs_trainval", "nearest_spr_outer_train"]].astype(int)
    diagnostics[DISTANT] = diagnostics.ge(3).all(axis=1)
    if not diagnostics["nearest_spr_outer_train"].ge(3).all() or diagnostics[DISTANT].sum() != 75:
        raise ValueError("Distance diagnostic does not match the prespecified 75-variant subset")
    return diagnostics.reset_index(drop=True), source


def load_frozen_features(npz_path, metadata_path, sequence_ids):
    """Load named fixed embeddings from a non-object NPZ in exact sequence order."""
    metadata = json.loads(metadata_path.read_text())
    if not isinstance(metadata, dict) or not metadata:
        raise ValueError("Frozen feature metadata must be a nonempty JSON object")
    json.dumps(metadata, allow_nan=False)
    npz_source = frozen.source(npz_path)
    if "ordered_sequence_ids" in metadata and metadata["ordered_sequence_ids"] != list(sequence_ids):
        raise ValueError("Frozen feature metadata IDs must match the exact primary cohort order")
    declared_output = metadata.get("output", {})
    if not isinstance(declared_output, dict):
        raise ValueError("Frozen feature metadata output must be an object")
    if "sha256" in declared_output and declared_output["sha256"] != npz_source["sha256"]:
        raise ValueError("Frozen feature NPZ does not match the metadata hash")
    with np.load(npz_path, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError("Frozen feature array names must be unique")
        if "sequence_ids" not in archive.files:
            raise ValueError("Frozen features require a sequence_ids array")
        ids = archive["sequence_ids"]
        if ids.dtype.kind != "U" or not np.array_equal(ids, np.asarray(sequence_ids, dtype=str)):
            raise ValueError("Frozen feature sequence IDs must match the exact primary cohort order")
        features = {}
        for name in sorted(set(archive.files) - {"sequence_ids"}):
            if not name.isidentifier() or not name.isascii():
                raise ValueError(f"Invalid frozen feature name: {name!r}")
            values = archive[name]
            if (values.ndim != 2 or values.shape[0] != len(ids) or values.shape[1] == 0
                    or values.dtype.kind not in "fiu"):
                raise ValueError(f"Frozen feature {name} must be a real numeric matrix with one row per sequence")
            features[name] = frozen.finite_array(values, f"Frozen feature {name}")
    if not features:
        raise ValueError("Frozen features require at least one named numeric matrix")
    return features, {"npz": npz_source, "metadata": frozen.source(metadata_path), "extraction": metadata}


def metric_difference(left, right):
    return None if left is None or right is None else left - right


def summarize(labels, scores, groups, comparisons):
    """Compute point metrics and paired percentile intervals for fixed OOF scores."""
    members = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    metrics = fold_metrics(labels, scores)
    pairs = {f"{left} minus {right}": (left, right) for left, right in comparisons}
    draws = {name: {metric: [] for metric in values} for name, values in metrics.items()}
    deltas = {pair: {metric: [] for metric in ("spearman", "rmse")} for pair in pairs}
    rng = np.random.default_rng(spr.SEED)
    for _ in range(spr.BOOTSTRAP_DRAWS):
        rows = np.concatenate([members[i] for i in rng.integers(len(members), size=len(members))])
        drawn = fold_metrics(labels[rows], {name: values[rows] for name, values in scores.items()})
        for name, values in drawn.items():
            for metric, value in values.items():
                if value is not None:
                    draws[name][metric].append(value)
        for pair, (left, right) in pairs.items():
            for metric in deltas[pair]:
                value = metric_difference(drawn[left][metric], drawn[right][metric])
                if value is not None:
                    deltas[pair][metric].append(value)

    def interval(values):
        return {"95_percent": np.percentile(values, [2.5, 97.5]).tolist() if values else None,
                "valid_draws": len(values)}

    return {
        "n": len(labels),
        "groups": len(members),
        "metrics": {
            name: {metric: {"estimate": value, **interval(draws[name][metric])} for metric, value in values.items()}
            for name, values in metrics.items()
        },
        "paired_differences": {
            pair: {
                metric: {"estimate": metric_difference(metrics[left][metric], metrics[right][metric]),
                         **interval(values)}
                for metric, values in deltas[pair].items()
            }
            for pair, (left, right) in pairs.items()
        },
    }


def code_sources():
    root = Path(__file__).resolve().parents[1]
    return {
        str(path): {**frozen.source(root / path), "source": str(path)}
        for path in (
            Path(__file__).resolve().relative_to(root),
            *(Path(module.__file__).resolve().relative_to(root) for module in (facs, spr, frozen)),
        )
    }


def run(data_dir, metrics_path, oof_path, distance_path, output, feature_npz=None, feature_metadata=None,
        feature_contrasts=()):
    benchmark, cohort = frozen.load_benchmark(metrics_path)
    primary, inputs, _ = frozen.load_primary(data_dir, benchmark)
    original, inputs["frozen_oof"] = frozen.load_oof(oof_path, (POSITION, facs.SCORE))
    folds, aligned = frozen.original_folds(primary, original, cohort)
    distances, inputs["distance_diagnostic"] = load_distances(distance_path, primary, aligned)
    inputs["frozen_metrics"] = frozen.source(metrics_path)
    labels, groups = primary[spr.TARGET].to_numpy(), primary["group"].to_numpy()
    label_gap = frozen.max_gap(labels, aligned[spr.TARGET], "OOF labels", frozen.OOF_TOLERANCE)
    if not aligned[facs.SCORE].between(0, 1).all():
        raise ValueError("Frozen FACS probability must be between 0 and 1")

    ids = primary["cdrh3"].tolist()
    features = {
        POSITION: facs.sequence_features(primary["cdrh3"], "position").astype(float),
        PAIRWISE: facs.sequence_features(primary["cdrh3"], "pairwise").astype(float),
    }
    if (feature_npz is None) != (feature_metadata is None):
        raise ValueError("Supply both --feature-npz and --feature-metadata")
    extra = {}
    if feature_npz is not None:
        extra, inputs["frozen_features"] = load_frozen_features(feature_npz, feature_metadata, ids)
        if set(extra) & (set(features) | set(primary) | set(aligned) | set(distances)):
            raise ValueError("Frozen feature names must not replace benchmark columns")
        features.update(extra)
    for left, right in feature_contrasts:
        if left not in extra or right not in extra or left == right:
            raise ValueError("Feature contrasts require two different supplied embedding names")

    scores, fit = nested_oof(features, labels, groups, folds)
    scores[facs.SCORE] = aligned[facs.SCORE].to_numpy()
    replay = {
        "tolerance": frozen.OOF_TOLERANCE,
        "max_label_gap": label_gap,
        "max_position_oof_gap": frozen.max_gap(
            aligned[POSITION], scores[POSITION], "Position ridge OOF replay", frozen.OOF_TOLERANCE
        ),
        "pooled_position_spearman": spr.spearman(labels, scores[POSITION]),
        "frozen_pooled_position_spearman": BASELINE_RHO,
    }
    frozen.max_gap(BASELINE_RHO, replay["pooled_position_spearman"], "Position ridge pooled Spearman", frozen.OOF_TOLERANCE)
    for (_, _, test, _), record, reference in zip(folds, fit, cohort["folds"]):
        if record["inner_folds"] != reference["inner_folds"]:
            raise ValueError(f"Inner group folds differ from the frozen benchmark in outer fold {record['fold']}")
        errors = record["tuning"][POSITION]["inner_pooled_rmse_by_alpha"]
        ref_errors = reference["tuning"][POSITION]["inner_pooled_rmse_by_alpha"]
        record["max_position_inner_rmse_gap"] = frozen.max_gap(
            [ref_errors[str(a)] for a in spr.ALPHAS], list(errors.values()),
            f"Position inner RMSE fold {record['fold']}", frozen.OOF_TOLERANCE,
        )
        record["metrics"] = fold_metrics(labels[test], {name: values[test] for name, values in scores.items()})

    comparisons = [(name, baseline) for name in (PAIRWISE, *extra)
                   for baseline in (POSITION, PAIRWISE, facs.SCORE) if name != baseline]
    comparisons.extend(feature_contrasts)
    mask = distances[DISTANT].to_numpy()
    predictions = primary.join(distances).assign(fold=aligned["fold"].to_numpy(), **scores)
    result = {
        "status": "retrospective_method_development_on_already_viewed_data",
        "target": spr.TARGET,
        "cohort": frozen.COHORT,
        "spr_audit": benchmark["spr_audit"],
        "alpha_grid": spr.ALPHAS,
        "alpha_selection": "Pooled inner GroupKFold RMSE, ties prefer larger alpha",
        "levenshtein_group_cutoff": spr.MAX_DISTANCE,
        "max_folds": spr.MAX_FOLDS,
        "fit_membership": "Join groups to oof_predictions.csv. Inner training uses outer train_groups minus inner val_groups.",
        "distant75_definition": "Frozen rows at Levenshtein distance >=3 from FACS train/val and the row's SPR outer training fold",
        "bootstrap": {
            "draws": spr.BOOTSTRAP_DRAWS, "seed": spr.SEED,
            "method": "Within each subset, resample original groups with replacement, pair all models on each draw, no refits",
        },
        "models": {
            **{name: {"feature_columns": matrix.shape[1], "preprocessing": (
                "Fixed sparse one-hot features and train-only intercept centering" if issparse(matrix)
                else "Frozen dense features and train-only StandardScaler and intercept centering"
            ), "solver": "Centered dual Gram matrix, positive-definite solve"} for name, matrix in features.items()},
            facs.SCORE: {
                "preprocessing": "Frozen FACS high probability from the original SPR OOF CSV",
                "original_model_source": benchmark["inputs"]["facs_model"],
                "solver": None,
            },
        },
        "replay": replay,
        "folds": fit,
        "subsets": {
            "pooled": summarize(labels, scores, groups, comparisons),
            DISTANT: summarize(labels[mask], {name: values[mask] for name, values in scores.items()}, groups[mask], comparisons),
        },
        "inputs": inputs,
        "code": code_sources(),
        "versions": {"python": platform.python_version(), **{name: version(name) for name in (
            "numpy", "pandas", "scikit-learn", "scipy", "rapidfuzz", "threadpoolctl",
        )}},
    }
    output.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output / "oof_predictions.csv", index=False)
    (output / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, default=Path("results/spr-benchmark/metrics.json"))
    parser.add_argument("--oof", type=Path, default=Path("results/spr-benchmark/oof_predictions.csv"))
    parser.add_argument("--distance-rows", type=Path, default=Path("results/spr-distance-diagnostic/distance_rows.csv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--feature-npz", type=Path)
    parser.add_argument("--feature-metadata", type=Path)
    parser.add_argument("--feature-contrast", nargs=2, action="append", default=None, metavar=("MODEL", "BASELINE"))
    args = parser.parse_args()
    result = run(args.data_dir, args.metrics, args.oof, args.distance_rows, args.output,
                 args.feature_npz, args.feature_metadata, args.feature_contrast or ())
    print(json.dumps({name: subset["metrics"] for name, subset in result["subsets"].items()}, indent=2))


if __name__ == "__main__":
    main()
