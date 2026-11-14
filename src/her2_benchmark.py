"""Benchmark classifiers for the published HER2 CDRH3 FACS splits."""

import argparse
import hashlib
import itertools
import json
import math
import platform
import re
from pathlib import Path
from zipfile import BadZipFile, ZipFile

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr, spearmanr
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from threadpoolctl import threadpool_limits

SPLIT_PATH = (
    "data/affinity_data/her2/her2_aff_large/processed/remove_overlap/"
    "clonotype_split/0.7_0.15_0.15/"
)
SPLIT_NAMES = ("train", "val", "test")
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
WT_CDRH3 = "WGGDGFYAMD"
LOOP_LENGTH = len(WT_CDRH3)
VH_PREFIX = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTI"
    "SADTSKNTAYLQMNSLRAEDTAVYYCSR"
)
VH_SUFFIX = "YWGQGTLVTVSS"
WT_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG"
    "TDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
SCORE = "p_facs_high"
MODEL_VERSION = 1
PAIR_LEFT, PAIR_RIGHT = np.triu_indices(LOOP_LENGTH, k=1)
RESIDUE_CODES = np.full(256, -1, dtype=np.int32)
RESIDUE_CODES[np.frombuffer(AMINO_ACIDS.encode("ascii"), dtype=np.uint8)] = np.arange(
    len(AMINO_ACIDS)
)


def read_csv(path, *, archive=None, columns=None):
    """Hash the uncompressed CSV bytes, regardless of input form."""
    source = str(path) if archive is None else f"{archive.filename}!{path}"
    stream = path.open("rb") if archive is None else archive.open(path)
    with stream:
        sha256 = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha256.update(block)
        stream.seek(0)
        frame = pd.read_csv(
            stream, dtype=str, keep_default_na=False,
            usecols=None if columns is None else lambda name: name in columns,
        )
    missing = set(columns or ()) - set(frame.columns)
    if missing:
        raise ValueError(f"{source}: missing columns {', '.join(sorted(missing))}")
    return frame, {"source": source, "sha256": sha256.hexdigest()}


def normalize_sequences(values, source, *, length=None):
    sequences = values.str.strip()
    valid = sequences.str.fullmatch(
        f"[{AMINO_ACIDS}]+", flags=re.IGNORECASE | re.ASCII, na=False
    )
    if not valid.all():
        examples = sequences[~valid].head(3).tolist()
        raise ValueError(f"{source}: missing sequence or unknown residues: {examples}")
    if length is not None and not sequences.str.len().eq(length).all():
        examples = sequences[sequences.str.len().ne(length)].head(3).tolist()
        raise ValueError(f"{source}: expected {length} residues: {examples}")
    return sequences.str.upper()


def normalize_classes(values, source):
    classes = values.str.strip().str.lower().replace({"medium": "mid"})
    if not classes.isin(["high", "mid", "low"]).all():
        raise ValueError(f"{source}: expected high, mid/medium or low")
    return classes


def load_split(path, *, archive=None):
    split, provenance = read_csv(
        path, archive=archive, columns=("seq", "label", "class", "cluster_centre")
    )
    source = provenance["source"]
    for column in ("seq", "cluster_centre"):
        split[column] = normalize_sequences(
            split[column], f"{source} {column}", length=LOOP_LENGTH
        )
    split["class"] = normalize_classes(split["class"], source)
    if not split["label"].isin(["0", "1"]).all():
        raise ValueError(f"{source}: expected binary labels 0 and 1")
    split["label"] = split["label"].astype(np.int8)
    if not split["label"].eq(split["class"].eq("high")).all():
        raise ValueError(f"{source}: binary label disagrees with published FACS class")
    if split["label"].nunique() != 2:
        raise ValueError(f"{source}: both binary classes are required")
    return split, provenance


def audit_splits(splits):
    counts, overlaps, sequences, centres = {}, {}, {}, {}
    for name, split in splits.items():
        sequences[name] = set(split["seq"])
        centres[name] = set(split["cluster_centre"])
        duplicates = len(split) - len(sequences[name])
        if duplicates:
            raise ValueError(f"{name}: {duplicates} duplicate CDRH3 rows")
        counts[name] = {
            "rows": len(split),
            "unique_cdrh3": len(sequences[name]),
            "duplicate_cdrh3_rows": duplicates,
            "unique_cluster_centres": len(centres[name]),
            "facs_class_counts": {
                label: int(split["class"].eq(label).sum()) for label in ("high", "mid", "low")
            },
            "positive_high": int(split["label"].sum()),
            "negative_mid_or_low": int(split["label"].eq(0).sum()),
        }
    for left, right in itertools.combinations(splits, 2):
        overlap = {
            "cdrh3": len(sequences[left] & sequences[right]),
            "cluster_centre": len(centres[left] & centres[right]),
            "left_cdrh3_right_centre": len(sequences[left] & centres[right]),
            "left_centre_right_cdrh3": len(centres[left] & sequences[right]),
        }
        if any(overlap.values()):
            raise ValueError(f"{left}/{right} split overlap: {overlap}")
        overlaps[f"{left}/{right}"] = overlap
    return {"counts": counts, "overlaps": overlaps}


def audit_facs(path, splits):
    facs, provenance = read_csv(path, columns=("cdrh3", "binding_class"))
    facs["cdrh3"] = normalize_sequences(facs["cdrh3"], "FACS cdrh3", length=LOOP_LENGTH)
    facs["binding_class"] = normalize_classes(facs["binding_class"], "FACS binding_class")
    if facs["cdrh3"].duplicated().any():
        raise ValueError("FACS mapping contains repeated CDRH3 sequences")
    class_by_sequence = facs.set_index("cdrh3")["binding_class"]
    matched = {}
    for name, split in splits.items():
        mapped = split["seq"].map(class_by_sequence)
        missing = int(mapped.isna().sum())
        disagreed = int((mapped.notna() & mapped.ne(split["class"])).sum())
        if missing or disagreed:
            raise ValueError(
                f"{name}: FACS mapping has {missing} missing and {disagreed} conflicting rows"
            )
        matched[name] = len(split)
    return provenance, {
        "rows": len(facs),
        "matched_rows_by_split": matched,
        "unused_mapping_rows": len(facs) - sum(matched.values()),
    }


def sequence_features(sequences, kind):
    """Build fixed sparse columns for amino acid identities and position pairs."""
    sequences = normalize_sequences(sequences, "CDRH3 features", length=LOOP_LENGTH)
    if kind == "dummy":
        return np.zeros((len(sequences), 1), dtype=np.float32)
    if kind not in ("position", "pairwise"):
        raise ValueError(f"Unknown feature set: {kind}")
    residues = RESIDUE_CODES[
        np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8)
    ].reshape(-1, LOOP_LENGTH)

    alphabet_size = len(AMINO_ACIDS)
    position_width = LOOP_LENGTH * alphabet_size
    active = residues + np.arange(LOOP_LENGTH, dtype=np.int32) * alphabet_size
    width = position_width
    if kind == "pairwise":
        pair_width = alphabet_size**2
        pair_offsets = position_width + np.arange(len(PAIR_LEFT), dtype=np.int32) * pair_width
        pair_codes = residues[:, PAIR_LEFT] * alphabet_size + residues[:, PAIR_RIGHT]
        pair_codes += pair_offsets
        active = np.concatenate((active, pair_codes), axis=1)
        width += len(PAIR_LEFT) * pair_width

    indices = active.ravel()
    indptr = np.arange(len(sequences) + 1, dtype=np.int32) * active.shape[1]
    return csr_matrix(
        (np.ones(indices.size, dtype=np.float32), indices, indptr),
        shape=(len(sequences), width),
    )


def top_metrics(labels, scores, fraction):
    # If the cutoff bisects a score tie, use expected positives within the tie.
    n_top = math.ceil(fraction * len(scores))
    cutoff = np.partition(scores, len(scores) - n_top)[-n_top]
    above, tied = scores > cutoff, scores == cutoff
    slots = n_top - above.sum()
    precision = (labels[above].sum() + slots * labels[tied].mean()) / n_top
    return {
        "n": n_top,
        "precision": float(precision),
        "enrichment": float(precision / labels.mean()),
    }


def classification_metrics(labels, scores):
    labels = np.asarray(labels)
    return {
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc_average_precision": float(average_precision_score(labels, scores)),
        "balanced_accuracy_at_0_5": float(balanced_accuracy_score(labels, scores > 0.5)),
        "brier": float(brier_score_loss(labels, scores)),
        "top_1_percent": top_metrics(labels, scores, 0.01),
        "top_5_percent": top_metrics(labels, scores, 0.05),
    }


def predict_scores(estimator, kind, sequences):
    return estimator.predict_proba(sequence_features(sequences, kind))[:, 1]


def prepare_spr(frame, split_sequences):
    """Apply assay and scaffold restrictions before aggregating quantified replicates."""
    frame = frame.copy()
    for column in ("cdrh3", "vh_sequence", "vl_sequence"):
        frame[column] = normalize_sequences(frame[column], f"SPR {column}")

    quantified = frame["binding_result"].str.strip().eq("quantified_binding")
    kd = pd.to_numeric(frame["kd_nM"], errors="coerce")
    valid_kd = np.isfinite(kd) & kd.gt(0)
    eligible = quantified & valid_kd
    correct_background = (
        frame["vl_sequence"].eq(WT_VL)
        & frame["vh_sequence"].eq(VH_PREFIX + frame["cdrh3"] + VH_SUFFIX)
    )
    is_wt = frame["cdrh3"].eq(WT_CDRH3)
    length_supported = frame["cdrh3"].str.len().eq(LOOP_LENGTH)
    audit = {
        "input_rows": len(frame),
        "non_quantified_rows": int((~quantified).sum()),
        "quantified_rows_with_missing_or_invalid_kd": int((quantified & ~valid_kd).sum()),
        "quantified_rows_with_finite_positive_kd": int(eligible.sum()),
        "excluded_non_trastuzumab_background_rows": int((eligible & ~correct_background).sum()),
    }
    eligible &= correct_background
    audit["excluded_wild_type_rows"] = int((eligible & is_wt).sum())
    eligible &= ~is_wt

    overlaps = {name: frame["cdrh3"].isin(values) for name, values in split_sequences.items()}
    audit["quantified_variant_overlap_rows_by_facs_split"] = {
        name: int((eligible & overlap).sum()) for name, overlap in overlaps.items()
    }
    eligible &= ~(overlaps["train"] | overlaps["val"])
    audit["excluded_unsupported_length_rows"] = int((eligible & ~length_supported).sum())
    eligible &= length_supported
    variants = (
        frame.loc[eligible, ["cdrh3"]]
        .assign(neg_log10_kd_nm=-np.log10(kd[eligible]))
        .groupby("cdrh3", as_index=False)
        .agg(
            neg_log10_kd_nm=("neg_log10_kd_nm", "median"),
            quantified_measurements=("neg_log10_kd_nm", "size"),
        )
    )
    variants["in_facs_test"] = variants["cdrh3"].isin(split_sequences["test"])
    audit["eligible_quantified_rows"] = int(eligible.sum())
    audit["eligible_unique_cdrh3"] = len(variants)
    audit["replicate_rows_collapsed"] = int(eligible.sum()) - len(variants)
    return variants, audit


def correlations(variants):
    scores, kd = variants[SCORE], variants["neg_log10_kd_nm"]
    result = {"n": len(variants), "spearman_rho": None, "pearson_r": None}
    if len(variants) < 3 or scores.nunique() < 2 or kd.nunique() < 2:
        return result
    result["spearman_rho"] = float(spearmanr(scores, kd)[0])
    result["pearson_r"] = float(pearsonr(scores, kd)[0])
    return result


def evaluate_spr(path, selected_model, splits):
    frame, provenance = read_csv(
        path, columns=("cdrh3", "vh_sequence", "vl_sequence", "kd_nM", "binding_result")
    )
    sequences = {name: set(split["seq"]) for name, split in splits.items()}
    variants, audit = prepare_spr(frame, sequences)
    variants[SCORE] = (
        predict_scores(selected_model["estimator"], selected_model["features"], variants["cdrh3"])
        if len(variants) else np.empty(0, dtype=float)
    )
    return variants, provenance, {
        "audit": audit,
        "primary": correlations(variants),
        "excluding_facs_test": correlations(variants[~variants["in_facs_test"]]),
    }


def render_report(results):
    audit = results["split_audit"]
    lines = [
        "# HER2 FACS supervised benchmark",
        "",
        "Target is the published high FACS class versus mid and low. Models use CDRH3 only.",
        f"Seed {results['seed']}. Selected {results['selected_model']} by validation average precision.",
        "",
        "| Split | Rows | CDRH3s | Clonotypes | High | Mid | Low |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for name, count in audit["counts"].items():
        classes = count["facs_class_counts"]
        lines.append(
            f"| {name} | {count['rows']:,} | {count['unique_cdrh3']:,} | "
            f"{count['unique_cluster_centres']:,} | {classes['high']:,} | "
            f"{classes['mid']:,} | {classes['low']:,} |"
        )
    lines.extend([
        "",
        "All split pairs had zero CDRH3, clonotype centre and sequence-to-centre overlap.",
        "",
    ])
    if results["facs_mapping"] is not None:
        matched = sum(results["facs_mapping"]["matched_rows_by_split"].values())
        lines.extend([f"The optional FACS mapping agreed on all {matched:,} split rows.", ""])
    lines.extend([
        "| Split | Model | ROC AUC | PR AUC (AP) | Balanced accuracy | Brier | "
        "Precision top 1% | Enrichment top 1% | Precision top 5% | Enrichment top 5% |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for split_name, models in results["classification"].items():
        for model_name, metrics in models.items():
            top1, top5 = metrics["top_1_percent"], metrics["top_5_percent"]
            lines.append(
                f"| {split_name} | {model_name} | {metrics['roc_auc']:.4f} | "
                f"{metrics['pr_auc_average_precision']:.4f} | "
                f"{metrics['balanced_accuracy_at_0_5']:.4f} | {metrics['brier']:.4f} | "
                f"{top1['precision']:.4f} | {top1['enrichment']:.4f} | "
                f"{top5['precision']:.4f} | {top5['enrichment']:.4f} |"
            )
    lines.extend([
        "",
        "Balanced accuracy uses a probability threshold of 0.5. Precision at a tied cutoff is "
        "the expected precision over tied rows. Enrichment is relative to split prevalence.",
        "",
        "## SPR evaluation",
        "",
    ])
    spr = results["spr"]
    if spr is None:
        lines.append("No SPR input was supplied.")
    else:
        lines.extend([
            "Correlations compare predicted high FACS probability to the median "
            "of −log10(KD in nM) for quantified replicates of each eligible CDRH3.",
            "Eligibility requires the trastuzumab VH framework and VL, a ten-residue CDRH3 "
            "distinct from wild type and no exact train or validation CDRH3 match.",
            "",
            "| Cohort | CDRH3s | Spearman ρ | Pearson r |",
            "| --- | ---: | ---: | ---: |",
        ])
        for name, label in (
            ("primary", "Eligible SPR"), ("excluding_facs_test", "Also excludes FACS test CDRH3s")
        ):
            values = spr[name]
            rho, pearson = values["spearman_rho"], values["pearson_r"]
            rho_text = "n/a" if rho is None else f"{rho:.4f}"
            pearson_text = "n/a" if pearson is None else f"{pearson:.4f}"
            lines.append(f"| {label} | {values['n']} | {rho_text} | {pearson_text} |")
        audit = spr["audit"]
        excluded_overlaps = audit["quantified_variant_overlap_rows_by_facs_split"]
        lines.extend([
            "",
            f"Of {audit['input_rows']} assay rows, "
            f"{audit['quantified_rows_with_finite_positive_kd']} had quantified, finite positive KD.",
            f"Exclusions from those rows were {audit['excluded_non_trastuzumab_background_rows']} "
            f"other backgrounds, {audit['excluded_wild_type_rows']} wild type, "
            f"{excluded_overlaps['train']} train matches, {excluded_overlaps['val']} validation "
            f"matches and {audit['excluded_unsupported_length_rows']} unsupported loop lengths.",
            f"Aggregation collapsed {audit['replicate_rows_collapsed']} replicate rows.",
            "These correlations describe the cohort with a reportable KD.",
        ])
    lines.extend(["", "## Input hashes", ""])
    for name, info in results["inputs"].items():
        lines.append(f"- {name} {info['sha256']} ({info['source']})")
    return "\n".join(lines) + "\n"


def save_model(artifact, path):
    model = artifact['estimator']
    payload = {'version': MODEL_VERSION, 'features': artifact['features']}
    if artifact['features'] == 'dummy':
        payload['prior'] = model.class_prior_.tolist()
    else:
        payload['coefficients'] = model.coef_.tolist()
        payload['intercept'] = model.intercept_.tolist()
    path.write_text(json.dumps(payload, allow_nan=False) + '\n')


def load_model(path):
    payload = json.loads(path.read_text())
    if payload['version'] != MODEL_VERSION:
        raise ValueError('Model artifact version is not supported')
    kind = payload['features']
    if kind == 'dummy':
        estimator = DummyClassifier(strategy='prior')
        estimator.classes_ = np.array([0,1]); estimator.n_classes_ = 2
        estimator.class_prior_ = np.array(payload['prior']); estimator.n_outputs_ = 1
    elif kind in ['position','pairwise']:
        estimator = LogisticRegression()
        estimator.classes_ = np.array([0,1])
        estimator.coef_ = np.asarray(payload['coefficients'], dtype=float)
        estimator.intercept_ = np.asarray(payload['intercept'], dtype=float)
        width = 200 if kind == 'position' else 18200
        if estimator.coef_.shape != (1,width) or estimator.intercept_.shape != (1,):
            raise ValueError('Invalid model parameter shape')
        if not np.isfinite(estimator.coef_).all() or not np.isfinite(estimator.intercept_).all():
            raise ValueError('Nonfinite model parameters')
        estimator.n_features_in_ = width
    else:
        raise ValueError('Unknown feature model')
    return {'version': MODEL_VERSION, 'features': kind, 'estimator': estimator}


def run_benchmark(args):
    if args.archive is not None:
        with ZipFile(args.archive) as archive:
            loaded = {
                name: load_split(f"{SPLIT_PATH}{name}.csv", archive=archive)
                for name in SPLIT_NAMES
            }
    else:
        loaded = {name: load_split(args.splits_dir / f"{name}.csv") for name in SPLIT_NAMES}
    splits = {name: data for name, (data, _) in loaded.items()}
    results = {
        "seed": args.seed,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "target": "high versus mid or low",
        "inputs": {name: info for name, (_, info) in loaded.items()},
        "split_audit": audit_splits(splits),
        "facs_mapping": None,
        "spr": None,
    }
    if args.facs is not None:
        results["inputs"]["facs"], results["facs_mapping"] = audit_facs(args.facs, splits)

    models = {
        "dummy": DummyClassifier(strategy="prior"),
        "position": LogisticRegression(C=1.0, solver="lbfgs", max_iter=500, random_state=args.seed),
    }
    if not args.skip_pairwise:
        models["pairwise"] = LogisticRegression(
            C=0.1, solver="lbfgs", max_iter=500, random_state=args.seed
        )
    results["model_settings"] = {name: model.get_params() for name, model in models.items()}
    train, val, test = (splits[name] for name in SPLIT_NAMES)
    val_scores = {}
    with threadpool_limits(limits=1):
        for kind, estimator in models.items():
            estimator.fit(sequence_features(train["seq"], kind), train["label"])
            if kind != "dummy" and estimator.n_iter_.max() >= estimator.max_iter:
                raise ValueError(f"{kind} logistic regression reached its iteration limit")
            val_scores[kind] = predict_scores(estimator, kind, val["seq"])

    val_metrics = {
        kind: classification_metrics(val["label"], scores) for kind, scores in val_scores.items()
    }
    selected = max(val_metrics, key=lambda kind: val_metrics[kind]["pr_auc_average_precision"])
    results["selection"] = {
        "metric": "validation_pr_auc_average_precision",
        "tie_break_order": list(models),
        "fit_split": "train",
    }
    results["selected_model"] = selected
    selected_model = {
        "version": MODEL_VERSION,
        "features": selected,
        "estimator": models[selected],
    }
    with threadpool_limits(limits=1):
        test_scores = {
            kind: predict_scores(estimator, kind, test["seq"]) for kind, estimator in models.items()
        }
        if args.spr is not None:
            spr_predictions, results["inputs"]["spr"], results["spr"] = evaluate_spr(
                args.spr, selected_model, splits
            )
    results["classification"] = {
        "validation": val_metrics,
        "test": {
            kind: classification_metrics(test["label"], scores) for kind, scores in test_scores.items()
        },
    }

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    results["outputs"] = {
        "metrics": "benchmark_metrics.json",
        "report": "benchmark_report.md",
        "model": "selected_model.json",
        "facs_predictions": "facs_predictions.csv",
    }
    facs_predictions = []
    for name, split, scores in (("val", val, val_scores), ("test", test, test_scores)):
        predictions = split[["seq", "class", "label", "cluster_centre"]].rename(
            columns={"seq": "cdrh3"}
        )
        predictions = predictions.assign(split=name, **{SCORE: scores[selected]})
        for kind, probabilities in scores.items():
            predictions[f"{SCORE}_{kind}"] = probabilities
        facs_predictions.append(predictions)
    pd.concat(facs_predictions, ignore_index=True).to_csv(out / "facs_predictions.csv", index=False)
    if args.spr is not None:
        spr_predictions.to_csv(out / "spr_predictions.csv", index=False)
        results["outputs"]["spr_predictions"] = "spr_predictions.csv"
    save_model(selected_model, out / "selected_model.json")
    (out / "benchmark_metrics.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    (out / "benchmark_report.md").write_text(render_report(results), encoding="utf-8")
    print(out / "benchmark_report.md")


def run_predict(args):
    artifact = load_model(args.model)
    if artifact["version"] != MODEL_VERSION:
        raise ValueError("Model artifact version is not supported")
    candidates, _ = read_csv(args.candidates)
    if args.sequence_column not in candidates:
        raise ValueError(f"Candidates CSV is missing column {args.sequence_column}")
    if SCORE in candidates:
        raise ValueError(f"Candidates CSV already contains the output column {SCORE}")
    if candidates.empty:
        raise ValueError("Candidates CSV has no rows")
    with threadpool_limits(limits=1):
        candidates[SCORE] = predict_scores(
            artifact["estimator"], artifact["features"], candidates[args.sequence_column]
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for column in candidates.select_dtypes(include=['object', 'str']).columns:
        candidates[column] = candidates[column].map(
            lambda value: "'" + value if isinstance(value, str) and value.lstrip().startswith(('=', '+', '-', '@')) else value
        )
    candidates.to_csv(args.output, index=False)
    print(args.output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    benchmark = commands.add_parser("benchmark", help="train on the published train split")
    source = benchmark.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive", type=Path, help="published Zenodo data.zip")
    source.add_argument(
        "--splits-dir", type=Path, help="directory containing train.csv, val.csv and test.csv"
    )
    benchmark.add_argument("--facs", type=Path, help="optional organizer FACS mapping to cross-check")
    benchmark.add_argument("--spr", type=Path, help="optional independent SPR CSV")
    benchmark.add_argument("--output-dir", type=Path, required=True)
    benchmark.add_argument("--seed", type=int, default=0)
    benchmark.add_argument("--skip-pairwise", action="store_true")
    benchmark.set_defaults(handler=run_benchmark)

    predict = commands.add_parser("predict", help="score candidate CDRH3s")
    predict.add_argument("--model", type=Path, required=True, help="selected_model.json")
    predict.add_argument("--candidates", type=Path, required=True, help="candidate CSV")
    predict.add_argument("--sequence-column", default="cdrh3")
    predict.add_argument("--output", type=Path, required=True)
    predict.set_defaults(handler=run_predict)

    args = parser.parse_args(argv)
    try:
        args.handler(args)
    except (ValueError, OSError, KeyError, BadZipFile) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
