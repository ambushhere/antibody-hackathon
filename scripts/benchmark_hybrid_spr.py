"""Evaluate three fixed blends of aligned, frozen outer-fold SPR predictions."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src import her2_interaction_benchmark as bench
from src import her2_spr_benchmark as spr

INPUT = Path("results/neural-all-benchmark/oof_predictions.csv")
INPUT_SHA256 = "95e8dfbbfeb2a7d6daef8cc4582252ef3a341c29f1869eb4503977ec8140150e"
BASELINES = (bench.POSITION, bench.PAIRWISE)
COMPONENTS = tuple(f"large_transformer_seed{seed}" for seed in (0, 1, 2))
RULES = ("transformer_seed_mean", "position_transformer_half", "pairwise_transformer_half")


def fixed_blends(frame):
    values = frame[[*BASELINES, *COMPONENTS]].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite component prediction")
    mean = values[:, 2:].mean(axis=1)
    return dict(zip(RULES, (mean, (values[:, 0] + mean) / 2, (values[:, 1] + mean) / 2)))


def run(path, output):
    if hashlib.sha256(path.read_bytes()).hexdigest() != INPUT_SHA256:
        raise ValueError("Frozen OOF fingerprint differs from protocol")
    frame = pd.read_csv(path)
    if len(frame) != 158 or frame.cdrh3.duplicated().any() or frame.groupby("group").fold.nunique().max() != 1:
        raise ValueError("Frozen cohort or fold membership mismatch")
    labels = frame[spr.TARGET].to_numpy()
    groups = frame.group.to_numpy()
    mask = frame.distant75.to_numpy()
    if mask.dtype.kind != "b" or mask.sum() != 75:
        raise ValueError("Frozen distant subset differs from protocol")
    scores = {name: frame[name].to_numpy() for name in BASELINES}
    scores.update(fixed_blends(frame))
    contrasts = [(rule, baseline) for rule in RULES for baseline in BASELINES]
    result = {
        "status": "retrospective diagnostic on previously viewed outcomes; no candidate promotion",
        "input": {"path": str(path), "sha256": INPUT_SHA256},
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "rules": dict(zip(RULES, ("Mean of seeds 0,1,2", "0.5 position + 0.5 seed mean", "0.5 pairwise + 0.5 seed mean"))),
        "bootstrap": {"draws": spr.BOOTSTRAP_DRAWS, "seed": spr.SEED,
                      "method": "Paired group resampling of fixed OOF predictions, no refits or multiplicity adjustment"},
        "subsets": {
            "pooled": bench.summarize(labels, scores, groups, contrasts),
            "distant75": bench.summarize(labels[mask], {k: v[mask] for k, v in scores.items()}, groups[mask], contrasts),
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    frame.assign(**fixed_blends(frame)).to_csv(output / "oof_predictions.csv", index=False)
    print(json.dumps({subset: {name: metrics["spearman"]["estimate"] for name, metrics in body["metrics"].items()}
                      for subset, body in result["subsets"].items()}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=Path("results/hybrid-spr"))
    args = parser.parse_args()
    run(args.input, args.output)
