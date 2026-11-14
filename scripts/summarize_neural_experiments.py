"""Render every prespecified neural arm, including failed transfer experiments."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    benchmark = json.loads((ROOT / "results/neural-all-benchmark/metrics.json").read_text())
    rows = []
    for path in sorted((ROOT / "results/neural-runs").glob("*/metadata.json")):
        run = json.loads(path.read_text())
        config = run["configuration"]
        family = "augmented" if run.get("augmentation") else config["size"]
        for model, metrics in run["models"].items():
            key = f"{family}_{model}_seed{config['seed']}"
            row = {"model": key, "family": family, "architecture": model, "seed": config["seed"],
                   "parameters": metrics["parameter_count"], "selected_epoch": metrics["best_epoch"],
                   "training_seconds": metrics["training_seconds"],
                   "facs_test_ap": metrics["facs_test"]["pr_auc_average_precision"],
                   "facs_test_roc_auc": metrics["facs_test"]["roc_auc"]}
            for subset in ("pooled", "distant75"):
                result = benchmark["subsets"][subset]["metrics"][key]
                row[f"{subset}_spearman"] = result["spearman"]["estimate"]
                row[f"{subset}_rmse"] = result["rmse"]["estimate"]
            rows.append(row)
    table = pd.DataFrame(rows).sort_values(["family", "architecture", "seed"])
    destination = ROOT / "results/neural-summary"
    destination.mkdir(exist_ok=True)
    table.to_csv(destination / "experiments.csv", index=False)
    colors = {"small": "#8c96a1", "large_mlp": "#ce8b36", "large_transformer": "#287b8e",
              "augmented": "#8462a5"}
    groups = [("small", table[table.family.eq("small")]),
              ("large_mlp", table[table.family.eq("large") & table.architecture.eq("mlp")]),
              ("large_transformer", table[table.family.eq("large") & table.architecture.eq("transformer")]),
              ("augmented", table[table.family.eq("augmented")])]
    labels = {"small": "Small pilot, 5 epochs", "large_mlp": "MLP, 5.13M parameters",
              "large_transformer": "Transformer, 4.75M parameters", "augmented": "Transformer + pseudolabels"}
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharex=True, sharey=True)
    for ax, subset, title in zip(axes, ("pooled", "distant75"),
                                 ("158 measured SPR variants", "75 distant SPR variants")):
        for family, data in groups:
            ax.scatter(data.facs_test_ap, data[f"{subset}_spearman"], c=colors[family],
                       s=70, label=labels[family], edgecolors="white", linewidths=0.8)
            for _, row in data.iterrows():
                ax.annotate(str(row.seed), (row.facs_test_ap, row[f"{subset}_spearman"]),
                            xytext=(4, 4), textcoords="offset points", fontsize=8)
        baseline = benchmark["subsets"][subset]["metrics"]["position_pairwise_ridge"]["spearman"]["estimate"]
        ax.axhline(baseline, color="#394954", ls="--", lw=1.2, label="SPR pairwise ridge")
        ax.set_title(title)
        ax.set_xlabel("FACS test average precision (higher is better)")
        ax.grid(alpha=.18)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("SPR rank correlation (higher is better)")
    axes[0].set_ylim(0, .72)
    axes[0].set_xlim(.915, 1.001)
    figure.suptitle("Better FACS classification did not establish better affinity ranking", fontsize=14)
    handles, names = axes[0].get_legend_handles_labels()
    figure.legend(handles, names, loc="lower center", ncol=3, frameon=False, fontsize=9)
    figure.text(.5, .13, "Numbers identify seeds. Frozen encoder + nested SPR ridge; same 158 rows and group folds.",
                ha="center", fontsize=9, color="#45515c")
    figure.tight_layout(rect=(0, .17, 1, .94))
    output = ROOT / "output/neural-benchmark"
    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(output / f"facs-vs-spr.{extension}", dpi=160, facecolor="white")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
