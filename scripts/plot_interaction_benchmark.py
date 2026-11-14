"""Render the frozen interaction benchmark, without selecting a winning model."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODELS = [
    ("p_facs_high", "FACS-скор", "#8493a4"),
    ("position_ridge", "Отдельные позиции", "#425b76"),
    ("position_pairwise_ridge", "Позиции + пары", "#236b9a"),
    ("cdrh3_mean", "AbLang2: вся петля", "#21806c"),
    ("wt_contact_mean", "AbLang2: контакты WT", "#a57c2d"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, default=Path("results/ablang-benchmark/metrics.json"))
    parser.add_argument("--output", type=Path, default=Path("output/interaction-benchmark"))
    args = parser.parse_args()
    data = json.loads(args.metrics.read_text())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, cohort, title in zip(axes, ("pooled", "distant75"),
                                 ("Все 158 вариантов", "75 вариантов дальше от обучения")):
        metrics = data["subsets"][cohort]["metrics"]
        for row, (name, label, color) in enumerate(MODELS):
            result = metrics[name]["spearman"]
            low, high = result["95_percent"]
            value = result["estimate"]
            ax.plot([low, high], [row, row], color=color, linewidth=3, alpha=.6)
            ax.scatter([value], [row], color=color, s=65, zorder=3)
            ax.text(value, row-.22, f"{value:.3f}", ha="center", color=color, fontsize=10)
        ax.set_title(title, fontsize=13, pad=15)
        ax.set_xlim(-.08, .85)
        ax.set_xlabel("Согласие с измеренным порядком (Spearman ρ)", fontsize=10)
        ax.grid(axis="x", color="#e3e8ed")
        ax.axvline(0, color="#9ca8b4", linestyle=":", linewidth=1)
        ax.set_yticks(range(len(MODELS)), [item[1] for item in MODELS])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
    axes[0].invert_yaxis()
    fig.suptitle("Трансформер пока не показал убедительного выигрыша над парами", fontsize=15, y=.98)
    fig.text(.02, .025, "Линии: 95% групповые bootstrap-интервалы. Интервалы разниц AbLang2 с парной моделью включают 0.\n"
             "Ретроспективный анализ уже просмотренных данных. Корреляция — не процент правильных ответов.", fontsize=9, color="#52616f")
    fig.tight_layout(rect=(0, .14, 1, .93))
    args.output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(args.output / f"comparison.{extension}", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
