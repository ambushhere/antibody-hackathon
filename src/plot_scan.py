"""График скана мутаций для слайда: тепловая карта позиция x аминокислота.

    PYTHONPATH=src .venv/bin/python src/plot_scan.py results/1MEL_A_scan.csv

Форма: знаковая величина по двум категориальным осям -> тепловая карта.
Цвет: расходящаяся палитра blue<->red с нейтральной серой серединой ровно
на нуле. Ноль здесь содержателен — это «так же, как у WT», поэтому
середина палитры привязана к нему (TwoSlopeNorm), а не к центру данных.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULTS = Path(__file__).resolve().parent.parent / "results"

# Аминокислоты сгруппированы по физико-химическим классам: соседние строки
# тогда сравнимы по смыслу, а не случайны как при алфавитном порядке.
AA_GROUPS = [
    ("отриц.", "DE"),
    ("полож.", "KRH"),
    ("полярные", "STNQ"),
    ("особые", "CGP"),
    ("гидрофоб.", "AVLIM"),
    ("аромат.", "FWY"),
]
AA_ORDER = [aa for _, block in AA_GROUPS for aa in block]

# Расходящаяся пара из палитры: blue <-> red, серая середина.
# Плечи заданы шагами одного тона с монотонной светлотой (проверяется ниже).
THEME = {
    "light": {
        "surface": "#fcfcfb",
        "ink": "#0b0b0b",
        "ink2": "#52514e",
        "muted": "#8a8984",
        "mid": "#f0efec",
        "neg": ["#0d366b", "#184f95", "#3987e5", "#9ec5f4", "#cde2fb"],
        "pos": ["#fbd5d4", "#f3a9a8", "#e34948", "#b02a2a", "#7a1a1a"],
    },
    "dark": {
        "surface": "#1a1a19",
        "ink": "#f0efec",
        "ink2": "#c3c2b7",
        "muted": "#8a8984",
        "mid": "#383835",
        "neg": ["#9ec5f4", "#6da7ec", "#3987e5", "#1c5cab", "#104281"],
        "pos": ["#5e1414", "#8c2020", "#e66767", "#f0a09f", "#f8cfce"],
    },
}


def _relative_luminance(hex_color: str) -> float:
    rgb = np.array(mcolors.to_rgb(hex_color))
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return float(lin @ (0.2126, 0.7152, 0.0722))


def build_cmap(theme: dict) -> mcolors.LinearSegmentedColormap:
    """Расходящаяся карта: одно тоновое плечо -> серый -> второе плечо.

    Проверяем монотонность светлоты внутри каждого плеча: для расходящейся
    палитры это замена проверки lightness band (та — про категориальные).
    """
    neg = list(theme["neg"])  # от крайнего полюса к середине
    pos = list(theme["pos"])  # от середины к крайнему полюсу
    mid_lum = _relative_luminance(theme["mid"])
    # Плечо ориентируем «полюс -> середина»: светлота обязана монотонно
    # двигаться и приходить к середине. Знак зависит от темы, поэтому
    # проверяем направление, а не возрастание.
    for arm, name in ((neg, "neg"), (pos[::-1], "pos")):
        lums = [_relative_luminance(c) for c in arm]
        direction = 1 if lums[-1] > lums[0] else -1
        if not all((b - a) * direction > 0 for a, b in zip(lums, lums[1:])):
            raise ValueError(f"светлота плеча {name} не монотонна: {lums}")
        if abs(lums[-1] - mid_lum) >= abs(lums[0] - mid_lum):
            raise ValueError(f"плечо {name} не приближается к середине: {lums} vs {mid_lum}")
    return mcolors.LinearSegmentedColormap.from_list(
        "esm_div", neg + [theme["mid"]] + pos, N=256
    )


def plot(csv_path: Path, mode: str = "light") -> Path:
    theme = THEME[mode]
    df = pd.read_csv(csv_path)
    if "cdr" not in df.columns:
        raise SystemExit("в CSV нет колонки cdr — перегенерируйте через src/mini_task.py")

    df["mut_aa"] = df["mutation"].str[-1]
    df["wt_aa"] = df["mutation"].str[0]
    positions = sorted(df["res_id"].unique())
    wt_of = dict(zip(df["res_id"], df["wt_aa"]))
    cdr_of = dict(zip(df["res_id"], df["cdr"]))

    grid = np.full((len(AA_ORDER), len(positions)), np.nan)
    row_of = {aa: i for i, aa in enumerate(AA_ORDER)}
    col_of = {p: j for j, p in enumerate(positions)}
    for r in df.itertuples():
        grid[row_of[r.mut_aa], col_of[r.res_id]] = r.esm_score

    vmax = max(0.2, float(np.nanmax(grid)))
    vmin = min(-0.2, float(np.nanmin(grid)))
    cmap = build_cmap(theme)
    cmap.set_bad(theme["mid"])  # WT-ячейки: «нет данных», это сам WT

    fig, ax = plt.subplots(figsize=(1.0 + 0.44 * len(positions), 7.4))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    mesh = ax.pcolormesh(
        np.arange(len(positions) + 1),
        np.arange(len(AA_ORDER) + 1),
        np.ma.masked_invalid(grid),
        cmap=cmap,
        norm=mcolors.TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax),
        edgecolors=theme["surface"],
        linewidth=1.5,  # зазор цвета поверхности между ячейками
    )

    # WT: точка в своей ячейке вместо пустоты
    for p in positions:
        aa = wt_of[p]
        if aa in row_of:
            ax.plot(col_of[p] + 0.5, row_of[aa] + 0.5, "o", ms=3.2,
                    color=theme["muted"], zorder=3)

    ax.set_xticks(np.arange(len(positions)) + 0.5)
    ax.set_xticklabels([f"{wt_of[p]}{p}" for p in positions], fontsize=9,
                       color=theme["ink2"], rotation=90)
    ax.set_yticks(np.arange(len(AA_ORDER)) + 0.5)
    ax.set_yticklabels(AA_ORDER, fontsize=9, color=theme["ink2"], family="monospace")
    ax.invert_yaxis()  # классы читаются сверху вниз, как в AA_GROUPS
    ax.tick_params(length=0, pad=4)
    for s in ax.spines.values():
        s.set_visible(False)

    # Разделители классов аминокислот и подписи классов справа
    pos_y = 0
    for label, block in AA_GROUPS:
        pos_y += len(block)
        if pos_y < len(AA_ORDER):
            ax.axhline(pos_y, color=theme["surface"], lw=3)
        ax.annotate(
            label,
            xy=(1.012, 1 - (pos_y - len(block) / 2) / len(AA_ORDER)),
            xycoords="axes fraction", fontsize=8, color=theme["muted"],
            va="center", ha="left", annotation_clip=False,
        )

    # Полосы CDR под осью X: какая позиция из какой петли
    start = 0
    for j, p in enumerate(positions + [None]):
        same = p is not None and cdr_of[p] == cdr_of[positions[start]]
        if not same:
            mid = (start + j) / 2
            ax.annotate(
                cdr_of[positions[start]].replace("CDR", "CDR-"),
                xy=(mid, -0.108), xycoords=("data", "axes fraction"),
                ha="center", va="top", fontsize=9.5, color=theme["ink2"],
                annotation_clip=False,
            )
            ax.plot([start + 0.12, j - 0.12], [-0.098, -0.098],
                    transform=ax.get_xaxis_transform(), clip_on=False,
                    color=theme["muted"], lw=1.2)
            start = j

    # Единственная прямая подпись: лучший вариант. Числа на каждой ячейке — анти-паттерн.
    best = df.loc[df["esm_score"].idxmax()]
    bx, by = col_of[best.res_id] + 0.5, row_of[best.mut_aa] + 0.5
    # Ось Y перевёрнута, поэтому «вверх» на картинке — это меньшее значение.
    dy = -2.8 if by > 3 else 2.8
    dx = 1.6 if bx < len(positions) - 4 else -1.6
    ax.annotate(
        f"{best.mutation}  {best.esm_score:+.2f}",
        xy=(bx, by), xytext=(bx + dx, by + dy),
        ha="left" if dx > 0 else "right",
        va="center",
        fontsize=10, color=theme["ink"], fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.32", fc=theme["surface"],
                  ec=theme["muted"], lw=0.6, alpha=0.95),
        arrowprops=dict(arrowstyle="-", color=theme["ink2"], lw=1.1,
                        shrinkA=2, shrinkB=4),
    )

    n_pos = int((df["esm_score"] > 0).sum())
    stem = csv_path.stem.replace("_scan", "")
    ax.set_title(
        f"{stem}: скан замен в остатках CDR на интерфейсе с антигеном",
        fontsize=13, color=theme["ink"], pad=26, loc="left", fontweight="bold",
    )
    ax.annotate(
        f"{len(df)} вариантов · {len(positions)} позиций · "
        f"{n_pos} с положительным скором · точка = остаток WT",
        xy=(0, 1.012), xycoords="axes fraction", fontsize=9.5, color=theme["ink2"],
    )

    cb = fig.colorbar(mesh, ax=ax, pad=0.16, fraction=0.028, extend="neither")
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, colors=theme["ink2"], labelsize=9)
    cb.set_label("ESM wt-marginal: log-odds замены против WT",
                 color=theme["ink2"], fontsize=9.5)
    ticks = sorted({round(vmin, 1), 0.0, round(vmax, 1)})
    cb.set_ticks(ticks)
    cb.set_ticklabels(["хуже WT" if t == ticks[0] else
                       ("как WT" if t == 0 else "лучше WT") for t in ticks])

    out = RESULTS / f"{stem}_heatmap_{mode}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor=theme["surface"])
    plt.close(fig)
    return out


if __name__ == "__main__":
    args = sys.argv[1:] or [str(RESULTS / "1N8Z_B_scan.csv"), str(RESULTS / "1MEL_A_scan.csv")]
    for a in args:
        for mode in ("light", "dark"):
            print("->", plot(Path(a), mode).relative_to(RESULTS.parent))
