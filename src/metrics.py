"""Оценка вариантов: ESM-2 wt-marginals + простые метрики developability.

ESM wt-marginals (Meier et al. 2021, ESM-1v): скор мутации =
log p(mut | контекст WT) - log p(wt | контекст WT) в этой позиции.
Один forward pass на всю последовательность даёт скоры для всех
мутаций сразу — на CPU это единственный практичный вариант.
Знак: >0 модель считает мутацию правдоподобнее исходного остатка.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np
import torch
from Bio.SeqUtils.ProtParam import ProteinAnalysis

AA20 = "ACDEFGHIKLMNPQRSTVWY"

# Остатки-факторы риска developability в CDR:
#   N-гликозилирование NxS/T, деамидирование NG/NS, изомеризация DG,
#   окисление M/W, свободный C.
RISK_MOTIFS = {
    "деамидирование": ("NG", "NS", "NT"),
    "изомеризация": ("DG", "DS"),
}


@dataclass
class MutationScore:
    """Скор одной точечной мутации."""

    position: int  # 0-based индекс в последовательности
    wt: str
    mut: str
    esm_score: float

    @property
    def label(self) -> str:
        return f"{self.wt}{self.position + 1}{self.mut}"


@functools.lru_cache(maxsize=1)
def load_esm(model_name: str = "esm2_t33_650M_UR50D"):
    """Загрузить ESM-2 один раз на процесс (веса ~2.4 ГБ)."""
    import esm

    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    return model, alphabet


def esm_log_probs(seq: str, model_name: str = "esm2_t33_650M_UR50D") -> np.ndarray:
    """Матрица log-вероятностей (len(seq), 20) по колонкам AA20."""
    model, alphabet = load_esm(model_name)
    _, _, tokens = alphabet.get_batch_converter()([("query", seq)])
    with torch.no_grad():
        logits = model(tokens)["logits"]
    log_probs = torch.log_softmax(logits, dim=-1)[0]
    # +1: токен BOS в начале
    idx = [alphabet.get_idx(aa) for aa in AA20]
    return log_probs[1 : len(seq) + 1, idx].numpy()


def score_mutations(
    seq: str, positions: list[int], model_name: str = "esm2_t33_650M_UR50D"
) -> list[MutationScore]:
    """Полный скан 19 замен в каждой из позиций, один forward pass."""
    lp = esm_log_probs(seq, model_name)
    col = {aa: i for i, aa in enumerate(AA20)}
    out: list[MutationScore] = []
    for pos in positions:
        wt = seq[pos]
        if wt not in col:
            continue
        for mut in AA20:
            if mut == wt:
                continue
            out.append(
                MutationScore(
                    position=pos,
                    wt=wt,
                    mut=mut,
                    esm_score=float(lp[pos, col[mut]] - lp[pos, col[wt]]),
                )
            )
    return sorted(out, key=lambda m: m.esm_score, reverse=True)


def sequence_metrics(seq: str) -> dict[str, float]:
    """Простые метрики developability по последовательности."""
    pa = ProteinAnalysis(seq)
    return {
        "length": len(seq),
        "pI": round(pa.isoelectric_point(), 2),
        "charge_pH7": round(pa.charge_at_pH(7.0), 2),
        "gravy": round(pa.gravy(), 3),  # >0 гидрофобна => риск агрегации
        "aromaticity": round(pa.aromaticity(), 3),
        "instability": round(pa.instability_index(), 1),  # >40 считается нестабильной
    }


def liability_flags(seq: str, region: tuple[int, int] | None = None) -> list[str]:
    """Мотивы-риски developability. region — полуинтервал индексов (lo, hi)."""
    lo, hi = region or (0, len(seq))
    window = seq[lo:hi]
    flags: list[str] = []
    for name, motifs in RISK_MOTIFS.items():
        for m in motifs:
            start = window.find(m)
            if start != -1:
                flags.append(f"{name} {m}@{lo + start + 1}")
    if "C" in window:
        flags.append(f"свободный C@{lo + window.find('C') + 1}")
    return flags
