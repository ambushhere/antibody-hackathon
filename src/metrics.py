"""Оценка вариантов: скоры ESM-2 плюс простые метрики developability.

Модель здесь ESM-2 (`esm2_t33_650M_UR50D`). Схема скоринга взята из работы
Meier et al. 2021, где её предложили для ESM-1v, но сама модель другая.

Две схемы для мутаций:
  wt-marginal: один forward pass по WT, скор = log p(mut) - log p(wt) в позиции.
    Дёшево, даёт все мутации сразу, единственный практичный вариант на CPU.
  masked-marginal: позиция маскируется перед предсказанием. По Meier et al.
    точнее примерно на 1% абсолютного Spearman, стоит L проходов. На GPU это
    секунды, поэтому при наличии CUDA или MPS лучше брать её.

Для сравнения ЦЕЛЫХ последовательностей разной длины есть `sequence_pll`,
и она нормирована по длине. Это не косметика: в раунде 2 конкурса Adaptyv
по EGFR организаторы не нормировали pseudo-likelihood, отбор уехал в сторону
коротких дизайнов около 50 остатков, и победитель провалился в промежуточном
лидерборде на 300-600 место. Ненормированный скор сравнивает длины, не качество.

ЧТО ЭТИ СКОРЫ ЗНАЧАТ. Насколько остаток типичен для этого места по корпусу
обучения. Это НЕ аффинность и НЕ ΔΔG. Для антител есть отдельная проблема:
последовательности в обучающих наборах почти германлайновые, а мутации,
улучшающие аффинность, как правило уходят ОТ германлайна (Olsen et al. 2024).
То есть модель систематически недооценивает именно то, что мы ищем. Поэтому
скор годится как слабый признак и как фильтр явно неправдоподобного, но не
как ранжирование по силе связывания.

Корреляция с аффинностью для ориентира: на DMS анти-лизоцима D44.1 Spearman
у ESM-2 около 0.264, у AntiFold около 0.418. Оба слабые, оба на одном антителе.
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


def pick_device() -> torch.device:
    """CUDA, затем MPS, иначе CPU. На VM это H100, на макбуке MPS."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


@functools.lru_cache(maxsize=1)
def load_esm(model_name: str = "esm2_t33_650M_UR50D"):
    """Загрузить ESM-2 один раз на процесс (веса ~2.4 ГБ) и положить на устройство."""
    import esm

    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    model.to(pick_device())
    return model, alphabet


def _tokens(seq: str, alphabet) -> torch.Tensor:
    _, _, tokens = alphabet.get_batch_converter()([("query", seq)])
    return tokens.to(pick_device())


def esm_log_probs(seq: str, model_name: str = "esm2_t33_650M_UR50D") -> np.ndarray:
    """Матрица log-вероятностей (len(seq), 20) по колонкам AA20, один forward pass."""
    model, alphabet = load_esm(model_name)
    with torch.no_grad():
        logits = model(_tokens(seq, alphabet))["logits"]
    log_probs = torch.log_softmax(logits, dim=-1)[0]
    # +1: токен BOS в начале
    idx = [alphabet.get_idx(aa) for aa in AA20]
    return log_probs[1 : len(seq) + 1, idx].float().cpu().numpy()


def esm_masked_log_probs(
    seq: str,
    positions: list[int] | None = None,
    model_name: str = "esm2_t33_650M_UR50D",
    batch_size: int = 16,
) -> dict[int, np.ndarray]:
    """Masked-marginal: по одному проходу на позицию, позиция замаскирована.

    Возвращает {позиция: вектор log-вероятностей по AA20}. Точнее wt-marginal
    (Meier et al. 2021), но стоит len(positions) проходов. Для VHH на H100
    это секунды, на CPU заметно дольше.
    """
    model, alphabet = load_esm(model_name)
    device = pick_device()
    idx = [alphabet.get_idx(aa) for aa in AA20]
    targets = list(range(len(seq))) if positions is None else sorted(set(positions))

    base = _tokens(seq, alphabet)
    out: dict[int, np.ndarray] = {}
    for start in range(0, len(targets), batch_size):
        chunk = targets[start : start + batch_size]
        batch = base.repeat(len(chunk), 1).clone()
        for row, pos in enumerate(chunk):
            batch[row, pos + 1] = alphabet.mask_idx  # +1: токен BOS
        with torch.no_grad():
            logits = model(batch)["logits"]
        lp = torch.log_softmax(logits, dim=-1)
        for row, pos in enumerate(chunk):
            out[pos] = lp[row, pos + 1, idx].float().cpu().numpy()
    return out


def score_mutations(
    seq: str,
    positions: list[int],
    model_name: str = "esm2_t33_650M_UR50D",
    scheme: str = "wt-marginal",
) -> list[MutationScore]:
    """Полный скан 19 замен в каждой из позиций.

    scheme: "wt-marginal" (один проход, дёшево) или "masked-marginal"
    (проход на позицию, точнее). Знак: >0 модель считает мутацию
    правдоподобнее исходного остатка. Это НЕ предсказание аффинности.
    """
    col = {aa: i for i, aa in enumerate(AA20)}

    if scheme == "masked-marginal":
        masked = esm_masked_log_probs(seq, positions, model_name)
        rows = {pos: masked[pos] for pos in masked}
    elif scheme == "wt-marginal":
        lp = esm_log_probs(seq, model_name)
        rows = {pos: lp[pos] for pos in positions}
    else:
        raise ValueError(f"неизвестная схема: {scheme}")

    out: list[MutationScore] = []
    for pos in sorted(rows):
        wt = seq[pos]
        if wt not in col:
            continue
        row = rows[pos]
        for mut in AA20:
            if mut == wt:
                continue
            out.append(
                MutationScore(
                    position=pos,
                    wt=wt,
                    mut=mut,
                    esm_score=float(row[col[mut]] - row[col[wt]]),
                )
            )
    return sorted(out, key=lambda m: m.esm_score, reverse=True)


def sequence_pll(
    seq: str, model_name: str = "esm2_t33_650M_UR50D", normalize: bool = True
) -> float:
    """Pseudo-log-likelihood последовательности, по умолчанию нормированный по длине.

    Для сравнения целых дизайнов между собой. `normalize=True` делит на длину,
    и это обязательно, когда длины разные: ненормированный скор монотонно падает
    с длиной, поэтому сравнивает длины, а не качество. Именно на этом в раунде 2
    конкурса Adaptyv отбор уехал к коротким дизайнам.

    Считается по masked-marginal, то есть len(seq) проходов. На CPU медленно.
    """
    masked = esm_masked_log_probs(seq, None, model_name)
    col = {aa: i for i, aa in enumerate(AA20)}
    total = 0.0
    counted = 0
    for pos, row in masked.items():
        aa = seq[pos]
        if aa not in col:
            continue  # нестандартный остаток пропускаем, иначе скор не определён
        total += float(row[col[aa]])
        counted += 1
    if counted == 0:
        return float("nan")
    return total / counted if normalize else total


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


def liability_flags(
    seq: str,
    region: tuple[int, int] | None = None,
    bonded_cys_idx: set[int] | None = None,
) -> list[str]:
    """Мотивы-риски developability. region — полуинтервал индексов (lo, hi).

    bonded_cys_idx — 0-based индексы цистеинов, уже связанных дисульфидом
    (из структуры). Они не помечаются как свободные: связанный цистеин
    нормален и часто структурно необходим.
    """
    lo, hi = region or (0, len(seq))
    lo, hi = max(0, lo), min(len(seq), hi)
    window = seq[lo:hi]
    bonded = bonded_cys_idx or set()
    flags: list[str] = []
    for name, motifs in RISK_MOTIFS.items():
        for m in motifs:
            start = window.find(m)
            if start != -1:
                flags.append(f"{name} {m}@{lo + start + 1}")
    for off, aa in enumerate(window):
        if aa == "C" and (lo + off) not in bonded:
            flags.append(f"свободный C@{lo + off + 1}")
            break
    return flags
