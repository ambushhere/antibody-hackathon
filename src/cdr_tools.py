"""Минимальные утилиты для работы с антителами: IMGT-нумерация, CDR, мутации.

Нумерация через anarcii (не требует HMMER, в отличие от классического ANARCI).
IMGT-границы CDR: CDR1 27-38, CDR2 56-65, CDR3 105-117.
"""

from __future__ import annotations

from dataclasses import dataclass

from anarcii import Anarcii

IMGT_CDR_RANGES = {
    "CDR1": (27, 38),
    "CDR2": (56, 65),
    "CDR3": (105, 117),
}


@dataclass
class NumberedChain:
    """Одна пронумерованная цепь антитела."""

    seq: str
    chain_type: str  # 'H', 'K', 'L'
    score: float
    numbering: list[tuple[tuple[int, str], str]]

    @property
    def cdrs(self) -> dict[str, str]:
        """CDR1-3 по границам IMGT, без gap-символов."""
        out: dict[str, str] = {}
        for name, (lo, hi) in IMGT_CDR_RANGES.items():
            out[name] = "".join(
                aa for (pos, _ins), aa in self.numbering if lo <= pos <= hi and aa != "-"
            )
        return out

    def position_to_index(self, pos: int, ins: str = " ") -> int:
        """Индекс остатка в исходной последовательности по IMGT-позиции."""
        idx = 0
        for (p, i), aa in self.numbering:
            if aa == "-":
                continue
            if p == pos and i == ins:
                return idx
            idx += 1
        raise KeyError(f"IMGT position {pos}{ins!r} not present in this chain")


def number_chains(seqs: list[str], mode: str = "accuracy") -> list[NumberedChain]:
    """Пронумеровать список последовательностей. Нераспознанные пропускаются."""
    model = Anarcii(seq_type="antibody", mode=mode, verbose=False)
    result = model.number(seqs)
    chains: list[NumberedChain] = []
    for seq, (_key, rec) in zip(seqs, result.items()):
        if rec.get("error") not in (None, "None") or not rec.get("numbering"):
            continue
        chains.append(
            NumberedChain(
                seq=seq,
                chain_type=str(rec["chain_type"]),
                score=float(rec["score"]),
                numbering=rec["numbering"],
            )
        )
    return chains


def mutate(seq: str, index: int, new_aa: str) -> str:
    """Точечная мутация по 0-based индексу в последовательности."""
    if not 0 <= index < len(seq):
        raise IndexError(f"index {index} outside sequence of length {len(seq)}")
    return seq[:index] + new_aa + seq[index + 1 :]


def mutate_imgt(chain: NumberedChain, pos: int, new_aa: str, ins: str = " ") -> str:
    """Точечная мутация по IMGT-позиции."""
    return mutate(chain.seq, chain.position_to_index(pos, ins), new_aa)


if __name__ == "__main__":
    VH = (
        "EVQLQQSGAELVRPGALVKLSCKASGFNIKDYYMHWVKQRPEQGLEWIGLIDPENGNTIYDPKFQGK"
        "ASITADTSSNTAYLQLSSLTSEDTAVYYCARDNSYYFDYWGQGTTLTVSS"
    )
    (chain,) = number_chains([VH])
    print(f"chain_type={chain.chain_type} score={chain.score:.1f}")
    for name, cdr in chain.cdrs.items():
        print(f"  {name}: {cdr}")
    variant = mutate_imgt(chain, 107, "A")
    print(f"mutant differs at: {[i for i, (a, b) in enumerate(zip(chain.seq, variant)) if a != b]}")
