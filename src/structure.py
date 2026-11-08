"""Загрузка комплексов антитело-антиген из RCSB PDB и разбор на цепи.

Что делает:
  - скачивает структуру (BinaryCIF, кэш в data/raw)
  - делит цепи на антитела (H/K/L по anarcii) и антиген
  - сопоставляет IMGT-позиции с номерами остатков в структуре
  - находит остатки на интерфейсе антитело-антиген по расстоянию
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import biotite.database.rcsb as rcsb
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np

from cdr_tools import IMGT_CDR_RANGES, NumberedChain, number_chains

DATA_RAW = Path(__file__).resolve().parent.parent / "data" / "raw"

# Остатки ближе этого расстояния (Å) считаем контактирующими.
# 4.5 Å — обычный порог для нековалентных контактов в анализе интерфейсов.
CONTACT_CUTOFF = 4.5


@dataclass
class Chain:
    """Одна полипептидная цепь структуры."""

    chain_id: str
    seq: str
    res_ids: np.ndarray  # номера остатков в структуре, параллельно seq
    numbered: NumberedChain | None = None  # None для антигена

    @property
    def is_antibody(self) -> bool:
        return self.numbered is not None

    @property
    def kind(self) -> str:
        if self.numbered is None:
            return "antigen"
        return {"H": "heavy", "K": "light", "L": "light"}.get(
            self.numbered.chain_type, self.numbered.chain_type
        )


@dataclass
class Complex:
    """Разобранный комплекс: цепи плюс сама структура."""

    pdb_id: str
    atoms: struc.AtomArray
    chains: dict[str, Chain] = field(default_factory=dict)

    @property
    def antibody_chains(self) -> dict[str, Chain]:
        return {c: ch for c, ch in self.chains.items() if ch.is_antibody}

    @property
    def antigen_chains(self) -> dict[str, Chain]:
        return {c: ch for c, ch in self.chains.items() if not ch.is_antibody}


def _chain_seq_and_ids(atoms: struc.AtomArray, chain_id: str) -> tuple[str, np.ndarray]:
    """Однобуквенная последовательность цепи и номера её остатков."""
    sub = atoms[atoms.chain_id == chain_id]
    res_ids, res_names = struc.get_residues(sub)
    seq = "".join(struc.info.one_letter_code(name) or "X" for name in res_names)
    return seq, res_ids


def load_complex(pdb_id: str, cache_dir: Path = DATA_RAW) -> Complex:
    """Скачать (или взять из кэша) структуру и разобрать её на цепи."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = rcsb.fetch(pdb_id, "bcif", str(cache_dir))
    atoms = pdbx.get_structure(pdbx.BinaryCIFFile.read(path), model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]

    chain_ids = list(struc.get_chains(atoms))
    seqs_ids = {c: _chain_seq_and_ids(atoms, c) for c in chain_ids}

    # Одним прогоном anarcii: распознанные цепи — антитела, остальные антиген.
    numbered = number_chains([seqs_ids[c][0] for c in chain_ids])
    by_seq = {n.seq: n for n in numbered}

    cx = Complex(pdb_id=pdb_id.upper(), atoms=atoms)
    for c in chain_ids:
        seq, res_ids = seqs_ids[c]
        cx.chains[c] = Chain(chain_id=c, seq=seq, res_ids=res_ids, numbered=by_seq.get(seq))
    return cx


def imgt_to_res_id(chain: Chain) -> dict[tuple[int, str], int]:
    """IMGT-позиция -> номер остатка в структуре.

    anarcii нумерует только вариабельный домен, поэтому отсчёт идёт
    от query_start, а gap-символы ('-') позиций в цепи не занимают.
    """
    if chain.numbered is None:
        raise ValueError(f"chain {chain.chain_id} is not an antibody chain")
    mapping: dict[tuple[int, str], int] = {}
    idx = chain.numbered.query_start
    for (pos, ins), aa in chain.numbered.numbering:
        if aa == "-":
            continue
        if idx < len(chain.res_ids):
            mapping[(pos, ins)] = int(chain.res_ids[idx])
        idx += 1
    return mapping


def cdr_res_ids(chain: Chain) -> dict[str, list[int]]:
    """Номера остатков структуры для CDR1-3 этой цепи."""
    mapping = imgt_to_res_id(chain)
    out: dict[str, list[int]] = {}
    for name, (lo, hi) in IMGT_CDR_RANGES.items():
        out[name] = [rid for (pos, _ins), rid in sorted(mapping.items()) if lo <= pos <= hi]
    return out


def interface_res_ids(
    cx: Complex, ab_chain_id: str, cutoff: float = CONTACT_CUTOFF
) -> set[int]:
    """Остатки цепи антитела, контактирующие с любой цепью антигена."""
    ab = cx.atoms[cx.atoms.chain_id == ab_chain_id]
    antigen_ids = list(cx.antigen_chains)
    if not antigen_ids:
        return set()
    ag = cx.atoms[np.isin(cx.atoms.chain_id, antigen_ids)]

    cell = struc.CellList(ag, cell_size=cutoff)
    contacts = cell.get_atoms(ab.coord, radius=cutoff)
    # get_atoms возвращает -1 там, где соседей нет
    touching = (contacts != -1).any(axis=1)
    return {int(r) for r in np.unique(ab.res_id[touching])}


def epitope_summary(cx: Complex, ab_chain_id: str) -> dict[str, list[int]]:
    """Для каждой CDR — какие её остатки реально трогают антиген."""
    chain = cx.chains[ab_chain_id]
    iface = interface_res_ids(cx, ab_chain_id)
    return {
        name: [rid for rid in rids if rid in iface]
        for name, rids in cdr_res_ids(chain).items()
    }
