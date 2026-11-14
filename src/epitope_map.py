"""Карта эпитопов на HER2: какие остатки мишени занимают известные антитела.

Считает контакты из координат PDB напрямую, без литературы. Для каждой структуры
берёт цепь антигена и цепи антитела, находит остатки антигена ближе CUTOFF к любому
атому антитела, переводит их в нумерацию UniProt P04626 и сопоставляет с доменами ECD.

Зачем: чтобы выбрать свободный эпитоп, надо знать занятые. Литературные списки
остатков расходятся по схемам нумерации, а координаты однозначны.

Запуск:
    python src/epitope_map.py
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import biotite.database.rcsb as rcsb
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np

DATA_RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
CONTACT_CUTOFF = 4.5

# Границы доменов ECD в нумерации предшественника UniProt P04626.
# Сигнальный пептид 1-22, ECD 23-652, TM 653-675.
# Зрелая нумерация = нумерация предшественника минус 22.
ECD_DOMAINS = [
    ("I", 23, 217),
    ("II", 218, 341),
    ("III", 342, 510),
    ("IV", 511, 652),
]

# Сайты N-гликозилирования в нумерации предшественника (UniProt P04626).
# Гликан физически закрывает участок, поэтому эпитоп рядом с ним рискованный.
GLYCOSYLATION_SITES = [68, 124, 187, 259, 530, 571, 629]


@dataclass
class Case:
    """Структура комплекса: что считаем и как называются цепи."""

    pdb_id: str
    antigen_chain: str
    antibody_chains: tuple[str, ...]
    label: str


CASES = [
    Case("1N8Z", "C", ("A", "B"), "трастузумаб (Fab)"),
    Case("1S78", "A", ("C", "D"), "пертузумаб (Fab), копия 1"),
    Case("1S78", "B", ("E", "F"), "пертузумаб (Fab), копия 2"),
    Case("6OGE", "A", ("B", "C"), "пертузумаб (крио-ЭМ)"),
    Case("6OGE", "A", ("D", "E"), "трастузумаб (крио-ЭМ)"),
    # Цепи определены по составу: у нанотел около 110-120 остатков, у ECD 400-600.
    Case("5MY6", "A", ("B",), "нанотело против HER2"),
    Case("9MTE", "C", ("A",), "нанотело 27A05"),
    Case("9MTE", "C", ("B",), "нанотело 47D05"),
    Case("9MTX", "B", ("A",), "нанотело 29E09"),
]


def load(pdb_id: str) -> struc.AtomArray:
    """Скачать структуру (с кэшем) и вернуть первую модель без воды и гетероатомов."""
    DATA_RAW.mkdir(parents=True, exist_ok=True)
    path = rcsb.fetch(pdb_id, "bcif", DATA_RAW)
    block = pdbx.BinaryCIFFile.read(path)
    atoms = pdbx.get_structure(block, model=1)
    return atoms[struc.filter_amino_acids(atoms) & ~atoms.hetero]


def chain_report(atoms: struc.AtomArray) -> list[tuple[str, int, int, int]]:
    """Для каждой цепи: id, число остатков, первый и последний номер остатка."""
    out = []
    for cid in struc.get_chains(atoms):
        sub = atoms[atoms.chain_id == cid]
        res = np.unique(sub.res_id)
        out.append((cid, len(res), int(res.min()), int(res.max())))
    return out


def contact_residues(
    atoms: struc.AtomArray, antigen_chain: str, antibody_chains: tuple[str, ...]
) -> list[tuple[int, str]]:
    """Остатки антигена ближе CONTACT_CUTOFF к любому атому антитела.

    Возвращает список (номер остатка, трёхбуквенный код), отсортированный по номеру.
    """
    ag = atoms[atoms.chain_id == antigen_chain]
    ab_mask = np.isin(atoms.chain_id, list(antibody_chains))
    ab = atoms[ab_mask]
    if ag.array_length() == 0 or ab.array_length() == 0:
        return []

    # Пространственный индекс по атомам антитела, затем запрос для каждого атома антигена.
    cell = struc.CellList(ab, cell_size=CONTACT_CUTOFF)
    hits = cell.get_atoms(ag.coord, radius=CONTACT_CUTOFF)
    touching = (hits != -1).any(axis=1)

    found = {}
    for res_id, res_name in zip(ag.res_id[touching], ag.res_name[touching]):
        found[int(res_id)] = str(res_name)
    return sorted(found.items())


def domain_of(res_id: int) -> str:
    """Домен ECD по номеру остатка в нумерации предшественника."""
    for name, lo, hi in ECD_DOMAINS:
        if lo <= res_id <= hi:
            return name
    return "вне ECD"


def summarize(res_ids: list[int]) -> dict[str, list[int]]:
    """Сгруппировать остатки по доменам ECD."""
    groups: dict[str, list[int]] = {}
    for r in res_ids:
        groups.setdefault(domain_of(r), []).append(r)
    return groups


def ranges(nums: list[int]) -> str:
    """Сжать список номеров в диапазоны: [1,2,3,7] -> '1-3, 7'."""
    if not nums:
        return "нет"
    nums = sorted(nums)
    out, start, prev = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        out.append(f"{start}-{prev}" if start != prev else f"{start}")
        start = prev = n
    out.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ", ".join(out)


def main() -> int:
    occupied: set[int] = set()

    for case in CASES:
        if case.antigen_chain is None:
            # Цепи неизвестны: печатаем состав, чтобы определить их глазами.
            print(f"\n=== {case.pdb_id}: {case.label} — состав цепей")
            try:
                atoms = load(case.pdb_id)
            except Exception as exc:  # структура может отсутствовать или не читаться
                print(f"  не удалось загрузить: {exc}")
                continue
            for cid, n, lo, hi in chain_report(atoms):
                print(f"  цепь {cid}: {n} остатков, номера {lo}-{hi}")
            continue

        print(f"\n=== {case.pdb_id}: {case.label}")
        try:
            atoms = load(case.pdb_id)
        except Exception as exc:
            print(f"  не удалось загрузить: {exc}")
            continue

        contacts = contact_residues(atoms, case.antigen_chain, case.antibody_chains)
        if not contacts:
            print(f"  контактов не найдено (проверить имена цепей)")
            for cid, n, lo, hi in chain_report(atoms):
                print(f"    цепь {cid}: {n} остатков, номера {lo}-{hi}")
            continue

        res_ids = [r for r, _ in contacts]
        occupied.update(res_ids)
        print(f"  остатков в контакте: {len(res_ids)}  (порог {CONTACT_CUTOFF} Å)")
        for dom, nums in summarize(res_ids).items():
            print(f"  домен {dom}: {ranges(nums)}")
        near_glycan = [g for g in GLYCOSYLATION_SITES if any(abs(g - r) <= 8 for r in res_ids)]
        if near_glycan:
            print(f"  рядом сайты гликозилирования: {near_glycan}")

    if occupied:
        print("\n=== свободные участки ECD (не касается ни одно из посчитанных антител)")
        print(f"занято всего остатков: {len(occupied)}")
        for name, lo, hi in ECD_DOMAINS:
            free = [r for r in range(lo, hi + 1) if r not in occupied]
            taken = [r for r in range(lo, hi + 1) if r in occupied]
            pct = 100 * len(taken) / (hi - lo + 1)
            print(f"\nдомен {name} ({lo}-{hi}, зрелые {lo-22}-{hi-22}): занято {pct:.1f}%")
            if taken:
                print(f"  занятые: {ranges(taken)}")
            else:
                print(f"  занятых нет, домен свободен целиком")

    print(
        "\nВАЖНО: нумерация остатков берётся из структуры PDB. Для HER2 она обычно "
        "совпадает с нумерацией предшественника UniProt P04626 (сигнальный пептид 1-22), "
        "но это надо сверять для каждой структуры: диапазоны цепей напечатаны выше.\n"
        "Отсутствие контакта в этих структурах не означает, что участок доступен in vivo: "
        "гликаны, мембрана и конформация тетеринга тоже закрывают поверхность."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
