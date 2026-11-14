"""Карта контактов антител с HER2 по первой модели каждой структуры PDB.

Геометрия: тяжёлые атомы аминокислот, порог 4.5 Å. Нумерация: авторские
остатки mmCIF → label_seq_id полимерной сущности → UniProt P04626 по SIFTS.
Только сопоставленные контакты используются в аннотациях и сводке.

Запуск: python src/epitope_map.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict
from dataclasses import dataclass
from functools import cache
from io import StringIO
from urllib.request import urlopen

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from biotite import InvalidFileError

CONTACT_CUTOFF = 4.5
ACCESSION = "P04626"
CANONICAL_LENGTH = 1255
CIF_URL = "https://files.rcsb.org/download/{}.cif"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}"
LABEL_FIELDS = ("label_asym_id", "label_entity_id", "label_seq_id")

# Границы доменов ECD в нумерации предшественника UniProt P04626.
ECD_DOMAINS = [
    ("I", 23, 217),
    ("II", 218, 341),
    ("III", 342, 510),
    ("IV", 511, 652),
]

GLYCOSYLATION_SITES = [68, 124, 187, 259, 530, 571, 629]

# (auth_seq_id, pdbx_PDB_ins_code), scoped to one author chain.
AuthorResidue = tuple[int | None, str]


@dataclass
class Case:
    """Структура комплекса и авторские идентификаторы цепей."""

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
    Case("5MY6", "A", ("B",), "нанотело против HER2"),
    Case("9MTE", "C", ("A",), "нанотело 27A05"),
    Case("9MTE", "C", ("B",), "нанотело 47D05"),
    Case("9MTX", "B", ("A",), "нанотело 29E09"),
]


@cache
def read_public(url: str) -> bytes:
    with urlopen(url, timeout=20) as response:
        return response.read()


def protein_atoms(cif: bytes) -> struc.AtomArray:
    """Read author identifiers and unmodified label fields from the same atom rows."""
    block = pdbx.CIFFile.read(StringIO(cif.decode("utf-8"))).block
    atom_site = block["atom_site"]
    if not {"auth_asym_id", "auth_seq_id"} <= set(atom_site):
        raise ValueError("В mmCIF отсутствуют авторские идентификаторы цепей или остатков")
    atoms = pdbx.get_structure(
        block,
        model=1,
        extra_fields=["auth_seq_id", *(field for field in LABEL_FIELDS if field in atom_site)],
        use_author_fields=True,
    )
    return atoms[struc.filter_amino_acids(atoms) & ~atoms.hetero]


def author_residue(number: str, insertion: str) -> AuthorResidue:
    return (None if number in {"", ".", "?"} else int(number), str(insertion))


def contact_residues(
    atoms: struc.AtomArray, antigen_chain: str, antibody_chains: tuple[str, ...]
) -> set[AuthorResidue]:
    """Find author residues with a heavy atom within the cutoff of the antibody."""
    heavy = ~np.isin(np.char.upper(atoms.element), ["H", "D", "T"])
    ag = atoms[heavy & (atoms.chain_id == antigen_chain)]
    ab = atoms[heavy & np.isin(atoms.chain_id, antibody_chains)]
    if not len(ag) or not len(ab):
        return set()

    cell = struc.CellList(ab, cell_size=CONTACT_CUTOFF)
    hits = cell.get_atoms(ag.coord, radius=CONTACT_CUTOFF)
    touching = (hits != -1).any(axis=1)
    return {author_residue(number, ins) for number, ins in zip(ag.auth_seq_id[touching], ag.ins_code[touching])}


def author_to_label(
    antigen: struc.AtomArray,
) -> tuple[str, str, dict[AuthorResidue, int]]:
    """Resolve one polymer instance and retain only unambiguous residue assignments."""
    if not set(LABEL_FIELDS) <= set(antigen.get_annotation_categories()):
        raise ValueError("В mmCIF отсутствуют идентификаторы label_asym_id, label_entity_id или label_seq_id")
    instances = set(zip(antigen.label_asym_id, antigen.label_entity_id))
    if len(instances) != 1:
        raise ValueError("Авторской цепи не соответствует ровно один экземпляр полимера mmCIF")
    label_chain, entity_id = map(str, instances.pop())
    if label_chain in {"", ".", "?"} or not entity_id.isdecimal() or int(entity_id) < 1:
        raise ValueError("В mmCIF отсутствует идентификатор экземпляра или сущности полимера")

    labels_by_author = defaultdict(set)
    authors_by_label = defaultdict(set)
    for number, ins, label_seq_id in set(zip(antigen.auth_seq_id, antigen.ins_code, antigen.label_seq_id)):
        author = author_residue(number, ins)
        label = int(label_seq_id) if label_seq_id.isdecimal() else None
        labels_by_author[author].add(label)
        authors_by_label[label].add(author)

    mapping = {}
    for author, labels in labels_by_author.items():
        if author[0] is None or len(labels) != 1:
            continue
        label = next(iter(labels))
        if label is not None and label > 0 and len(authors_by_label[label]) == 1:
            mapping[author] = label
    return label_chain, entity_id, mapping


def sifts_mapping(record: dict, case: Case, label_chain: str, entity_id: str) -> dict[int, int]:
    """Expand verified RCSB SIFTS regions from entity label positions to P04626."""
    identifiers = record["rcsb_polymer_entity_container_identifiers"]
    if (
        record["rcsb_id"] != f"{case.pdb_id}_{entity_id}"
        or identifiers["entity_id"] != entity_id
        or label_chain not in identifiers["asym_ids"]
        or case.antigen_chain not in identifiers["auth_asym_ids"]
    ):
        raise ValueError("RCSB и mmCIF расходятся в идентификаторах полимера")

    alignments = [
        alignment for alignment in record.get("rcsb_polymer_entity_align", [])
        if alignment["provenance_source"] == "SIFTS"
        and alignment["reference_database_name"] == "UniProt"
        and alignment["reference_database_accession"] == ACCESSION
    ]
    if len(alignments) != 1 or not alignments[0].get("aligned_regions"):
        raise ValueError(f"Для сущности нет однозначного выравнивания SIFTS к {ACCESSION}")

    mapping: dict[int, int] = {}
    seen: set[int] = set()
    for region in alignments[0]["aligned_regions"]:
        label_start = region["entity_beg_seq_id"]
        canonical_start = region["ref_beg_seq_id"]
        length = region["length"]
        if (
            any(type(value) is not int or value < 1 for value in (label_start, canonical_start, length))
            or canonical_start + length - 1 > CANONICAL_LENGTH
        ):
            raise ValueError("Некорректный диапазон SIFTS")
        for offset in range(length):
            label_pos = label_start + offset
            canonical_pos = canonical_start + offset
            if label_pos in mapping or canonical_pos in seen:
                raise ValueError("Перекрывающиеся диапазоны SIFTS")
            mapping[label_pos] = canonical_pos
            seen.add(canonical_pos)
    return mapping


def domain_of(res_id: int) -> str:
    for name, lo, hi in ECD_DOMAINS:
        if lo <= res_id <= hi:
            return name
    return "вне ECD"


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


def author_list(residues: set[AuthorResidue]) -> str:
    numbered = sorted((number, ins) for number, ins in residues if number is not None)
    return ", ".join(
        [f"{number}{ins}" for number, ins in numbered]
        + [f"?{ins}" for ins in sorted(ins for number, ins in residues if number is None)]
    )


def mapped_contacts(
    atoms: struc.AtomArray, case: Case, contacts: set[AuthorResidue]
) -> dict[AuthorResidue, int]:
    antigen = atoms[atoms.chain_id == case.antigen_chain]
    label_chain, entity_id, by_author = author_to_label(antigen)
    url = ENTITY_URL.format(case.pdb_id, entity_id)
    payload = read_public(url)
    print(f"  RCSB SIFTS: {url} SHA-256 {hashlib.sha256(payload).hexdigest()}")
    by_label = sifts_mapping(json.loads(payload), case, label_chain, entity_id)
    print(f"  mmCIF: auth_asym_id={case.antigen_chain}, label_asym_id={label_chain}, entity_id={entity_id}")
    return {
        author: by_label[label]
        for author, label in by_author.items()
        if author in contacts and label in by_label
    }


def main() -> int:
    observed: set[int] = set()
    complete = 0
    for case in CASES:
        print(f"\n=== {case.pdb_id}: {case.label}")
        url = CIF_URL.format(case.pdb_id)
        try:
            cif = read_public(url)
            print(f"  mmCIF: {url} SHA-256 {hashlib.sha256(cif).hexdigest()}")
            atoms = protein_atoms(cif)
        except (OSError, ValueError, KeyError, InvalidFileError) as exc:
            print(f"  координаты недоступны: {exc}")
            continue

        missing_chains = {case.antigen_chain, *case.antibody_chains} - set(atoms.chain_id)
        if missing_chains:
            print(f"  отсутствуют авторские цепи: {', '.join(sorted(missing_chains))}")
            continue
        contacts = contact_residues(atoms, case.antigen_chain, case.antibody_chains)
        print(
            f"  авторские контакты цепи {case.antigen_chain} ({len(contacts)}, "
            f"порог {CONTACT_CUTOFF} Å): {author_list(contacts) or 'нет'}"
        )
        try:
            mapped = mapped_contacts(atoms, case, contacts)
        except (OSError, KeyError, TypeError, ValueError) as exc:
            print(f"  сопоставление недоступно: {exc}")
            mapped = {}
        print(f"  покрытие контактов {ACCESSION}: {len(mapped)}/{len(contacts)}")
        unmapped = contacts - mapped.keys()
        if unmapped:
            print(f"  авторские контакты без сопоставления: {author_list(unmapped)}")
        elif mapped:
            complete += 1

        canonical = sorted(set(mapped.values()))
        observed.update(canonical)
        groups: dict[str, list[int]] = {}
        for pos in canonical:
            groups.setdefault(domain_of(pos), []).append(pos)
        for domain, positions in groups.items():
            print(f"  {ACCESSION}, домен {domain}: {ranges(positions)}")
        near_glycan = [g for g in GLYCOSYLATION_SITES if any(abs(g - r) <= 8 for r in canonical)]
        if near_glycan:
            print(f"  сайты гликозилирования {ACCESSION} в пределах ±8 остатков по последовательности: {near_glycan}")

    print(f"\n=== Сводка наблюдаемых контактов {ACCESSION}")
    print(f"Полное сопоставление контактов: {complete}/{len(CASES)} случаев")
    print(f"Уникальных сопоставленных остатков: {len(observed)}")
    for name, lo, hi in ECD_DOMAINS:
        positions = [pos for pos in observed if lo <= pos <= hi]
        pct = 100 * len(positions) / (hi - lo + 1)
        print(f"  домен {name} ({lo}-{hi}): {len(positions)} наблюдаемых контактов ({pct:.1f}% позиций); {ranges(positions)}")
    print("Доступность поверхности требует отдельной оценки.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
