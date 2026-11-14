"""Offline checks for author, mmCIF label and canonical HER2 numbering."""

import json

import pytest

from src import epitope_map as epitope


ATOM_SITE = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""


def cif_with_atoms(*rows):
    return (ATOM_SITE + "\n".join(rows) + "\n").encode()


def sifts_record(pdb_id, label_chain, author_chain, entity, start=1, reference=23, length=607):
    return {
        "rcsb_id": f"{pdb_id}_{entity}",
        "rcsb_polymer_entity_container_identifiers": {
            "entity_id": entity,
            "asym_ids": [label_chain],
            "auth_asym_ids": [author_chain],
        },
        "rcsb_polymer_entity_align": [{
            "provenance_source": "SIFTS",
            "reference_database_name": "UniProt",
            "reference_database_accession": "P04626",
            "aligned_regions": [{
                "entity_beg_seq_id": start,
                "ref_beg_seq_id": reference,
                "length": length,
            }],
        }],
    }


@pytest.mark.parametrize(
    "pdb_id,chain,entity,author_num,length,expected",
    [("1N8Z", "C", "3", 557, 607, 579), ("6OGE", "A", "1", 579, 622, 579)],
)
def test_structure_specific_numbering(monkeypatch, pdb_id, chain, entity, author_num, length, expected):
    # Both entity sequences begin at P04626 23, but their author numbering differs.
    case = epitope.Case(pdb_id, chain, ("Z",), "trastuzumab")
    atoms = epitope.protein_atoms(cif_with_atoms(
        f"ATOM 1 C CA . ALA {chain} {entity} 557 ? 0 0 0 {author_num} ALA {chain} CA 1",
    ))
    record = sifts_record(pdb_id, chain, chain, entity, length=length)

    def read_public(url):
        assert url == epitope.ENTITY_URL.format(pdb_id, entity)
        return json.dumps(record).encode()

    monkeypatch.setattr(epitope, "read_public", read_public)
    assert epitope.mapped_contacts(atoms, case, {(author_num, "")}) == {(author_num, ""): expected}
    mapping = epitope.sifts_mapping(record, case, chain, entity)
    assert mapping[1] == 23
    assert mapping[length] == 22 + length
    assert length + 1 not in mapping


def test_insertion_codes_and_heavy_atoms(monkeypatch):
    atoms = epitope.protein_atoms(cif_with_atoms(
        "ATOM 1 C CA . ALA X 3 2 ? 0 0 0 40 ALA C CA 1",
        "ATOM 2 C CA . GLY X 3 3 A 1 0 0 40 GLY C CA 1",
        "ATOM 3 H H  . ALA X 3 4 ? 1 0 0 41 ALA C H  1",
        "ATOM 4 C CA . ALA X 3 4 ? 9 0 0 41 ALA C CA 1",
        "ATOM 5 C CA . ALA Y 1 1 ? 0 1 0 1 ALA A CA 1",
        "ATOM 6 D D  . ALA Y 1 1 ? 9 1 0 1 ALA A D  1",
        "ATOM 7 C CA . ALA X 3 5 ? 0 0 0 42 ALA C CA 2",
    ))
    case = epitope.Case("TEST", "C", ("A",), "test")
    record = sifts_record("TEST", "X", "C", "3", start=2, reference=580, length=3)
    monkeypatch.setattr(epitope, "read_public", lambda url: json.dumps(record).encode())

    contacts = epitope.contact_residues(atoms, case.antigen_chain, case.antibody_chains)
    assert contacts == {(40, ""), (40, "A")}
    assert epitope.mapped_contacts(atoms, case, contacts) == {(40, ""): 580, (40, "A"): 581}


@pytest.mark.parametrize("sifts_available", [True, False])
def test_missing_or_ambiguous_mapping_stays_out_of_aggregate(monkeypatch, capsys, sifts_available):
    case = epitope.Case("TEST", "C", ("A",), "test")
    cif = cif_with_atoms(
        "ATOM 1 C CA . ALA X 3 2 ? 0 0 0 550 ALA C CA 1",
        "ATOM 2 C CA . ALA X 3 3 A 1 0 0 550 ALA C CA 1",
        "ATOM 3 C CA . ALA X 3 4 A 2 0 0 550 ALA C CA 1",
        "ATOM 4 C CA . ALA X 3 . ? 0 0 0 570 ALA C CA 1",
        "ATOM 5 C CA . ALA X 3 5 ? 0 0 0 620 ALA C CA 1",
        "ATOM 6 C CA . ALA Y 1 1 ? 0 1 0 1 ALA A CA 1",
        "ATOM 7 C CA . ALA X 3 7 ? 0 0 0 ? ALA C CA 1",
    )
    record = sifts_record("TEST", "X", "C", "3", start=2, reference=580, length=3)

    def read_public(url):
        if url == epitope.CIF_URL.format(case.pdb_id):
            return cif
        assert url == epitope.ENTITY_URL.format(case.pdb_id, "3")
        if not sifts_available:
            raise OSError("source unavailable")
        return json.dumps(record).encode()

    monkeypatch.setattr(epitope, "read_public", read_public)
    monkeypatch.setattr(epitope, "CASES", [case])
    assert epitope.main() == 0
    output = capsys.readouterr().out
    assert "550, 550A, 570, 620, ?" in output
    assert "Полное сопоставление контактов: 0/1" in output
    summary = output.split("=== Сводка наблюдаемых контактов P04626")[1]
    if sifts_available:
        assert "покрытие контактов P04626: 1/5" in output
        assert "авторские контакты без сопоставления: 550A, 570, 620, ?" in output
        assert "Уникальных сопоставленных остатков: 1" in summary
        assert "; 580" in summary
    else:
        assert "покрытие контактов P04626: 0/5" in output
        assert "авторские контакты без сопоставления: 550, 550A, 570, 620, ?" in output
        assert "Уникальных сопоставленных остатков: 0" in summary
    assert all(str(pos) not in summary for pos in (550, 570, 620))


def test_ambiguous_sifts_regions_are_rejected():
    case = epitope.Case("TEST", "C", ("A",), "test")
    record = sifts_record("TEST", "X", "C", "3", length=2)
    record["rcsb_polymer_entity_align"][0]["aligned_regions"].append({
        "entity_beg_seq_id": 2,
        "ref_beg_seq_id": 30,
        "length": 2,
    })
    with pytest.raises(ValueError, match="Перекрывающиеся диапазоны SIFTS"):
        epitope.sifts_mapping(record, case, "X", "3")
