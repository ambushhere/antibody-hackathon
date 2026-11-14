import pytest

from scripts import thermal_proxy_pilot as pilot


def pdb_atom(serial, atom, amino_acid, chain, number, insertion=" ", alt=" ", occupancy=1):
    return (
        f"ATOM  {serial:5d} {atom:^4s}{alt}{amino_acid} {chain}{number:4d}{insertion}   "
        f"{serial:8.3f}{0:8.3f}{0:8.3f}{occupancy:6.2f}{10:6.2f}          {atom[0]:>2}"
    )


def test_prepare_fv_preserves_insertions_and_disambiguates_chain_numbering(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot, "DOMAINS", {"A": "AH", "B": "Y"})
    source, destination = (tmp_path / filename for filename in ("in.pdb", "out.pdb"))
    lines = []
    for chain, number, insertion, aa in (("A", 5, " ", "ALA"), ("A", 5, "A", "HIS"), ("B", 5, " ", "TYR")):
        for atom in pilot.BACKBONE:
            if aa == "HIS" and atom == "CA":
                lines.append(pdb_atom(len(lines) + 1, atom, aa, chain, number, insertion, "B", 0.5))
                lines.append(pdb_atom(len(lines) + 1, atom, aa, chain, number, insertion, "A", 0.5))
            else:
                lines.append(pdb_atom(len(lines) + 1, atom, aa, chain, number, insertion))
    source.write_text("\n".join(lines + ["END", ""]))

    mapping, coordinates = pilot.prepare_fv(source, destination)
    assert [(r.source_chain, r.source_resnum, r.source_icode, r.model_resnum) for r in mapping] == [
        ("A", 5, "", 1), ("A", 5, "A", 2), ("B", 5, "", 3),
    ]
    assert mapping[1].alternate_backbone == "CA:A(0.50)"
    assert coordinates["A", 5, "A"] == (7.0, 0.0, 0.0)
    output = pilot.PDBParser(PERMISSIVE=False, QUIET=True).get_structure("fv", destination)
    assert [(r.parent.id, r.id[1:]) for r in output.get_residues()] == [
        ("A", (1, " ")), ("A", (2, " ")), ("B", (3, " ")),
    ]
    assert output[0]["A"][2]["CA"].altloc == " "


def test_prepare_fv_rejects_missing_backbone(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot, "DOMAINS", {"A": "A"})
    source = tmp_path / "in.pdb"
    source.write_text(pdb_atom(1, "CA", "ALA", "A", 1) + "\nEND\n")
    with pytest.raises(ValueError, match="incomplete backbone"):
        pilot.prepare_fv(source, tmp_path / "out.pdb")


def test_lookup_requires_source_chain_wildtype_and_blank_insertion():
    mapping = [
        pilot.Mapping("A", 91, "", "H", 91, ""),
        pilot.Mapping("B", 91, "", "Y", 198, ""),
        pilot.Mapping("A", 92, "A", "R", 92, ""),
    ]
    assert pilot.lookup_mutations(mapping, (("A", 91, "H", "Y"),)) == [(mapping[0], "Y")]
    for mutation in (("B", 91, "H", "Y"), ("A", 92, "R", "G"), ("A", 91, "H", "H")):
        with pytest.raises(ValueError, match="Invalid substitution"):
            pilot.lookup_mutations(mapping, [mutation])
    with pytest.raises(ValueError, match="Repeated position"):
        pilot.lookup_mutations(mapping, [("A", 91, "H", "Y"), ("A", 91, "H", "A")])


def test_check_input_rejects_position_order_even_with_matching_wt():
    mapping = [
        pilot.Mapping("A", i, "", wt, i, "") for i, wt in enumerate(pilot.VL, 1)
    ] + [
        pilot.Mapping("B", i, "", wt, len(pilot.VL) + i, "") for i, wt in enumerate(pilot.VH, 1)
    ]
    pdb = {
        "num_of_chains": 2, "seq": pilot.VL + pilot.VH,
        "seq_chain_A": pilot.VL, "seq_chain_B": pilot.VH,
        "resn_list_A": [str(r.model_resnum) for r in mapping if r.source_chain == "A"],
        "resn_list_B": [str(r.model_resnum) for r in mapping if r.source_chain == "B"],
    }
    pilot.check_input(pdb, mapping)
    pdb["resn_list_B"] = pdb["resn_list_B"][::-1]
    with pytest.raises(ValueError, match="differs for chain B"):
        pilot.check_input(pdb, mapping)
