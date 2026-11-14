"""Local regression checks for target sequence, numbering and matched templates."""

import hashlib
import io
import json
from pathlib import Path

import gemmi
import pytest
from Bio.PDB import MMCIFParser

from scripts import prepare_her2_context as context

CONFIG = Path(__file__).resolve().parents[1] / "configs/her2-context"
MANIFEST = json.loads((CONFIG / "manifest.json").read_text())
CANONICAL = context.canonical_sequence((CONFIG / "P04626.fasta").read_bytes())
EXPECTED = {
    "split": {"C": (511, 602), "D": (613, 629)},
    "domain-iv": {"T": (511, 652)},
    "ecd": {"T": (23, 652)},
}
DEPOSITED_MISSING = set(range(102, 111)) | set(range(303, 306)) | set(range(361, 365)) | set(range(581, 591))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def test_mapping_and_design_hotspots():
    assert context.target_sequence(CANONICAL, (1, 1)) == "T"
    assert context.target_sequence(CANONICAL, (581, 590)) == "GVKPDLSYMP"
    assert context.target_sequence(CANONICAL, (608, 630)) == "CTHSCVDLDDKGCPAEQRASPLT"
    mapping = MANIFEST["sifts"]
    assert mapping["rcsb"] == {"entity_beg_seq_id": 1, "length": 607, "ref_beg_seq_id": 23}
    assert (mapping["pdbe"]["unp_start"], mapping["pdbe"]["unp_end"]) == (23, 629)
    audit = MANIFEST["audit"]
    assert set(audit["rcsb_unobserved_deposited_mature_positions"]) == DEPOSITED_MISSING
    assert audit["rcsb_unobserved_deposited_canonical_positions"] == [r + 22 for r in sorted(DEPOSITED_MISSING)]
    assert audit["footprint_canonical_positions"] == [
        579, 580, 582, 583, 591, 592, 593, 594, 595, 601, 613, 615, 618, 620, 622, 624, 625, 626, 627,
    ]
    assert [(r["amino_acid"], r["canonical_position"], r["domain_iv_query_position"])
            for r in audit["design_hotspots"]] == [("K", 615, 105), ("F", 595, 85), ("I", 613, 103), ("P", 594, 84)]


def test_file_hashes():
    assert sha256(CANONICAL.encode()) == MANIFEST["canonical_sequence_sha256"]
    for filename, checksum in MANIFEST["artifact_sha256"].items():
        assert sha256((CONFIG / filename).read_bytes()) == checksum


@pytest.mark.parametrize("arm", EXPECTED)
def test_polymer_and_atom_indices(arm):
    info = MANIFEST["arms"][arm]
    template = (CONFIG / info["template"]).read_text()
    block = gemmi.cif.read_string(template).sole_block()
    author = MMCIFParser(QUIET=True).get_structure(arm, io.StringIO(template))[0]
    label = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(
        arm, io.StringIO(template)
    )[0]
    expected_chains = EXPECTED[arm]
    asym = {row["id"]: row["entity_id"] for row in context.rows(block, "_struct_asym.")}
    assert set(asym) == {chain.id for chain in author} == {chain.id for chain in label} == set(expected_chains)
    assert len(asym) == len({row["id"] for row in context.rows(block, "_entity.")})
    assert {row["id"] for row in context.rows(block, "_entity.")} == set(asym.values())
    components = {row["id"]: row["type"] for row in context.rows(block, "_chem_comp.")}
    poly_seq = context.rows(block, "_entity_poly_seq.")
    for chain, (first, last) in expected_chains.items():
        sequence = CANONICAL[first - 1:last]
        mature = set(range(first - 22, last - 21))
        deposited_gaps = mature & DEPOSITED_MISSING
        tail = {position for position in mature if position > 607}
        resolved = mature - deposited_gaps - tail
        metadata = info["target_chains"][chain]
        assert metadata["query_position_to_canonical"] == list(range(first, last + 1))
        assert metadata["length"] == len(sequence)
        assert metadata["sequence_sha256"] == sha256(sequence.encode())
        assert set(metadata["unresolved_deposited_mature"]) == deposited_gaps
        assert set(metadata["beyond_deposit_mature"]) == tail
        assert metadata["cysteine_canonical_positions"] == [
            i for i in range(first, last + 1) if CANONICAL[i - 1] == "C"
        ]
        polymer, = [row for row in context.rows(block, "_entity_poly.") if row["entity_id"] == asym[chain]]
        assert polymer["pdbx_strand_id"] == chain
        assert polymer["pdbx_seq_one_letter_code_can"] == sequence
        assert polymer["type"] == "polypeptide(L)"
        entity_poly_seq = [row for row in poly_seq if row["entity_id"] == asym[chain]]
        assert [int(row["num"]) for row in entity_poly_seq] == list(range(1, len(sequence) + 1))
        assert all("peptide" in components[row["mon_id"]] for row in entity_poly_seq)
        assert "".join(context.residue_aa(row["mon_id"]) for row in entity_poly_seq) == sequence
        scheme = [row for row in context.rows(block, "_pdbx_poly_seq_scheme.") if row["asym_id"] == chain]
        assert [int(row["pdb_seq_num"]) for row in scheme] == sorted(mature)
        assert [int(row["seq_id"]) for row in scheme] == list(range(1, len(sequence) + 1))
        assert {int(row["pdb_seq_num"]) for row in scheme if row["auth_seq_num"] is None} == deposited_gaps | tail
        unobserved = [row for row in context.rows(block, "_pdbx_unobs_or_zero_occ_residues.")
                      if row["auth_asym_id"] == chain]
        assert {int(row["auth_seq_id"]) for row in unobserved} == deposited_gaps
        assert {int(row["label_seq_id"]) for row in unobserved} == {r + 23 - first for r in deposited_gaps}
        assert {residue.id[1] for residue in author[chain]} == resolved
        assert {residue.id[1] for residue in label[chain]} == {r + 23 - first for r in resolved}
        atoms = context.protein_atoms(block, chain)
        assert context.residue_numbers(atoms) == resolved
        assert all(atom["label_entity_id"] == asym[chain] and atom["label_atom_id"] != "OXT"
                   and int(atom["label_seq_id"]) == int(atom["auth_seq_id"]) + 23 - first
                   and atom["label_asym_id"] == chain for atom in atoms)
    assert len(context.rows(block, "_atom_site.")) == sum(len(context.protein_atoms(block, c)) for c in asym)


def test_coordinate_identity_and_cross_gap_disulfide():
    split, domain, ecd = (gemmi.cif.read(str(CONFIG / f"{name}-target.cif")).sole_block()
                          for name in EXPECTED)
    split_atoms = context.atom_coordinates(context.protein_atoms(split, "C") + context.protein_atoms(split, "D"))
    domain_atoms = context.atom_coordinates(context.protein_atoms(domain, "T"))
    ecd_atoms = context.atom_coordinates(context.protein_atoms(ecd, "T"))
    assert len(split_atoms) == len(domain_atoms) == 836
    assert len(ecd_atoms) == 4490
    assert split_atoms == domain_atoms == {key: ecd_atoms[key] for key in domain_atoms}

    def disulfides(block):
        return {(int(row["ptnr1_auth_seq_id"]), int(row["ptnr2_auth_seq_id"]))
                for row in context.rows(block, "_struct_conn.")}
    assert (578, 601) in disulfides(domain) & disulfides(ecd)
    assert (578, 601) not in disulfides(split)


def test_matched_antibody_pairs_and_controls():
    expected_panel = {
        "control_wt", "control_facs_low", "control_spr_strong", "control_spr_weak",
        "denovo_fullseq_k1_1", "opt_WTAVGFFAMD", "opt_YTGDGFPAMN",
    }
    panel = {p["query_id"]: p for p in MANIFEST["panel"]}
    assert set(panel) == expected_panel
    assert len(panel["control_wt"]["all_wt_kd_nM"]) == 8
    assert panel["control_facs_low"]["cdrh3"] == "WSMFGHYMFA"
    assert "kd_nM" not in panel["control_facs_low"]
    assert panel["control_spr_strong"]["kd_nM"] == 0.02343889
    assert panel["control_spr_weak"]["kd_nM"] == 13863.16
    assert panel["denovo_fullseq_k1_1"]["historical_footprint_recall"] == pytest.approx(13 / 19)
    for arm, expected_chains in EXPECTED.items():
        queries = json.loads((CONFIG / MANIFEST["arms"][arm]["queries"]).read_text())["queries"]
        assert [query["query_id"] for query in queries] == list(panel)
        for query in queries:
            entry = panel[query["query_id"]]
            assert set(query) == {"query_id", "chains"}
            assert [chain["chain_id"] for chain in query["chains"]] == [*expected_chains, "H", "L"]
            assert all(set(chain) == {"chain_id", "sequence", "msa_path"} and chain["msa_path"] is None
                       for chain in query["chains"])
            sequences = context.chain_sequences(query)
            for chain, (start, end) in expected_chains.items():
                assert sequences[chain] == CANONICAL[start - 1:end]
            assert sequences["H"] == entry["vh_sequence"]
            assert sequences["L"] == entry["vl_sequence"]
            assert sha256(sequences["H"].encode()) == entry["vh_sha256"]
            assert sha256(sequences["L"].encode()) == entry["vl_sha256"]
