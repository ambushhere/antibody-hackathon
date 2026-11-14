import json

import gemmi

from scripts.prepare_her2_template_controls import (
    CONFIG, MASK, PANEL, POSITIONS, QUERY_FILE, SAMPLES, SEEDS, TEMPLATES,
    frozen_context, rows, sha256,
)

DISULFIDES = [[511, 520], [515, 528], [531, 540], [544, 560],
             [563, 576], [567, 584], [587, 596], [600, 623]]


def test_frozen_queries_hashes_sequences_and_pose_budget():
    _, panel, expected, sources = frozen_context()
    manifest = json.loads((CONFIG / "manifest.json").read_bytes())
    assert all(manifest["sources"][path] == digest for path, digest in sources.items())
    assert manifest["panel"] == panel
    assert json.loads((CONFIG / QUERY_FILE).read_bytes()) == expected
    assert set(manifest["artifact_sha256"]) == {QUERY_FILE, *TEMPLATES.values()}
    for path, digest in manifest["artifact_sha256"].items():
        assert sha256((CONFIG / path).read_bytes()) == digest
    assert manifest["reference"]["common_canonical_positions"] == [p for p in POSITIONS if p not in MASK]
    assert manifest["run"] == dict(model="protenix-v2", msa=False, seeds=list(SEEDS),
                                  samples=list(SAMPLES), expected_poses=48)
    assert [q["query_id"] for q in expected["queries"]] == list(PANEL)
    assert all(c["msa_path"] is None for q in expected["queries"] for c in q["chains"])


def test_template_mask_only_changes_atom_rows():
    manifest = json.loads((CONFIG / "manifest.json").read_bytes())
    blocks = {arm: gemmi.cif.read(str(CONFIG / name)).sole_block() for arm, name in TEMPLATES.items()}
    complete, masked = blocks.values()
    categories = complete.get_mmcif_category_names()
    assert categories == masked.get_mmcif_category_names()
    for category in categories:
        if category != "_atom_site.":
            assert complete.get_mmcif_category(category) == masked.get_mmcif_category(category)
    original = rows(complete, "_atom_site.")
    assert rows(masked, "_atom_site.") == [a for a in original if int(a["auth_seq_id"]) not in MASK]
    for arm, block in blocks.items():
        atoms = rows(block, "_atom_site.")
        assert len(atoms) == manifest["arms"][arm]["atom_count"] == (912 if arm == "complete" else 836)
        assert set((a["group_PDB"], a["label_asym_id"], a["auth_asym_id"], a["label_entity_id"])
                   for a in atoms) == {("ATOM", "A", "A", "1")}
        positions = [int(a["auth_seq_id"]) for a in atoms if a["label_atom_id"] == "CA"]
        assert positions == [p for p in POSITIONS if arm == "complete" or p not in MASK]
        assert all(int(a["label_seq_id"]) + 510 == int(a["auth_seq_id"]) for a in atoms)
        assert rows(block, "_struct_asym.") == [{"id": "A", "entity_id": "1"}]
        scheme = rows(block, "_pdbx_poly_seq_scheme.")
        assert [int(r["seq_id"]) for r in scheme] == list(range(1, 120))
        assert [int(r["auth_seq_num"]) for r in scheme] == list(POSITIONS)
        records = rows(block, "_struct_conn.")
        assert all(r["conn_type_id"] == "disulf" for r in records)
        assert [[int(r[f"ptnr{i}_auth_seq_id"]) for i in (1, 2)] for r in records] == DISULFIDES
        assert all(int(r[f"ptnr{i}_label_seq_id"]) + 510 == int(r[f"ptnr{i}_auth_seq_id"])
                   for r in records for i in (1, 2))
        assert manifest["target"]["disulfides_crossing_crop_boundary_p04626"] == [[626, 634]]
