import json

import gemmi
import pytest

from scripts import audit_template_controls as audit
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


def make_receipt(campaign, arm, seed=None):
    """Build the actual Foundry receipt layout with placeholder pose payloads."""
    arm_folder = campaign / arm
    folder = arm_folder / f"seed-{seed}" if seed is not None else arm_folder
    seeds = [seed] if seed is not None else list(SEEDS)
    config = {
        "missing_msas": "fail", "num_diffusion_samples": 2, "query_failures": "fail",
        "seeds": seeds, "template_paths": [f"/staged/{TEMPLATES[arm]}"],
        "weights": "/app/weights/protenix-v2.pt", "write_detailed_scores": True,
    }
    predictions = []
    for run_seed in seeds:
        for i, query in enumerate(PANEL):
            for sample in SAMPLES:
                paths = audit.artifact_paths(query, run_seed, sample)
                predictions.append(dict(query_id=query, query_index=i, seed=run_seed, sample=sample, **paths))
                for path in paths.values():
                    file = folder / path
                    file.parent.mkdir(parents=True, exist_ok=True)
                    file.touch()
    run = {
        "task": "predict", "schema_version": audit.SCHEMA,
        "provenance": {
            "cli_contract_version": audit.SCHEMA, "module": {"name": "apheris-protenix", "version": "1.3.0"},
            "inputs": [{"path": f"/staged/{QUERY_FILE}", "tier": "none"}], "resolved_config": config,
        },
        "result": {"predictions": predictions},
    }
    (folder / "result.json").write_text(json.dumps(run))
    (folder / "config.json").write_text(json.dumps({"schema_version": audit.SCHEMA, **config}))
    count, atoms = ("119", "912") if arm == "complete" else ("109", "836")
    folder.with_suffix(".log").write_text(
        "Assigned templates to 4 entities across 4 queries\n"
        "Read 1 chain with a total of 119 residues from `_pdbx_poly_seq_scheme` table.\n"
        f"Read {atoms} atoms from `_atom_site` table.\n"
        f"Identified {count} resolved residue(s) from `_atom_site` table.\n"
        f"Seed set to {seeds[0]}\n"
        "Loaded constraint feature: #atom contact:0 #contact:0 #pocket:0\n"
    )
    return folder, run


def test_actual_layout_all_48_samples(tmp_path):
    for arm in TEMPLATES:
        for seed in SEEDS:
            make_receipt(tmp_path, arm, seed)
    poses, runs = audit.collect_runs(tmp_path)
    assert set(poses) == set(audit.EXPECTED)
    assert len(runs) == 6
    assert sum(run["pose_count"] for run in runs.values()) == 48


def test_combined_msa_layout_covers_only_complete_mask(tmp_path):
    make_receipt(tmp_path, "complete")
    poses, runs = audit.collect_runs(tmp_path, arms=("complete",))
    assert set(poses) == {key for key in audit.EXPECTED if key[0] == "complete"}
    assert runs["complete"]["pose_count"] == 24
    assert runs["complete"]["lightning_initialization_seed"] == 42
    make_receipt(tmp_path, "masked", 42)
    with pytest.raises(ValueError, match="Unexpected Protenix run receipt location"):
        audit.collect_runs(tmp_path, arms=("complete",))


@pytest.mark.parametrize("bad", ["missing", "duplicate", "incorrect", "traversal", "unreceipted", "mask"])
def test_incorrect_coverage_or_provenance_rejected(tmp_path, bad):
    folder, run = make_receipt(tmp_path, "complete", 42)
    predictions = run["result"]["predictions"]
    if bad == "missing":
        predictions.pop()
    elif bad == "duplicate":
        predictions.append(predictions[0])
    elif bad == "incorrect":
        predictions[0]["sample"] = 3
    elif bad == "traversal":
        predictions[0]["path"] = "../other_model.cif"
    elif bad == "unreceipted":
        (folder / "other_model.cif").touch()
    else:
        log = folder.with_suffix(".log")
        log.write_text(log.read_text().replace("119 resolved", "109 resolved"))
    (folder / "result.json").write_text(json.dumps(run))
    with pytest.raises(ValueError, match="prediction receipt|unexpected path|artifacts|template, seed"):
        audit.collect_runs(tmp_path)


def test_seed42_is_incomplete_and_keeps_all_pose_rows(tmp_path, monkeypatch):
    for arm in TEMPLATES:
        make_receipt(tmp_path, arm, 42)
    fake_input = ({"data_origin": "pinned"}, None, None, {"control_wt": None}, None, None, None, None, None)
    monkeypatch.setattr(audit, "validate_inputs", lambda *_: fake_input)
    monkeypatch.setattr(audit, "read_reference",
                        lambda *_: (None, dict.fromkeys(audit.COMMON), None, None,
                                    list(range(19)), audit.COMMON, "reference-hash"))
    monkeypatch.setattr(audit, "ecd_reference", lambda *_: (None, None, []))
    monkeypatch.setattr(audit, "score_pose",
                        lambda key, *_: make_score(*key, pair_mean=0.5, interface_rmsd=10))
    report = audit.audit(tmp_path, None, None)
    assert report["status"] == "incomplete"
    assert report["audited_pose_count"] == len(report["predictions"]) == 16
    assert report["expected_pose_count"] == 48
    assert len(report["missing"]) == 32
    audit.write_audit(report, tmp_path / "report")
    assert len((tmp_path / "report/scores.csv").read_text().splitlines()) == 17


def make_score(arm, query, seed, sample, pair_mean, interface_rmsd, max_pair=0.99):
    return {
        "arm": arm, "query_id": query, "seed": seed, "sample": sample,
        "mean_antibody_target_pair_iptm": pair_mean, "max_antibody_target_pair_iptm": max_pair,
        "target_fit_ca_rmsd_A": 2, "ref_interface_ca_rmsd_in_target_frame_A": interface_rmsd,
        "wt_contact_count": 13, "wt_footprint_recall": 13 / 19, "wt_footprint_jaccard": 0.3,
        "heavy_atom_clash_pairs_below_1p5A": 0,
        "ecd_6oge_outside_crop_heavy_atom_clash_pairs_below_1p5A": 0,
        "operational_conditions": {}, "operational_pose_pass": False,
    }


def test_top1_within_each_seed_is_by_mean_pair_iptm_and_lower_sample_ties():
    prefix = ("complete", "control_wt")
    samples = [
        make_score(*prefix, 42, 1, 0.60, 8, max_pair=0.60),
        make_score(*prefix, 42, 2, 0.55, 2, max_pair=0.95),
        make_score(*prefix, 43, 2, 0.75, 1),
        make_score(*prefix, 43, 1, 0.75, 7),
    ]
    seed42, seed43 = audit.rank_samples(samples)
    for selection in [seed42, seed43]:
        assert selection["top1"]["sample"] == 1
        assert selection["best_reference_geometry_diagnostic"]["sample"] == 2
        assert selection["best_reference_geometry_diagnostic"]["confidence_rank_in_seed"] == 2


def test_operational_rule_uses_two_passing_selected_wt_seeds_and_all_four_conditions():
    prefix = ("complete", "control_wt")
    samples = [make_score(*prefix, seed, 1, 0.5, 4, max_pair=0.6) for seed in SEEDS]
    samples[-1]["heavy_atom_clash_pairs_below_1p5A"] = 1
    for row in samples:
        row["operational_conditions"] = audit.operational_conditions(row)
        row["operational_pose_pass"] = all(row["operational_conditions"].values())
    assert all(samples[0]["operational_conditions"].values())
    assert sum(samples[-1]["operational_conditions"].values()) == 3
    groups = audit.rank_samples(samples)
    assert audit.operational_rule(groups, complete=False)["arms"]["complete"]["qualifies"] is None
    complete = audit.operational_rule(groups, complete=True)["arms"]["complete"]
    assert complete == {"wt_top1_seeds_observed": [42, 43, 44],
                        "passing_wt_top1_seeds": [42, 43], "qualifies": True}


def test_target_msa_can_only_add_alignment_path_to_target(tmp_path):
    _, _, expected, _ = frozen_context()
    target = expected["queries"][0]["chains"][0]["sequence"]
    alignment_path = "/staged/target.a3m"
    queries = {"queries": [
        {**q, "chains": [{**c, "msa_path": alignment_path if c["chain_id"] == "T" else None}
                         for c in q["chains"]]} for q in expected["queries"]
    ]}
    sequence_rows = [target] + [target[:i] + "-" + target[i + 1:] for i in range(17)]
    sequence_rows += [target] * 25
    a3m = tmp_path / "target.a3m"
    a3m.write_text("".join(f">{i}\n{seq}\n" for i, seq in enumerate(sequence_rows)))
    query_file = tmp_path / QUERY_FILE
    query_file.write_text(json.dumps(queries))
    manifest = {
        "template_sha256": "same-template",
        "a3m_sha256": audit.file_hash(a3m), "queries_sha256": audit.file_hash(query_file),
        "total_rows_including_query": 43, "unique_aligned_sequences": 18, "query_length": 119,
        "query_sha256": sha256(target.encode()),
    }
    manifest_file = tmp_path / "manifest.json"
    manifest_file.write_text(json.dumps(manifest))
    staged, path = audit.validate_target_msa(tmp_path, expected, "same-template")
    assert path == alignment_path
    assert set(staged) == {"complete-target.cif", "target.a3m", "queries.json", "manifest.json"}
    queries["queries"][1]["chains"][1]["msa_path"] = alignment_path
    query_file.write_text(json.dumps(queries))
    manifest["queries_sha256"] = audit.file_hash(query_file)
    manifest_file.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="change more than the target alignment"):
        audit.validate_target_msa(tmp_path, expected, "same-template")
