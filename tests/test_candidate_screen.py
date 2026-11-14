import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from biotite.structure import AtomArray

from scripts import audit_candidate_screen as audit
from scripts import prepare_candidate_screen as prepare


def dump(path, value):
    path.write_text(json.dumps(value))


def test_frozen_panel_exact_coverage_without_ranking():
    manifest, canonical, queries, hashes = audit.validate_frozen(prepare.CONFIG)
    panel = manifest["panel"]
    assert manifest["panel_counts"] == {"control": 4, "optimization": 10, "denovo": 43}
    assert list(queries)[:4] == list(prepare.CONTROL_IDS)
    assert len({p["paired_sequence_sha256"] for p in panel}) == len(queries) == 57
    assert all(tuple(chains) == ("T", "H", "L") and chains["T"] == canonical[510:629]
               and 0 < len(chains["H"]) + len(chains["L"]) <= 250 for chains in queries.values())
    assert all(c["msa_path"] is None
               for q in audit.read_json(prepare.CONFIG / prepare.QUERY_FILE)["queries"] for c in q["chains"])
    with (prepare.ROOT / prepare.DENOVO).open(newline="") as stream:
        expected = {r["design_id"]: r for r in csv.DictReader(stream)
                    if r["local_only_eligibility"] == "eligible_local_only"}
    de_novo = {p["generation"]["design_id"]: p for p in panel if p["track"] == "denovo"}
    assert de_novo.keys() == expected.keys()
    for design_id, row in de_novo.items():
        source = expected[design_id]
        generation = row["generation"]
        assert row["query_id"] == f"{source['arm']}_{generation['index']}"
        assert generation["structure_sha256"] == source["structure_sha256"]
        assert row["official_adaptyv_novelty"] == "unverified_database_unavailable"
    with (prepare.ROOT / prepare.OPTIMIZATION).open(newline="") as stream:
        assert [p["query_id"] for p in panel if p["track"] == "optimization"] == [
            f"opt_{row['cdrh3']}" for row in csv.DictReader(stream)]
    baseline = audit.read_json(prepare.CONTROLS / "manifest.json")
    assert hashes[prepare.TEMPLATE] == baseline["artifact_sha256"][prepare.TEMPLATE]


def test_prepare_refuses_to_replace_a_changed_frozen_input(tmp_path):
    prepare.prepare(tmp_path)
    (tmp_path / "queries.json").write_text("{}")
    with pytest.raises(ValueError, match="Refusing to replace frozen input"):
        prepare.prepare(tmp_path)
    assert (tmp_path / "queries.json").read_text() == "{}"


def receipt(campaign, ids, prediction_ids, *, config_dir=prepare.CONFIG, staged_dir=Path("/staged"), target_msa=False):
    manifest, _, _, hashes = audit.validate_frozen(config_dir, target_msa=target_msa)
    campaign.mkdir()
    staged_paths = {key: str(staged_dir / key) for key in hashes}
    (campaign / "input.sha256").write_text("".join(f"{digest}  {staged_paths[k]}\n" for k, digest in hashes.items()))
    (campaign / "image-id.txt").write_text(prepare.IMAGE_ID + "\n")
    config = {
        "missing_msas": "fail", "num_diffusion_samples": 1, "query_failures": "fail",
        "seeds": [42], "template_paths": [staged_paths[prepare.TEMPLATE]],
        "weights": "/app/weights/protenix-v2.pt", "write_detailed_scores": True,
    }
    predictions = []
    for q in prediction_ids:
        paths = audit.artifact_paths(q, 42, 1)
        predictions.append(dict(query_id=q, query_index=ids.index(q), seed=42, sample=1, **paths))
        for relative in paths.values():
            path = campaign / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    run = {
        "task": "predict", "schema_version": audit.SCHEMA,
        "provenance": {
            "cli_contract_version": audit.SCHEMA, "module": {"name": "apheris-protenix", "version": "1.3.0"},
            "inputs": [{"path": staged_paths[prepare.QUERY_FILE], "tier": "none"}],
            "resolved_config": config, "started_at": "2026-11-14T14:45:00Z", "ended_at": "2026-11-14T14:50:00Z",
        },
        "result": {"predictions": predictions},
    }
    dump(campaign / "result.json", run)
    dump(campaign / "config.json", {"schema_version": audit.SCHEMA, **config})
    (campaign / "run.log").write_text(
        f"Assigned templates to {len(manifest['panel'])} entities across {len(ids)} queries\n"
        "Read 1 chain with a total of 119 residues from `_pdbx_poly_seq_scheme` table.\n"
        "Read 912 atoms from `_atom_site` table.\n"
        "Identified 119 resolved residue(s) from the template\n"
        "Seed set to 42\nLoaded constraint feature: #atom contact:0 #contact:0 #pocket:0\n"
    )
    return staged_paths, run


@pytest.fixture
def ids():
    return list(audit.validate_frozen(prepare.CONFIG)[2])


def test_receipt_exactly_matches_frozen_keys_and_staged_paths(tmp_path, ids):
    campaign = tmp_path / "campaign"
    paths, _ = receipt(campaign, ids, ids[::-1])
    _, _, _, hashes = audit.validate_frozen(prepare.CONFIG)
    staged, staged_paths, image = audit.validate_staged(campaign, hashes)
    files, run = audit.validate_run(campaign, ids, staged_paths)
    assert set(files) == set(ids)
    assert run["receipted_pose_count"] == len(ids) == 57
    assert staged == hashes and staged_paths == paths and image == prepare.IMAGE_ID


@pytest.mark.parametrize("bad", ["duplicate", "extra", "traversal", "sample", "index",
                                  "raw", "another_receipt", "settings", "template_log"])
def test_ambiguous_attribution_rejected(tmp_path, ids, bad):
    campaign = tmp_path / "campaign"
    paths, run = receipt(campaign, ids, ids[:1])
    rows = run["result"]["predictions"]
    if bad == "duplicate":
        rows.append(rows[0].copy())
    elif bad == "extra":
        rows[0]["query_id"] = "unknown_query"
    elif bad == "traversal":
        rows[0]["path"] = "../example.cif"
    elif bad == "sample":
        rows[0]["sample"] = 2
    elif bad == "index":
        rows[0]["query_index"] = 7
    elif bad == "raw":
        (campaign / "unreceipted_model.cif").touch()
    elif bad == "another_receipt":
        (campaign / ids[0] / "result.json").touch()
    elif bad == "settings":
        run["provenance"]["resolved_config"]["weights"] = "/app/weights/protenix-v1.pt"
    else:
        path = campaign / "run.log"
        path.write_text(path.read_text().replace("119 resolved", "109 resolved"))
    dump(campaign / "result.json", run)
    with pytest.raises(ValueError):
        audit.validate_run(campaign, ids, paths)


def test_missing_and_invalid_poses_remain_unknown_and_zero_is_observed(tmp_path, ids, monkeypatch):
    campaign = tmp_path / "campaign"
    receipt(campaign, ids, ids[:3])
    (campaign / audit.artifact_paths(ids[0], 42, 1)["scores_path"]).unlink()
    monkeypatch.setattr(audit, "validate_file", lambda _, digest: digest)
    monkeypatch.setattr(audit, "read_1n8z", lambda *_: (
        None, dict.fromkeys(audit.COMMON), None, {}, list(range(19)), audit.COMMON, "1n8z",
    ))
    monkeypatch.setattr(audit.gemmi.cif, "read", lambda *_: SimpleNamespace(sole_block=lambda: None))
    monkeypatch.setattr(audit, "ecd_reference", lambda *_: (None, None, []))

    def score(files, expected, entry, *_):
        if entry["query_id"] == ids[1]:
            raise ValueError("Wrong submitted sequence")
        return {audit.CROP_CLASHES: 0, "wt_contact_count": 0}

    monkeypatch.setattr(audit, "score_pose", score)
    report = audit.audit(campaign, Path("1n8z"), Path("6oge"), control_audit=tmp_path / "absent")
    assert report["status"] == "incomplete"
    assert report["audited_pose_count"] == 1
    assert report["unknown_pose_count"] == 56
    missing, invalid, known, unreceipted = report["predictions"][:4]
    assert "Missing receipted artifacts" in missing["reason"]
    assert missing["row_artifact_sha256"]["scores_path"] is None
    assert invalid["reason"] == "Wrong submitted sequence"
    assert "no prediction receipt" in unreceipted["reason"]
    assert [r["wt_contact_count"] for r in (missing, invalid, known, unreceipted)] == [None, None, 0, None]
    assert all(r["descriptive_pareto_front"] is None for r in report["predictions"])
    audit.write_audit(report, tmp_path / "output")
    with (tmp_path / "output/scores.csv").open(newline="") as stream:
        assert len(list(csv.DictReader(stream))) == 57


def test_directional_predicted_interface_pae_and_no_contact_null():
    mapped = AtomArray(3)
    mapped.chain_id = ["T", "H", "L"]
    mapped.atom_name = ["CA"] * 3
    mapped.res_id = [579, 3, 2]
    pae = np.arange(9).reshape(3, 3)
    result = audit.predicted_interface_pae(pae, mapped, {
        "target_contacts_p04626": [579], "antibody_contact_residues": ["H:3", "L:2"],
    })
    assert result == {
        "predicted_interface_pae_ab_rows_target_cols_A": 4.5,
        "predicted_interface_pae_target_rows_ab_cols_A": 1.5,
    }
    empty = audit.predicted_interface_pae(pae, mapped, {
        "target_contacts_p04626": [], "antibody_contact_residues": [],
    })
    assert list(empty.values()) == [None, None]


def test_descriptive_front_uses_track_specific_supported_axes():
    base = {"track": "denovo", "status": "audited", audit.MEAN_IPTM: 0.3,
            "wt_contact_count": 13, "target_fit_ca_rmsd_A": 1.0,
            audit.CROP_CLASHES: 0, audit.ECD_CLASHES: 0, audit.INTERFACE_RMSD: None,
            "descriptive_pareto_front": None}
    rows = [
        {**base, "query_id": "best"}, {**base, "query_id": "tie"},
        {**base, "query_id": "dominated", audit.CROP_CLASHES: 1},
        {**base, "query_id": "tradeoff", audit.MEAN_IPTM: 0.4, "wt_contact_count": 12},
        {**base, "query_id": "unknown", audit.CROP_CLASHES: None},
        {**base, "query_id": "needs_mapping", "track": "optimization"},
    ]
    fronts = audit.descriptive_pareto(rows, supported=True)
    assert fronts["denovo"]["front"] == ["best", "tie", "tradeoff"]
    assert fronts["optimization"]["front"] is None
    assert rows[2]["descriptive_pareto_front"] is False
    assert all(row["descriptive_pareto_front"] is None for row in rows[-2:])


def test_failed_complete_arm_retains_unrun_poses_and_uses_top_1(tmp_path):
    _, _, _, hashes = audit.validate_frozen(prepare.CONFIG)
    predictions = [
        {"arm": "complete", "query_id": q, "seed": seed, "sample": sample,
         audit.MEAN_IPTM: 0.4 if sample == 1 else 0.3,
         "wt_contact_count": 14 if seed == 42 or sample == 2 else 12,
         "target_fit_ca_rmsd_A": 1.0, audit.INTERFACE_RMSD: 3.0, audit.CROP_CLASHES: 0}
        for q in prepare.CONTROL_IDS for seed in audit.SEEDS for sample in audit.SAMPLES
    ]
    baseline = {
        "status": "incomplete", "msa_policy": "none", "protenix_image_id": prepare.IMAGE_ID,
        "inputs_sha256": {prepare.TEMPLATE: hashes[prepare.TEMPLATE]},
        "predictions": predictions,
        "operational_rule": {
            "source": prepare.ADMISSION, "status": "incomplete", "minimum_passing_seeds": 2,
            "arms": {"complete": {
                "wt_top1_seeds_observed": [42, 43, 44],
                "passing_wt_top1_seeds": [42], "qualifies": None,
            }},
        },
    }
    control = tmp_path / "baseline.json"
    dump(control, baseline)
    report = audit.audit(tmp_path / "no-campaign", Path("1n8z"), Path("6oge"), control_audit=control)
    assert report["status"] == "not_run_admission_failed"
    assert report["baseline_admission"]["passing_wt_top1_seeds"] == [42]
    assert report["audited_pose_count"] == 0 and report["unknown_pose_count"] == 57
    assert all(r[audit.CROP_CLASHES] is None for r in report["predictions"])
    baseline["msa_policy"] = "target_only"
    dump(control, baseline)
    assert audit.control_admission(control, hashes[prepare.TEMPLATE])["status"] == "unknown"


def test_denovo_does_not_map_to_the_wt_antibody(monkeypatch):
    manifest, _, _, _ = audit.validate_frozen(prepare.CONFIG)
    panel = {row["query_id"]: row for row in manifest["panel"]}
    entry = next(row for row in panel.values() if row["track"] == "denovo")
    files = {name: Path(name) for name in audit.SUFFIXES}
    atoms = AtomArray(3)
    atoms.atom_name = ["CA"] * 3
    atoms.chain_id = list(prepare.ROLES)
    monkeypatch.setattr(audit, "validate_prediction", lambda *_args, **_kw: (
        atoms, dict(zip(prepare.ROLES, prepare.ROLES, strict=True))))
    monkeypatch.setattr(audit, "read_json", lambda *_: {
        "num_chains": 3, "num_atoms": 3, "num_tokens": 3, "num_residues": 3,
        "interface_predicted_tm_score": 0.5, "clash_flag": False,
        "chain_chain_interface_predicted_tm_score": [[1, 0.4, 0.2], [0.4, 1, 0.1], [0.2, 0.1, 1]],
    })
    monkeypatch.setattr(audit, "protenix_pae", lambda *_: np.zeros((3, 3)))
    monkeypatch.setattr(audit, "overlay_clashes", lambda *_: {audit.ECD_CLASHES: 0})

    def geometry(_atoms, _roles, _reference, _mapping, matched):
        assert matched is False
        return atoms, {"target_contacts_p04626": [], "antibody_contact_residues": [],
                       "wt_contact_count": 0, audit.CROP_CLASHES: 0}

    monkeypatch.setattr(audit, "geometry_row", geometry)
    values = audit.score_pose(files, {}, entry, panel, [None, None, None, None, None, []], None)
    assert values[audit.INTERFACE_RMSD] is None
    assert values["ref_interface_pae_ab_rows_target_cols_A"] is None
    assert values["wt_contact_count"] == values[audit.CROP_CLASHES] == values[audit.ECD_CLASHES] == 0
    assert values[audit.MEAN_IPTM] == pytest.approx(0.3)
