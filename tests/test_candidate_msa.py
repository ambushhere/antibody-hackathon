import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import audit_candidate_screen as audit
from scripts import prepare_candidate_msa as msa
from scripts import prepare_candidate_screen as parent
from test_candidate_screen import receipt


def test_frozen_msa_is_exact_parent_panel_with_one_changed_field_per_query():
    original, canonical = parent.panel_bundle()
    variant, variant_canonical = msa.msa_bundle()
    assert canonical == variant_canonical
    assert variant[parent.TEMPLATE] == original[parent.TEMPLATE]
    assert variant[msa.ALIGNMENT] == (msa.MSA_CONFIG / msa.ALIGNMENT).read_bytes()
    assert all((msa.CONFIG / name).read_bytes() == payload for name, payload in variant.items())
    manifest = json.loads(variant["manifest.json"])
    assert manifest["panel"] == json.loads(original["manifest.json"])["panel"]
    assert manifest["panel_counts"] == {"control": 4, "optimization": 10, "denovo": 43}
    assert manifest["parent_input_sha256"] == {k: parent.sha256(v) for k, v in original.items()}
    assert manifest["source_sha256"][msa.FALLBACK] == parent.sha256((parent.ROOT / msa.FALLBACK).read_bytes())
    assert (manifest["run"]["seed"], manifest["run"]["sample"]) == (42, 1)
    queries = json.loads(variant[parent.QUERY_FILE])
    assert len(queries["queries"]) == 57
    for query in queries["queries"]:
        target, heavy, light = query["chains"]
        assert target["chain_id"] == "T"
        assert target["msa_path"] == str(msa.STAGED_DIR / msa.ALIGNMENT)
        assert heavy["msa_path"] is light["msa_path"] is None
        target["msa_path"] = None
    assert queries == json.loads(original[parent.QUERY_FILE])


def test_prepare_refuses_to_change_parent_input():
    with pytest.raises(ValueError, match="Refusing to replace frozen input"):
        msa.prepare(parent.CONFIG)
    assert {name: audit.file_hash(parent.CONFIG / name) for name in msa.PARENT_HASHES} == msa.PARENT_HASHES


@pytest.fixture
def complete_control():
    manifest, _, _, _ = audit.validate_frozen(msa.CONFIG, target_msa=True)
    hashes = manifest["target_msa_control_input_sha256"]
    rows = []
    for query in parent.CONTROL_IDS:
        for seed in audit.SEEDS:
            for sample in audit.SAMPLES:
                # Exact confidence tie at seed 42, higher confidence for sample 1 otherwise.
                confidence = 0.4 if seed == 42 or sample == 1 else 0.3
                contacts = 14 if seed in (42, 43) or sample == 2 else 12
                row = {
                    "arm": "complete", "query_id": query, "seed": seed, "sample": sample,
                    "msa_policy": "target_only", "confidence_chain_order": list(audit.ROLES),
                    "cif_chain_roles": dict(zip(audit.ROLES, audit.ROLES, strict=True)),
                    "row_artifact_sha256": dict.fromkeys(audit.SUFFIXES, "a" * 64),
                    "iptm_H_T": confidence, "iptm_L_T": confidence, audit.MEAN_IPTM: confidence,
                    "wt_contact_count": contacts, "wt_footprint_recall": contacts / 19,
                    "wt_footprint_jaccard": contacts / 19, "target_fit_ca_rmsd_A": 2.0,
                    audit.INTERFACE_RMSD: 4.0, audit.CROP_CLASHES: 0, audit.ECD_CLASHES: 0,
                }
                conditions = audit.operational_conditions(row)
                row.update(operational_conditions=conditions, operational_pose_pass=all(conditions.values()))
                rows.append(row)
    groups = audit.rank_samples(rows)
    return {
        "status": "complete", "msa_policy": "target_only", "protenix_image_id": parent.IMAGE_ID,
        "expected_pose_count": 24, "audited_pose_count": 24, "missing": [], "inputs_sha256": hashes,
        "predictions": rows, "selection": groups,
        "selection_policy": {
            "group_by": ["arm", "query_id", "seed"], "score": "mean of VH–HER2 and VL–HER2 pair ipTM",
            "matrix_order": list(audit.ROLES), "tie_break": "lower sample index",
        },
        "operational_rule": audit.operational_rule(groups, complete=True, eligible_arms=("complete",)),
    }


def write_control(tmp_path, control):
    path = tmp_path / "control.json"
    path.write_text(json.dumps(control))
    return path


def test_complete_control_requires_two_selected_wt_passes(tmp_path, complete_control):
    path = write_control(tmp_path, complete_control)
    result = audit.target_msa_admission(path, complete_control["inputs_sha256"])
    assert result["status"] == "qualified"
    assert result["audit_sha256"] == audit.file_hash(path)
    assert result["passing_wt_top1_seeds"] == [42, 43]
    assert [row["sample"] for row in result["wt_top1"]] == [1, 1, 1]
    assert result["wt_top1"][-1]["operational_pose_pass"] is False


@pytest.mark.parametrize("failure", ["partial", "duplicate", "other_arm", "misranked", "mean", "rule", "alignment",
                                     "pose_hash", "chain_mapping", "wrong_policy"])
def test_unverified_control_cannot_admit(tmp_path, complete_control, failure):
    expected = complete_control["inputs_sha256"].copy()
    control = complete_control
    rows = control["predictions"]
    if failure == "partial":
        rows.pop()  # Non-WT control, so all WT passes still appear to be observed.
    elif failure == "duplicate":
        rows[-1] = copy.deepcopy(rows[-2])
    elif failure == "other_arm":
        rows[-1]["arm"] = "masked"
    elif failure == "misranked":
        control["selection"][0]["top1"]["sample"] = 2
    elif failure == "mean":
        rows[0]["iptm_H_T"] = 0.9
    elif failure == "rule":
        control["operational_rule"]["minimum_passing_seeds"] = 1
    elif failure == "alignment":
        control["inputs_sha256"][msa.ALIGNMENT] = "b" * 64
    elif failure == "pose_hash":
        rows[-1]["row_artifact_sha256"].pop("path")
    elif failure == "chain_mapping":
        rows[-1]["confidence_chain_order"] = ["H", "L", "T"]
    else:
        control["msa_policy"] = "none"
    result = audit.target_msa_admission(write_control(tmp_path, control), expected)
    assert result["status"] == "unknown"
    assert "passing_wt_top1_seeds" not in result


def test_failed_control_keeps_every_candidate_unknown(tmp_path, complete_control):
    rows = complete_control["predictions"]
    wt = next(r for r in rows if r["query_id"] == "control_wt" and r["seed"] == 43 and r["sample"] == 1)
    wt[audit.CROP_CLASHES] = 1
    wt["operational_conditions"] = audit.operational_conditions(wt)
    wt["operational_pose_pass"] = False
    complete_control["selection"] = audit.rank_samples(rows)
    complete_control["operational_rule"] = audit.operational_rule(
        complete_control["selection"], complete=True, eligible_arms=("complete",),
    )
    report = audit.audit(tmp_path / "absent", Path("1n8z"), Path("6oge"), target_msa=True,
                         control_audit=write_control(tmp_path, complete_control))
    assert report["status"] == "not_run_admission_failed"
    assert report["baseline_admission"]["passing_wt_top1_seeds"] == [42]
    assert report["run"] is None
    assert report["audited_pose_count"] == 0 and report["unknown_pose_count"] == 57
    assert report["track_counts"] == {"control": {"unknown": 4}, "optimization": {"unknown": 10},
                                      "denovo": {"unknown": 43}}
    assert all(r["descriptive_pareto_front"] is None for r in report["predictions"])
    audit.write_audit(report, tmp_path / "output")
    with (tmp_path / "output/scores.csv").open(newline="") as stream:
        assert len(list(csv.DictReader(stream))) == 57


def test_pending_control_never_scores_present_candidate_files(tmp_path, monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("Unverified controls cannot admit campaign receipt or pose scoring")

    monkeypatch.setattr(audit, "validate_staged", unexpected)
    monkeypatch.setattr(audit, "score_pose", unexpected)
    report = audit.audit(tmp_path, Path("1n8z"), Path("6oge"), target_msa=True, control_audit=tmp_path / "absent")
    assert report["status"] == "unverified"
    assert report["unknown_pose_count"] == 57


@pytest.fixture
def candidate_receipt(tmp_path):
    _, _, queries, hashes = audit.validate_frozen(msa.CONFIG, target_msa=True)
    ids = list(queries)
    campaign = tmp_path / "campaign"
    receipt(campaign, ids, ids[:2], config_dir=msa.CONFIG, staged_dir=msa.STAGED_DIR, target_msa=True)
    return campaign, ids, hashes


@pytest.mark.parametrize("input_name", [msa.ALIGNMENT, parent.QUERY_FILE, parent.TEMPLATE, "manifest.json"])
def test_receipted_inputs_require_exact_target_msa_staging_paths(candidate_receipt, input_name):
    campaign, _, hashes = candidate_receipt
    assert audit.validate_staged(campaign, hashes, target_msa=True)[0] == hashes
    inputs = campaign / "input.sha256"
    path = str(msa.STAGED_DIR / input_name)
    inputs.write_text(inputs.read_text().replace(path, f"/some/other/location/{input_name}"))
    # Matching digests with altered paths remain acceptable in the original mode.
    assert audit.validate_staged(campaign, hashes)[0] == hashes
    with pytest.raises(ValueError, match="paths differ"):
        audit.validate_staged(campaign, hashes, target_msa=True)


@pytest.mark.parametrize("bad_marker", [None, "control_digest", "manifest_digest", "timestamp"])
def test_optional_receipts_bind_the_candidate_run_to_admission(candidate_receipt, complete_control, tmp_path, bad_marker):
    campaign, ids, hashes = candidate_receipt
    control_path = write_control(tmp_path, complete_control)
    admission = audit.target_msa_admission(control_path, complete_control["inputs_sha256"])
    _, paths, _ = audit.validate_staged(campaign, hashes, target_msa=True)
    _, run = audit.validate_run(campaign, ids, paths)
    assert audit.validate_msa_run_markers(campaign, run, admission, hashes) == {}
    record = {
        "admitted": True, "selected_arm": "complete_target_msa", "expected_poses": 57,
        "control_audit_sha256": "wrong" if bad_marker == "control_digest" else admission["audit_sha256"],
        "input_manifest_sha256": "wrong" if bad_marker == "manifest_digest" else hashes["manifest.json"],
    }
    (campaign / "admission.json").write_text(json.dumps(record))
    (campaign / "started.txt").write_text("2026-11-14T14:44:59Z\n")
    (campaign / "completed.txt").write_text("2026-11-14T14:49:59Z\n" if bad_marker == "timestamp"
                                           else "2026-11-14T14:50:01Z\n")
    if bad_marker is None:
        assert audit.validate_msa_run_markers(campaign, run, admission, hashes) == {
            name: audit.file_hash(campaign / name) for name in ("admission.json", "started.txt", "completed.txt")
        }
    else:
        with pytest.raises(ValueError, match="does not (match|bound)"):
            audit.validate_msa_run_markers(campaign, run, admission, hashes)


def test_qualified_mode_scans_only_receipted_poses_and_preserves_unknowns(
    candidate_receipt, complete_control, tmp_path, monkeypatch,
):
    campaign, ids, _ = candidate_receipt
    (campaign / audit.artifact_paths(ids[0], 42, 1)["scores_path"]).unlink()
    monkeypatch.setattr(audit, "validate_file", lambda _, digest: digest)
    monkeypatch.setattr(audit, "read_1n8z", lambda *_: (
        None, dict.fromkeys(audit.COMMON), None, {}, list(range(19)), audit.COMMON, "1n8z",
    ))
    monkeypatch.setattr(audit.gemmi.cif, "read", lambda *_: SimpleNamespace(sole_block=lambda: None))
    monkeypatch.setattr(audit, "ecd_reference", lambda *_: (None, None, []))
    monkeypatch.setattr(audit, "score_pose", lambda *_: {audit.CROP_CLASHES: 0})
    report = audit.audit(campaign, Path("1n8z"), Path("6oge"), target_msa=True,
                         control_audit=write_control(tmp_path, complete_control))
    assert report["status"] == "incomplete"
    assert report["audited_pose_count"] == 1 and report["unknown_pose_count"] == 56
    assert "Missing receipted artifacts" in report["predictions"][0]["reason"]
    assert report["predictions"][1][audit.CROP_CLASHES] == 0
    assert "no prediction receipt" in report["predictions"][2]["reason"]
    assert report["candidate_promotion"] == "requires_replicates"
    assert all(front["front"] is None for front in report["descriptive_pareto"].values())
