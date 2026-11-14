"""Audit the fixed 57-query HER2 candidate screen on CPU."""

import argparse
import csv
import json
import re
from collections import Counter
from datetime import datetime
from decimal import Decimal
from pathlib import Path

import gemmi
import numpy as np
from biotite import InvalidFileError

from scripts.audit_context_panel import (
    REFERENCE_CHAINS, antibody_mapping, exact_coverage, file_hash, footprint_pae,
    geometry_row, protenix_pae, score_pairs, validate_file, validate_prediction,
)
from scripts.audit_context_panel import read_reference as read_1n8z
from scripts.audit_template_controls import (
    COMMON, SCHEMA, SUFFIXES, artifact_paths, ecd_reference, operational_conditions, operational_rule,
    overlay_clashes, rank_samples,
)
from scripts import prepare_candidate_msa as msa
from scripts.prepare_candidate_screen import (
    ADMISSION, CONFIG, CONTROL_IDS, IMAGE_ID, POSITIONS, QUERY_FILE, ROLES, ROOT, SAMPLE, SEED,
    SOURCE_SHA256, TEMPLATE, chain_sequences, panel_bundle, require, sha256,
)
from scripts.prepare_her2_template_controls import CONTEXT, SAMPLES, SEEDS

CAMPAIGN = ROOT / "results/campaign/candidate-screen"
CONTROL_AUDIT = ROOT / "results/template-controls/audit.json"
OUTPUT = ROOT / "results/candidate-screen"
REFERENCE = Path("/tmp/her2-workers/template-controls-1n8z-source.cif")
TEMPLATE_SOURCE = Path("/tmp/her2-workers/template-controls-6oge-source.cif")
MSA_CAMPAIGN = ROOT / "results/campaign/candidate-msa-screen"
MSA_CONTROL_AUDIT = ROOT / "results/template-controls/target-msa/audit.json"
MSA_OUTPUT = ROOT / "results/candidate-msa-screen"

CROP_CLASHES = "heavy_atom_clash_pairs_below_1p5A"
ECD_CLASHES = "ecd_6oge_outside_crop_heavy_atom_clash_pairs_below_1p5A"
INTERFACE_RMSD = "ref_interface_ca_rmsd_in_target_frame_A"
MEAN_IPTM = "mean_antibody_target_pair_iptm"
REFERENCE_GEOMETRY = (
    "antibody_mapped_ca_count", "antibody_ca_rmsd_in_target_frame_A",
    "antibody_centroid_displacement_A", "antibody_ref_interface_ca_count",
    "antibody_ref_interface_ca_rmsd_in_target_frame_A", "ref_interface_ca_count", INTERFACE_RMSD,
    "ref_interface_pae_ab_rows_target_cols_A", "ref_interface_pae_target_rows_ab_cols_A",
)
POSE_FIELDS = (
    "cif_chain_roles", "confidence_chain_order", "overall_iptm", "provider_clash_flag", MEAN_IPTM,
    "vh_vl_iptm", "max_antibody_target_pair_iptm", "iptm_H_T", "iptm_L_T",
    "target_contacts_p04626", "target_contacts_outside_common_p04626", "wt_positions_recovered_p04626",
    "wt_contact_count", "wt_footprint_recall", "wt_footprint_precision", "wt_footprint_jaccard",
    "antibody_contact_residues", "heavy_atom_contact_pairs", CROP_CLASHES, "minimum_target_distance_A",
    "contact_cutoff_A", "target_fit_ca_count", "target_fit_ca_rmsd_A", *REFERENCE_GEOMETRY,
    "ecd_6oge_common_target_ca_rmsd_A", ECD_CLASHES, "ecd_6oge_outside_crop_clash_positions_p04626",
    "ecd_6oge_outside_crop_minimum_distance_A", "predicted_interface_pae_ab_rows_target_cols_A",
    "predicted_interface_pae_target_rows_ab_cols_A",
)
PARETO_AXES = {
    MEAN_IPTM: 1, "wt_contact_count": 1, "target_fit_ca_rmsd_A": -1,
    CROP_CLASHES: -1, ECD_CLASHES: -1,
}
# Expected read/parse failures become unknown rows. Programming errors still surface.
INPUT_ERRORS = (OSError, ValueError, KeyError, TypeError, InvalidFileError, EOFError, RuntimeError)


def read_json(path):
    return json.loads(path.read_bytes())


def validate_frozen(config, *, target_msa=False):
    payloads, canonical = msa.msa_bundle() if target_msa else panel_bundle()
    hashes = {name: validate_file(config / name, sha256(payload)) for name, payload in payloads.items()}
    manifest = json.loads(payloads["manifest.json"])
    queries = {q["query_id"]: chain_sequences(q) for q in json.loads(payloads[QUERY_FILE])["queries"]}
    return manifest, canonical, queries, hashes


def campaign_file(campaign, path):
    file = campaign / path
    require(file.resolve().is_relative_to(campaign.resolve()), f"Artifact escapes its run directory: {path}")
    return file


def validate_staged(campaign, expected, *, target_msa=False):
    staged = {}
    paths = {}
    for line in campaign_file(campaign, "input.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        name = Path(path).name
        require(name not in staged and Path(path).is_absolute(), f"Invalid staged input: {name}")
        staged[name], paths[name] = digest, path
    require(staged == expected, "Staged input hashes differ from frozen files")
    if target_msa:
        # ASVS 5.3.2: Match the exact VM paths, including the target alignment embedded in the queries.
        require(paths == {name: str(msa.STAGED_DIR / name) for name in expected},
                "Staged input paths differ from the frozen target MSA directory")
    image = campaign_file(campaign, "image-id.txt").read_text().strip()
    require(image == IMAGE_ID, "Unexpected Protenix v2 image ID")
    return staged, paths, image


def validate_run_log(log_path, query_count):
    log = log_path.read_text()
    assigned = set(re.findall(r"Assigned templates to (\d+) entities across (\d+) queries", log))
    polymer = set(re.findall(r"Read (\d+) chains? with a total of (\d+) residues from `_pdbx_poly_seq_scheme`", log))
    resolved = set(re.findall(r"Identified (\d+) resolved residue\(s\) from", log))
    atoms = set(re.findall(r"Read (\d+) atoms from `_atom_site` table", log))
    seeds = set(re.findall(r"Seed set to (\d+)", log))
    constraints = set(re.findall(r"Loaded constraint feature: #atom contact:(\d+) #contact:(\d+) #pocket:(\d+)", log))
    count = str(query_count)
    require(assigned == {(count, count)} and polymer == {("1", "119")} and resolved == {"119"}
            and atoms == {"912"} and seeds == {str(SEED)} and constraints == {("0", "0", "0")},
            "Run log does not confirm the complete template, template assignments, seed and constraints")


def validate_run(campaign, ids, staged_paths):
    result_path = campaign_file(campaign, "result.json")
    exact_coverage(campaign.rglob("result.json"), [result_path], "Single-folder result receipt")
    config_path = campaign_file(campaign, "config.json")
    log_path = campaign_file(campaign, "run.log")
    run = read_json(result_path)
    provenance = run["provenance"]
    config = {
        "missing_msas": "fail", "num_diffusion_samples": 1, "query_failures": "fail",
        "seeds": [SEED], "template_paths": [staged_paths[TEMPLATE]],
        "weights": "/app/weights/protenix-v2.pt", "write_detailed_scores": True,
    }
    require(run["task"] == "predict" and run["schema_version"] == provenance["cli_contract_version"] == SCHEMA
            and provenance["module"] == {"name": "apheris-protenix", "version": "1.3.0"}
            and provenance["inputs"] == [{"path": staged_paths[QUERY_FILE], "tier": "none"}]
            and provenance["resolved_config"] == config
            and read_json(config_path) == {"schema_version": SCHEMA, **config},
            "Foundry receipt/config disagrees with the pinned model, runtime settings or inputs")
    started, ended = [datetime.fromisoformat(provenance[k]) for k in ("started_at", "ended_at")]
    require(started.tzinfo is not None and ended.tzinfo is not None and started <= ended,
            "Missing timezone or inverted run timestamps")
    validate_run_log(log_path, len(ids))

    indexes = dict(zip(ids, range(len(ids)), strict=True))
    result = run["result"]
    require(set(result) == {"predictions"} and isinstance(result["predictions"], list),
            "Unexpected Foundry prediction receipt schema")
    received = {}
    for prediction in result["predictions"]:
        query = prediction["query_id"]
        require(query in indexes and query not in received, f"Unexpected or duplicate query in receipt: {query}")
        require((prediction["query_index"], prediction["seed"], prediction["sample"])
                == (indexes[query], SEED, SAMPLE), f"Wrong query index, seed or sample in receipt: {query}")
        paths = artifact_paths(query, SEED, SAMPLE)
        require(all(prediction[key] == path for key, path in paths.items()),
                f"Unexpected prediction artifact path: {query}")
        received[query] = {key: campaign_file(campaign, path) for key, path in paths.items()}
    for key, suffix in SUFFIXES.items():
        known = {files[key] for files in received.values()}
        unreceipted = sorted(p.relative_to(campaign).as_posix()
                            for p in campaign.rglob(f"*{suffix}") if p not in known)
        require(not unreceipted, f"Unreceipted {key} artifacts: {unreceipted}")

    return received, {
        "receipt_sha256": {p.name: file_hash(p) for p in (result_path, config_path, log_path)},
        "started_at": provenance["started_at"], "ended_at": provenance["ended_at"],
        "config": config, "receipted_pose_count": len(received),
    }


def control_admission(path, template_sha256):
    if not path.is_file():
        return {"status": "unknown", "reason": "baseline control audit unavailable"}
    report = {"audit_sha256": file_hash(path), "arm": "complete", "msa_policy": "none"}
    try:
        baseline = read_json(path)
        rule = baseline["operational_rule"]
        arm = rule["arms"]["complete"]
        observed, passing = arm["wt_top1_seeds_observed"], arm["passing_wt_top1_seeds"]
        # Check that the copied audit reports the same complete no-MSA template and admission rule.
        require(baseline["msa_policy"] == "none"
                and baseline["protenix_image_id"] == IMAGE_ID
                and baseline["inputs_sha256"][TEMPLATE] == template_sha256
                and rule["status"] in ("evaluated", "incomplete") and rule["minimum_passing_seeds"] == 2
                and rule["source"] == ADMISSION,
                "Baseline control audit uses a different setup or rule")
        complete = [row for row in baseline["predictions"] if row["arm"] == "complete"]
        exact_coverage(
            [(r["query_id"], r["seed"], r["sample"]) for r in complete],
            [(query, seed, sample) for query in CONTROL_IDS for seed in SEEDS for sample in SAMPLES],
            "Complete no-MSA baseline pose receipts",
        )
        wt = [row for row in complete if row["query_id"] == "control_wt"]
        top = [min((r for r in wt if r["seed"] == seed), key=lambda r: (-r[MEAN_IPTM], r["sample"]))
               for seed in SEEDS]
        verified_passing = [row["seed"] for row in top if all(operational_conditions(row).values())]
        qualifies = len(verified_passing) >= 2
        require(sorted(observed) == list(SEEDS) and sorted(passing) == verified_passing
                and arm["qualifies"] is (qualifies if rule["status"] == "evaluated" else None),
                "Inconsistent complete/no-MSA WT top-1 coverage or decision")
        return {**report, "status": "qualified" if qualifies else "failed",
                "wt_top1_seeds_observed": observed, "passing_wt_top1_seeds": passing,
                "wt_top1": [
                    {k: row[k] for k in ("seed", "sample", MEAN_IPTM, "wt_contact_count",
                                         "target_fit_ca_rmsd_A", INTERFACE_RMSD, CROP_CLASHES)}
                    for row in top
                ]}
    except INPUT_ERRORS as error:
        return {**report, "status": "unknown", "reason": str(error)}


def require_control_coverage(control, expected_inputs):
    require(control["msa_policy"] == "target_only" and control["protenix_image_id"] == IMAGE_ID
            and control["inputs_sha256"] == expected_inputs,
            "Control audit uses a different MSA, template, query panel or image")
    require(control["status"] == "complete" and control["missing"] == []
            and control["expected_pose_count"] == control["audited_pose_count"] == 24,
            "Target MSA admission requires all 24 control poses to be audited")
    rows = control["predictions"]
    exact_coverage(
        [(row["arm"], row["query_id"], row["seed"], row["sample"]) for row in rows],
        [("complete", q, seed, sample) for q in CONTROL_IDS for seed in SEEDS for sample in SAMPLES],
        "Complete target MSA control pose receipts",
    )
    require(all(type(row["seed"]) is int and type(row["sample"]) is int for row in rows),
            "Control audit seed and sample identifiers must be integers")
    return rows


def verify_control_pose(row):
    require(row["msa_policy"] == "target_only" and row["confidence_chain_order"] == list(ROLES)
            and row["cif_chain_roles"] == dict(zip(ROLES, ROLES, strict=True)),
            "Control pose does not have the expected MSA and chain mapping")
    require(set(row["row_artifact_sha256"]) == set(SUFFIXES)
            and all(re.fullmatch(r"[0-9a-f]{64}", digest) for digest in row["row_artifact_sha256"].values()),
            "Control pose is missing verified artifact hashes")
    pairs = [row[k] for k in ("iptm_H_T", "iptm_L_T")]
    require(all(type(ip) in (int, float) and 0 <= ip <= 1 for ip in pairs)
            and row[MEAN_IPTM] == float(sum(Decimal(str(ip)) for ip in pairs) / 2),
            "Control pose mean pair ipTM disagrees with the chain pair scores")
    require(type(row["wt_contact_count"]) is int and 0 <= row["wt_contact_count"] <= 19
            and type(row[CROP_CLASHES]) is int and row[CROP_CLASHES] >= 0
            and all(type(row[key]) in (int, float) and 0 <= row[key] < float("inf")
                    for key in ("target_fit_ca_rmsd_A", INTERFACE_RMSD)),
            "Control pose has invalid gate metrics")
    conditions = operational_conditions(row)
    require(row["operational_conditions"] == conditions
            and row["operational_pose_pass"] is all(conditions.values()),
            "Control pose gate does not match its metrics")


def verified_msa_selection(control, rows):
    selection = control["selection_policy"]
    require(selection["group_by"] == ["arm", "query_id", "seed"]
            and selection["score"] == "mean of VH–HER2 and VL–HER2 pair ipTM"
            and selection["matrix_order"] == list(ROLES) and selection["tie_break"] == "lower sample index",
            "Target MSA control uses a different top-1 selection rule")
    for row in rows:
        verify_control_pose(row)
    groups = rank_samples([row.copy() for row in rows])
    require(control["selection"] == groups, "Control selection disagrees with the frozen ranking rule")
    verified = operational_rule(groups, complete=True, eligible_arms=("complete",))
    require(control["operational_rule"] == verified, "Target MSA control disagrees with the frozen reproduction rule")
    return groups, verified["arms"]["complete"]


def target_msa_admission(path, expected_inputs):
    report = {"arm": "complete", "msa_policy": "target_only"}
    if not path.is_file():
        return {**report, "status": "unknown", "reason": "Target MSA control audit unavailable"}
    try:
        payload = path.read_bytes()
        report["audit_sha256"] = sha256(payload)
        control = json.loads(payload)
        rows = require_control_coverage(control, expected_inputs)
        groups, complete = verified_msa_selection(control, rows)
        return {**report, "status": "qualified" if complete["qualifies"] else "failed",
                "audited_control_pose_count": len(rows),
                "wt_top1_seeds_observed": complete["wt_top1_seeds_observed"],
                "passing_wt_top1_seeds": complete["passing_wt_top1_seeds"],
                "wt_top1": [{"seed": g["seed"], **g["top1"]} for g in groups if g["query_id"] == "control_wt"]}
    except INPUT_ERRORS as error:
        return {**report, "status": "unknown", "reason": str(error)}


def validate_msa_run_markers(campaign, run, admission, hashes):
    markers = {}
    path = campaign_file(campaign, "admission.json")
    if path.exists():
        payload = path.read_bytes()
        receipt = json.loads(payload)
        require(receipt["admitted"] is True and receipt["selected_arm"] == "complete_target_msa"
                and receipt["expected_poses"] == 57
                and receipt["control_audit_sha256"] == admission["audit_sha256"]
                and receipt["input_manifest_sha256"] == hashes["manifest.json"],
                "Run admission receipt does not match the qualified control and candidate manifest")
        markers[path.name] = sha256(payload)
    for marker, bound, start in (("started.txt", "started_at", True), ("completed.txt", "ended_at", False)):
        path = campaign_file(campaign, marker)
        if path.exists():
            payload = path.read_bytes()
            timestamp = datetime.fromisoformat(payload.decode().strip())
            receipt_time = datetime.fromisoformat(run[bound])
            require(timestamp.tzinfo is not None
                    and (timestamp <= receipt_time if start else timestamp >= receipt_time.replace(microsecond=0)),
                    f"{marker} does not bound the Foundry run timestamp")
            markers[path.name] = sha256(payload)
    return markers


def predicted_interface_pae(pae, mapped, geometry):
    tokens = mapped[mapped.atom_name == "CA"]
    target = np.flatnonzero((tokens.chain_id == "T")
                            & np.isin(tokens.res_id, geometry["target_contacts_p04626"]))
    sites = set(geometry["antibody_contact_residues"])
    antibody = [i for i, (chain, index) in enumerate(zip(tokens.chain_id, tokens.res_id))
                if f"{chain}:{index}" in sites]
    require(len(target) == len(geometry["target_contacts_p04626"]) and len(antibody) == len(sites),
            "Predicted interface PAE residue mapping is incomplete")
    return {
        "predicted_interface_pae_ab_rows_target_cols_A": (
            float(pae[np.ix_(antibody, target)].mean()) if len(target) and antibody else None),
        "predicted_interface_pae_target_rows_ab_cols_A": (
            float(pae[np.ix_(target, antibody)].mean()) if len(target) and antibody else None),
    }


def score_pose(files, expected, entry, panel, reference, ecd):
    atoms, roles = validate_prediction(files["path"], expected, protenix=True)
    require(tuple(roles) == tuple(roles.values()) == ROLES, "Unexpected Protenix confidence chain order")
    summary = read_json(files["scores_path"])
    tokens = int(np.count_nonzero(atoms.atom_name == "CA"))
    require((summary["num_chains"], summary["num_atoms"], summary["num_tokens"], summary["num_residues"])
            == (len(ROLES), len(atoms), tokens, tokens), "Structure and confidence counts disagree")
    overall = summary["interface_predicted_tm_score"]
    provider_clash = summary["clash_flag"]
    require(type(overall) in (float, int) and 0 <= overall <= 1 and isinstance(provider_clash, bool),
            "Invalid global confidence or provider clash flag")
    matrix = np.asarray(summary["chain_chain_interface_predicted_tm_score"], dtype=float)
    pairs = score_pairs(ROLES, matrix)
    pae = protenix_pae(files["detailed_scores_path"], atoms, tokens)

    matched = entry["track"] != "denovo"
    require(not matched or antibody_mapping(entry["query_id"], panel), "WT scaffold mapping is unavailable")
    mapped, geometry = geometry_row(atoms, roles, reference, {"T": list(POSITIONS)}, matched)
    common_contact_count = len(set(geometry["target_contacts_p04626"]) & set(reference[5]))
    reference_pae = footprint_pae(pae, mapped, roles, reference[4], reference[3]) if matched else {}
    return {
        **dict.fromkeys(REFERENCE_GEOMETRY),
        "cif_chain_roles": roles, "confidence_chain_order": list(ROLES),
        "overall_iptm": overall, "provider_clash_flag": provider_clash,
        MEAN_IPTM: (pairs["iptm_H_T"] + pairs["iptm_L_T"]) / 2,
        "wt_footprint_precision": (
            geometry["wt_contact_count"] / common_contact_count if common_contact_count else None),
        **pairs, **geometry, **reference_pae, **predicted_interface_pae(pae, mapped, geometry),
        **overlay_clashes(mapped, ecd),
    }


def pose_row(entry, files, queries, panel, reference, ecd, campaign, run_issue):
    query = entry["query_id"]
    row = {
        "query_id": query, "track": entry["track"], "seed": SEED, "sample": SAMPLE,
        "cdrh3": entry["cdrh3"], "paired_sequence_sha256": entry["paired_sequence_sha256"],
        "paired_length": entry["paired_length"], "source": entry["source"],
        "source_sha256": entry["source_sha256"], "source_record_sha256": entry["source_record_sha256"],
        "generation": entry.get("generation"), "shortlist_pick": entry.get("shortlist_pick"),
        "local_eligibility": entry.get("local_eligibility"),
        "official_validation": entry.get("official_adaptyv_novelty", entry.get("official_validation")),
        "control_spr_kd_nM": entry.get("kd_nM"), "control_facs_class": entry.get("facs_class"),
        "status": "unknown", "reason": run_issue,
        "wt_antibody_mapping": "inapplicable" if entry["track"] == "denovo" else "substitution_only",
        "wt_interface_rmsd_null_reason": (
            "de novo scaffold has no WT antibody residue correspondence" if entry["track"] == "denovo"
            else "pose unavailable"
        ),
        "row_artifact_paths": None, "row_artifact_sha256": None,
        **dict.fromkeys(POSE_FIELDS), "descriptive_pareto_front": None,
    }
    if files is None:
        row["reason"] = run_issue or "expected query has no prediction receipt"
        return row
    row["row_artifact_paths"] = {k: p.relative_to(campaign).as_posix() for k, p in files.items()}
    try:
        row["row_artifact_sha256"] = {k: file_hash(p) if p.is_file() else None for k, p in files.items()}
        missing = [key for key, digest in row["row_artifact_sha256"].items() if digest is None]
        require(not missing, f"Missing receipted artifacts: {', '.join(missing)}")
        if run_issue is not None:
            return row
        row.update(score_pose(files, queries[query], entry, panel, reference, ecd))
        row["status"], row["reason"] = "audited", None
        if entry["track"] != "denovo":
            row["wt_interface_rmsd_null_reason"] = None
    except INPUT_ERRORS as error:
        row["reason"] = str(error)
    return row


def descriptive_pareto(rows, supported):
    fronts = {}
    for track in ("optimization", "denovo"):
        axes = dict(PARETO_AXES)
        if track == "optimization":
            axes[INTERFACE_RMSD] = -1
        eligible = [row for row in rows if row["track"] == track and row["status"] == "audited"
                    and all(row[k] is not None for k in axes)]
        fronts[track] = {"axes": {k: "maximize" if sign > 0 else "minimize" for k, sign in axes.items()},
                         "supported_count": len(eligible), "front": None}
        if not supported or not eligible:
            continue
        values = np.asarray([[row[k] * sign for k, sign in axes.items()] for row in eligible])
        for i, row in enumerate(eligible):
            dominates = (values >= values[i]).all(axis=1) & (values > values[i]).any(axis=1)
            row["descriptive_pareto_front"] = not bool(dominates.any())
        fronts[track]["front"] = [row["query_id"] for row in eligible if row["descriptive_pareto_front"]]
    return fronts


def audit(campaign, reference_path, template_source, config=None, control_audit=None, *, target_msa=False):
    config = config or (msa.CONFIG if target_msa else CONFIG)
    control_audit = control_audit or (MSA_CONTROL_AUDIT if target_msa else CONTROL_AUDIT)
    manifest, canonical, queries, hashes = validate_frozen(config, target_msa=target_msa)
    panel = {row["query_id"]: row for row in manifest["panel"]}
    admission = (target_msa_admission(control_audit, manifest["target_msa_control_input_sha256"])
                 if target_msa else control_admission(control_audit, hashes[TEMPLATE]))

    received, run, run_issue, staged = {}, None, None, None
    if target_msa and admission["status"] != "qualified":
        reason = admission.get("reason", "Complete target MSA WT control failed the admission rule.")
        run_issue = (f"{reason} Raw campaign artifacts are ineligible for scoring." if campaign.is_dir()
                     else f"{reason} No candidate campaign supplied.")
    elif not campaign.is_dir():
        run_issue = (
            "Complete no-MSA WT control arm failed the admission rule. Candidate poses were not run."
            if admission["status"] == "failed" else "Raw campaign directory unavailable."
        )
    else:
        try:
            staged, paths, image = validate_staged(campaign, hashes, target_msa=target_msa)
            candidate_files, candidate_run = validate_run(campaign, list(queries), paths)
            if target_msa:
                candidate_run["staged_input_paths"] = paths
                candidate_run["receipt_sha256"].update({
                    name: file_hash(campaign_file(campaign, name)) for name in ("input.sha256", "image-id.txt")
                })
                candidate_run["optional_receipt_sha256"] = validate_msa_run_markers(
                    campaign, candidate_run, admission, hashes)
            received, run = candidate_files, {**candidate_run, "image_id": image}
        except INPUT_ERRORS as error:
            run_issue = str(error)

    reference, ecd, references = None, None, None
    if received:
        try:
            validate_file(template_source, SOURCE_SHA256)
            context = read_json(CONTEXT / "manifest.json")
            reference = read_1n8z(
                reference_path, context, canonical, queries["control_wt"], {"split": {"T": COMMON}},
            )
            require(set(reference[1]) & set(POSITIONS) == set(COMMON) and len(reference[4]) == 19,
                    "Unexpected shared target or WT footprint coverage")
            ecd = ecd_reference(gemmi.cif.read(str(template_source)).sole_block(), canonical)
            references = {
                "1n8z_sha256": reference[-1], "1n8z_antibody_role_to_chain": REFERENCE_CHAINS,
                "target_fit_positions_p04626": COMMON, "wt_footprint_positions_p04626": sorted(reference[4]),
                "wt_antibody_interface_positions": reference[3], "6oge_sha256": SOURCE_SHA256,
                "6oge_outside_crop_heavy_atoms": len(ecd[2]),
            }
        except INPUT_ERRORS as error:
            run_issue = f"Reference geometry unavailable. {error}"

    rows = [pose_row(entry, received.get(query), queries, panel, reference, ecd, campaign, run_issue)
            for query, entry in panel.items()]
    known = sum(row["status"] == "audited" for row in rows)
    status = "complete" if known == len(rows) else "incomplete" if run else "unverified"
    if not campaign.is_dir() and admission["status"] == "failed":
        status = "not_run_admission_failed"
    elif target_msa and not campaign.is_dir():
        status = "awaiting_candidate_artifacts" if admission["status"] == "qualified" else "not_run_admission_pending"
    track_counts = {track: dict(Counter(row["status"] for row in rows if row["track"] == track))
                    for track in manifest["panel_counts"]}
    return {
        "status": status, "expected_pose_count": len(rows), "audited_pose_count": known,
        "unknown_pose_count": len(rows) - known, "track_counts": track_counts,
        "run_issue": run_issue, "baseline_admission": admission, "run": run,
        "frozen_input_sha256": hashes, "staged_input_sha256": staged,
        "source_sha256": manifest["source_sha256"], "reference": references,
        "reference_data_origin": manifest["data_origin"],
        "replicates_required": True,
        "descriptive_pareto": descriptive_pareto(rows, status == "complete" and admission["status"] == "qualified"),
        "predictions": rows,
        **({"msa_policy": "target_only", "candidate_promotion": "requires_replicates"} if target_msa else {}),
    }


def result_markdown(report):
    counts = report["track_counts"]
    admission = report["baseline_admission"]
    target_msa = report.get("msa_policy") == "target_only"
    title = "HER2 candidate refold screen with target MSA" if target_msa else "HER2 candidate refold screen results"
    setup = "complete target MSA" if target_msa else "complete, no MSA"
    lines = [
        f"# {title}", "",
        *([f"Audit status is `{report['status']}`.", ""] if target_msa else []),
        f"Audited {report['audited_pose_count']} of {report['expected_pose_count']} expected poses. "
        f"{report['unknown_pose_count']} remain unknown. "
        f"The {setup} control admission is **{admission['status']}**.", "",
        "| Track | Audited | Unknown |", "| --- | ---: | ---: |",
        *(f"| {label} | {counts[track].get('audited', 0)} | {counts[track].get('unknown', 0)} |"
          for track, label in (("control", "Controls"), ("optimization", "Optimization"), ("denovo", "De novo"))),
        "",
    ]
    if report["run_issue"]:
        lines.extend([report["run_issue"], ""])
    if admission.get("passing_wt_top1_seeds") is not None:
        seeds = ", ".join(map(str, admission["passing_wt_top1_seeds"])) or "none"
        lines.extend([f"Passing baseline WT top-pose seeds: {seeds}.", ""])
    controls = [row for row in report["predictions"] if row["track"] == "control"]
    lines.extend(["## Matched controls", "",
                  "| Control | Historical SPR KD (nM) | Mean pair ipTM | WT contacts (of 19) | "
                  "Target fit (Å) | WT interface (Å) | Crop clashes | ECD clashes |",
                  "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"])
    for row in controls:
        keys = ("control_spr_kd_nM", MEAN_IPTM, "wt_contact_count", "target_fit_ca_rmsd_A",
                INTERFACE_RMSD, CROP_CLASHES, ECD_CLASHES)
        values = ["—" if row[k] is None else f"{row[k]:.5g}" for k in keys]
        lines.append("| " + " | ".join([row["query_id"], *values]) + " |")
    lines.extend(["", "## Descriptive Pareto fronts", ""])
    if any(info["front"] is None for info in report["descriptive_pareto"].values()):
        lines.append("Fronts require a complete, admitted panel and known values for the listed axes.")
    for track, info in report["descriptive_pareto"].items():
        if info["front"] is not None:
            lines.append(f"{track.capitalize()} ({len(info['front'])} of {info['supported_count']}): "
                         + ", ".join(f"`{query}`" for query in info["front"]) + ".")
    lines.extend([
        "",
        "The protocol specifies one model draw per query at seed 42 and requires replicates before "
        "advancing candidates. Affinity, competitive binding and melting temperature require "
        "experimental measurements. De novo eligibility is local and Adaptyv validation is pending.",
    ])
    return "\n".join(lines) + "\n"


def write_audit(report, output):
    output.mkdir(parents=True, exist_ok=True)
    (output / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    rows = report["predictions"]
    with (output / "scores.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(",", ":")) if isinstance(v, (dict, list)) else v
                             for k, v in row.items()})
    (output / "report.md").write_text(result_markdown(report))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-msa", action="store_true", help="Audit the separate target-only MSA experiment")
    parser.add_argument("--check-admission", action="store_true", help="Verify target MSA controls without scoring poses")
    parser.add_argument("--campaign", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--template-source", type=Path, default=TEMPLATE_SOURCE)
    parser.add_argument("--control-audit", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check_admission:
        if not args.target_msa:
            parser.error("--check-admission requires --target-msa")
        manifest, _, _, _ = validate_frozen(args.config or msa.CONFIG, target_msa=True)
        admission = target_msa_admission(args.control_audit or MSA_CONTROL_AUDIT,
                                         manifest["target_msa_control_input_sha256"])
        print(json.dumps(admission, indent=2, allow_nan=False))
        raise SystemExit(0 if admission["status"] == "qualified" else 2)
    campaign = args.campaign or (MSA_CAMPAIGN if args.target_msa else CAMPAIGN)
    output = args.output or (MSA_OUTPUT if args.target_msa else OUTPUT)
    require(not args.target_msa or output.resolve() != OUTPUT.resolve(),
            "Target MSA audit cannot write to the no-MSA results directory")
    report = audit(campaign, args.reference, args.template_source, args.config, args.control_audit,
                   target_msa=args.target_msa)
    write_audit(report, output)
    print(f"{report['status']}: {report['audited_pose_count']}/{report['expected_pose_count']} poses audited")
    raise SystemExit(0 if report["status"] == "complete" and report["baseline_admission"]["status"] == "qualified"
                     else 2)


if __name__ == "__main__":
    main()
