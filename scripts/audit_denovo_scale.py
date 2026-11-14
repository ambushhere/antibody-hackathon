"""Audit two completed paired-Fv RFD4 arms using local sequence gates only.

This audit records provider provenance and independently numbers the explicit
H and L structure chains. It does not query the Adaptyv novelty database and
does not predict binding, competition, or thermal stability.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path, PurePosixPath
from typing import Any

from biotite.structure import get_residues
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.sequence import ProteinSequence
from immunum import Annotator
from rapidfuzz.distance import Levenshtein

from .audit_denovo_identity import numbered_identity
from src.her2_candidates import VH_PREFIX, VH_SUFFIX, WT_CDRH3, TRASTUZUMAB_VL


ARM_NAMES = ("denovo-scale-antibody-k1", "denovo-scale-her2novel-k1")
EXPECTED_SEEDS = {
    "denovo-scale-antibody-k1": 20261118,
    "denovo-scale-her2novel-k1": 20261119,
}
EXPECTED_CHECKPOINTS = {
    "denovo-scale-antibody-k1": "antibody-2026-10-30",
    "denovo-scale-her2novel-k1": "her2-novel-2026-11-09",
}
EXPECTED_MODULE = ("apheris-rfd4-proteina", "0.3.0")
EXPECTED_REPLICATES = 32
EXPECTED_TSR_K = 1.0
MIN_NUMBERING_CONFIDENCE = 0.1
MAX_PAIRED_LENGTH = 250
CANONICAL_CYSTEINE_POSITIONS = {"23", "104"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def _declared_design_paths(result: dict[str, Any] | None) -> list[str]:
    values = (((result or {}).get("result") or {}).get("designs") or [])
    return [value for value in values if isinstance(value, str)]


def _safe_design_path(design_root: Path, relative: str) -> Path | None:
    path = PurePosixPath(relative)
    if path.is_absolute() or ".." in path.parts:
        return None
    return design_root.joinpath(*path.parts)


def audit_arm_provenance(
    arm_root: Path,
    *,
    full_sequence_input: Path,
    arm_name: str,
) -> dict[str, Any]:
    """Check local copies of provider run manifests and artifact coverage."""
    design_root = arm_root / "design"
    config_path = design_root / "config.json"
    result_path = design_root / "result.json"
    config = read_json(config_path)
    result = read_json(result_path)
    result_provenance = (result or {}).get("provenance") or {}
    resolved_config = result_provenance.get("resolved_config") or {}
    provider_inputs = result_provenance.get("inputs") or []
    declared = _declared_design_paths(result)
    actual_paths = sorted(
        path.relative_to(design_root).as_posix()
        for path in design_root.glob("generation/inference_samples/*/*_generated.cif")
    )
    issues: list[str] = []

    config_hash = sha256_file(config_path) if config_path.is_file() else None
    result_hash = sha256_file(result_path) if result_path.is_file() else None
    input_actual_hash = sha256_file(full_sequence_input) if full_sequence_input.is_file() else None
    input_manifest_path = arm_root / "input.sha256"
    input_manifest_text = input_manifest_path.read_text(errors="replace") if input_manifest_path.is_file() else ""
    input_recorded_hash = input_manifest_text.split()[0] if input_manifest_text.split() else None
    parameters = parse_key_values(arm_root / "parameters.txt")
    image_tag = (arm_root / "image.txt").read_text(errors="replace").strip() if (arm_root / "image.txt").is_file() else ""
    image_id_text = (arm_root / "image-id.txt").read_text(errors="replace") if (arm_root / "image-id.txt").is_file() else ""
    image_digest_match = re.search(r"sha256:([0-9a-fA-F]{64})", image_id_text)
    image_digest = image_digest_match.group(1).lower() if image_digest_match else None
    checkpoint = parameters.get("tag")

    provenance_started = result_provenance.get("started_at")
    provenance_ended = result_provenance.get("ended_at")
    provider_completion = bool(
        result
        and result.get("task") == "design"
        and (result.get("result") or {}).get("artifacts") == "generation"
        and provenance_started
        and provenance_ended
        and len(declared) == EXPECTED_REPLICATES
    )
    if not provider_completion:
        issues.append("provider_completion_not_proven")

    if config is None:
        issues.append("config_missing_or_invalid")
    elif {key: value for key, value in config.items() if key != "schema_version"} != resolved_config:
        issues.append("resolved_config_mismatch")
    if config is not None:
        if config.get("num_replicates") != EXPECTED_REPLICATES:
            issues.append("configured_replicate_count_mismatch")
        try:
            tsr_k = float(config.get("sequence_tsr_k", math.nan))
        except (TypeError, ValueError):
            tsr_k = math.nan
        if tsr_k != EXPECTED_TSR_K:
            issues.append("sequence_tsr_k_mismatch")
        if config.get("seed") != EXPECTED_SEEDS[arm_name]:
            issues.append("config_seed_mismatch")
    if parameters.get("n") != str(EXPECTED_REPLICATES):
        issues.append("parameters_replicate_count_mismatch")
    if parameters.get("seed") != str(EXPECTED_SEEDS[arm_name]):
        issues.append("parameters_seed_mismatch")
    try:
        if float(parameters.get("sequence_tsr_k", "nan")) != EXPECTED_TSR_K:
            issues.append("parameters_sequence_tsr_k_mismatch")
    except ValueError:
        issues.append("parameters_sequence_tsr_k_mismatch")
    if config is not None and resolved_config and config.get("seed") != resolved_config.get("seed"):
        issues.append("resolved_seed_mismatch")

    if not input_actual_hash:
        issues.append("full_sequence_input_missing")
    if not input_recorded_hash or input_recorded_hash != input_actual_hash:
        issues.append("input_hash_mismatch")
    if not provider_inputs or not any(
        isinstance(item, dict) and item.get("path", "").endswith("full_sequence_design.json")
        for item in provider_inputs
    ):
        issues.append("provider_input_not_proven")

    if not image_digest:
        issues.append("image_digest_missing")
    if not image_tag or not checkpoint or not image_tag.endswith(checkpoint):
        issues.append("image_tag_checkpoint_mismatch")
    if checkpoint != EXPECTED_CHECKPOINTS[arm_name]:
        issues.append("checkpoint_tag_mismatch")
    module = result_provenance.get("module") or {}
    if (module.get("name"), module.get("version")) != EXPECTED_MODULE:
        issues.append("generation_module_mismatch")

    actual_set = set(actual_paths)
    declared_set = set(declared)
    safe_declared = all(_safe_design_path(design_root, item) is not None for item in declared)
    declared_files_exist = safe_declared and all(
        (_safe_design_path(design_root, item) or Path("/nonexistent")).is_file()
        for item in declared
    )
    coverage_exact = (
        safe_declared
        and len(declared) == EXPECTED_REPLICATES
        and len(declared_set) == EXPECTED_REPLICATES
        and len(actual_paths) == EXPECTED_REPLICATES
        and actual_set == declared_set
        and declared_files_exist
    )
    if not coverage_exact:
        issues.append("design_artifact_coverage_mismatch")
    return {
        "arm": arm_name,
        "provider_completion_status": "completed" if provider_completion and coverage_exact else "unproven",
        "provenance_status": "pass" if not issues else "fail",
        "provenance_issues": sorted(set(issues)),
        "provider_started_at": provenance_started,
        "provider_ended_at": provenance_ended,
        "generation_model": "/".join(
            str(value) for value in (module.get("name"), module.get("version")) if value
        ) or None,
        "checkpoint": checkpoint,
        "generation_seed": (config or {}).get("seed"),
        "sequence_tsr_k": (config or {}).get("sequence_tsr_k"),
        "requested_design_count": (config or {}).get("num_replicates"),
        "declared_design_count": len(declared),
        "actual_cif_count": len(actual_paths),
        "coverage_exact": coverage_exact,
        "full_sequence_input_sha256": input_actual_hash,
        "manifest_input_sha256": input_recorded_hash,
        "settings_sha256": config_hash,
        "provider_result_sha256": result_hash,
        "image_tag": image_tag or None,
        "image_sha256": image_digest,
        "declared_design_paths": sorted(declared),
        "actual_design_paths": actual_paths,
    }


def make_annotator() -> Annotator:
    return Annotator(chains=["H", "K", "L"], scheme="imgt", min_confidence=MIN_NUMBERING_CONFIDENCE)


def _number_reference_pair(
    sequences: dict[str, str], annotator: Annotator, label: str
) -> dict[str, Any]:
    numbered: dict[str, Any] = {}
    expected = {"H": {"H"}, "L": {"K", "L"}}
    for chain_id in ("H", "L"):
        if chain_id not in sequences:
            raise ValueError(f"{label}: missing {chain_id} reference chain")
        result = annotator.number(sequences[chain_id])
        if result.error:
            raise ValueError(f"{label} {chain_id}: numbering failed: {result.error}")
        if result.chain not in expected[chain_id]:
            raise ValueError(f"{label} {chain_id}: immunum mapped to {result.chain!r}")
        if result.confidence is None or result.confidence < MIN_NUMBERING_CONFIDENCE:
            raise ValueError(f"{label} {chain_id}: confidence below threshold")
        numbered[chain_id] = result
    return numbered


def load_references(reference_path: Path, annotator: Annotator) -> dict[str, Any]:
    trastuzumab_sequences = {
        "H": VH_PREFIX + WT_CDRH3 + VH_SUFFIX,
        "L": TRASTUZUMAB_VL,
    }
    trastuzumab = _number_reference_pair(trastuzumab_sequences, annotator, "trastuzumab")
    reference_doc = json.loads(reference_path.read_text())
    adalimumab_sequences: dict[str, str] = {}
    for entity in reference_doc.get("entities", []):
        description = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "").lower()
        sequence = entity.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        if "heavy chain" in description:
            adalimumab_sequences["H"] = sequence
        elif "light chain" in description:
            adalimumab_sequences["L"] = sequence
    adalimumab = _number_reference_pair(adalimumab_sequences, annotator, "adalimumab")
    segments: dict[str, dict[str, str]] = {"trastuzumab": {}, "adalimumab": {}}
    for name, sequences in (("trastuzumab", trastuzumab_sequences), ("adalimumab", adalimumab_sequences)):
        for chain_id, sequence in sequences.items():
            segment = annotator.segment(sequence)
            if segment.error or not isinstance(segment.cdr3, str) or not segment.cdr3:
                raise ValueError(f"{name} {chain_id}: immunum CDR3 segmentation failed")
            segments[name][chain_id] = segment.cdr3
    return {
        "sequences": {"trastuzumab": trastuzumab_sequences, "adalimumab": adalimumab_sequences},
        "numbered": {"trastuzumab": trastuzumab, "adalimumab": adalimumab},
        "cdr3": segments,
    }


def extract_structure_sequences(path: Path) -> dict[str, str]:
    """Read only explicit antibody H/L chains; never infer from target numbering."""
    atoms = get_structure(CIFFile.read(path), model=1)
    chain_ids = set(atoms.chain_id.tolist())
    sequences: dict[str, str] = {}
    for chain_id in ("H", "L"):
        if chain_id not in chain_ids:
            continue
        _, residue_names = get_residues(atoms[atoms.chain_id == chain_id])
        sequence = "".join(ProteinSequence.convert_letter_3to1(name) for name in residue_names)
        if sequence:
            sequences[chain_id] = sequence
    return sequences


def number_pair(sequences: dict[str, str], annotator: Annotator) -> dict[str, Any]:
    """Number explicit H and L chains independently and report mapping failures."""
    chain_results: dict[str, Any] = {}
    failures: list[str] = []
    expected = {"H": {"H"}, "L": {"K", "L"}}
    for chain_id in ("H", "L"):
        sequence = sequences.get(chain_id)
        if not sequence:
            failures.append(f"{chain_id}:missing_chain")
            chain_results[chain_id] = None
            continue
        try:
            result = annotator.number(sequence)
        except Exception as exc:  # retain the row if one chain cannot be numbered
            failures.append(f"{chain_id}:numbering_exception:{type(exc).__name__}")
            chain_results[chain_id] = None
            continue
        chain_results[chain_id] = result
        if result.error:
            failures.append(f"{chain_id}:numbering_error:{result.error}")
        elif result.chain not in expected[chain_id]:
            failures.append(f"{chain_id}:mapped_as_{result.chain}")
        elif result.confidence is None or result.confidence < MIN_NUMBERING_CONFIDENCE:
            failures.append(f"{chain_id}:confidence_below_0.1")
    return {
        "numbering_status": "pass" if not failures else "fail",
        "numbering_failure_reason": ";".join(failures),
        "numbered_chains": chain_results,
    }


def cdr3_identity(left: str | None, right: str | None) -> float | None:
    """Return direct positional identity only when IMGT CDR3 lengths agree."""
    if not left or not right or len(left) != len(right):
        return None
    return sum(a == b for a, b in zip(left, right)) / len(left)


def assess_sequence_gates(
    *,
    numbering_ok: bool,
    paired_length: int | None,
    canonical_cysteines_only: bool | None,
    trastuzumab_identity: float | None,
    adalimumab_identity: float | None,
    cdrh3_edit_distance: int | None,
) -> tuple[str, list[str]]:
    """Apply local sequence filters without claiming official database novelty."""
    failures: list[str] = []
    unknown: list[str] = []
    if paired_length is None:
        unknown.append("paired_length_unavailable")
    elif paired_length > MAX_PAIRED_LENGTH:
        failures.append("paired_length_over_250")
    if canonical_cysteines_only is False:
        failures.append("noncanonical_cysteine_count_or_position")
    elif canonical_cysteines_only is None:
        unknown.append("cysteine_numbering_unavailable")
    for name, identity in (("trastuzumab", trastuzumab_identity), ("adalimumab", adalimumab_identity)):
        if identity is None:
            unknown.append(f"{name}_identity_unavailable")
        elif identity >= 0.70:
            failures.append(f"{name}_paired_identity_not_below_70_percent")
    if cdrh3_edit_distance is None:
        unknown.append("trastuzumab_cdrh3_distance_unavailable")
    elif cdrh3_edit_distance < 2:
        failures.append("fewer_than_two_trastuzumab_cdrh3_changes")

    if failures:
        return "fail", failures + unknown
    if not numbering_ok or unknown:
        return "unknown", unknown or ["numbering_failed"]
    return "pass", []


def _cysteine_summary(sequence: str | None, numbered: Any) -> tuple[int | None, str | None, bool | None]:
    if sequence is None:
        return None, None, None
    count = sequence.count("C")
    if numbered is None or numbered.error or not numbered.numbering:
        return count, None, None
    positions = sorted(str(position) for position, residue in numbered.numbering.items() if residue == "C")
    canonical = count == 2 and set(positions) == CANONICAL_CYSTEINE_POSITIONS
    return count, ";".join(positions), canonical


def _identity_to_reference(candidate_numbered: dict[str, Any], reference_numbered: dict[str, Any]) -> float:
    counts = [
        numbered_identity(candidate_numbered[chain].numbering, reference_numbered[chain].numbering)
        for chain in ("H", "L")
    ]
    denominator = sum(total for _, total in counts)
    return sum(matches for matches, _ in counts) / denominator if denominator else math.nan


def _base_row(arm_name: str, design_name: str, structure_path: Path | None, provenance: dict[str, Any]) -> dict[str, Any]:
    return {
        "design_id": f"{arm_name}/{design_name}",
        "arm": arm_name,
        "design": design_name,
        "structure_path": str(structure_path) if structure_path else None,
        "structure_sha256": sha256_file(structure_path) if structure_path and structure_path.is_file() else None,
        "generation_model": provenance.get("generation_model"),
        "checkpoint": provenance.get("checkpoint"),
        "generation_seed": provenance.get("generation_seed"),
        "sequence_tsr_k": provenance.get("sequence_tsr_k"),
        "provider_completion_status": provenance.get("provider_completion_status"),
        "provenance_status": provenance.get("provenance_status"),
        "provenance_issues": ";".join(provenance.get("provenance_issues") or []),
        "full_sequence_input_sha256": provenance.get("full_sequence_input_sha256"),
        "manifest_input_sha256": provenance.get("manifest_input_sha256"),
        "settings_sha256": provenance.get("settings_sha256"),
        "image_sha256": provenance.get("image_sha256"),
        "provider_result_sha256": provenance.get("provider_result_sha256"),
        "numbering_status": "not_attempted",
        "numbering_failure_reason": "",
        "vh_sequence": None,
        "vl_sequence": None,
        "vh_length": None,
        "vl_length": None,
        "paired_length": None,
        "vh_numbered_chain_type": None,
        "vl_numbered_chain_type": None,
        "vh_confidence": None,
        "vl_confidence": None,
        "vh_cysteine_count": None,
        "vl_cysteine_count": None,
        "vh_cysteine_imgt_positions": None,
        "vl_cysteine_imgt_positions": None,
        "canonical_cysteines_only": None,
        "trastuzumab_paired_imgt_identity": None,
        "adalimumab_paired_imgt_identity": None,
        "trastuzumab_h_cdr3": None,
        "candidate_h_cdr3": None,
        "trastuzumab_h_cdr3_identity": None,
        "adalimumab_h_cdr3": None,
        "adalimumab_h_cdr3_identity": None,
        "trastuzumab_h_cdr3_edit_distance": None,
        "sequence_gate_status": "unknown",
        "sequence_gate_reasons": "",
        "official_adaptyv_novelty_status": "unverified_database_unavailable",
        "local_only_eligibility": "undetermined_structure",
        "is_duplicate_sequence": None,
        "duplicate_of": None,
        "nearest_eligible_full_paired_edit_distance": None,
    }


def audit_design(
    arm_name: str,
    design_name: str,
    structure_path: Path | None,
    provenance: dict[str, Any],
    annotator: Annotator,
    references: dict[str, Any],
) -> dict[str, Any]:
    row = _base_row(arm_name, design_name, structure_path, provenance)
    if structure_path is None or not structure_path.is_file():
        row["numbering_status"] = "fail"
        row["numbering_failure_reason"] = "structure_missing"
        row["sequence_gate_status"] = "unknown"
        row["sequence_gate_reasons"] = "structure_missing"
        row["local_only_eligibility"] = "undetermined_structure"
        return row
    try:
        sequences = extract_structure_sequences(structure_path)
    except Exception as exc:
        row["numbering_status"] = "fail"
        row["numbering_failure_reason"] = f"structure_parse_failed:{type(exc).__name__}:{exc}"
        row["sequence_gate_status"] = "unknown"
        row["sequence_gate_reasons"] = "structure_parse_failed"
        row["local_only_eligibility"] = "undetermined_structure"
        return row

    vh, vl = sequences.get("H"), sequences.get("L")
    row.update({
        "vh_sequence": vh,
        "vl_sequence": vl,
        "vh_length": len(vh) if vh else None,
        "vl_length": len(vl) if vl else None,
        "paired_length": len(vh) + len(vl) if vh and vl else None,
    })
    pair = number_pair(sequences, annotator)
    numbered = pair["numbered_chains"]
    row["numbering_status"] = pair["numbering_status"]
    row["numbering_failure_reason"] = pair["numbering_failure_reason"]
    for chain_id, prefix in (("H", "vh"), ("L", "vl")):
        result = numbered[chain_id]
        row[f"{prefix}_numbered_chain_type"] = result.chain if result is not None else None
        row[f"{prefix}_confidence"] = result.confidence if result is not None else None
        count, positions, canonical = _cysteine_summary(sequences.get(chain_id), result)
        row[f"{prefix}_cysteine_count"] = count
        row[f"{prefix}_cysteine_imgt_positions"] = positions
        row[f"{prefix}_canonical_cysteines_only"] = canonical
    cys_values = [row["vh_canonical_cysteines_only"], row["vl_canonical_cysteines_only"]]
    canonical_pair = None if any(value is None for value in cys_values) else all(cys_values)
    row["canonical_cysteines_only"] = canonical_pair

    identity_t, identity_a = None, None
    cdr3_h = None
    cdr3_edit = None
    cdr3_identity_t = None
    cdr3_identity_a = None
    cdr3_segment_error = None
    if pair["numbering_status"] == "pass":
        identity_t = _identity_to_reference(numbered, references["numbered"]["trastuzumab"])
        identity_a = _identity_to_reference(numbered, references["numbered"]["adalimumab"])
        try:
            h_segment = annotator.segment(vh)
            if not h_segment.error and isinstance(h_segment.cdr3, str) and h_segment.cdr3:
                cdr3_h = h_segment.cdr3
                tr_cdr3 = references["cdr3"]["trastuzumab"]["H"]
                ad_cdr3 = references["cdr3"]["adalimumab"]["H"]
                cdr3_edit = int(Levenshtein.distance(cdr3_h, tr_cdr3))
                cdr3_identity_t = cdr3_identity(cdr3_h, tr_cdr3)
                cdr3_identity_a = cdr3_identity(cdr3_h, ad_cdr3)
            else:
                cdr3_segment_error = "immunum_returned_no_h_cdr3"
        except Exception as exc:  # CDR3 metrics are optional; retain the numbered design row.
            cdr3_segment_error = f"immunum_segment_exception:{type(exc).__name__}"
    row.update({
        "trastuzumab_paired_imgt_identity": identity_t,
        "adalimumab_paired_imgt_identity": identity_a,
        "trastuzumab_h_cdr3": references["cdr3"]["trastuzumab"]["H"],
        "candidate_h_cdr3": cdr3_h,
        "trastuzumab_h_cdr3_identity": cdr3_identity_t,
        "adalimumab_h_cdr3": references["cdr3"]["adalimumab"]["H"],
        "adalimumab_h_cdr3_identity": cdr3_identity_a,
        "trastuzumab_h_cdr3_edit_distance": cdr3_edit,
    })
    status, reasons = assess_sequence_gates(
        numbering_ok=pair["numbering_status"] == "pass",
        paired_length=row["paired_length"],
        canonical_cysteines_only=canonical_pair,
        trastuzumab_identity=identity_t,
        adalimumab_identity=identity_a,
        cdrh3_edit_distance=cdr3_edit,
    )
    if cdr3_segment_error:
        reasons.append(cdr3_segment_error)
    row["sequence_gate_status"] = status
    row["sequence_gate_reasons"] = ";".join(reasons)
    if status == "fail":
        row["local_only_eligibility"] = "ineligible_sequence_gate"
    elif status == "unknown":
        row["local_only_eligibility"] = (
            "undetermined_numbering"
            if pair["numbering_status"] != "pass"
            else "undetermined_sequence_audit"
        )
    elif provenance.get("provenance_status") != "pass":
        row["local_only_eligibility"] = "undetermined_provenance"
    else:
        row["local_only_eligibility"] = "eligible_local_only"
    return row


def paired_edit_distance(left: tuple[str, str], right: tuple[str, str]) -> int:
    """Sum VH and VL Levenshtein distances, keeping the chain boundary explicit."""
    return int(Levenshtein.distance(left[0], right[0]) + Levenshtein.distance(left[1], right[1]))


def _sequence_key(row: dict[str, Any]) -> tuple[str, str] | None:
    vh, vl = row.get("vh_sequence"), row.get("vl_sequence")
    return (vh, vl) if isinstance(vh, str) and isinstance(vl, str) and vh and vl else None


def deduplicate_rows(rows: list[dict[str, Any]]) -> None:
    """Mark first occurrence across both arms while preserving every design row."""
    first_seen: dict[tuple[str, str], str] = {}
    for row in rows:
        key = _sequence_key(row)
        if key is None:
            row["is_duplicate_sequence"] = None
            row["duplicate_of"] = None
        elif key in first_seen:
            row["is_duplicate_sequence"] = True
            row["duplicate_of"] = first_seen[key]
        else:
            row["is_duplicate_sequence"] = False
            row["duplicate_of"] = None
            first_seen[key] = row["design_id"]


def _distance_summary(eligible_rows: list[dict[str, Any]]) -> dict[str, Any]:
    unique: dict[tuple[str, str], dict[str, Any]] = {}
    for row in eligible_rows:
        key = _sequence_key(row)
        if key is not None:
            unique.setdefault(key, row)
    keys = list(unique)
    distances: list[int] = []
    nearest: dict[tuple[str, str], int] = {}
    for i, left in enumerate(keys):
        for right in keys[i + 1 :]:
            distance = paired_edit_distance(left, right)
            distances.append(distance)
            nearest[left] = min(nearest.get(left, distance), distance)
            nearest[right] = min(nearest.get(right, distance), distance)
    for key, row in unique.items():
        row["nearest_eligible_full_paired_edit_distance"] = nearest.get(key)
    if not distances:
        return {
            "n_unique_eligible_sequences": len(keys),
            "pairwise_distance_method": "VH Levenshtein + VL Levenshtein",
            "n_pairwise_distances": 0,
            "minimum": None,
            "median": None,
            "mean": None,
            "maximum": None,
        }
    ordered = sorted(distances)
    middle = len(ordered) // 2
    median = ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2
    return {
        "n_unique_eligible_sequences": len(keys),
        "pairwise_distance_method": "VH Levenshtein + VL Levenshtein",
        "n_pairwise_distances": len(distances),
        "minimum": min(distances),
        "median": median,
        "mean": sum(distances) / len(distances),
        "maximum": max(distances),
    }


def _finite_or_none(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def run_audit(
    *,
    campaign_root: Path,
    full_sequence_input: Path,
    reference_path: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    annotator = make_annotator()
    references = load_references(reference_path, annotator)
    arm_provenance = {
        arm: audit_arm_provenance(
            campaign_root / arm,
            full_sequence_input=full_sequence_input,
            arm_name=arm,
        )
        for arm in ARM_NAMES
    }
    rows: list[dict[str, Any]] = []
    for arm in ARM_NAMES:
        provenance = arm_provenance[arm]
        arm_root = campaign_root / arm
        design_root = arm_root / "design"
        declared = provenance["declared_design_paths"]
        design_names: set[str] = set()
        for relative in declared:
            path = _safe_design_path(design_root, relative)
            if path is not None:
                design_names.add(path.name)
        for actual in provenance["actual_design_paths"]:
            design_names.add(Path(actual).name)
        # Emit declared names and discovered extra files, with explicit missing rows.
        paths_by_name = {
            Path(relative).name: _safe_design_path(design_root, relative)
            for relative in declared
            if _safe_design_path(design_root, relative) is not None
        }
        for relative in provenance["actual_design_paths"]:
            paths_by_name.setdefault(Path(relative).name, design_root / relative)
        for design_name in sorted(design_names):
            rows.append(audit_design(
                arm,
                design_name,
                paths_by_name.get(design_name),
                provenance,
                annotator,
                references,
            ))

    deduplicate_rows(rows)
    eligible = [row for row in rows if row["local_only_eligibility"] == "eligible_local_only"]
    eligible_keys = {_sequence_key(row) for row in eligible}
    eligible_keys.discard(None)
    unique_pairs = {_sequence_key(row) for row in rows}
    unique_pairs.discard(None)
    sequence_counts_by_arm = {
        arm: {
            "design_rows": sum(row["arm"] == arm for row in rows),
            "unique_full_paired_sequences": len({
                _sequence_key(row) for row in rows if row["arm"] == arm and _sequence_key(row) is not None
            }),
        }
        for arm in ARM_NAMES
    }
    eligibility_counts_by_arm = {
        arm: {
            status: sum(row["arm"] == arm and row["local_only_eligibility"] == status for row in rows)
            for status in (
                "eligible_local_only",
                "ineligible_sequence_gate",
                "undetermined_numbering",
                "undetermined_structure",
                "undetermined_sequence_audit",
                "undetermined_provenance",
            )
        }
        for arm in ARM_NAMES
    }
    diversity_by_arm = {
        arm: _distance_summary([row for row in eligible if row["arm"] == arm])
        for arm in ARM_NAMES
    }
    summary = {
        "audit_scope": "local-only sequence and provider provenance audit; no binding or thermal prediction",
        "official_adaptyv_novelty_status": "unverified_database_unavailable",
        "expected_designs_per_arm": EXPECTED_REPLICATES,
        "design_rows": len(rows),
        "unique_paired_sequences_across_arms": len(unique_pairs),
        "paired_sequence_counts_by_arm": sequence_counts_by_arm,
        "duplicate_design_rows": sum(row["is_duplicate_sequence"] is True for row in rows),
        "locally_eligible_design_rows_with_provenance": len(eligible),
        "unique_locally_eligible_molecules_after_cross_arm_deduplication": len(eligible_keys),
        "local_only_eligibility_counts_by_arm": eligibility_counts_by_arm,
        "local_only_eligibility_counts": {
            label: sum(row["local_only_eligibility"] == label for row in rows)
            for label in (
                "eligible_local_only",
                "ineligible_sequence_gate",
                "undetermined_numbering",
                "undetermined_structure",
                "undetermined_sequence_audit",
                "undetermined_provenance",
            )
        },
        "arms": [arm_provenance[arm] for arm in ARM_NAMES],
        "diversification_by_arm": diversity_by_arm,
        "diversification_of_unique_locally_eligible_molecules": _distance_summary(eligible),
        "interpretation_limits": [
            "The checkpoint arms differ in checkpoint image and seed, so the comparison is not a matched causal checkpoint experiment.",
            "RFD4 structural outputs and local sequence filters do not establish binding, trastuzumab epitope competition, thermal stability, or official Adaptyv novelty.",
            "No SPR score was transferred to these de novo designs.",
        ],
    }
    return rows, summary


def _csv_value(value: Any) -> Any:
    return _finite_or_none(value)


def write_outputs(rows: list[dict[str, Any]], summary: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "designs.csv"
    columns = list(rows[0]) if rows else []
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, default=Path("results/campaign/denovo-scale"))
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/campaign/campaign-20261114/denovo-input/full_sequence_design.json"),
    )
    parser.add_argument("--reference", type=Path, default=Path("configs/adalimumab_reference.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/denovo-scale-audit"))
    args = parser.parse_args()
    rows, summary = run_audit(
        campaign_root=args.campaign_root,
        full_sequence_input=args.input,
        reference_path=args.reference,
    )
    write_outputs(rows, summary, args.output_dir)
    print(json.dumps({
        "design_rows": summary["design_rows"],
        "unique_paired_sequences": summary["unique_paired_sequences_across_arms"],
        "locally_eligible_rows": summary["locally_eligible_design_rows_with_provenance"],
        "unique_eligible_molecules": summary["unique_locally_eligible_molecules_after_cross_arm_deduplication"],
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
