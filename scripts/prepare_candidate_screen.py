"""Freeze the matched HER2 candidate refold panel from the three audited sources."""

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

from scripts.audit_context_panel import IMAGE_ID, antibody_mapping
from scripts.prepare_her2_template_controls import (
    CONFIG as CONTROLS, PANEL as CONTROL_IDS, POSITIONS, ROOT, SOURCE_SHA256, SOURCE_URL,
    chain_sequences, frozen_context, json_text, require, sha256,
)

CONFIG = ROOT / "configs/her2-candidate-screen"
QUERY_FILE = "queries.json"
TEMPLATE = "complete-target.cif"
SEED, SAMPLE = 42, 1
ROLES = ("T", "H", "L")
AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
COUNTS = {"control": 4, "optimization": 10, "denovo": 43}
DENOVO = "results/denovo-scale-audit/designs.csv"
DENOVO_SUMMARY = "results/denovo-scale-audit/summary.json"
OPTIMIZATION = "results/spr-expanded-ranking/shortlist.csv"
ADMISSION = "docs/research/structure-admission-rule.md"
SOURCE_HASHES = {
    DENOVO: "dd5012da6b57418cc8b43e67aad638d241992829c9feb94e522d2d2f4d869d36",
    DENOVO_SUMMARY: "0fa9de02b1ec0c342734b66ccf6bde3fbf7b8bcd205dd0197d19ba81ad745310",
    OPTIMIZATION: "8b79fc635a3036c202806ac24acedd5ab976ebe0245aaaeb3d071f1448486fe6",
    "configs/her2-template-controls/queries.json":
        "481a0a4c4de707fb6769d0ab6dfd70ec0a8a3ba09e981824719f3611f3794140",
    "configs/her2-template-controls/manifest.json":
        "bb4269e094b0b1e8b5061656e2b86b94fddb78e38f3f6fb30014452bd1e6ace2",
    ADMISSION: "8c6f342357599f1d147b75e1f19860aa37c4bff8f1240c9c3d2705913a6aaa17",
}
ARM_COUNTS = {"denovo-scale-antibody-k1": 18, "denovo-scale-her2novel-k1": 25}
GENERATION_FIELDS = (
    "generation_model", "checkpoint", "generation_seed", "sequence_tsr_k",
    "provider_completion_status", "full_sequence_input_sha256", "manifest_input_sha256",
    "settings_sha256", "image_sha256", "provider_result_sha256",
)


def read_sources():
    payloads = {}
    for name, digest in SOURCE_HASHES.items():
        payload = (ROOT / name).read_bytes()
        require(sha256(payload) == digest, f"Frozen source SHA-256 changed: {name}")
        payloads[name] = payload
    return payloads


def source_record(source, row):
    return {"source": source, "source_sha256": SOURCE_HASHES[source],
            "source_record_sha256": sha256(json_text(row).encode())}


def candidate(query_id, track, row, source, **metadata):
    return {
        "query_id": query_id, "track": track, **source_record(source, row),
        "vh_sequence": row["vh_sequence"], "vl_sequence": row["vl_sequence"], **metadata,
    }


def denovo_candidates(rows, summary):
    arms = {arm["arm"]: arm for arm in summary["arms"]}
    require(len(summary["arms"]) == len(ARM_COUNTS) and set(arms) == set(ARM_COUNTS)
            and summary["design_rows"] == len(rows) == 64
            and summary["unique_paired_sequences_across_arms"] == 64
            and summary["unique_locally_eligible_molecules_after_cross_arm_deduplication"] == 43
            and len({(r["vh_sequence"], r["vl_sequence"]) for r in rows}) == 64
            and len({r["design_id"] for r in rows}) == 64
            and Counter(r["arm"] for r in rows) == dict.fromkeys(ARM_COUNTS, 32),
            "Unexpected coverage in the pre-existing 64-design VM batch")
    selected = []
    for row in rows:
        if row["local_only_eligibility"] != "eligible_local_only":
            continue
        match = re.fullmatch(r"her2_denovo_full_sequence_(\d+)_generated\.cif", row["design"])
        require(match is not None, f"Unexpected generated design name: {row['design']}")
        index = int(match[1])
        arm = arms[row["arm"]]
        require(0 <= index < 32 and row["design_id"] == f"{arm['arm']}/{row['design']}"
                and row["is_duplicate_sequence"] == "False" and not row["duplicate_of"]
                and row["provenance_status"] == "pass" and not row["provenance_issues"]
                and row["numbering_status"] == row["sequence_gate_status"] == "pass"
                and row["official_adaptyv_novelty_status"] == "unverified_database_unavailable"
                and row["canonical_cysteines_only"] == "True"
                and row["provider_completion_status"] == "completed"
                and row["paired_length"] == str(len(row["vh_sequence"]) + len(row["vl_sequence"]))
                and all(row[k] == str(arm[k]) for k in GENERATION_FIELDS)
                and re.fullmatch(r"[0-9a-f]{64}", row["structure_sha256"]) is not None,
                f"Unexpected eligibility or provenance: {row['design_id']}")
        selected.append(candidate(
            f"{arm['arm']}_{index}", "denovo", row, DENOVO,
            cdrh3=row["candidate_h_cdr3"], local_eligibility=row["local_only_eligibility"],
            official_adaptyv_novelty=row["official_adaptyv_novelty_status"],
            generation={
                "origin": "pre-existing team VM RFD4 scale runs",
                "arm": arm["arm"], "index": index, "design_id": row["design_id"],
                "structure_path": row["structure_path"], "structure_sha256": row["structure_sha256"],
                "provider_started_at": arm["provider_started_at"],
                "provider_ended_at": arm["provider_ended_at"],
                **{k: arm[k] for k in GENERATION_FIELDS},
            },
        ))
    require(Counter(row["generation"]["arm"] for row in selected) == ARM_COUNTS,
            "Exact locally eligible de novo coverage required")
    return sorted(selected, key=lambda row: (row["generation"]["arm"], row["generation"]["index"]))


def optimization_candidates(rows):
    require([int(row["shortlist_pick"]) for row in rows] == list(range(1, 11)),
            "Exact frozen optimization shortlist coverage required")
    candidates = []
    for row in rows:
        loop = row["cdrh3"]
        require(len(loop) == 10 and set(loop) <= AA and row["strict_local_gate"] == "True"
                and row["official_validation"] == "pending"
                and row["paired_length"] == str(len(row["vh_sequence"]) + len(row["vl_sequence"])),
                f"Unexpected shortlist sequence or gate status: {loop}")
        candidates.append(candidate(
            f"opt_{loop}", "optimization", row, OPTIMIZATION, cdrh3=loop,
            shortlist_pick=int(row["shortlist_pick"]), origin=row["origin"],
            parent_cdrh3=row["parent_cdrh3"], shortlist_seed=int(row["seed"]),
            local_eligibility="eligible_local_only", official_validation=row["official_validation"],
        ))
    return candidates


def validate_panel(panel, queries, canonical):
    require(Counter(row["track"] for row in panel) == COUNTS, "Expected exactly 4 + 10 + 43 queries")
    ids = [row["query_id"] for row in panel]
    require(len(ids) == len(set(ids)) and ids == [q["query_id"] for q in queries],
            "Duplicate panel ID or mismatched query coverage")
    pairs = {(row["vh_sequence"], row["vl_sequence"]) for row in panel}
    require(len(pairs) == len(panel), "Duplicate full paired antibody sequence")
    indexed = {row["query_id"]: row for row in panel}
    for row, query in zip(panel, queries, strict=True):
        query_id = row["query_id"]
        require(re.fullmatch(r"[A-Za-z0-9_-]+", query_id) is not None, f"Unsafe query ID: {query_id}")
        seqs = chain_sequences(query)
        require(tuple(seqs) == ROLES and len(query["chains"]) == len(ROLES)
                and all(chain["msa_path"] is None for chain in query["chains"])
                and seqs["T"] == canonical[POSITIONS.start - 1:POSITIONS.stop - 1]
                and (seqs["H"], seqs["L"]) == (row["vh_sequence"], row["vl_sequence"])
                and all(seq and set(seq) <= AA for seq in seqs.values())
                and len(seqs["H"]) + len(seqs["L"]) <= 250,
                f"Invalid sequence, role, MSA or paired length: {query_id}")
        if row["track"] != "denovo":
            require(antibody_mapping(query_id, indexed), f"Missing WT scaffold mapping: {query_id}")


def panel_bundle():
    payloads = read_sources()
    canonical, controls, expected, context_sources = frozen_context()
    control_path = str((CONTROLS / QUERY_FILE).relative_to(ROOT))
    control_queries = json.loads(payloads[control_path])
    control_manifest = json.loads(payloads[str((CONTROLS / "manifest.json").relative_to(ROOT))])
    require(control_queries == expected and control_manifest["panel"] == controls
            and tuple(row["query_id"] for row in controls) == CONTROL_IDS,
            "Control queries or metadata changed")
    template = (CONTROLS / TEMPLATE).read_bytes()
    require(sha256(template) == control_manifest["artifact_sha256"][TEMPLATE]
            and control_manifest["arms"]["complete"] == {
                "template": TEMPLATE, "query": QUERY_FILE, "masked_canonical_positions": [],
                "coordinate_residues": len(POSITIONS), "atom_count": 912,
            },
            "Complete target-only template disagrees with the control manifest")
    control_panel = [
        {**row, "track": "control", **source_record(control_path, query)}
        for row, query in zip(controls, control_queries["queries"], strict=True)
    ]
    optimization = optimization_candidates(list(csv.DictReader(payloads[OPTIMIZATION].decode().splitlines())))
    denovo = denovo_candidates(list(csv.DictReader(payloads[DENOVO].decode().splitlines())),
                               json.loads(payloads[DENOVO_SUMMARY]))
    panel = control_panel + optimization + denovo
    queries = control_queries["queries"] + [
        {"query_id": row["query_id"], "chains": [
            {"chain_id": role, "sequence": sequence, "msa_path": None}
            for role, sequence in zip(ROLES, (canonical[510:629], row["vh_sequence"], row["vl_sequence"]), strict=True)
        ]}
        for row in optimization + denovo
    ]
    validate_panel(panel, queries, canonical)
    for row in panel:
        row["vh_sha256"] = sha256(row["vh_sequence"].encode())
        row["vl_sha256"] = sha256(row["vl_sequence"].encode())
        row["paired_sequence_sha256"] = sha256(f"{row['vh_sequence']}\n{row['vl_sequence']}".encode())
        row["paired_length"] = len(row["vh_sequence"]) + len(row["vl_sequence"])
    artifacts = {QUERY_FILE: json_text({"queries": queries}).encode(), TEMPLATE: template}
    manifest = {
        "schema_version": 1,
        "source_sha256": {
            **SOURCE_HASHES, **context_sources, SOURCE_URL: SOURCE_SHA256,
            str((CONTROLS / TEMPLATE).relative_to(ROOT)): sha256(template),
            control_manifest["reference"]["url"]: control_manifest["reference"]["sha256"],
        },
        "target": control_manifest["target"],
        "reference": control_manifest["reference"],
        "data_origin": control_manifest["data_origin"],
        "panel_counts": COUNTS,
        "pre_existing_generation_count": 64,
        "panel": panel,
        "run": {
            "model": "protenix-v2", "image_id": IMAGE_ID, "msas": "none", "template": TEMPLATE,
            "seed": SEED, "diffusion_samples": 1, "sample": SAMPLE, "expected_poses": len(panel),
            "admission_rule": ADMISSION, "required_control_arm": "complete_no_msa",
        },
        "artifact_sha256": {name: sha256(payload) for name, payload in artifacts.items()},
    }
    return {**artifacts, "manifest.json": json_text(manifest).encode()}, canonical


def prepare(output):
    artifacts, _ = panel_bundle()
    output.mkdir(parents=True, exist_ok=True)
    for name, payload in artifacts.items():
        path = output / name
        require(not path.exists() or path.read_bytes() == payload, f"Refusing to replace frozen input {path}")
    for name, payload in artifacts.items():
        path = output / name
        if not path.exists():
            path.write_bytes(payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CONFIG)
    args = parser.parse_args()
    prepare(args.output)


if __name__ == "__main__":
    main()
