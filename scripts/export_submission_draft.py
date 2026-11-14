"""Export the research pool in the organizer's CSV schema without submitting it."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from immunum import Annotator

from scripts.prepare_candidate_screen import panel_bundle

FIELDS = ("sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description")
AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def validate_pair(heavy, light, annotator):
    if not heavy or not light or set(heavy + light) - AA:
        raise ValueError("Both variable domains must contain only the 20 standard amino acids")
    if len(heavy) + len(light) > 250:
        raise ValueError("Paired variable domains exceed the combined total length limit")
    checks = {}
    for role, sequence, allowed in (("heavy", heavy, {"H"}), ("light", light, {"K", "L"})):
        if sequence.count("C") % 2:
            raise ValueError(f"Odd cysteine count in {role} domain")
        numbered = annotator.number(sequence)
        if (numbered.error or numbered.chain not in allowed or numbered.confidence is None
                or not math.isfinite(numbered.confidence) or numbered.confidence < 0.1):
            raise ValueError(f"Incorrect or uncertain immunum {role} domain assignment")
        checks[role] = {"length": len(sequence), "cysteines": sequence.count("C"),
                        "numbered_chain": numbered.chain, "confidence": numbered.confidence}
    return checks


def export(output, selected_ids=None):
    bundle, _ = panel_bundle()
    source = json.loads(bundle["manifest.json"])
    candidates = [r for r in source["panel"] if r["track"] != "control"]
    if selected_ids is not None:
        unknown = set(selected_ids) - {r["query_id"] for r in candidates}
        if unknown or len(selected_ids) != len(set(selected_ids)):
            raise ValueError(f"Unknown or duplicate candidate IDs: {sorted(unknown)}")
        candidates = [r for r in candidates if r["query_id"] in selected_ids]
    if not candidates:
        raise ValueError("No research candidates selected")
    annotator = Annotator(chains=["H", "K", "L"], scheme="imgt", min_confidence=0.1)
    rows, evidence = [], []
    seen = set()
    for candidate in candidates:
        heavy, light = candidate["vh_sequence"], candidate["vl_sequence"]
        if (heavy, light) in seen:
            raise ValueError("Duplicate full paired sequence")
        seen.add((heavy, light))
        checks = validate_pair(heavy, light, annotator)
        description = (
            f"DRAFT research candidate; intended track {candidate['track']}; CDRH3 {candidate['cdrh3']}; "
            + ("trastuzumab framework, statistical proposal and SPR model ranking; "
               if candidate["track"] == "optimization" else
               f"RFD4-Proteina {candidate['generation']['checkpoint']}, paired VH/VL design; ")
            + "local sequence checks passed; official novelty pending; binding, epitope and thermal results unmeasured"
        )
        rows.append(dict(zip(FIELDS, (candidate["query_id"], "paired", heavy, light, description))))
        evidence.append({"sequence_id": candidate["query_id"], "track": candidate["track"],
                         "source": candidate["source"], "source_sha256": candidate["source_sha256"],
                         "paired_sequence_sha256": candidate["paired_sequence_sha256"],
                         "checks": checks, "official_abguard_status": "not_run_by_this_exporter",
                         "official_novelty_status": "pending", "laboratory_status": "not_measured"})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "status": "draft_research_pool_not_submitted_not_final_slot_selection",
        "file": output.name, "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "source_panel_sha256": hashlib.sha256(bundle["manifest.json"]).hexdigest(),
        "rows": len(rows), "kind": "paired",
        "counts_by_track": {t: sum(r["track"] == t for r in evidence) for t in ("optimization", "denovo")},
        "validation": "local format and immunum checks, not the official abguard gate",
        "official_track_assignment": "abguard determines the track; descriptions record intent only",
        "final_ranked_submission_limit": 20,
        "within_final_submission_row_limit": len(rows) <= 20,
        "team_filename": "team-3 is the VM team label; confirm the actual team name before manual submission",
        "designs": evidence,
    }
    output.with_suffix(".validation.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("output/submission-draft/team-3-draft.csv"))
    parser.add_argument("--ids", nargs="+", help="Explicit subset of frozen query IDs; default exports all 53 candidates")
    args = parser.parse_args()
    report = export(args.output, args.ids)
    print(json.dumps({k: report[k] for k in ("status", "rows", "counts_by_track", "sha256")}, indent=2))
