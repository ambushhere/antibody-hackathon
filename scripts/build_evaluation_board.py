"""Assemble candidate evidence without converting model confidence into lab success."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def source(path):
    return {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def build(output):
    ranking_path = ROOT / "results/spr-expanded-ranking/shortlist.csv"
    context_path = ROOT / "results/context-panel/scores.csv"
    replica_path = ROOT / "results/structural-replicates/scores.csv"
    initial_path = ROOT / "results/experiment-20261114/structure_comparison.csv"
    paths = [ranking_path, context_path, initial_path]
    if replica_path.exists():
        paths.append(replica_path)
    ranking = pd.read_csv(ranking_path)
    if ranking.cdrh3.duplicated().any():
        raise ValueError("Shortlist contains duplicate loop identifiers")
    frames = []
    for path in paths[1:]:
        frame = pd.read_csv(path)
        if path == replica_path:
            for name, expected in (("model", "opendde-abag"), ("target_context", "domain-iv")):
                if name in frame and not frame[name].eq(expected).all():
                    raise ValueError(f"Unexpected replica {name}")
                frame[name] = expected
        elif path == initial_path:
            frame["target_context"] = "initial_split"
            frame["wt_contact_count"] = float("nan")
            frame["geometry_audit_status"] = "confidence-only legacy table; footprint not included"
            frame["model"] = frame["model"].replace({"protenix_v2": "protenix-v2"})
        if frame[["query_id", "model", "target_context", "seed"]].isna().any().any():
            raise ValueError(f"Missing structural identifiers: {path}")
        frames.append(frame)
    structures = pd.concat(frames, ignore_index=True)
    required = {"query_id", "model", "target_context", "seed", "wt_contact_count"}
    if not required <= set(structures):
        raise ValueError("Structural audit schema missing required identifiers or geometry")
    key = ["query_id", "model", "target_context", "seed"]
    if structures.duplicated(key).any():
        raise ValueError("Structural audit has duplicate query/model/context/seed rows")
    controls = structures[structures.query_id.str.startswith("control_")].copy()
    candidates = []
    for _, row in ranking.sort_values("shortlist_pick").iterrows():
        # Match optimisation only in its fixed VH/VL scaffold. Never transfer this score to de novo.
        matched = structures[structures.query_id.eq("opt_" + row.cdrh3)]
        candidates.append({
            "candidate_id": "opt_" + row.cdrh3, "track": "optimization", "cdrh3": row.cdrh3,
            "frozen_shortlist_order": int(row.shortlist_pick),
            "local_sequence_gate": "pass" if str(row.strict_local_gate).lower() == "true" else "fail",
            "official_validation": "pending",
            "affinity_evidence": "retrospective SPR position-ridge hypothesis in fixed scaffold",
            "predicted_neg_log10_kd_nm": row.predicted_neg_log10_kd_nm,
            "spr_fit_sensitivity_min": row.spr_sensitivity_min,
            "spr_fit_sensitivity_max": row.spr_sensitivity_max,
            "structural_runs": len(matched),
            "structure_evidence": "control qualification pending or failed; inspect per-run evidence" if len(matched) else "not evaluated",
            "thermal_evidence": "no qualified candidate Tm prediction or measurement",
            "lab_status": "unmeasured", "decision": "provisional research priority, not submission approval",
        })
    for name in sorted(set(structures.loc[structures.query_id.str.startswith("denovo_"), "query_id"])):
        matched = structures[structures.query_id.eq(name)]
        candidates.append({"candidate_id": name, "track": "de_novo", "cdrh3": "",
                           "frozen_shortlist_order": None, "local_sequence_gate": "see design-specific validation",
                           "official_validation": "pending", "affinity_evidence": "fixed-scaffold SPR model not applicable",
                           "predicted_neg_log10_kd_nm": None, "spr_fit_sensitivity_min": None, "spr_fit_sensitivity_max": None,
                           "structural_runs": len(matched), "structure_evidence": "inspect controls and all replicas",
                           "thermal_evidence": "no qualified candidate Tm prediction or measurement",
                           "lab_status": "unmeasured", "decision": "exploratory design"})
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(candidates).to_csv(output / "candidates.csv", index=False)
    controls.to_csv(output / "controls.csv", index=False)
    structures.to_csv(output / "structures.csv", index=False)
    board = {"schema_version": 1, "candidate_count": len(candidates), "structural_rows": len(structures),
             "controls": len(controls), "sources": [source(path) for path in paths],
             "scoring_policy": "Separate evidence axes. No weighted universal score, inferred binding probability or automatic submission.",
             "interpretation": "Failure of a control removes support from that setup; it is not proof that every candidate is a nonbinder.",
             "neural_decision": "No affinity model promotion from FACS-only or synthetic transfer experiments.",
             "replica_audit_included": replica_path.exists(), "builder": source(Path(__file__).resolve())}
    (output / "manifest.json").write_text(json.dumps(board, indent=2) + "\n")
    lines = ["# Candidate evaluation board", "", board["scoring_policy"], "", board["interpretation"], "",
             "The frozen optimisation shortlist is unchanged. All candidates remain experimentally unmeasured.", "",
             "| Candidate | Track | Local gate | Structural runs | Status |",
             "| --- | --- | --- | ---: | --- |"]
    for row in candidates:
        lines.append(f"| {row['candidate_id']} | {row['track']} | {row['local_sequence_gate']} | {row['structural_runs']} | {row['decision']} |")
    lines.extend(["", "## WT control diagnostics", "",
                  "Contact recovery is a diagnostic, not full pose qualification. Read target distortion and interface RMSD alongside it.", "",
                  "| Model | Target context | Seed | WT contacts out of 19 |",
                  "| --- | --- | ---: | ---: |"])
    for _, row in controls[controls.query_id.eq("control_wt")].sort_values(["model", "target_context", "seed"]).iterrows():
        contacts = int(row.wt_contact_count) if pd.notna(row.wt_contact_count) else "not in this audit"
        lines.append(f"| {row.model} | {row.target_context} | {row.seed} | {contacts} |")
    lines.extend(["", "Full numeric geometry and confidence are in `structures.csv`, measured controls in `controls.csv` and separate candidate evidence fields in `candidates.csv`.",
                  "The file manifest records source hashes. Rebuild after a new audited model arm; raw provider outputs do not enter directly."])
    (output / "report.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: board[k] for k in ("candidate_count", "structural_rows", "replica_audit_included")}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/evaluation-board")
    build(parser.parse_args().output)


if __name__ == "__main__":
    main()
