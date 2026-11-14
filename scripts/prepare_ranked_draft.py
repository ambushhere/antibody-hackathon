"""Prepare a ranked 30-row review draft without submitting designs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import tempfile
from collections import Counter
from pathlib import Path

from rapidfuzz.distance import Levenshtein


FIELDS = (
    "sequence_id",
    "kind",
    "heavy_sequence",
    "light_sequence",
    "design_description",
    "track",
)
ANNOTATED_FIELDS = FIELDS[:-1]
EXPECTED_CLASS = {"optimization": "novel", "denovo": "de_novo"}
OPTIMIZATION_COUNT = 10
DENOVO_COUNT = 20
OUTPUT_COUNT = OPTIMIZATION_COUNT + DENOVO_COUNT


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _read_inputs(input_csv: Path, annotation_json: Path):
    csv_bytes = input_csv.read_bytes()
    with io.TextIOWrapper(io.BytesIO(csv_bytes), encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError(f"Input CSV must have exactly these columns in order: {','.join(FIELDS)}")
        rows = list(reader)
    if not rows:
        raise ValueError("Input CSV is empty")

    input_ids = [row["sequence_id"] for row in rows]
    if any(type(value) is not str for row in rows for value in row.values()):
        raise ValueError("Input CSV contains a missing field or extra column value")
    if any(not isinstance(sequence_id, str) or not sequence_id for sequence_id in input_ids):
        raise ValueError("Input CSV contains an empty or invalid sequence_id")
    if len(input_ids) != len(set(input_ids)):
        raise ValueError("Duplicate sequence_id in input CSV")

    annotation_bytes = annotation_json.read_bytes()
    annotations = json.loads(annotation_bytes)
    if not isinstance(annotations, list):
        raise ValueError("Annotation JSON must contain a list of rows")
    annotation_ids = [row.get("sequence_id") if isinstance(row, dict) else None for row in annotations]
    if any(not isinstance(sequence_id, str) or not sequence_id for sequence_id in annotation_ids):
        raise ValueError("Annotation JSON contains a row without sequence_id")
    if len(annotation_ids) != len(set(annotation_ids)):
        raise ValueError("Duplicate sequence_id in annotations")
    annotation_by_id = dict(zip(annotation_ids, annotations))
    missing = sorted(set(input_ids) - set(annotation_by_id))
    if missing:
        raise ValueError(f"Annotation coverage is incomplete for input IDs: {missing}")

    return rows, annotation_by_id, csv_bytes, annotation_bytes, set(annotation_ids) - set(input_ids)


def _validate_evidence(row: dict, annotation: dict) -> None:
    for field in ANNOTATED_FIELDS:
        if field not in annotation or annotation[field] != row[field]:
            raise ValueError(f"Annotation input mismatch for {row['sequence_id']} field {field}")

    for field in ("is_valid", "is_novel_enough"):
        value = annotation.get(field)
        if type(value) is not bool:
            raise ValueError(f"Annotation {field} must be a boolean for {row['sequence_id']}")
        if value is not True:
            raise ValueError(f"Accepted input row is not marked {field}=true: {row['sequence_id']}")

    if not isinstance(annotation.get("novelty_class"), str):
        raise ValueError(f"Annotation novelty_class must be a string for {row['sequence_id']}")

    for field in ("n_high_severity", "cdrh3_min_distance"):
        value = annotation.get(field)
        if type(value) is not int or value < 0:
            raise ValueError(f"Annotation {field} must be a non-negative integer for {row['sequence_id']}")


def _denovo_pair_distance(candidate: dict, selected: list[dict]) -> int:
    if not selected:
        return 0
    row = candidate["row"]
    return min(
        Levenshtein.distance(row["heavy_sequence"], other["heavy_sequence"])
        + Levenshtein.distance(row["light_sequence"], other["light_sequence"])
        for other in selected
    )


def _rank_denovo(candidates: list[dict]) -> list[dict]:
    remaining = list(candidates)
    ranked = []
    selected = []
    while remaining:
        for candidate in remaining:
            candidate["min_pair_edit_distance_to_selected_denovo"] = _denovo_pair_distance(candidate, selected)
        remaining.sort(
            key=lambda candidate: (
                candidate["annotation"]["n_high_severity"],
                -candidate["min_pair_edit_distance_to_selected_denovo"],
                -candidate["annotation"]["cdrh3_min_distance"],
                candidate["row"]["sequence_id"],
            )
        )
        winner = remaining.pop(0)
        winner["track_rank"] = len(ranked) + 1
        ranked.append(winner)
        selected.append(winner["row"])
    return ranked


def _csv_bytes(rows: list[dict]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _publish_new_files(output_csv: Path, csv_bytes: bytes, selection_json: Path, report_bytes: bytes) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if output_csv.exists() or selection_json.exists():
        raise FileExistsError(f"Refusing to overwrite existing output or selection report: {output_csv}")

    staged = []
    published = []
    try:
        for data in (csv_bytes, report_bytes):
            with tempfile.NamedTemporaryFile(dir=output_csv.parent, prefix=".ranked-draft-", delete=False) as stream:
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
                staged.append(Path(stream.name))
        for temp_path, final_path in zip(staged, (output_csv, selection_json)):
            os.link(temp_path, final_path)
            published.append(final_path)
    except Exception:
        for path in published:
            path.unlink(missing_ok=True)
        raise
    finally:
        for path in staged:
            path.unlink(missing_ok=True)


def prepare_ranked_draft(
    input_csv: Path,
    annotation_json: Path,
    output_csv: Path,
) -> dict:
    """Write a six-column 30-row review draft and its selection diagnostics."""
    input_csv = Path(input_csv)
    annotation_json = Path(annotation_json)
    output_csv = Path(output_csv)
    selection_json = output_csv.with_suffix(".selection.json")
    if output_csv.exists() or selection_json.exists():
        raise FileExistsError(f"Refusing to overwrite existing output or selection report: {output_csv}")

    rows, annotations, csv_input_bytes, annotation_bytes, annotation_only_ids = _read_inputs(
        input_csv, annotation_json
    )
    for row in rows:
        if row["track"] not in EXPECTED_CLASS:
            raise ValueError(f"Invalid track {row['track']!r} for {row['sequence_id']}")
    declared_optimization_count = sum(row["track"] == "optimization" for row in rows)
    if declared_optimization_count != OPTIMIZATION_COUNT:
        raise ValueError(
            f"Expected exactly {OPTIMIZATION_COUNT} declared optimization rows, "
            f"found {declared_optimization_count}"
        )

    optimization = []
    denovo = []
    excluded = [
        {"sequence_id": sequence_id, "reason": "annotation_only_not_in_input_csv"}
        for sequence_id in sorted(annotation_only_ids)
    ]
    for input_order, row in enumerate(rows, start=1):
        track = row["track"]
        if not row["kind"] or not row["heavy_sequence"] or not row["light_sequence"]:
            raise ValueError(f"Missing kind or paired sequence for {row['sequence_id']}")
        annotation = annotations[row["sequence_id"]]
        _validate_evidence(row, annotation)
        candidate = {
            "row": row,
            "annotation": annotation,
            "input_order": input_order,
            "expected_novelty_class": EXPECTED_CLASS[track],
        }
        if annotation["novelty_class"] != EXPECTED_CLASS[track]:
            excluded.append({
                "sequence_id": row["sequence_id"],
                "track": track,
                "reason": "track_novelty_class_mismatch",
                "observed_novelty_class": annotation["novelty_class"],
                "expected_novelty_class": EXPECTED_CLASS[track],
                "n_high_severity": annotation["n_high_severity"],
                "cdrh3_min_distance": annotation["cdrh3_min_distance"],
            })
            continue
        target = optimization if track == "optimization" else denovo
        candidate["track_rank"] = len(target) + 1
        target.append(candidate)

    if len(optimization) != OPTIMIZATION_COUNT:
        raise ValueError(f"Expected exactly {OPTIMIZATION_COUNT} eligible optimization rows, found {len(optimization)}")

    ranked_denovo = _rank_denovo(denovo)
    if len(ranked_denovo) < DENOVO_COUNT:
        raise ValueError(
            f"Need at least {DENOVO_COUNT} eligible de novo rows after class/track exclusions, "
            f"found {len(ranked_denovo)}"
        )
    chosen_denovo = ranked_denovo[:DENOVO_COUNT]
    for candidate in ranked_denovo[DENOVO_COUNT:]:
        excluded.append({
            "sequence_id": candidate["row"]["sequence_id"],
            "track": "denovo",
            "reason": "ranked_out_after_top_20_de_novo",
            "de_novo_selection_rank": candidate["track_rank"],
            "n_high_severity": candidate["annotation"]["n_high_severity"],
            "min_pair_edit_distance_to_selected_denovo": candidate[
                "min_pair_edit_distance_to_selected_denovo"
            ],
            "cdrh3_min_distance": candidate["annotation"]["cdrh3_min_distance"],
        })

    ordered = []
    opt_index = 0
    denovo_index = 0
    while opt_index < len(optimization):
        ordered.extend(optimization[opt_index:opt_index + 3])
        opt_index += 3
        if denovo_index < len(chosen_denovo):
            ordered.append(chosen_denovo[denovo_index])
            denovo_index += 1
    ordered.extend(chosen_denovo[denovo_index:])

    output_rows = []
    design_rows = []
    for output_rank, candidate in enumerate(ordered, start=1):
        row = candidate["row"]
        annotation = candidate["annotation"]
        basis = (
            "optimization existing frozen FACS/SPR shortlist"
            if row["track"] == "optimization"
            else "denovo motif/diversity screen only"
        )
        output_row = dict(row)
        output_row["design_description"] = (
            f"priority={output_rank:02d}; selection_basis={basis}; {row['design_description']}"
        )
        output_rows.append(output_row)
        design_rows.append({
            "sequence_id": row["sequence_id"],
            "rank": output_rank,
            "track_rank": candidate["track_rank"],
            "input_order": candidate["input_order"],
            "track": row["track"],
            "novelty_class": annotation["novelty_class"],
            "selection_method": (
                "preserved_frozen_optimization_order"
                if row["track"] == "optimization"
                else "greedy_motif_and_pairwise_diversity"
            ),
            "n_high_severity": annotation["n_high_severity"],
            "min_pair_edit_distance_to_selected_denovo": candidate.get(
                "min_pair_edit_distance_to_selected_denovo"
            ),
            "cdrh3_min_distance": annotation["cdrh3_min_distance"],
            "selection_basis": basis,
        })

    output_csv_bytes = _csv_bytes(output_rows)
    exclusion_counts = dict(sorted(Counter(record["reason"] for record in excluded).items()))
    report = {
        "status": "review_draft_not_submitted",
        "file": output_csv.name,
        "sha256": _sha256(output_csv_bytes),
        "input_hashes": {
            str(input_csv): _sha256(csv_input_bytes),
            str(annotation_json): _sha256(annotation_bytes),
        },
        "policy": {
            "requested_rows": OUTPUT_COUNT,
            "minimum_rows": 28,
            "suggested_rows": 30,
            "track_counts": {"optimization": OPTIMIZATION_COUNT, "denovo": DENOVO_COUNT},
            "optimization": "preserve all eligible rows in original CSV order",
            "denovo": [
                "greedy selection over eligible de novo rows",
                "lowest n_high_severity motif occurrence count",
                "highest minimum pair edit distance to selected de novo rows, Lev VH plus Lev VL, zero when empty",
                "highest cdrh3_min_distance",
                "lexical sequence_id",
            ],
            "interleaving": "three optimization rows then one de novo row while optimization rows remain, then remaining de novo rows",
            "score_scope": "sequence warnings and diversity only; not a calibrated biological score",
            "description_prefix": "priority=NN; selection_basis=...; original description follows unchanged",
        },
        "row_count": len(output_rows),
        "input_row_count": len(rows),
        "eligible_row_count": len(optimization) + len(denovo),
        "eligible_counts_by_track": {
            "optimization": len(optimization),
            "denovo": len(denovo),
        },
        "counts_by_track": dict(Counter(row["track"] for row in output_rows)),
        "no_lab_evidence": {
            "binding": "not_measured",
            "epitope_competition": "not_measured",
            "thermal_stability": "not_measured",
        },
        "annotation_only_ids": sorted(annotation_only_ids),
        "exclusion_counts": exclusion_counts,
        "exclusions": excluded,
        "design_rows": design_rows,
    }
    report_bytes = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    _publish_new_files(output_csv, output_csv_bytes, selection_json, report_bytes)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("output/submission-draft/team-3-abguard-pass.csv"),
        help="Six-column CSV with explicit tracks",
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        default=Path("results/abguard-validation/raw/annotated.json"),
        help="Saved per-row abguard annotations",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/submission-draft/team-3-ranked-draft30.csv"),
    )
    args = parser.parse_args()
    report = prepare_ranked_draft(args.input, args.annotations, args.output)
    print(json.dumps({key: report[key] for key in ("status", "row_count", "counts_by_track", "sha256")}, indent=2))


if __name__ == "__main__":
    main()
