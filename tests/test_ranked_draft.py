import csv
import json

import pytest

from scripts.prepare_ranked_draft import FIELDS, prepare_ranked_draft


def _candidate(sequence_id, track, input_order, *, novelty_class=None, n_high=0, distance=5):
    heavy = "A" * 40 + chr(ord("A") + (input_order % 20)) + "G" * 39
    light = "C" * 22 + chr(ord("A") + ((input_order * 3) % 20)) + "S" * 21
    return {
        "sequence_id": sequence_id,
        "kind": "paired",
        "heavy_sequence": heavy,
        "light_sequence": light,
        "design_description": f"original description for {sequence_id}",
        "track": track,
    }, {
        "sequence_id": sequence_id,
        "kind": "paired",
        "heavy_sequence": heavy,
        "light_sequence": light,
        "design_description": f"original description for {sequence_id}",
        "is_valid": True,
        "is_novel_enough": True,
        "novelty_class": novelty_class or ("novel" if track == "optimization" else "de_novo"),
        "n_high_severity": n_high,
        "cdrh3_min_distance": distance,
    }


def _write_inputs(tmp_path, *, denovo_count=22, mismatch_id="dn_005"):
    source_rows = []
    annotations = []
    for index in range(10):
        row, annotation = _candidate(f"opt_{index:02d}", "optimization", index, n_high=index % 3)
        source_rows.append(row)
        annotations.append(annotation)
    for index in range(denovo_count):
        sequence_id = f"dn_{index:03d}"
        row, annotation = _candidate(
            sequence_id,
            "denovo",
            20 + index,
            novelty_class="novel" if sequence_id == mismatch_id else "de_novo",
            n_high=index % 4,
            distance=(index * 7) % 13,
        )
        source_rows.append(row)
        annotations.append(annotation)

    input_path = tmp_path / "input.csv"
    with input_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(source_rows)
    # This extra row models the saved annotation for a rejected raw candidate.
    _, extra = _candidate("annotation_only_rejected", "denovo", 99, novelty_class="de_novo")
    annotation_path = tmp_path / "annotated.json"
    annotation_path.write_text(json.dumps(annotations + [extra]), encoding="utf-8")
    return input_path, annotation_path, source_rows


def test_builds_exact_interleaved_30_row_review_draft_and_selection_report(tmp_path):
    input_path, annotation_path, source_rows = _write_inputs(tmp_path)
    output = tmp_path / "team-3-ranked-draft30.csv"

    report = prepare_ranked_draft(input_path, annotation_path, output)

    with output.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        assert tuple(reader.fieldnames) == FIELDS
    assert len(rows) == report["row_count"] == 30
    assert report["counts_by_track"] == {"optimization": 10, "denovo": 20}
    assert [row["sequence_id"] for row in rows if row["track"] == "optimization"] == [
        row["sequence_id"] for row in source_rows if row["track"] == "optimization"
    ]
    assert [row["track"] for row in rows[:4]] == ["optimization"] * 3 + ["denovo"]
    assert all(row["kind"] == "paired" for row in rows)

    source_by_id = {row["sequence_id"]: row for row in source_rows}
    for rank, row in enumerate(rows, start=1):
        assert row["heavy_sequence"] == source_by_id[row["sequence_id"]]["heavy_sequence"]
        assert row["light_sequence"] == source_by_id[row["sequence_id"]]["light_sequence"]
        priority, basis, original = row["design_description"].split("; ", 2)
        expected_basis = (
            "selection_basis=optimization existing frozen FACS/SPR shortlist"
            if row["track"] == "optimization"
            else "selection_basis=denovo motif/diversity screen only"
        )
        assert priority == f"priority={rank:02d}"
        assert basis == expected_basis
        assert original == source_by_id[row["sequence_id"]]["design_description"]

    report_path = output.with_suffix(".selection.json")
    saved_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert saved_report == report
    assert report["status"] == "review_draft_not_submitted"
    assert report["no_lab_evidence"] == {
        "binding": "not_measured",
        "epitope_competition": "not_measured",
        "thermal_stability": "not_measured",
    }
    assert report["annotation_only_ids"] == ["annotation_only_rejected"]
    assert report["exclusion_counts"] == {
        "annotation_only_not_in_input_csv": 1,
        "ranked_out_after_top_20_de_novo": 1,
        "track_novelty_class_mismatch": 1,
    }
    mismatch = next(item for item in report["exclusions"] if item["reason"] == "track_novelty_class_mismatch")
    assert mismatch["sequence_id"] == "dn_005"
    assert {row["rank"] for row in report["design_rows"]} == set(range(1, 31))
    denovo_diagnostics = sorted(
        (row for row in report["design_rows"] if row["track"] == "denovo"),
        key=lambda row: row["track_rank"],
    )
    assert denovo_diagnostics[0]["min_pair_edit_distance_to_selected_denovo"] == 0
    assert denovo_diagnostics[0]["sequence_id"] == "dn_020"
    assert len(denovo_diagnostics) == 20


def test_selection_is_deterministic(tmp_path):
    input_path, annotation_path, _ = _write_inputs(tmp_path)
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"

    prepare_ranked_draft(input_path, annotation_path, first)
    prepare_ranked_draft(input_path, annotation_path, second)

    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("tamper", ["boolean", "sequence", "missing_annotation"])
def test_rejects_tampered_or_incomplete_annotation_before_writing(tmp_path, tamper):
    input_path, annotation_path, _ = _write_inputs(tmp_path)
    annotations = json.loads(annotation_path.read_text(encoding="utf-8"))
    if tamper == "boolean":
        annotations[0]["is_novel_enough"] = 1
    elif tamper == "sequence":
        annotations[0]["heavy_sequence"] += "A"
    else:
        annotations.pop(0)
    annotation_path.write_text(json.dumps(annotations), encoding="utf-8")
    output = tmp_path / "should-not-exist.csv"

    with pytest.raises(ValueError):
        prepare_ranked_draft(input_path, annotation_path, output)

    assert not output.exists()
    assert not output.with_suffix(".selection.json").exists()


def test_rejects_invalid_track_before_writing(tmp_path):
    input_path, annotation_path, source_rows = _write_inputs(tmp_path)
    source_rows[0]["track"] = "Optimization"
    with input_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(source_rows)
    output = tmp_path / "invalid-track.csv"

    with pytest.raises(ValueError, match="Invalid track"):
        prepare_ranked_draft(input_path, annotation_path, output)
    assert not output.exists()


def test_fails_when_class_mismatch_leaves_too_few_denovo_rows(tmp_path):
    input_path, annotation_path, _ = _write_inputs(tmp_path, denovo_count=20)
    output = tmp_path / "insufficient.csv"

    with pytest.raises(ValueError, match="Need at least 20 eligible de novo rows"):
        prepare_ranked_draft(input_path, annotation_path, output)
    assert not output.exists()


def test_refuses_to_overwrite_existing_output(tmp_path):
    input_path, annotation_path, _ = _write_inputs(tmp_path)
    output = tmp_path / "existing.csv"
    output.write_text("keep this file\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        prepare_ranked_draft(input_path, annotation_path, output)

    assert output.read_text(encoding="utf-8") == "keep this file\n"
    assert not output.with_suffix(".selection.json").exists()
