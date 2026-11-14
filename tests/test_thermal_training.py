import csv
from pathlib import Path

import pytest

from scripts.benchmark_thermal_training import (
    TABLES,
    WT_VH,
    WT_VL,
    Observation,
    duplicate_inventory,
    link_absci,
    parse_number,
    read_csv,
    read_table,
    similarity,
    target_reason,
)


def test_numeric_measurement_uses_original_column(tmp_path):
    stem = "shanehsazzadeh2023unlocking_DSF1"
    path = tmp_path / f"{stem}.csv"
    path.write_text(
        f"heavy,light,DSF - IP #1 (oC),fitness\n{WT_VH},{WT_VL},79.52\u200b,999\n",
        encoding="utf-8",
    )
    count, rows, invalid = read_table(path, TABLES[stem])
    assert count == 1 and invalid == []
    assert rows[0].value == 79.52
    assert rows[0].table.endpoint == "first_inflection"
    assert rows[0].target_reason == "known study lineage"
    assert parse_number("NaN") is None
    assert parse_number("Inf") is None
    assert parse_number(">57") is None


def test_family_holdout_uses_lineage_and_both_sequences():
    assert target_reason(("AAAA", "CCCC"), "anti-Grb2 from hu4D5")
    assert target_reason((WT_VH, WT_VL), "") == "exact WT pair"
    assert target_reason(("AAAA", WT_VL), "") == "conservative WT similarity flag"
    assert target_reason((WT_VH.replace("WGGDGFYAMD", "TGGDASYAMD"), "AAAA"), "")
    assert target_reason(("AAAA", "CCCC"), "") is None
    assert similarity(("AC", "FF"), ("ACC", "FF")) == 0.8


def test_duplicate_inventory_preserves_transitions_and_flags_cross_source():
    pair = ("AC", "DE")
    rows = [
        Observation("shanehsazzadeh2023unlocking_DSF1.csv", 2, pair, 71.0, "nanoDSF", "target"),
        Observation("shanehsazzadeh2023unlocking_DSF1.csv", 3, pair, 71.0, "nanoDSF", "target"),
        Observation("shanehsazzadeh2023unlocking_DSF2.csv", 4, pair, 81.0, "nanoDSF", "target"),
        Observation("jain2017biophysical_Tm.csv", 2, pair, 74.0, "dye DSF", None),
        Observation("jain2017biophysical_Tm.csv", 3, pair, 75.0, "dye DSF", None),
    ]
    duplicates, unique, conflicts, shared_pairs = duplicate_inventory(rows)
    assert len(duplicates) == 1 and duplicates[0]["cross_source"]
    assert shared_pairs == {pair}
    assert unique == [rows[0], rows[2]]
    assert conflicts == [[rows[3].ref, rows[4].ref]]


def test_author_join_rejects_label_disagreement(tmp_path):
    path = tmp_path / "absci.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Candidate", "HCDR3", "DLS - Cumulant PDI", "DSF - IP #1 (oC)", "DSF - IP #2 (oC)"])
        writer.writerow(["Trastuzumab", "SRWGGDGFYAMDY", 0.11, 79.52, 79.52])
    row = Observation(
        "shanehsazzadeh2023unlocking_DSF1.csv", 2, (WT_VH, WT_VL), 79.53, "nanoDSF", "target"
    )
    with pytest.raises(ValueError, match="does not match"):
        link_absci(path, [row])


@pytest.mark.parametrize("content", ["heavy,heavy,light\nAC,AC,QQ\n", "heavy,light\nAC\n"])
def test_malformed_csv_is_rejected(tmp_path: Path, content: str):
    path = tmp_path / "bad.csv"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError):
        list(read_csv(path, {"heavy", "light"}))
