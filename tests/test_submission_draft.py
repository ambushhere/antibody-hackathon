import csv
from types import SimpleNamespace

import pytest

from scripts.export_submission_draft import FIELDS, export, validate_pair


class Numberer:
    def number(self, sequence):
        return SimpleNamespace(error=None, chain="H" if sequence.startswith("H") else "K", confidence=0.1)


def test_pair_rejects_invalid_alphabet_odd_cysteine_and_long_construct():
    assert validate_pair("HCC", "LCC", Numberer())["heavy"]["confidence"] == 0.1
    for heavy, light in (("HCX", "LCC"), ("HC", "LCC"), ("H" * 249, "LCC")):
        with pytest.raises(ValueError):
            validate_pair(heavy, light, Numberer())
    with pytest.raises(ValueError, match="heavy domain"):
        validate_pair("LCC", "LCC", Numberer())


def test_export_exact_schema_and_source_sequences_without_controls(tmp_path):
    output = tmp_path / "draft.csv"
    report = export(output)
    with output.open() as stream:
        reader = csv.DictReader(stream)
        assert tuple(reader.fieldnames) == FIELDS
        rows = list(reader)
    assert len(rows) == report["rows"] == 53
    assert report["counts_by_track"] == {"optimization": 10, "denovo": 43}
    assert len({r["sequence_id"] for r in rows}) == 53
    assert all(r["kind"] == "paired" and r["light_sequence"] for r in rows)
    assert all(not r["sequence_id"].startswith("control_") for r in rows)
    assert all(r["design_description"].startswith("DRAFT") for r in rows)
    assert all(r["official_abguard_status"] == "not_run_by_this_exporter" for r in report["designs"])


def test_unknown_or_duplicate_subset_is_rejected_before_writing(tmp_path):
    output = tmp_path / "draft.csv"
    with pytest.raises(ValueError, match="Unknown or duplicate"):
        export(output, ["does_not_exist"])
    assert not output.exists()
