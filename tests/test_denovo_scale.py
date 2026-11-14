from types import SimpleNamespace

from scripts.audit_denovo_scale import (
    _distance_summary,
    assess_sequence_gates,
    cdr3_identity,
    deduplicate_rows,
    number_pair,
)


class FakeAnnotator:
    def __init__(self, results):
        self.results = results
        self.seen = []

    def number(self, sequence):
        self.seen.append(sequence)
        return self.results[sequence]


def numbered(chain, confidence=0.1, error=None):
    return SimpleNamespace(chain=chain, confidence=confidence, error=error, numbering={})


def test_numbers_explicit_heavy_and_light_chains_independently_at_boundary():
    annotator = FakeAnnotator({"heavy": numbered("H"), "light": numbered("K")})
    result = number_pair({"H": "heavy", "L": "light"}, annotator)

    assert result["numbering_status"] == "pass"
    assert annotator.seen == ["heavy", "light"]
    assert result["numbered_chains"]["H"].chain == "H"
    assert result["numbered_chains"]["L"].chain == "K"


def test_chain_mapping_failure_is_kept_as_numbering_failure():
    annotator = FakeAnnotator({"heavy": numbered("K"), "light": numbered("H")})
    result = number_pair({"H": "heavy", "L": "light"}, annotator)

    assert result["numbering_status"] == "fail"
    assert "H:mapped_as_K" in result["numbering_failure_reason"]
    assert "L:mapped_as_H" in result["numbering_failure_reason"]


def test_strict_filter_boundaries_are_inclusive_length_and_strict_identity():
    base = {
        "numbering_ok": True,
        "paired_length": 250,
        "canonical_cysteines_only": True,
        "trastuzumab_identity": 0.699,
        "adalimumab_identity": 0.699,
        "cdrh3_edit_distance": 2,
    }
    assert assess_sequence_gates(**base)[0] == "pass"
    assert "paired_length_over_250" in assess_sequence_gates(**{**base, "paired_length": 251})[1]
    assert "trastuzumab_paired_identity_not_below_70_percent" in assess_sequence_gates(
        **{**base, "trastuzumab_identity": 0.70}
    )[1]
    assert "adalimumab_paired_identity_not_below_70_percent" in assess_sequence_gates(
        **{**base, "adalimumab_identity": 0.70}
    )[1]
    assert "fewer_than_two_trastuzumab_cdrh3_changes" in assess_sequence_gates(
        **{**base, "cdrh3_edit_distance": 1}
    )[1]
    assert "noncanonical_cysteine_count_or_position" in assess_sequence_gates(
        **{**base, "canonical_cysteines_only": False}
    )[1]


def test_cdr3_identity_is_blank_when_lengths_do_not_align_positionally():
    assert cdr3_identity("ACD", "AED") == 2 / 3
    assert cdr3_identity("ACD", "ACDE") is None


def test_duplicate_full_pairs_across_arms_count_once_for_diversity():
    rows = [
        {"design_id": "arm-a/0", "vh_sequence": "ACD", "vl_sequence": "EFG"},
        {"design_id": "arm-b/0", "vh_sequence": "ACD", "vl_sequence": "EFG"},
        {"design_id": "arm-b/1", "vh_sequence": "ACE", "vl_sequence": "EFG"},
    ]

    deduplicate_rows(rows)
    summary = _distance_summary(rows)

    assert [row["is_duplicate_sequence"] for row in rows] == [False, True, False]
    assert rows[1]["duplicate_of"] == "arm-a/0"
    assert summary["n_unique_eligible_sequences"] == 2
    assert summary["n_pairwise_distances"] == 1
    assert summary["minimum"] == 1
