"""Check the frozen feature contract with the upstream tokenizer and synthetic states."""

import numpy as np
import pytest

from scripts import extract_ablang_features as extraction

VARIANT = "ARDYEQYAFD"
SECOND = "VGGAGFYAMD"


@pytest.fixture
def tokenizer():
    module = pytest.importorskip("ablang2.models.ablang2.tokenizers")
    return module.ABtokenizer()


def spr_row(sequence_id, heavy, light):
    return {"cdrh3": sequence_id, "vh_sequence": heavy, "vl_sequence": light}


def expected_heavy(sequence_id):
    return extraction.VH_PREFIX + sequence_id + extraction.VH_SUFFIX


def primary_row(sequence_id):
    return {"cdrh3": sequence_id, "cohort": extraction.PRIMARY}


def test_validates_single_distinct_pair_and_preserves_oof_order():
    wt = spr_row(extraction.WT_LOOP, extraction.WT_VH, extraction.WT_VL)
    first = spr_row(VARIANT, expected_heavy(VARIANT), extraction.WT_VL)
    second = spr_row(SECOND, expected_heavy(SECOND), extraction.WT_VL)
    oof = [primary_row(SECOND), {"cohort": "broader", "cdrh3": SECOND}, primary_row(VARIANT)]
    assert extraction.validate_pairs(oof, [wt, first, first, second], expected_count=2) == [
        SECOND, VARIANT
    ]
    with pytest.raises(ValueError, match="unique primary"):
        extraction.validate_pairs([primary_row(VARIANT)] * 2, [wt, first], expected_count=2)


@pytest.mark.parametrize("include_correct_pair", [False, True])
def test_rejects_swapped_or_ambiguous_spr_pair(include_correct_pair):
    rows = [
        spr_row(extraction.WT_LOOP, extraction.WT_VH, extraction.WT_VL),
        spr_row(VARIANT, extraction.WT_VL, expected_heavy(VARIANT)),
    ]
    if include_correct_pair:
        rows.append(spr_row(VARIANT, expected_heavy(VARIANT), extraction.WT_VL))
    with pytest.raises(ValueError, match="ambiguous SPR VH/VL"):
        extraction.validate_pairs([primary_row(VARIANT)], rows, expected_count=1)


def test_official_tokenizer_agrees_with_absolute_pair_and_window_offsets(tokenizer):
    tokens = extraction.tokenize(tokenizer, [VARIANT, SECOND])
    assert tuple(tokens.shape) == (2, 232)
    assert tokenizer.decode(tokens[0][:1]) == "<"
    assert tokenizer.decode(tokens[0][1:121]) == expected_heavy(VARIANT)
    assert tokenizer.decode(tokens[0][121:124]) == ">|<"
    assert tokenizer.decode(tokens[0][124:231]) == extraction.WT_VL
    assert tokenizer.decode(tokens[0][231:]) == ">"
    assert tokenizer.decode(tokens[0][98:110]) == "R" + VARIANT + "Y"
    assert tokenizer.decode(tokens[1][99:109]) == SECOND


@pytest.mark.parametrize("offset", [-1, 1])
def test_token_validation_rejects_shifted_window(monkeypatch, offset, tokenizer):
    monkeypatch.setattr(extraction, "H3_TOKEN_INDICES", tuple(range(99 + offset, 109 + offset)))
    with pytest.raises(ValueError, match="token indexing"):
        extraction.tokenize(tokenizer, [VARIANT])


def test_token_validation_rejects_swapped_chain_order(monkeypatch, tokenizer):
    monkeypatch.setattr(
        extraction, "format_pairs", lambda ids: [
            f"<{extraction.WT_VL}>|<{expected_heavy(sequence_id)}>" for sequence_id in ids
        ]
    )
    with pytest.raises(ValueError, match="token indexing"):
        extraction.tokenize(tokenizer, [VARIANT])


def test_pooling_uses_exact_token_positions_and_same_hidden_states_for_both_arms():
    hidden = np.broadcast_to(np.arange(232, dtype=np.float32)[None, :, None], (2, 232, 480))
    hidden = hidden.copy()
    hidden[1] += 10
    primary, contact = extraction.pool_embeddings(hidden, 2, [0, 4, 5, 6])
    np.testing.assert_array_equal(primary, np.repeat([[103.5], [113.5]], 480, axis=1))
    np.testing.assert_array_equal(contact, np.repeat([[102.75], [112.75]], 480, axis=1))
    assert primary.dtype == contact.dtype == np.float32


@pytest.mark.parametrize(
    "hidden, message",
    [
        (np.zeros((1, 231, 480), dtype=np.float32), "shape"),
        (np.zeros((1, 232, 479), dtype=np.float32), "shape"),
        (np.zeros((1, 232, 480), dtype=np.float64), "finite float32"),
        (np.full((1, 232, 480), np.nan, dtype=np.float32), "finite float32"),
        (np.full((1, 232, 480), np.inf, dtype=np.float32), "finite float32"),
    ],
)
def test_rejects_malformed_embeddings(hidden, message):
    with pytest.raises(ValueError, match=message):
        extraction.pool_embeddings(hidden, 1, [0, 4, 5, 6])


@pytest.mark.parametrize("indices", [[], [-1], [10], [0, 0]])
def test_rejects_empty_or_invalid_contact_mask(indices):
    hidden = np.zeros((1, 232, 480), dtype=np.float32)
    with pytest.raises(ValueError, match="contact mask"):
        extraction.pool_embeddings(hidden, 1, indices)


def test_output_round_trips_order_and_has_matching_digest(tmp_path):
    ids = [SECOND, VARIANT]
    primary = np.zeros((2, 480), dtype=np.float32)
    contact = np.ones_like(primary)
    metadata = {}
    output = tmp_path / "complete"
    extraction.write_output(output, ids, primary, contact, metadata)
    with np.load(output / "features.npz", allow_pickle=False) as archive:
        assert archive.files == ["sequence_ids", "cdrh3_mean", "wt_contact_mean"]
        assert archive["sequence_ids"].tolist() == ids
        np.testing.assert_array_equal(archive["cdrh3_mean"], primary)
        np.testing.assert_array_equal(archive["wt_contact_mean"], contact)
    assert metadata["output"]["sha256"] == extraction.digest((output / "features.npz").read_bytes())
