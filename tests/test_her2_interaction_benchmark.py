import json
from hashlib import sha256

import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from src import her2_benchmark as facs
from src import her2_interaction_benchmark as interaction
from src import her2_spr_benchmark as spr


@pytest.mark.parametrize("kind", ["position", "pairwise"])
def test_sparse_dual_matches_dense_primal(kind):
    features = facs.sequence_features(pd.Series([
        "AAAAAAAAAA", "CCAAAAAAAA", "ADAAAAAAAA", "DEAAAAAAAA", "EEAAAAAAAA", "AAAAAAAACD",
    ]), kind).astype(float)
    # All omitted columns are zero across both train and test, so the small
    # reference fits the same primal problem.
    dense = features[:, np.unique(features.indices)].toarray()
    train, test = np.array([0, 2, 3, 4]), np.array([1, 5])
    y = np.array([-3., -2., 1., 0.])
    actual = interaction.ridge_predict(features, train, test, y, spr.ALPHAS)
    expected = [Ridge(alpha=a, solver="cholesky").fit(dense[train], y).predict(dense[test]) for a in spr.ALPHAS]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


def test_dense_preprocessing_is_fit_only_on_training_rows():
    features = np.array([[1., 10., 3.], [2., 8., 3.], [4., 4., 3.], [5., -9., 3.], [100., 100., -100.]])
    labels = np.array([1., 4., -5.])
    train, test = np.array([0, 1, 2]), np.array([3, 4])
    expected = make_pipeline(StandardScaler(), Ridge(alpha=1)).fit(features[train], labels).predict(features[test])
    actual = interaction.ridge_predict(features, train, test, labels, (1.,))[0]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
    features[4] = -1e8
    modified = interaction.ridge_predict(features, train, test, labels, (1.,))[0]
    assert actual[0] == pytest.approx(modified[0], abs=1e-12)


def test_nested_fit_uses_group_disjoint_rows_and_ignores_outer_test_labels(monkeypatch):
    groups = np.repeat(np.arange(8), [3, 2, 2, 3, 1, 1, 1, 2])
    labels = np.sin(np.arange(len(groups)))
    features = np.column_stack((groups, np.arange(len(groups)))).astype(float)
    folds = [(i, train, test, None) for i, (train, test) in enumerate(spr.group_splits(groups), start=1)]
    calls = []
    predict = interaction.ridge_predict

    def record(matrix, train, test, fit_labels, alphas):
        assert not set(groups[train]) & set(groups[test])
        calls.append((train, test))
        return predict(matrix, train, test, fit_labels, alphas)

    monkeypatch.setattr(interaction, "ridge_predict", record)
    scores, provenance = interaction.nested_oof({"dense": features}, labels, groups, folds)
    offset = 0
    for (_, outer_train, outer_test, _), fold in zip(folds, provenance):
        n_inner = len(fold["inner_folds"])
        inner = calls[offset:offset + n_inner]
        np.testing.assert_array_equal(np.sort(np.concatenate([val for _, val in inner])), outer_train)
        for (fit, val), record in zip(inner, fold["inner_folds"]):
            np.testing.assert_array_equal(np.sort(np.r_[fit, val]), outer_train)
            assert record["val_groups"] == np.unique(groups[val]).tolist()
        fit, test = calls[offset + n_inner]
        np.testing.assert_array_equal(fit, outer_train)
        np.testing.assert_array_equal(test, outer_test)
        offset += n_inner + 1
    assert offset == len(calls)

    held_out = folds[0][2]
    modified_labels = labels.copy()
    modified_labels[held_out] += 1000
    modified, modified_provenance = interaction.nested_oof({"dense": features}, modified_labels, groups, folds)
    np.testing.assert_array_equal(scores["dense"][held_out], modified["dense"][held_out])
    assert modified_provenance[0] == provenance[0]


@pytest.mark.parametrize(("contents", "message"), [
    ({"embedding": [[1.], [2.]]}, "require a sequence_ids"),
    ({"sequence_ids": ["AAA"], "embedding": [[1.]]}, "exact primary cohort order"),
    ({"sequence_ids": ["CCC", "AAA"], "embedding": [[1.], [2.]]}, "exact primary cohort order"),
    ({"sequence_ids": ["AAA", "AAA"], "embedding": [[1.], [2.]]}, "exact primary cohort order"),
    ({"sequence_ids": ["AAA", "CCC"], "embedding": [[1.], [np.nan]]}, "finite numbers"),
    ({"sequence_ids": ["AAA", "CCC"], "embedding": [[1.], [np.inf]]}, "finite numbers"),
    ({"sequence_ids": ["AAA", "CCC"], "embedding": [1., 2.]}, "real numeric matrix"),
])
def test_frozen_features_reject_invalid_archives(tmp_path, contents, message):
    npz, metadata = tmp_path / "features.npz", tmp_path / "metadata.json"
    np.savez(npz, **contents)
    metadata.write_text(json.dumps({"extractor": "test"}))
    with pytest.raises(ValueError, match=message):
        interaction.load_frozen_features(npz, metadata, ["AAA", "CCC"])


def test_frozen_features_read_named_arrays_and_record_hashes(tmp_path):
    npz, metadata = tmp_path / "features.npz", tmp_path / "metadata.json"
    values = np.array([[0.5, 2.], [1., 3.]])
    np.savez(npz, sequence_ids=["AAA", "CCC"], cdrh3_mean=values, contact_mean=values[:, :1])
    extraction = {"extractor": "test", "ordered_sequence_ids": ["AAA", "CCC"],
                  "output": {"sha256": sha256(npz.read_bytes()).hexdigest()}}
    metadata.write_text(json.dumps(extraction))
    features, provenance = interaction.load_frozen_features(npz, metadata, ["AAA", "CCC"])
    np.testing.assert_array_equal(features["cdrh3_mean"], values)
    np.testing.assert_array_equal(features["contact_mean"], values[:, :1])
    assert provenance["extraction"] == extraction
    assert provenance["npz"]["sha256"] == extraction["output"]["sha256"]
    assert len(provenance["metadata"]["sha256"]) == 64


@pytest.mark.parametrize(("extraction", "message"), [
    ({"ordered_sequence_ids": ["CCC", "AAA"]}, "metadata IDs"),
    ({"output": {"sha256": "incorrect"}}, "metadata hash"),
])
def test_frozen_features_reject_mismatched_metadata(tmp_path, extraction, message):
    npz, metadata = tmp_path / "features.npz", tmp_path / "metadata.json"
    np.savez(npz, sequence_ids=["AAA", "CCC"], embedding=[[1.], [2.]])
    metadata.write_text(json.dumps(extraction))
    with pytest.raises(ValueError, match=message):
        interaction.load_frozen_features(npz, metadata, ["AAA", "CCC"])


def test_paired_group_bootstrap_with_identical_predictions(monkeypatch):
    monkeypatch.setattr(spr, "BOOTSTRAP_DRAWS", 30)
    labels = np.arange(6.)
    groups = np.repeat([0, 1, 2], 2)
    result = interaction.summarize(labels, {"pair": labels, "position": labels.copy()}, groups, [("pair", "position")])
    for metric in ("spearman", "rmse"):
        assert result["paired_differences"]["pair minus position"][metric] == {
            "estimate": 0., "95_percent": [0., 0.], "valid_draws": 30,
        }


def test_cli_accepts_explicit_feature_contrast(monkeypatch):
    import sys
    captured = []
    monkeypatch.setattr(sys, "argv", ["benchmark", "--data-dir", "data", "--output", "out",
                                    "--feature-contrast", "contact", "mean"])
    monkeypatch.setattr(interaction, "run", lambda *args: captured.append(args) or {"subsets": {}})
    interaction.main()
    assert captured[0][-1] == [["contact", "mean"]]
