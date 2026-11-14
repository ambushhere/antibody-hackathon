import csv
import json
import unittest
from itertools import combinations

import numpy as np
import pandas as pd
import pytest
from rapidfuzz.distance import Levenshtein

from scripts import rank_spr_candidates as rank
from src import her2_benchmark as facs
from src import her2_spr_benchmark as spr


def completed_benchmark():
    sequences = [residue * facs.LOOP_LENGTH for residue in facs.AMINO_ACIDS[:8]]
    variants = pd.DataFrame({
        "cdrh3": sequences, facs.SCORE: np.linspace(0.1, 0.9, len(sequences)),
        spr.TARGET: [-1.0, -1.2, -2.0, -0.3, 1.0, 0.2, 0.9, -0.5],
    })
    oof, folds, issue = spr.cross_validate(variants)
    assert issue is None
    primary = oof[["cdrh3", "group", spr.TARGET]].copy()
    cohort = {"n": len(primary), "groups": primary["group"].nunique(), "folds": folds}
    return primary, oof.assign(cohort=rank.COHORT), cohort


@pytest.fixture
def replay_files(tmp_path):
    primary, _, cohort = completed_benchmark()
    inputs = {}
    for name, sequences in {
        "train": ["  aaaaaaaacc  ", "YYYYYYYYYY"],
        "val": ["WWWWWWWWWW", "VVVVVVVVVV"],
        "test": ["AAAAAAAAAA", "TTTTTTTTTT"],
    }.items():
        path = tmp_path / f"{name}.csv"
        pd.DataFrame({
            "seq": sequences, "cluster_centre": sequences,
            "class": ["high", "low"], "label": [1, 0],
        }).to_csv(path, index=False)
        inputs[name] = rank.source(path)

    seqs = primary["cdrh3"]
    assay_path = tmp_path / "spr.csv"
    pd.DataFrame({
        "cdrh3": seqs, "vh_sequence": facs.VH_PREFIX + seqs + facs.VH_SUFFIX,
        "vl_sequence": facs.WT_VL,
        "kd_nM": np.power(10.0, -primary[spr.TARGET]), "binding_result": "quantified_binding",
    }).to_csv(assay_path, index=False)
    inputs["spr"] = rank.source(assay_path)
    model_path = tmp_path / "facs.json"
    model_path.write_text(json.dumps({
        "version": facs.MODEL_VERSION, "features": "position",
        "coefficients": [np.linspace(-0.2, 0.2, len(rank.FEATURE_ORDER)).tolist()],
        "intercept": [0.0],
    }))
    inputs["facs_model"] = rank.source(model_path)
    model = facs.load_model(model_path)
    primary[facs.SCORE] = facs.predict_scores(model["estimator"], model["features"], seqs)
    oof, folds, issue = spr.cross_validate(primary)
    oof_path = tmp_path / "oof.csv"
    oof.assign(cohort=rank.COHORT).to_csv(oof_path, index=False)
    cohort.update(folds=folds, split_issue=issue, metrics={
        score: {"spearman": spr.spearman(oof[spr.TARGET], oof[score])}
        for score in ("position_ridge", facs.SCORE)
    })
    metrics_path = tmp_path / "metrics.json"
    metrics_path.write_text(json.dumps({
        "seed": spr.SEED, "versions": {}, "alpha_grid": spr.ALPHAS,
        "max_folds": spr.MAX_FOLDS, "levenshtein_cutoff": spr.MAX_DISTANCE,
        "cohorts": {rank.COHORT: cohort}, "inputs": inputs,
    }))

    seqs = pd.Series(["MMMMMMMMMM", "AAAAAAAAMM"])
    candidates = pd.DataFrame({
        "cdrh3": seqs, "vh_sequence": facs.VH_PREFIX + seqs + facs.VH_SUFFIX,
        "vl_sequence": facs.WT_VL, "structural_status": "pending", "thermal_status": "pending",
        "strict_local_gate": True, "numbering": None, "cdrh3_by_scheme": None,
        facs.SCORE: facs.predict_scores(model["estimator"], model["features"], seqs),
        "origin": ["=1+1", "ordinary text"],
    })
    preflight_path = tmp_path / "preflight.json"
    old_path = tmp_path / "old.csv"
    preflight_path.write_text(json.dumps(candidates.to_dict("records")))
    candidates.head(1).to_csv(old_path, index=False)
    return {
        "data_dir": tmp_path, "preflight_path": preflight_path, "oof_path": oof_path,
        "metrics_path": metrics_path, "old_path": old_path, "output": tmp_path / "output",
    }


@pytest.mark.parametrize("column", [spr.TARGET, "position_ridge"])
@pytest.mark.parametrize("invalid", ["", "NaN", "inf"])
def test_run_rejects_incomplete_oof(replay_files, column, invalid):
    path = replay_files["oof_path"]
    oof = pd.read_csv(path, dtype=str)
    oof.loc[0, column] = invalid
    oof.to_csv(path, index=False)
    with pytest.raises(ValueError, match=f"OOF {column} must contain finite numbers"):
        rank.run(**replay_files)
    assert not replay_files["output"].exists()


@pytest.mark.parametrize(("column", "error"), [
    ("position_ridge", "Benchmark replay position_ridge differs"),
    (spr.TARGET, f"Benchmark replay {spr.TARGET} differs"),
])
def test_run_rejects_changed_oof(replay_files, column, error):
    path = replay_files["oof_path"]
    oof = pd.read_csv(path)
    oof.loc[0, column] += 0.01
    oof.to_csv(path, index=False)
    with pytest.raises(ValueError, match=error):
        rank.run(**replay_files)


@pytest.mark.parametrize("row_index", [0, 1])
def test_run_rejects_changed_preflight_scores(replay_files, row_index):
    path = replay_files["preflight_path"]
    candidates = json.loads(path.read_text())
    candidates[row_index][facs.SCORE] += 0.01
    path.write_text(json.dumps(candidates))
    with pytest.raises(ValueError, match="Preflight FACS scores differs"):
        rank.run(**replay_files)
    assert not replay_files["output"].exists()


@pytest.mark.parametrize("invalid", ["", "NaN", "inf", "0.9"])
def test_run_rejects_invalid_old_scores(replay_files, invalid):
    path = replay_files["old_path"]
    old = pd.read_csv(path, dtype=str)
    old.loc[0, facs.SCORE] = invalid
    old.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Old shortlist FACS scores"):
        rank.run(**replay_files)


def test_run_rejects_mismatched_facs_model(replay_files):
    model_path = replay_files["data_dir"] / "facs.json"
    model_path.write_text(model_path.read_text() + "\n")
    with pytest.raises(ValueError, match="FACS model does not match"):
        rank.run(**replay_files)


def test_run_uses_verified_model_scores_and_escapes_csv_text(replay_files):
    path = replay_files["preflight_path"]
    candidates = json.loads(path.read_text())
    verified_score = candidates[0][facs.SCORE]
    candidates[0][facs.SCORE] += rank.FACS_TOLERANCE / 2
    path.write_text(json.dumps(candidates))
    rank.run(**replay_files)
    for filename in ["candidates.csv", "shortlist.csv"]:
        with (replay_files["output"] / filename).open(newline="") as file:
            rows = {row["cdrh3"]: row for row in csv.DictReader(file)}
        row = rows[candidates[0]["cdrh3"]]
        assert row["origin"] == "'=1+1"
        assert float(row[facs.SCORE]) == verified_score
        assert float(row[rank.SPR_SCORE]) < 0


class FitTests(unittest.TestCase):
    def test_oof_predictions_and_labels_are_never_fit_inputs(self):
        primary, oof, cohort = completed_benchmark()
        folds, _ = rank.original_folds(primary, oof, cohort)
        final, outer, errors = rank.fit_models(primary, folds)
        altered = oof.copy()
        altered[list(spr.MODEL_NAMES) + [spr.TARGET]] = 500.0
        altered_folds, _ = rank.original_folds(primary, altered, cohort)
        altered_final, altered_outer, altered_errors = rank.fit_models(primary, altered_folds)
        self.assertEqual(errors, altered_errors)
        self.assertEqual(rank.export_parameters(final), rank.export_parameters(altered_final))
        for left, right in zip(outer, altered_outer, strict=True):
            self.assertEqual(rank.export_parameters(left), rank.export_parameters(right))

    def test_outer_refit_excludes_its_held_out_raw_assay_labels(self):
        primary, oof, cohort = completed_benchmark()
        folds, _ = rank.original_folds(primary, oof, cohort)
        final, outer, _ = rank.fit_models(primary, folds)
        held_out = folds[0][2]
        changed = primary.copy()
        changed.loc[held_out, spr.TARGET] += 500.0
        changed_final, changed_outer, _ = rank.fit_models(changed, folds)
        self.assertEqual(rank.export_parameters(outer[0]), rank.export_parameters(changed_outer[0]))
        self.assertEqual(outer[0].alpha, cohort["folds"][0]["tuning"]["position_ridge"]["alpha"])
        self.assertNotEqual(rank.export_parameters(final), rank.export_parameters(changed_final))
        self.assertNotEqual(rank.export_parameters(outer[1]), rank.export_parameters(changed_outer[1]))

    def test_json_parameter_predictions_match_ridge(self):
        primary, oof, cohort = completed_benchmark()
        folds, _ = rank.original_folds(primary, oof, cohort)
        final, outer, _ = rank.fit_models(primary, folds)
        sequences = pd.Series(["ACDEFGHIKL", "YWWWWWWWWW", "GCGCGCGCGC"])
        for model in [final, *outer]:
            parameters = json.loads(json.dumps({
                "feature_order": rank.FEATURE_ORDER, **rank.export_parameters(model),
            }, allow_nan=False))
            np.testing.assert_allclose(
                rank.predict_from_parameters(parameters, sequences),
                model.predict(rank.position_features(sequences)),
                rtol=0, atol=1e-12,
            )
            parameters["feature_order"].reverse()
            with self.assertRaisesRegex(ValueError, "feature ordering"):
                rank.predict_from_parameters(parameters, sequences)

    def test_mismatched_oof_fold_is_rejected(self):
        primary, oof, cohort = completed_benchmark()
        oof.loc[oof["fold"].eq(1), "fold"] = 2
        with self.assertRaisesRegex(ValueError, "Outer fold 1"):
            rank.original_folds(primary, oof, cohort)


class ShortlistTests(unittest.TestCase):
    def test_pareto_fronts_keep_ties_and_require_a_strict_improvement(self):
        scores = np.array([[1, 4], [3, 2], [1, 4], [1, 3], [0, 2], [0, 0]])
        np.testing.assert_array_equal(rank.pareto_fronts(scores), [1, 1, 1, 2, 3, 4])

    def test_shortlist_enforces_edit_distance_and_keeps_both_axis_extremes(self):
        candidates = pd.DataFrame({
            "cdrh3": ["AAAAAAAAAA", "AAAAAAAAAC", "CCCCCCCCCC", "GGGGGGGGGG"],
            rank.SPR_SCORE: [4, 3, 2, 1], facs.SCORE: [0.1, 0.9, 1.0, 0.8],
            "pareto_front": [1, 1, 1, 2],
        })
        picked = rank.exploratory_shortlist(candidates)["cdrh3"].tolist()
        self.assertEqual(picked, ["AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG"])
        self.assertTrue(all(Levenshtein.distance(a, b) >= 2 for a, b in combinations(picked, 2)))


if __name__ == "__main__":
    unittest.main()
