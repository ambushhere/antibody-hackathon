import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from src import her2_benchmark as benchmark


class FeatureTests(unittest.TestCase):
    def test_position_and_pairwise_columns(self):
        sequences = pd.Series([" acaaaaaaaa ", "CAAAAAAAAA"])
        position = benchmark.sequence_features(sequences, "position")
        pairwise = benchmark.sequence_features(sequences, "pairwise")
        self.assertEqual(position.shape, (2, 200))
        self.assertEqual(pairwise.shape, (2, 18_200))
        np.testing.assert_array_equal(position.getnnz(axis=1), [10, 10])
        np.testing.assert_array_equal(pairwise.getnnz(axis=1), [55, 55])
        np.testing.assert_array_equal(pairwise[:, :200].toarray(), position.toarray())
        self.assertEqual(position[0, 0], 1)
        self.assertEqual(position[0, 21], 1)
        self.assertEqual(pairwise[0, 201], 1)
        self.assertEqual(pairwise[1, 220], 1)
        self.assertEqual(pairwise[0, 220], 0)
        self.assertEqual(pairwise[0, 600], 1)
        self.assertEqual(pairwise[1, 620], 1)

    def test_rejects_unknown_residue_and_unsupported_length(self):
        for model in ("position", "pairwise", "dummy"):
            with self.subTest(model=model):
                with self.assertRaisesRegex(ValueError, "unknown residues"):
                    benchmark.sequence_features(pd.Series(["ACDEFGHIKX"]), model)
                with self.assertRaisesRegex(ValueError, "expected 10 residues"):
                    benchmark.sequence_features(pd.Series(["ACDEFGHIK"]), model)

    def test_tied_scores_do_not_use_label_order(self):
        labels = np.array([1, 0, 1, 0, 0])
        result = benchmark.top_metrics(labels, np.full(5, 0.4), 0.01)
        self.assertEqual(result, {"n": 1, "precision": 0.4, "enrichment": 1.0})
        result = benchmark.top_metrics(labels, np.array([1.0, 0.6, 0.6, 0.6, 0.2]), 0.4)
        self.assertAlmostEqual(result["precision"], 2 / 3)


def split_frame(sequences, centres=None):
    return pd.DataFrame({
        "seq": sequences,
        "label": [i % 2 for i in range(len(sequences))],
        "class": ["mid" if i % 2 == 0 else "high" for i in range(len(sequences))],
        "cluster_centre": sequences if centres is None else centres,
    })


class LeakageTests(unittest.TestCase):
    def test_rejects_sequence_and_clonotype_overlap(self):
        first = split_frame(["AAAAAAAAAA"], ["CCCCCCCCCC"])
        examples = (
            ("cdrh3", ["AAAAAAAAAA"], ["DDDDDDDDDD"]),
            ("cluster_centre", ["DDDDDDDDDD"], ["CCCCCCCCCC"]),
            ("left_centre_right_cdrh3", ["CCCCCCCCCC"], ["DDDDDDDDDD"]),
        )
        for conflicting_field, sequences, centres in examples:
            with self.subTest(field=conflicting_field):
                second = split_frame(sequences, centres)
                with self.assertRaisesRegex(ValueError, f"overlap.*'{conflicting_field}': 1"):
                    benchmark.audit_splits({"train": first, "val": second})

    def test_rejects_duplicate_cdrh3_within_split(self):
        repeated = split_frame(["AAAAAAAAAA", "AAAAAAAAAA"])
        with self.assertRaisesRegex(ValueError, "duplicate CDRH3"):
            benchmark.audit_splits({"train": repeated})

    def test_normalizes_sequences_before_overlap_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            split_frame(["AAAAAAAAAA", "CCCCCCCCCC"]).to_csv(path / "train.csv", index=False)
            split_frame([" aaaaaaaaaa ", "DDDDDDDDDD"]).to_csv(path / "val.csv", index=False)
            splits = {
                name: benchmark.load_split(path / f"{name}.csv")[0] for name in ("train", "val")
            }
        with self.assertRaisesRegex(ValueError, "split overlap"):
            benchmark.audit_splits(splits)

    def test_rejects_facs_mapping_conflict(self):
        split = split_frame(["AAAAAAAAAA", "CCCCCCCCCC"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "facs.csv"
            pd.DataFrame({
                "cdrh3": ["AAAAAAAAAA", "CCCCCCCCCC"],
                "binding_class": ["high", "medium"],
            }).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "2 conflicting"):
                benchmark.audit_facs(path, {"train": split})

    def test_spr_filters_before_aggregation(self):
        def measurement(sequence, kd, status="quantified_binding", light_chain=benchmark.WT_VL):
            return {
                "cdrh3": sequence,
                "vh_sequence": benchmark.VH_PREFIX + sequence + benchmark.VH_SUFFIX,
                "vl_sequence": light_chain,
                "kd_nM": kd,
                "binding_result": status,
            }

        kept = "ACDEFGHIKL"
        facs_test = "AAAAAAAAAA"
        rows = [
            measurement(kept, "1"),
            measurement(kept, "100"),
            measurement(kept, "10000", light_chain="CCCCCCCCCC"),
            measurement(kept, "10000", status="inconclusive"),
            measurement(facs_test, "4"),
            measurement("CCCCCCCCCC", "1"),
            measurement("DDDDDDDDDD", "1"),
            measurement(benchmark.WT_CDRH3, "1"),
            measurement("EEEEEEEEEE", "", status="no_reportable_kd"),
            measurement("FFFFFFFFFF", ""),
            measurement("GGGGGGGGGG", "0"),
            measurement("HHHHHHHHHH", "inf"),
            measurement("ACDEFGHIK", "1"),
        ]
        variants, audit = benchmark.prepare_spr(
            pd.DataFrame(rows),
            {"train": {"CCCCCCCCCC"}, "val": {"DDDDDDDDDD"}, "test": {facs_test}},
        )
        variants = variants.set_index("cdrh3")
        self.assertEqual(set(variants.index), {kept, facs_test})
        self.assertEqual(variants.loc[kept, "quantified_measurements"], 2)
        self.assertAlmostEqual(variants.loc[kept, "neg_log10_kd_nm"], -1.0)
        self.assertTrue(variants.loc[facs_test, "in_facs_test"])
        self.assertEqual(audit["quantified_rows_with_missing_or_invalid_kd"], 3)
        self.assertEqual(audit["excluded_non_trastuzumab_background_rows"], 1)
        self.assertEqual(audit["excluded_wild_type_rows"], 1)
        self.assertEqual(audit["excluded_unsupported_length_rows"], 1)
        self.assertEqual(audit["replicate_rows_collapsed"], 1)
        self.assertEqual(
            audit["quantified_variant_overlap_rows_by_facs_split"], {"train": 1, "val": 1, "test": 1}
        )


class CommandTests(unittest.TestCase):
    def test_archive_to_model_and_candidate_scores_with_validation_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            archive = path / "data.zip"
            splits = {
                "train": ["AAAAAAAAAA", "CAAAAAAAAA", "AAAAAAAACA", "CAAAAAAACA"],
                "val": ["AAAAAAAADA", "CAAAAAAADA"],
                "test": ["CAAAAAAAEA", "AAAAAAAAEA"],
            }
            with ZipFile(archive, "w") as output:
                for name, sequences in splits.items():
                    output.writestr(
                        f"{benchmark.SPLIT_PATH}{name}.csv",
                        split_frame(sequences).to_csv(index=False),
                    )
            output_dir = path / "results"
            script = Path(benchmark.__file__)
            subprocess.run([
                sys.executable, str(script), "benchmark", "--archive", str(archive),
                "--output-dir", str(output_dir),
            ], check=True, capture_output=True, text=True)
            results = json.loads((output_dir / "benchmark_metrics.json").read_text())
            # Position and pairwise tie on validation. Dummy scores higher on test.
            self.assertEqual(results["selected_model"], "position")
            self.assertEqual(results["classification"]["validation"]["position"]["roc_auc"], 1.0)
            self.assertEqual(results["classification"]["test"]["position"]["roc_auc"], 0.0)
            self.assertEqual(results["classification"]["test"]["dummy"]["roc_auc"], 0.5)
            with ZipFile(archive) as source:
                test_bytes = source.read(f"{benchmark.SPLIT_PATH}test.csv")
            self.assertEqual(results["inputs"]["test"]["sha256"], hashlib.sha256(test_bytes).hexdigest())

            candidates = path / "candidates.csv"
            scored = path / "candidate_predictions.csv"
            pd.DataFrame({"id": ["0001", "0002"], "loop": splits["test"]}).to_csv(
                candidates, index=False
            )
            subprocess.run([
                sys.executable, str(script), "predict",
                "--model", str(output_dir / results["outputs"]["model"]),
                "--candidates", str(candidates), "--sequence-column", "loop",
                "--output", str(scored),
            ], check=True, capture_output=True, text=True)
            predictions = pd.read_csv(scored, dtype={"id": str})
            benchmark_predictions = pd.read_csv(output_dir / results["outputs"]["facs_predictions"])
            self.assertEqual(predictions["id"].tolist(), ["0001", "0002"])
            np.testing.assert_allclose(
                predictions[benchmark.SCORE],
                benchmark_predictions.loc[benchmark_predictions["split"].eq("test"), benchmark.SCORE],
            )


if __name__ == "__main__":
    unittest.main()


def test_json_model_roundtrip(tmp_path):
    import numpy as np
    import pandas as pd
    from sklearn.linear_model import LogisticRegression
    seqs=pd.Series(['AAAAAAAAAA','WWWWWWWWWW','AAAAAAAAGA','WWWWWWWWGW'])
    estimator=LogisticRegression().fit(benchmark.sequence_features(seqs,'position'),[0,1,0,1])
    path=tmp_path/'model.json'
    benchmark.save_model({'features':'position','estimator':estimator},path)
    restored=benchmark.load_model(path)
    np.testing.assert_allclose(benchmark.predict_scores(estimator,'position',seqs),benchmark.predict_scores(restored['estimator'],'position',seqs),rtol=1e-6)
