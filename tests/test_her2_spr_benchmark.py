import unittest

import numpy as np
import pandas as pd

from src import her2_benchmark as facs
from src import her2_spr_benchmark as spr


def measurement(sequence, kd, status="quantified_binding", *, light_chain=facs.WT_VL):
    return {
        "cdrh3": sequence, "vh_sequence": facs.VH_PREFIX + sequence + facs.VH_SUFFIX,
        "vl_sequence": light_chain, "kd_nM": kd, "binding_result": status,
    }


def variants_for_cv():
    sequences = [residue * 10 for residue in facs.AMINO_ACIDS[:8]]
    sequences += [sequence[:8] + "WY" for sequence in sequences]
    rng = np.random.default_rng(4)
    return pd.DataFrame({
        "cdrh3": sequences, facs.SCORE: rng.random(len(sequences)),
        spr.TARGET: rng.normal(size=len(sequences)),
    })


class AssayTests(unittest.TestCase):
    def test_aggregates_only_valid_replicates_and_keeps_missing_statuses(self):
        variant = "AAAAAAAAAA"
        rows = [
            measurement(variant, "0.1"),
            measurement(variant.lower(), "10"),
            measurement(variant, "0.001", light_chain="AAA"),
            measurement(variant, "0.001", "not_reported"),
            measurement(variant, "", "not_reported"),
            measurement(variant, "", "no_reportable_kd"),
            measurement(variant, "", "inconclusive"),
            measurement(variant, "0"),
            measurement(variant, "inf"),
            measurement(variant, ""),
            measurement(facs.WT_CDRH3, "0.15"),
            measurement(facs.WT_CDRH3, "0.18"),
            measurement("AAAAAAAAA", "0.4"),
        ]
        variants, audit = spr.prepare_spr(
            pd.DataFrame(rows), {"train": set(), "val": {variant}, "test": set()},
        )
        self.assertEqual(variants["cdrh3"].tolist(), [variant])
        self.assertEqual(variants[spr.TARGET].tolist(), [0.0])
        self.assertEqual(variants["quantified_measurements"].tolist(), [2])
        self.assertTrue(variants["in_facs_val"].item())
        self.assertEqual(audit["replicate_rows_collapsed"], 1)
        self.assertEqual(audit["excluded_other_background_rows"], 1)
        self.assertEqual(audit["excluded_other_length_rows_by_length"], {9: 1})
        self.assertEqual(audit["wt_reference_kd_nm"], [0.15, 0.18])
        counts = audit["kd_state_by_binding_result"]
        self.assertEqual(counts["quantified_binding"], {"missing": 1, "invalid": 2, "finite_positive": 6})
        self.assertEqual([counts[status]["missing"] for status in (
            "not_reported", "no_reportable_kd", "inconclusive",
        )], [1, 1, 1])
        self.assertEqual(counts["not_reported"]["finite_positive"], 1)


class CrossValidationTests(unittest.TestCase):
    def test_grouping_handles_indels_and_transitive_links(self):
        groups = spr.sequence_groups([
            "AAAAAAAAAA", "CCAAAAAAAA", "CCCCAAAAAA", "WWWWWWWWWW",
            "ACDEFGHIKL", "CDEFGHIKLA",
        ])
        self.assertEqual(len(set(groups[:3])), 1)
        self.assertNotIn(groups[3], set(groups[:3]))
        self.assertEqual(groups[4], groups[5])
        self.assertEqual(len(set(groups)), 3)

    def test_outer_and_inner_folds_keep_each_group_intact(self):
        oof, folds, issue = spr.cross_validate(variants_for_cv())
        self.assertIsNone(issue)
        self.assertEqual(len(folds), 5)
        self.assertEqual(oof["group"].nunique(), 8)
        self.assertTrue(oof.groupby("group")["fold"].nunique().eq(1).all())
        all_groups = set(oof["group"])
        for fold in folds:
            outer_test = set(fold["test_groups"])
            self.assertEqual(outer_test, set(oof.loc[oof["fold"].eq(fold["fold"]), "group"]))
            validated = set()
            for inner in fold["inner_folds"]:
                inner_val = set(inner["val_groups"])
                self.assertFalse(outer_test & inner_val)
                self.assertFalse(validated & inner_val)
                self.assertEqual(inner["train_n"] + inner["val_n"], fold["train_n"])
                validated |= inner_val
            self.assertEqual(validated, all_groups - outer_test)

    def test_outer_test_labels_do_not_affect_its_tuning_or_predictions(self):
        original = variants_for_cv()
        before, folds_before, _ = spr.cross_validate(original)
        held_out = before["fold"].eq(1)
        self.assertTrue(held_out.any())
        changed = original.copy()
        changed.loc[held_out, spr.TARGET] += 1000
        after, folds_after, _ = spr.cross_validate(changed)
        self.assertEqual(folds_before[0]["tuning"], folds_after[0]["tuning"])
        models = list(spr.MODEL_NAMES)
        pd.testing.assert_frame_equal(before.loc[held_out, models], after.loc[held_out, models])
        self.assertNotEqual(folds_before[1]["tuning"], folds_after[1]["tuning"])

    def test_two_groups_report_infeasible_nested_cv(self):
        variants = variants_for_cv().iloc[:2]
        oof, folds, issue = spr.cross_validate(variants)
        self.assertEqual(oof["group"].nunique(), 2)
        self.assertEqual(folds, [])
        self.assertIn("three nonempty sequence groups", issue)
        self.assertNotIn("position_ridge", oof)


if __name__ == "__main__":
    unittest.main()
