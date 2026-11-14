import csv
import itertools
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from src import her2_candidates as her2


PARENT = "AGGAGFYAAD"
OTHER_PARENT = "VSGNGFYASD"
HELD_OUT = "TTGSGFYTMD"


def write_table(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(fields)
        writer.writerows(rows)


def edit_distance(left, right):
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        current = [i]
        for j, b in enumerate(right, 1):
            current.append(min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (a != b)))
        previous = current
    return previous[-1]


class CandidateTests(unittest.TestCase):
    def test_exact_novelty_at_all_relevant_lengths(self):
        query = "AGGAGFYAAD"
        for known in (query, "AGGAGFYAAC", "AGAGFYAAD", "AGGAGCFYAAD"):
            with self.subTest(known=known):
                self.assertTrue(her2.NoveltyIndex([known]).has_neighbor(query))
        for known in ("AGAGFYAD", "ACGAGFYAAC", "AGGAGCCFYAAD", "GAGGAGFYAADAC"):
            with self.subTest(known=known):
                self.assertFalse(her2.NoveltyIndex([known]).has_neighbor(query))

        # Unpositioned deletion signatures would wrongly reject this two-edit pair.
        self.assertFalse(her2.NoveltyIndex(["CAGAGFYAAD"]).has_neighbor(query))
        sequences = ["".join(aa) for size in range(1, 5) for aa in itertools.product("AC", repeat=size)]
        for known in sequences:
            index = her2.NoveltyIndex([known])
            for sequence in sequences:
                self.assertEqual(index.has_neighbor(sequence), edit_distance(sequence, known) < 2)

    def test_screen_rejects_duplicates_invalid_cysteine_and_wt_neighbors(self):
        known = her2.NoveltyIndex(["AGAGFYAAD"])
        seen = {OTHER_PARENT}
        for rejected in (
            PARENT, OTHER_PARENT, "AAGAGFYAAX", "aAGAGFYAAD",
            "CAGAGFYAAD", "AAGAGFYAA", "AAGDGFYAMD", her2.WT_CDRH3,
        ):
            with self.subTest(sequence=rejected):
                self.assertIsNone(her2.screen_candidate(rejected, known, seen))
        self.assertEqual(her2.screen_candidate("AAADGFYAMD", known, seen), 3)

    def test_liability_motifs_are_flags(self):
        self.assertEqual(her2.liability_motifs("AAPAAAYAAD"), "")
        self.assertEqual(set(her2.liability_motifs("NGTDGMWAAD").split("|")), {
            "n_linked_sequon", "deamidation_motif", "isomerization_motif",
            "oxidation_susceptible_residue",
        })

    def test_generation_is_deterministic_bounded_and_traceable(self):
        known = her2.NoveltyIndex([PARENT, OTHER_PARENT, HELD_OUT, "AGAGFYAAD", "AVGGAGFYAAA"])
        first = her2.generate_candidates([PARENT, OTHER_PARENT], known, pool_size=60, seed=101)
        second = her2.generate_candidates([OTHER_PARENT, PARENT, PARENT], known, pool_size=60, seed=101)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 60)
        self.assertEqual(len({row["cdrh3"] for row in first}), len(first))
        self.assertEqual({row["origin"] for row in first}, set(her2.ORIGINS))
        for row in first:
            sequence = row["cdrh3"]
            self.assertTrue(her2.is_candidate_sequence(sequence))
            self.assertEqual(row["wt_hamming_distance"], her2.hamming_distance(sequence, her2.WT_CDRH3))
            self.assertGreaterEqual(row["wt_hamming_distance"], 3)
            self.assertTrue(all(edit_distance(sequence, entry) >= 2 for entry in known.sequences))
            self.assertEqual(row["seed"], 101)
            self.assertEqual(row["structural_status"], "missing")
            self.assertEqual(row["thermal_status"], "missing")
            self.assertEqual(row["abguard_status"], "unverified")
            self.assertEqual(row["de_novo_novelty_status"], "unverified")
            self.assertEqual(row["submission_ready"], "false")
            if row["origin"] == "train_high_mutation":
                self.assertIn(row["parent_cdrh3"], [PARENT, OTHER_PARENT])
                self.assertGreaterEqual(edit_distance(sequence, row["parent_cdrh3"]), 2)
            elif row["origin"] == "wt_mutation":
                self.assertEqual(row["parent_cdrh3"], her2.WT_CDRH3)
                self.assertIn(row["wt_hamming_distance"], (3, 4))
            else:
                self.assertEqual(row["parent_cdrh3"], "")
        bounded = her2.generate_candidates([PARENT], known, pool_size=60, max_attempts=2)
        self.assertLessEqual(len(bounded), 2)

    def test_greedy_selection_uses_stable_ties_and_skips_missing_scores(self):
        near_first = "AGGAGFYAAE"
        missing = "KARTSYEFMD"
        pool = [PARENT, near_first, OTHER_PARENT, HELD_OUT, missing]
        chosen = her2.greedy_select(
            pool, {PARENT: 2, near_first: 2, OTHER_PARENT: 3, HELD_OUT: 1},
            count=5, min_hamming=2,
        )
        self.assertEqual(chosen, [OTHER_PARENT, PARENT, HELD_OUT])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            her2.greedy_select([PARENT, PARENT], {}, count=1, min_hamming=2)

    def test_template_uses_fixed_framework_and_light_chain(self):
        row = {
            "cdrh3": PARENT,
            "vh_sequence": her2.VH_PREFIX + PARENT + her2.VH_SUFFIX,
            "vl_sequence": her2.TRASTUZUMAB_VL,
        }
        template = her2.TrastuzumabTemplate.from_facs_row(row)
        self.assertIsNotNone(template)
        chains = template.chains(OTHER_PARENT)
        self.assertEqual(chains["vh_sequence"], row["vh_sequence"].replace(PARENT, OTHER_PARENT))
        self.assertEqual(chains["vl_sequence"], row["vl_sequence"])
        for change in (
            {"cdrh3": OTHER_PARENT},
            {"vh_sequence": "A" + row["vh_sequence"][1:]},
            {"vl_sequence": "A" + row["vl_sequence"][1:]},
        ):
            self.assertIsNone(her2.TrastuzumabTemplate.from_facs_row({**row, **change}))

    def test_cli_uses_held_out_rows_only_for_novelty(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            facs = root / "facs.csv"
            spr = root / "spr.csv"
            train = root / "train.csv"
            val = root / "val.csv"
            test = root / "test.csv"
            archive = root / "facs.csv.zip"
            pool = root / "pool.csv"
            scores = root / "scores.csv"
            shortlist = root / "shortlist.csv"
            facs_fields = ("cdrh3", "binding_class", "vh_sequence", "vl_sequence")
            train_fields = ("seq", "class", "label", "edit_distance", "cluster_centre")
            write_table(train, train_fields, [
                (PARENT, "high", 1, 3, PARENT),
                (OTHER_PARENT, "high", 1, 4, OTHER_PARENT),
                ("AAGAGFYAAD", "mid", 0, 4, PARENT),
                ("CAGAGFYAAD", "high", 1, 4, PARENT),
            ])
            shorter = "AGAGFYAAD"
            write_table(spr, ("cdrh3", "kd_nM"), [(shorter, ""), (her2.WT_CDRH3, "0.152")])
            command = [
                "generate", "--facs", str(archive), "--spr", str(spr), "--train", str(train),
                "--output", str(pool), "--pool-size", "60", "--export-vh-vl",
            ]
            for label in ("high", "low"):
                write_table(val, train_fields, [(HELD_OUT, label, 1, 4, HELD_OUT)])
                write_table(test, train_fields, [("CTGSGFYTMD", label, 1, 4, HELD_OUT)])
                write_table(facs, facs_fields, [
                    (sequence, binding_class, her2.VH_PREFIX + sequence + her2.VH_SUFFIX,
                     her2.TRASTUZUMAB_VL)
                    for sequence, binding_class in (
                        (PARENT, "high"), (OTHER_PARENT, "high"), (HELD_OUT, label), ("CTGSGFYTMD", label)
                    )
                ])
                with ZipFile(archive, "w") as file:
                    file.write(facs, "facs.csv")
                    file.writestr("__MACOSX/._facs.csv", "metadata")
                self.assertEqual(her2.main(command), 0)
                if label == "high":
                    original = pool.read_bytes()
                else:
                    self.assertEqual(pool.read_bytes(), original)

            index, _ = her2.load_corpus(archive, spr)
            self.assertEqual(
                index.sequences,
                {PARENT, OTHER_PARENT, HELD_OUT, "CTGSGFYTMD", shorter, her2.WT_CDRH3},
            )
            self.assertEqual(her2.load_high_parents(train), {PARENT, OTHER_PARENT})
            rows = list(her2.csv_rows(pool, her2.CANDIDATE_COLUMNS))
            self.assertEqual(len(rows), 60)
            self.assertTrue(all(not index.has_neighbor(row["cdrh3"]) for row in rows))

            write_table(scores, ("cdrh3", "predicted_score"), [
                (rows[0]["cdrh3"], ""), (rows[1]["cdrh3"], -0.6), (rows[2]["cdrh3"], 0.9)
            ])
            self.assertEqual(her2.main([
                "select", "--candidates", str(pool), "--scores", str(scores), "--output",
                str(shortlist), "--count", "3", "--min-hamming", "0",
            ]), 0)
            selected = list(her2.csv_rows(shortlist, ("cdrh3", "predicted_score", *her2.CHAIN_COLUMNS)))
            self.assertEqual([row["cdrh3"] for row in selected], [rows[2]["cdrh3"], rows[1]["cdrh3"]])
            self.assertEqual([row["predicted_score"] for row in selected], ["0.9", "-0.6"])
            self.assertEqual([row["structural_status"] for row in selected], ["missing", "missing"])
            self.assertEqual([row["thermal_status"] for row in selected], ["missing", "missing"])


if __name__ == "__main__":
    unittest.main()
