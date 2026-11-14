# HER2 SPR candidate ranking

Scored all 124 strict local passes from 743 existing candidates. This is an exploratory ranking for review. The SPR axis predicts −log10(KD in nM) from CDRH3 positions. The frozen FACS high class score remains a separate axis. Candidate KD, superiority to trastuzumab and thermal improvement remain experimentally unverified. The benchmark and fit condition on variants with reportable KD.

## Fit and replay

The completed benchmark measured primary cohort OOF Spearman ρ 0.613 for position ridge and 0.347 for FACS. The post-benchmark fit used 158 measured primary variants in 115 sequence groups. Full-cohort grouped tuning selected alpha **10**. Its pooled tuning RMSE is a selection diagnostic.

| Alpha | Pooled tuning RMSE |
| ---: | ---: |
| 0.1 | 1.082603 |
| 1 | 0.753108 |
| 10 | 0.745426 |
| 100 | 0.884024 |

The sensitivity columns give the minimum, maximum and range across the five original outer training refits at their recorded alphas. They describe sensitivity to training composition and regularization. They carry no calibrated coverage.

Across candidates the median range is 0.543 and the maximum is 1.640 score units. Nearest primary SPR training distance is 3–7 edits. All screened rows retain the source VH and VL, official validation and thermal status.

Maximum source mismatch is 0 score units between replayed held-out predictions and the original OOF position ridge scores. The rebuilt assay labels differ from the exported OOF target by at most 0. The maximum preflight versus old top 20 FACS difference is 0.

## Rank disagreements

Candidate SPR and FACS ranks have Spearman ρ -0.517. Their individual top 20 sets share 0 rows. The old top 20 shares 20 with FACS top 20 and 0 with SPR top 20. Equal axis scores share their best rank. Top 20 boundaries break ties by CDRH3.

| CDRH3 | FACS rank | SPR rank | FACS rank minus SPR rank |
| --- | ---: | ---: | ---: |
| WMKPAFYSMA | 122 | 2 | +120 |
| WSNVGFYTMS | 119 | 6 | +113 |
| WHQSGFYVME | 114 | 3 | +111 |
| WRSDGSYAMP | 120 | 12 | +108 |
| WGAPGFYLMT | 108 | 1 | +107 |

Pareto front sizes, listed as front (count), are 1 (9), 2 (13), 3 (11), 4 (11), 5 (9), 6 (9), 7 (11), 8 (9), 9 (8), 10 (8), 11 (7), 12 (4), 13 (5), 14 (4), 15 (4), 16 (1), 17 (1). Front 1 contains all nondominated rows when maximizing both axes. Later fronts apply the same rule after prior fronts are removed.

## Exploratory shortlist

Selected 10 sequences with a realized minimum pairwise CDRH3 Levenshtein distance of 4. The review list takes up to 10 rows by Pareto front, alternates SPR and FACS priority within the remaining eligible front and requires distance ≥2 from every earlier choice. Score ties use lexical CDRH3 order.

Overlap with the old top 20 is 4 (YTGDGFPAMN, WYGPGFPVMD, WTAVGFFAMD, YTADGVYHLK). Front 1 overlaps with the old top 20 in 3 rows.

| Pick | CDRH3 | Front | FACS score | SPR score | Fold spread | Nearest SPR edits | Old rank | Structural flag |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | WGAPGFYLMT | 1 | 0.681096 | -0.534 | 1.559 | 3 | — | — |
| 2 | YTGDGFPAMN | 1 | 0.999699 | -1.426 | 0.515 | 3 | 1 | — |
| 3 | WEKGGFYRMD | 1 | 0.962774 | -0.612 | 0.769 | 4 | — | — |
| 4 | WYGPGFPVMD | 1 | 0.998657 | -0.978 | 1.079 | 4 | 5 | — |
| 5 | FTSDRFYAMD | 1 | 0.964204 | -0.725 | 0.558 | 4 | — | — |
| 6 | WTAVGFFAMD | 1 | 0.996654 | -0.926 | 0.353 | 4 | 18 | — |
| 7 | YQADGFYHMD | 1 | 0.973365 | -0.760 | 0.732 | 3 | — | — |
| 8 | WVGGWFYLMD | 1 | 0.979902 | -0.780 | 1.508 | 3 | — | — |
| 9 | WSGRRFYHMD | 1 | 0.978392 | -0.769 | 0.544 | 4 | — | — |
| 10 | YTADGVYHLK | 2 | 0.999294 | -1.793 | 0.062 | 6 | 2 | Review |

Structural review flag from the old shortlist: YTADGVYHLK. Low predicted VH/VL confidence calls for review. Experimental rejection requires experimental evidence.

## Reproduction

Run from the repository root with the existing environment. The fit and shortlist are deterministic. Seed 0 identifies the completed benchmark's bootstrap.

```sh
source .venv/bin/activate
PYTHONDONTWRITEBYTECODE=1 python -m scripts.rank_spr_candidates \
  --data-dir "${HER2_DATA:?Set directory containing the organizer CSVs}" \
  --preflight results/screen-inputs/expanded-preflight.json \
  --oof results/spr-benchmark/oof_predictions.csv \
  --metrics results/spr-benchmark/metrics.json \
  --facs-model results/experiment-20261114/selected_model.json \
  --old-top20 results/experiment-20261114/provisional_top20.csv \
  --output "${HER2_REPLAY_OUT:?Set a new output directory}"
```

`model.json` records input paths and SHA-256 hashes, versions, feature order, numeric parameters, all fit CDRH3 IDs and original outer test IDs. `candidates.csv` includes both axes, ranks, fronts, fold predictions and flags. `shortlist.csv` contains the selected rows in pick order.
