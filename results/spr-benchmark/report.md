# HER2 SPR affinity benchmark

The target is −log10(KD in nM), aggregated by the median across quantified replicates of each CDRH3. Rows must have finite positive KD, the exact trastuzumab VH framework and VL and a ten-residue loop. Wild-type reference rows only define the enrichment thresholds.

## Assay audit

| Binding result | Missing KD | Invalid KD | Finite positive KD |
| --- | ---: | ---: | ---: |
| inconclusive | 104 | 0 | 0 |
| no_reportable_kd | 460 | 0 | 0 |
| not_reported | 14 | 0 | 0 |
| quantified_binding | 0 | 0 | 190 |

From 190 quantified rows, exclusions were 0 other frameworks, 8 wild-type references and 3 unsupported lengths (3 of length 9). Aggregation collapsed 0 variant replicate rows.

Training and validation FACS overlaps are excluded from the primary cohort by exact CDRH3 match. FACS was fit on train and selected on validation. The broader cohort includes those overlaps. The primary cohort can still contain near neighbors of the FACS training data.

The 8 wild-type measurements span 0.02851802 to 0.178 nM. Top-k enrichment uses two thresholds, strictly below the best wild type and at or below the weakest wild type. The observed reference range is descriptive and supplies no statistical equivalence margin. All results condition on reportable KD. The missing statuses supply no quantitative labels.

Each cohort is grouped separately by connected components at Levenshtein distance ≤2. Outer GroupKFold produces pooled out-of-fold predictions. In each outer training set, inner GroupKFold minimizes pooled RMSE over ridge alpha 0.1, 1.0, 10.0, 100.0. Both splitters use at most 5 folds.

Position ridge uses a fixed amino-acid one-hot encoding. The scalar baseline fits a linear calibration of the frozen FACS probability against SPR labels. Residual ridge adds a one-hot ridge prediction of that calibration's residuals. Both calibration and residual fitting are repeated inside each training fold. The FACS probability is evaluated by rank.

Intervals are 95% percentile intervals from 2,000 draws of whole sequence groups with replacement (seed 0). They summarize uncertainty in the fixed out-of-fold predictions. Paired differences use the same draws. The mean baseline ties within each fold and its pooled rank can reflect differences between fold means.

## Primary novelty cohort

158 unique CDRH3s in 115 groups. Group size × count is 1 × 110, 2 × 3, 7 × 1, 35 × 1. FACS overlaps are train 0, val 0, test 6.

| Outer fold | Train sequences / groups | Test sequences / groups | Inner folds | Position alpha | Residual alpha |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 123 / 114 | 35 / 1 | 5 | 10 | 10 |
| 2 | 127 / 90 | 31 / 25 | 5 | 10 | 10 |
| 3 | 127 / 85 | 31 / 30 | 5 | 10 | 10 |
| 4 | 127 / 85 | 31 / 30 | 5 | 1 | 1 |
| 5 | 128 / 86 | 30 / 29 | 5 | 10 | 10 |

| Model | Spearman ρ | 95% interval | RMSE | 95% interval |
| --- | ---: | --- | ---: | --- |
| Training mean | -0.252 | [-0.412, 0.158] | 0.978 | [0.794, 1.101] |
| Frozen FACS probability | 0.347 | [0.061, 0.504] | — | — |
| FACS scalar calibration | 0.017 | [-0.113, 0.253] | 0.969 | [0.802, 1.084] |
| Position ridge | 0.613 | [0.448, 0.704] | 0.741 | [0.627, 0.820] |
| FACS residual ridge | 0.587 | [0.416, 0.691] | 0.768 | [0.663, 0.845] |

- Position ridge minus frozen facs probability Δρ 0.265 [0.110, 0.515] (2,000 defined draws).

- FACS residual ridge minus frozen facs probability Δρ 0.240 [0.093, 0.471] (2,000 defined draws).

- FACS residual ridge minus facs scalar calibration Δρ 0.570 [0.244, 0.720] (2,000 defined draws).

Removing the largest group from evaluation leaves 123 predictions. Spearman ρ is 0.556 for position ridge, 0.525 for residual ridge and 0.198 for frozen FACS. This check uses the original out-of-fold predictions.

1/158 variants beat the best wild type and 4/158 match or beat the weakest.

| Model | Top-k | Below best WT precision / enrichment | At or below weakest WT precision / enrichment |
| --- | ---: | ---: | ---: |
| Training mean | 5 | 0.000 / 0.000 | 0.000 / 0.000 |
| Training mean | 10 | 0.000 / 0.000 | 0.000 / 0.000 |
| Frozen FACS probability | 5 | 0.000 / 0.000 | 0.200 / 7.900 |
| Frozen FACS probability | 10 | 0.000 / 0.000 | 0.200 / 7.900 |
| FACS scalar calibration | 5 | 0.000 / 0.000 | 0.000 / 0.000 |
| FACS scalar calibration | 10 | 0.000 / 0.000 | 0.000 / 0.000 |
| Position ridge | 5 | 0.200 / 31.600 | 0.200 / 7.900 |
| Position ridge | 10 | 0.100 / 15.800 | 0.200 / 7.900 |
| FACS residual ridge | 5 | 0.200 / 31.600 | 0.200 / 7.900 |
| FACS residual ridge | 10 | 0.100 / 15.800 | 0.200 / 7.900 |

## Broader cohort with FACS overlaps

179 unique CDRH3s in 126 groups. Group size × count is 1 × 121, 2 × 3, 7 × 1, 45 × 1. FACS overlaps are train 14, val 7, test 6.

| Outer fold | Train sequences / groups | Test sequences / groups | Inner folds | Position alpha | Residual alpha |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 134 / 125 | 45 / 1 | 5 | 10 | 10 |
| 2 | 145 / 98 | 34 / 28 | 5 | 10 | 10 |
| 3 | 145 / 93 | 34 / 33 | 5 | 10 | 10 |
| 4 | 146 / 94 | 33 / 32 | 5 | 10 | 10 |
| 5 | 146 / 94 | 33 / 32 | 5 | 10 | 10 |

| Model | Spearman ρ | 95% interval | RMSE | 95% interval |
| --- | ---: | --- | ---: | --- |
| Training mean | -0.332 | [-0.446, -0.027] | 1.062 | [0.835, 1.227] |
| Frozen FACS probability | 0.287 | [-0.045, 0.457] | — | — |
| FACS scalar calibration | -0.177 | [-0.382, 0.173] | 1.070 | [0.843, 1.235] |
| Position ridge | 0.682 | [0.559, 0.749] | 0.790 | [0.618, 0.900] |
| FACS residual ridge | 0.672 | [0.541, 0.742] | 0.804 | [0.640, 0.910] |

- Position ridge minus frozen facs probability Δρ 0.395 [0.238, 0.699] (2,000 defined draws).

- FACS residual ridge minus frozen facs probability Δρ 0.385 [0.231, 0.678] (2,000 defined draws).

- FACS residual ridge minus facs scalar calibration Δρ 0.849 [0.443, 1.080] (2,000 defined draws).

Removing the largest group from evaluation leaves 134 predictions. Spearman ρ is 0.644 for position ridge, 0.630 for residual ridge and 0.079 for frozen FACS. This check uses the original out-of-fold predictions.

1/179 variants beat the best wild type and 7/179 match or beat the weakest.

| Model | Top-k | Below best WT precision / enrichment | At or below weakest WT precision / enrichment |
| --- | ---: | ---: | ---: |
| Training mean | 5 | 0.000 / 0.000 | 0.000 / 0.000 |
| Training mean | 10 | 0.000 / 0.000 | 0.000 / 0.000 |
| Frozen FACS probability | 5 | 0.000 / 0.000 | 0.200 / 5.114 |
| Frozen FACS probability | 10 | 0.000 / 0.000 | 0.200 / 5.114 |
| FACS scalar calibration | 5 | 0.000 / 0.000 | 0.000 / 0.000 |
| FACS scalar calibration | 10 | 0.000 / 0.000 | 0.000 / 0.000 |
| Position ridge | 5 | 0.000 / 0.000 | 0.000 / 0.000 |
| Position ridge | 10 | 0.000 / 0.000 | 0.100 / 2.557 |
| FACS residual ridge | 5 | 0.000 / 0.000 | 0.200 / 5.114 |
| FACS residual ridge | 10 | 0.100 / 17.900 | 0.400 / 10.229 |

Top-k ties contribute expected precision. Enrichment uses each cohort's prevalence. With very few positives these values are descriptive. The JSON includes conditional enrichment intervals and valid draw counts, since draws without positives have undefined enrichment.

## Provenance

- python `3.13.18`
- numpy `2.5.4`
- pandas `3.1.1`
- scikit-learn `1.9.1`
- scipy `1.18.1`
- rapidfuzz `3.14.6`
- threadpoolctl `3.8.0`
- train `b639077ff116daa8b8c3a2122b6ed9291361892469f6237963967ca0772cb262` (`/tmp/her2-benchmark.NjMJ71m1/train.csv`)
- val `ea5e2e0908f018f579e476cef4452f9de11ba9144513bffa783588d49123f5d1` (`/tmp/her2-benchmark.NjMJ71m1/val.csv`)
- test `b350c10c7df4eb5e1b966f76909a805fc7d8011486f1ddb5dd1d7a6e19556a06` (`/tmp/her2-benchmark.NjMJ71m1/test.csv`)
- facs `7ff5a95bcc1942bb6e34b1e4ed6c02a9c177ee30ca160440f8abe82ed787735e` (`/tmp/her2-benchmark.NjMJ71m1/facs.csv`)
- spr `4d6f53a6183f377670e9577b3233a7c4c14888965ac63824863219d13424fa9b` (`/tmp/her2-benchmark.NjMJ71m1/spr.csv`)
- facs_model `4067c8e66e703a4db828cec37086f3d157fcad38c172ab265f1aebffe9f3a487` (`results/experiment-20261114/selected_model.json`)
