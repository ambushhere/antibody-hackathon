# HER2 SPR interaction experiment

14 November 2026. On the 158-variant primary cohort, adding pair indicators increased out-of-fold Spearman from 0.612740329 to 0.648471485 and reduced log-KD RMSE from 0.7406 to 0.6869. The paired 95% intervals for both changes span zero in the full cohort and in the distant 75. The advantage over position ridge is inconclusive. The previously reported pairwise values of 0.648 overall and 0.472 on the distant subset reproduce to three decimals under the fixed protocol. [Protocol](interaction-protocol.md) · [Metrics](../../results/interaction-benchmark/metrics.json).

## Pooled and distant results

| Evaluation set | Model | Spearman ρ (95% interval) | Log-KD RMSE (95% interval) |
| --- | --- | ---: | ---: |
| Pooled, 158 variants and 115 groups | Position ridge | 0.613 [0.448, 0.704] | 0.741 [0.627, 0.820] |
| Pooled | Position plus pair ridge | 0.648 [0.465, 0.747] | 0.687 [0.602, 0.769] |
| Pooled | Frozen FACS high probability | 0.347 [0.061, 0.504] | — |
| Distant, 75 variants and 67 groups | Position ridge | 0.425 [0.203, 0.608] | 0.625 [0.519, 0.730] |
| Distant | Position plus pair ridge | 0.472 [0.267, 0.635] | 0.579 [0.483, 0.682] |
| Distant | Frozen FACS high probability | 0.197 [−0.038, 0.421] | — |

| Evaluation set | Paired contrast | Δ Spearman ρ (95% interval) | Δ log-KD RMSE (95% interval) |
| --- | --- | ---: | ---: |
| Pooled | Position plus pair minus position | +0.036 [−0.042, +0.090] | −0.054 [−0.100, +0.017] |
| Distant | Position plus pair minus position | +0.047 [−0.073, +0.161] | −0.046 [−0.105, +0.019] |
| Pooled | Position plus pair minus frozen FACS | +0.301 [+0.159, +0.530] | — |
| Distant | Position plus pair minus frozen FACS | +0.275 [−0.040, +0.549] | — |

Positive Δρ and negative ΔRMSE favor the pair model. Frozen FACS is a probability on a different scale and is evaluated by rank. For each evaluation set, 2,000 draws with seed 0 resampled original sequence groups with replacement, using the same rows across models on each draw. On the distant subset each group contributes only its eligible rows. All reported intervals had 2,000 defined draws. Percentile intervals describe variation in the fixed held-out predictions under group resampling. Refitting variability and a new assay population are outside their scope.

## Outer folds

| Fold | Fit rows / groups | Held-out rows / groups | Position α | Pair α | Position ρ | Pair ρ | FACS ρ | Position RMSE | Pair RMSE |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 123 / 114 | 35 / 1 | 10 | 10 | 0.734 | 0.785 | 0.456 | 0.863 | 0.710 |
| 2 | 127 / 90 | 31 / 25 | 10 | 1 | 0.716 | 0.638 | −0.063 | 0.702 | 0.665 |
| 3 | 127 / 85 | 31 / 30 | 10 | 1 | 0.425 | 0.476 | 0.114 | 0.589 | 0.618 |
| 4 | 127 / 85 | 31 / 30 | 1 | 1 | 0.630 | 0.667 | 0.467 | 0.927 | 0.832 |
| 5 | 128 / 86 | 30 / 29 | 10 | 0.1 | 0.531 | 0.463 | 0.172 | 0.515 | 0.574 |

The pair model improved each metric in three of the five outer folds, with different folds improving on the two metrics. Fold 1 contains a single 35-variant sequence group.

## Fit and replay

The frozen primary set uses quantified positive SPR KD in the exact trastuzumab VH and VL background and ten-residue variant loops. The target is the median of replicate −log10(KD in nM) per CDRH3. WT and exact FACS train and validation overlaps are excluded. Sequence groups are connected components at Levenshtein distance ≤2. The previously viewed distant 75 use the exact rows in the original distance artifact with distance ≥3 from both FACS train/validation and the row’s outer SPR training set. The same 158 out-of-fold predictions supply both evaluations. [Frozen benchmark](../../results/spr-benchmark/report.md) · [Distance rows](../../results/spr-distance-diagnostic/distance_rows.csv).

Both fits have an unpenalized intercept. Position ridge uses 200 one-hot columns and the pair model adds 18,000 one-hot columns across all 45 position pairs, using the existing FACS encoder. Each outer training set independently selects alpha over 0.1, 1, 10 and 100 by pooled inner GroupKFold RMSE, with ties favoring the larger alpha. The CPU solver forms centered Gram matrices from sparse products and solves a system whose width is the number of training rows. Dense primal reference tests compare both encodings. The replay reproduced the frozen outer and inner groups and all five position alphas. The maximum discrepancies were 0 for labels, 2.89 × 10⁻¹⁵ for the 158 held-out position predictions and 6.44 × 10⁻¹⁵ for the inner position RMSE grid, against a tolerance of 10⁻¹⁰. [Evaluator](../../src/her2_interaction_benchmark.py) · [Tests](../../tests/test_her2_interaction_benchmark.py).

This is retrospective method development on already inspected data conditioned on reportable KD. The pair coefficients represent sequence associations in this scaffold. The experiment provides no blind evidence about new designs, unquantified variants, physical contacts, epitope retention or thermal stability. This run covers the one-hot interaction comparison. AbLang2 pooling and the WT contact ablation require a separate evaluation.

## Reproduction

```sh
source .venv/bin/activate
PYTHONDONTWRITEBYTECODE=1 python -m src.her2_interaction_benchmark \
  --data-dir /tmp/her2-benchmark.NjMJ71m1 \
  --output results/interaction-benchmark
```

The CSV contains every sequence, target, original group and fold, distance measurements and all three held-out scores. The JSON records input and code SHA-256 hashes, package versions, fit group membership, inner validation groups, all alpha scores and uncertainty. Outer training group IDs and inner validation group IDs resolve to exact fit rows in the CSV. [OOF predictions](../../results/interaction-benchmark/oof_predictions.csv) · [Fit and metrics](../../results/interaction-benchmark/metrics.json).

The evaluator can also accept `--feature-npz` and `--feature-metadata` with a separate `--output` directory. The NPZ must contain a Unicode `sequence_ids` vector in exact primary-cohort order and named finite real feature matrices of shape `(158, d)`. Each matrix gets a `StandardScaler` fitted on the corresponding training rows at every inner and outer fit. The JSON records the supplied extraction metadata and both file hashes. When the prespecified AbLang2 matrices are available, `--feature-contrast wt_contact_mean cdrh3_mean` reports the contact-mask minus unmasked paired comparison.
