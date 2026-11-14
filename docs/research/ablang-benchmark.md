# Frozen AbLang2 SPR comparison

14 November 2026. Neither frozen representation establishes an improvement over the pairwise ridge baseline on this retrospective dataset. Candidate priorities are unchanged. The protocol was recorded before extraction; the 158-variant cohort, five outer group folds and 75-variant distant diagnostic were already inspected in previous experiments.

| Model | Pooled Spearman, 158 variants | Distant Spearman, 75 variants |
| --- | ---: | ---: |
| Frozen FACS probability | 0.3473306 | 0.1970312 |
| Position ridge | 0.6127403 | 0.4247601 |
| Position plus pair ridge | 0.6484715 | 0.4717739 |
| AbLang2 CDRH3 mean plus ridge | 0.6051904 | 0.5098544 |
| AbLang2 WT-contact mean plus ridge | 0.6564109 | 0.4848468 |

The primary CDRH3 mean versus pairwise difference is -0.0433 (95% paired group interval -0.1388 to 0.0733) overall and +0.0381 (-0.1384 to 0.2283) on distant variants. WT-contact mean versus pairwise is +0.0079 (-0.0971 to 0.1009) overall and +0.0131 (-0.1796 to 0.2065) on distant variants. All four intervals include zero.

The structural pooling ablation changes the point estimate in opposite directions across cohorts. Compared with the full CDRH3 mean, its difference is +0.0512 (-0.0391 to 0.1170) overall and -0.0250 (-0.1495 to 0.0992) on distant variants. These results do not establish either a benefit or an absence of benefit. The fixed-prediction bootstrap does not include all training uncertainty or adjustment for repeated method development on the same dataset.

## Method and verification

Both representations use the same frozen AbLang2 paired VH/VL model, with 480 final hidden-state features. The primary averages all ten library-window residues. The geometric ablation averages WT-contact positions 1, 5, 6 and 7, identified from deposited 1N8Z without SPR labels. It changes pooling, not attention weights. The encoder does not receive a HER2 sequence. Its pretraining overlap with antibody families is not ruled out.

Each ridge head fits the same measured SPR target and uses train-only StandardScaler inside each inner and outer fold. Alpha is selected from 0.1, 1, 10 and 100 by inner grouped RMSE. Neither pooled nor distant test performance selects scaling, alpha or the contact mask. Both heads are reported. No synthetic labels or fine-tuning were used.

The root reproduced every value in both 158 by 480 feature arrays exactly on the same CPU environment and independently reproduced the pairwise OOF predictions. The feature extraction's 18 tests passed with the actual upstream tokenizer. The main test environment skips four tokenizer-specific cases when optional AbLang2 is absent. A CLI append-default bug discovered in integration was fixed and covered by a regression test before the final benchmark replay.

[Feature source and mapping](ablang-feature-protocol.md), [numeric metrics and paired intervals](../../results/ablang-benchmark/metrics.json), [OOF predictions](../../results/ablang-benchmark/oof_predictions.csv), [comparison figure](../../output/interaction-benchmark/comparison.png).

```bash
source .venv/bin/activate
PYTHONPATH=.:src python -m src.her2_interaction_benchmark \
  --data-dir "$HER2_DATA" \
  --feature-npz results/ablang-features/features.npz \
  --feature-metadata results/ablang-features/metadata.json \
  --feature-contrast wt_contact_mean cdrh3_mean \
  --output "$HER2_REPLAY_OUT"
python scripts/plot_interaction_benchmark.py
```
