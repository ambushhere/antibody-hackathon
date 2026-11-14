# Fixed SPR ensemble results

14 November 2026. The three fixed rules in the [protocol](hybrid-spr-protocol.md) were evaluated without fitting blend weights or selecting a seed. Averaging the frozen neural heads reduced seed sensitivity but did not exceed direct pairwise SPR ridge on the full cohort. Both half-weight blends gave higher distant-subset point estimates than their respective direct ridge component. The paired intervals include zero, so candidate priorities remain unchanged.

| Rule | Spearman, all 158 | RMSE, all 158 | Spearman, distant 75 | RMSE, distant 75 |
| --- | ---: | ---: | ---: | ---: |
| Position ridge | 0.613 | 0.741 | 0.425 | 0.625 |
| Position plus pair ridge | 0.648 | 0.687 | 0.472 | 0.579 |
| Mean of three Transformer heads | 0.483 | 0.839 | 0.502 | 0.647 |
| Half position, half Transformer mean | 0.601 | 0.759 | 0.505 | 0.587 |
| Half pairwise, half Transformer mean | 0.627 | 0.731 | 0.527 | 0.570 |

The position blend minus position ridge distant Spearman difference is +0.080 (95% paired group interval −0.017 to +0.188). The pairwise blend minus pairwise ridge difference is +0.055 (−0.053 to +0.166). On the full cohort these differences are −0.011 (−0.063 to +0.055) and −0.022 (−0.070 to +0.044). Intervals condition on fixed OOF predictions and are not adjusted for the multiple rules and subsets.

A possible explanation is complementary sequence information far from training data. This diagnostic does not distinguish that explanation from noise or cohort composition. It is retrospective development on already inspected outcomes. The external SPR check independently evaluates the previously fixed sequence models; this blend has not been externally evaluated. No thermal, epitope or laboratory conclusion follows.

[All metrics and contrasts](../../results/hybrid-spr/metrics.json) · [Per-row predictions](../../results/hybrid-spr/oof_predictions.csv)

```bash
source .venv/bin/activate
PYTHONPATH=.:src python -m scripts.benchmark_hybrid_spr
PYTHONPATH=.:src python -m pytest -q tests/test_hybrid_spr.py
```
