# Fixed hybrid SPR diagnostic

14 November 2026. This is retrospective method development on the already inspected 158-row grouped SPR cohort, not a new independent validation. Hypothesis: direct position effects and the frozen FACS-trained Transformer representations may contain complementary affinity information even when the latter alone rank less well.

Freeze three additional prediction rules before computing their results:

1. Arithmetic mean of the three real-only large Transformer SPR heads (seeds 0, 1 and 2).
2. Equal average of that ensemble and the position SPR ridge.
3. Equal average of that ensemble and the position-plus-pair SPR ridge.

All components already predict the same negative-log10 KD in nM using matching outer group folds. No seed selection, learned weights, rescaling, clipping or external label tuning. No new fitting. Retain all three rules regardless of result. Check row uniqueness, exact frozen file fingerprint and finite predictions. Report Spearman, RMSE and paired group bootstrap intervals on all 158 rows and the unchanged 75 distant rows using the existing bootstrap seed and draw count. Compare every new rule with both direct SPR baselines. Keep the existing candidate ranking unchanged; any apparent advantage requires external or prospective validation, and multiple-comparison intervals remain descriptive.

Input: `results/neural-all-benchmark/oof_predictions.csv` at commit `ed0f9431db129c9740ad81a1fe4bea8bf86aa49a`. All individual results were already viewed before this proposal. The external Absci check remains a separate experiment of the previously fixed sequence models and does not select these blend weights.
