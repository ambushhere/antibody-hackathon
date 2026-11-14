# Fourth Transformer seed result

A fixed fourth initialization completed on the H100 using the unchanged 4,748,035-parameter architecture and 15 epochs. Best epoch 15 was chosen only by FACS validation AP. FACS test AP is 0.992967. The existing primary 158 SPR variants and distant subset of 75 were reused for frozen-feature transfer evaluation, not an independent test.

| Method | SPR Spearman, 158 | Distant Spearman, 75 |
| --- | ---: | ---: |
| position_ridge | 0.6127 | 0.4248 |
| position_pairwise_ridge | 0.6485 | 0.4718 |
| facs_transformer | 0.3578 | 0.3004 |

Transformer minus pairwise ridge Spearman difference is −0.291, paired group interval [−0.423, −0.081] on 158 variants; on 75 it is −0.171 [−0.357, 0.042]. These fixed-OOF intervals omit training uncertainty and do not correct for repeated method development. No candidate priorities change.

[Protocol](neural-seed3-protocol.md), [training metadata](../../results/neural-runs/large-seed3/metadata.json), [evaluation](../../results/neural-seed3-benchmark/metrics.json), [OOF predictions](../../results/neural-seed3-benchmark/oof_predictions.csv), [selected-weight manifest](../../results/archives/neural-seed3-weights-manifest.json).

This result extends the previous eleven-model comparison by one run. It does not supersede or rewrite that frozen report.
