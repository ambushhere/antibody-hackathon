# SPR adaptation protocol

Fixed on 14 November 2026 before the SPR adaptation run. This is one retrospective
experiment on the previously examined 158 SPR variants with reportable KD. It uses
the five existing outer sequence group folds, their existing inner group folds
and the frozen distant subset of 75 variants. The target is −log10(KD in nM).
The cohort is reused for method development and provides no independent
replication or basis for candidate promotion.

Both arms start from the archived FACS large Transformer seed 0 `best.pt` selected
by FACS validation AP. Its SHA-256 is
`54c6dd9de3ded6f1112e7c1337b89b48be21c5348a97d6e96124e866bcafdd69`.
Use the existing large architecture (width 256, six blocks, eight heads and
feedforward width 1024) and replace the FACS classifier with a newly initialized
scalar linear head on the normalized CLS representation.

| Arm | Trainable parameters | AdamW learning rate |
| --- | --- | --- |
| Frozen encoder | Scalar head | 0.001 |
| Final block adaptation | Same head, final Transformer block and final encoder LayerNorm | 0.001 for the head, 0.00001 for the block and norm |

Seed 0 resets for every inner and outer fit. All encoder layers stay in evaluation
mode, including the trainable final block and norm, so dropout is disabled in both
arms. Use float32, full batch mean squared error and AdamW weight decay 0.01 for
each parameter group. Compute the mean and population standard deviation of the
target on the current training rows only. Train on that scale and invert the
transform before computing errors on validation or held-out rows.

For each arm and outer fold, train each inner fit through 100 epochs and record
its validation predictions after 20, 50 and 100 epochs. Choose the epoch by
minimum pooled inner RMSE on the original target scale, with ties going to the
earlier epoch. Reset from the seed 0 checkpoint and a fresh seed 0 head, then fit
on the entire outer training fold for that epoch count. Generate one held-out
prediction per variant per arm. Outer held-out targets are used for evaluation.

Report pooled and distant Spearman and RMSE alongside the archived aligned
position and pairwise ridge OOF predictions. Report paired differences between
the two neural arms and from each neural arm to both ridge baselines. Use 2,000
paired percentile bootstrap draws of sequence groups within each subset, with
the existing helper and its seed. Intervals describe variation in the fixed
OOF predictions under group resampling.

```bash
source .venv/bin/activate
python -m scripts.benchmark_spr_finetune \
  --data-dir "$SPR_DATA_DIR" \
  --checkpoint "$SPR_SEED0_CHECKPOINT" \
  --output results/spr-finetune-seed0 \
  --device cuda
```

Set `SPR_DATA_DIR` to the directory with uncompressed `train.csv`, `val.csv`,
`test.csv` and `spr.csv`. Set `SPR_SEED0_CHECKPOINT` to the extracted
`neural-runs/large-seed0/transformer/best.pt` from
`results/archives/neural-large-seed0-weights.tar.gz`. Use a new output directory.
The run writes OOF predictions, metrics with provenance and a Markdown report.
