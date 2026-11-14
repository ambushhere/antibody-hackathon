# Fourth neural seed replication

Fixed before launch on 14 November 2026 at the user's request for a short final GPU experiment. This is a fourth random-initialization replication of the existing large FACS Transformer, not a new independent dataset or end-to-end affinity fine-tune. Earlier seeds and outcomes have already been inspected.

Use the existing unchanged trainer, large Transformer only, seed 3, 15 epochs, CUDA. Keep published FACS splits, input fingerprints, batch size 2048, learning rate 0.0003, weight decay 0.01, dropout 0.1 and checkpoint selection by validation high-class AP. No synthetic augmentation. Export frozen features for the same 158 SPR variants.

Report FACS test AP and frozen-representation nested grouped SPR Spearman for all 158 and the existing distant subset of 75. Use the existing folds and alpha grid, with direct position and pairwise ridge baselines and paired group intervals. Preserve every result regardless of direction. No selection among seeds, no change to the final candidate priorities and no claim of independent confirmation from reuse of these outcomes.

Launch only if GPU is free, under the root GPU lock. Preserve selected weights, metadata, per-epoch history, frozen features and evaluation artifacts. Stop after this one 15-epoch run; no adaptive architecture or hyperparameter sweep. Existing completed experiments and slides remain deliverable if this replication cannot finish.
