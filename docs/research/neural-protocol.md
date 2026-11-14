# Neural FACS transfer and synthetic augmentation protocol

14 November 2026, recorded before neural training. This extends retrospective method development on previously inspected FACS and SPR benchmarks. It is not a fresh blind test. Do not change candidate rankings or claim thermal improvement from these experiments.

## Hypotheses and fixed design

Test whether bidirectional attention learns useful conditional residue effects beyond a similarly sized MLP, whether additional capacity helps and whether teacher pseudolabels improve generalisation. Inputs are only the ten variable amino acids in the fixed trastuzumab background. Framework and HER2 are constant; the model is not a general antibody affinity predictor.

Use the published remove_overlap clonotype train/validation/test partitions, with exact sequence and cluster overlap checks, source fingerprints and class mapping low/mid/high. Optimise three-class cross entropy on training rows. Select the best epoch by validation high-class average precision, ties choose the earlier epoch. Evaluate the FACS test once per prespecified arm after selection. Do not use SPR outcomes for encoder or checkpoint selection.

Small pilot, seed 0, five epochs: Transformer with residue dimension 64, two encoder layers, four heads, FFN width 128 and CLS output 64. MLP uses one-hot 200 to 256 to 64 to 3. Both have dropout 0.1. These verify the implementation and establish a small-budget reference.

Large comparison, seeds 0, 1 and 2, fifteen epochs per arm: Transformer with dimension 256, six layers, eight heads, FFN width 1024 and CLS output 256; MLP uses one-hot 200 to 2048 to 2048 to 256 to 3. Both use dropout 0.1, AdamW learning rate 0.0003, weight decay 0.01 and real-data batch size 2048. Report actual parameter counts. Architectures have comparable rather than identical capacity. Larger versus smaller comparisons also differ in training budget; only matched-budget Transformer versus MLP tests the architecture contrast. Use full float32. Do not tune the protocol after seeing test or SPR scores.

The Transformer uses ten residue tokens plus a learned CLS token, learned absolute positions, bidirectional attention without a causal mask, GELU and a three-class head. Attention coefficients are diagnostic, not evidence of molecular contacts. Save checkpoints after each epoch and preserve immutable run folders. On this shared H100, root controls sequential admission and checks other GPU processes. Export results before VM shutdown.

## Transfer evaluation

Freeze the selected encoders. Export latent vectors for the existing 158 primary SPR identifiers in their exact frozen order. All primary rows must remain absent from FACS training and validation. Fit ridge heads with the existing nested sequence-group folds, train-fold-only scaling and alpha grid. Preserve the 75 distant diagnostic rows. Report all prespecified arms and seeds, rather than selecting the best seed. The position and pairwise baselines remain unchanged. Compare paired OOF predictions with sequence-group bootstrap intervals, which do not capture all training uncertainty or repeated development.

Report FACS AP and ROC AUC, SPR pooled Spearman and log-KD RMSE, distant-subset metrics, per-seed variability and compute. A better FACS classifier alone is insufficient evidence for affinity selection. Neither dataset supplies new Tm measurements. No synthetic predictions enter the held-out evaluation labels.

## Separate synthetic arm

After the real-only large MLP seed 0 checkpoint is fixed, generate 100000 unique ten-residue proposals from real FACS training parents. Change uniformly one, two or three distinct positions to other amino acids. Reject exact membership in all real FACS splits, all primary SPR identifiers, WT and previously proposed rows; reject cysteine-containing proposals for this pilot. Validation/test sequences are used only for exclusion, never their labels. Save parent and mutation provenance and bounded generation counts.

The frozen real-only MLP teacher assigns three-class probabilities. Save these explicitly as pseudolabels with teacher and input hashes, not measured binding, KD or Tm. Do not confidence-filter proposals based on test results. Train a separate large Transformer for fifteen epochs, seeds 0, 1 and 2, with loss real CE + 0.2 synthetic soft-target CE. Use one synthetic batch per real batch, cycling independently, without replacing real examples. Keep real sample order and initialisation paired with the baseline where feasible, record any RNG differences from additional dropout. There are more forward passes per step, so this is an augmentation comparison, not a compute-matched architecture experiment.

The immutable real-only large Transformer is the control. Evaluate all arms on real held-out FACS and the existing real SPR benchmark. If augmentation only improves agreement with its teacher, that is not evidence of better biology. Report negative findings and do not generate synthetic experimental truth.
