# Neural training and synthetic augmentation results

14 November 2026. Eleven prespecified models were trained on the team H100. Large networks improved FACS classification, but their frozen representations did not establish an improvement over direct SPR ridge for affinity ranking. The synthetic arm did not consistently improve SPR transfer. Existing candidate priorities are unchanged.

| Model | FACS test AP | SPR Spearman, 158 variants | Distant SPR Spearman, 75 variants |
| --- | ---: | ---: | ---: |
| Historical pairwise FACS logistic | 0.9737 | 0.3473, direct FACS probability | 0.1970 |
| Direct SPR position ridge | not applicable | 0.6127 | 0.4248 |
| Direct SPR position plus pair ridge | not applicable | 0.6485 | 0.4718 |
| Small Transformer, seed 0 | 0.9212 | 0.3420 | 0.4667 |
| Small MLP, seed 0 | 0.9330 | 0.2919 | 0.2623 |
| Large MLP, seeds 0 / 1 / 2 | 0.9879 / 0.9881 / 0.9884 | 0.2646 / 0.0992 / 0.1470 | 0.1936 / 0.1049 / 0.0648 |
| Large Transformer, seeds 0 / 1 / 2 | 0.9924 / 0.9923 / 0.9930 | 0.3846 / 0.4916 / 0.4066 | 0.3591 / 0.5653 / 0.4141 |
| Transformer plus pseudolabels, seeds 0 / 1 / 2 | 0.9930 / 0.9935 / 0.9941 | 0.3986 / 0.3845 / 0.4579 | 0.4675 / 0.3868 / 0.5049 |

Every neural SPR entry uses a frozen encoder and the existing nested group-fold ridge head. The FACS probability row is a different, direct score. The direct SPR ridge baselines use measured SPR training rows and fixed sequence features. All heads share the same held-out rows, groups and alpha grid. Small pilots used five epochs and large comparisons used fifteen, so size and budget are confounded across those two tiers. Within each large architecture comparison, data, epochs, optimisation settings and seed numbers match.

## Interpretation

The large Transformer has 4,748,035 parameters and the large MLP has 5,133,315. The Transformer gives higher FACS test AP at each matched seed. On pooled SPR, Transformer minus MLP Spearman differences are +0.120, +0.392 and +0.260; the seed 0 paired interval includes zero. None of the three Transformer heads exceeds the pairwise SPR baseline's point estimate on the full cohort. On the distant subset, seed 1 exceeds the pairwise point estimate while seeds 0 and 2 do not. Picking that seed from the SPR result would reuse the evaluation for model selection, so all three remain reported.

Synthetic augmentation improves FACS test AP slightly at each seed. Its pooled SPR differences from the corresponding real-only Transformer are +0.014, -0.107 and +0.051. The paired 95% group intervals are [-0.144, 0.139], [-0.231, -0.020] and [-0.085, 0.202]. Distant-subset differences are +0.108, -0.179 and +0.091, with intervals [-0.095, 0.324], [-0.338, -0.003] and [-0.097, 0.273]. The seed 1 deterioration appears in both diagnostics; these intervals are not adjusted for multiple comparisons. There is no consistent affinity-transfer benefit across seeds.

A plausible explanation is that learning FACS categories emphasises features that do not preserve quantitative affinity information within those categories or across assays. This is an interpretation, not an established causal mechanism. The experiment tested frozen representation transfer, not end-to-end SPR fine-tuning or every possible neural architecture. The small retrospective SPR dataset and repeated method development limit the conclusion.

## Data and safeguards

The published FACS train/validation/test sizes are 367,044 / 78,651 / 78,651, with zero exact sequence and cluster-centre overlaps. The primary SPR cohort has zero exact matches to FACS training or validation and six matches to FACS test. SPR values were not read by the encoder training or checkpoint selector. Best epochs were selected only by validation high-class AP. Every prespecified seed was evaluated; no new candidate score was promoted from these results.

The synthetic set contains 100,000 unique proposals: 33,333 / 33,334 / 33,333 have one / two / three substitutions from a real training parent. Accepted mutation counts are balanced in randomised blocks. There are no exact matches to any real FACS split, the 158 primary SPR IDs or WT; no duplicate or cysteine-containing proposals. The large real-only MLP seed 0 generated soft FACS probabilities. The student used real CE plus 0.2 soft-target CE, preserving every real batch. Synthetic labels are not KD, Tm or new laboratory data. Additional dropout draws mean optimisation noise is not identical to the baseline even though initialisation and real batch order are paired.

## Reproduction and artifacts

The [protocol](neural-protocol.md) was committed before training. The VM ran PyTorch 2.15.0+cu132 on H100 80GB, with float32, deterministic algorithms and four CPU threads per process. The admitted trainer source was frozen and its SHA-256 matches every real run. Run metadata records the exact package versions, input hashes, training history, selected checkpoint and source bytes. The summed recorded fit durations are about 17 minutes, excluding loading, transfer evaluation and artifact export; this is not a provider billing measure.

- [All experiment rows](../../results/neural-summary/experiments.csv)
- [All-model SPR benchmark and paired intervals](../../results/neural-all-benchmark/metrics.json)
- [Per-sequence out-of-fold predictions](../../results/neural-all-benchmark/oof_predictions.csv)
- [Training runs and metadata](../../results/neural-runs/)
- [Synthetic data provenance](../../results/neural-synthetic/metadata.json)
- [Selected weight archive manifest](../../results/archives/neural-weights-manifest.json)
- [Figure](../../output/neural-benchmark/facs-vs-spr.png)

```bash
source .venv/bin/activate
PYTHONPATH=.:src python -m scripts.train_her2_neural \
  --data-dir "$HER2_DATA" --output "$NEW_RUN" \
  --size large --epochs 15 --seed 0 --models mlp transformer --device cuda
PYTHONPATH=.:src python -m scripts.generate_her2_synthetic \
  --data-dir "$HER2_DATA" --teacher-run "$NEW_RUN" \
  --output "$NEW_SYNTHETIC" --count 100000 --seed 0 --device cuda
PYTHONPATH=.:src python -m scripts.train_her2_augmented \
  --data-dir "$HER2_DATA" --baseline "$NEW_RUN" \
  --synthetic "$NEW_SYNTHETIC" --output "$NEW_AUGMENTED" --device cuda
PYTHONPATH=.:src python -m scripts.summarize_neural_experiments
```

Use fresh output directories. Repeat the large real and augmented arms at seeds 1 and 2 for the full protocol, using the same seed 0 synthetic teacher dataset. `scripts.collect_neural_features` combines every declared run for the existing `src.her2_interaction_benchmark` without choosing a winner. The committed `results/neural-all-features/metadata.json` lists the exact seven input run directories and eleven representations. Restore selected model tensors by verifying each archive and member hash from the manifest and extracting under `results/`; the archive members start with `neural-runs/`. Use `torch.load(..., weights_only=True)` for these state dictionaries.
