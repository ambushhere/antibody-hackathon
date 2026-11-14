# Experiment checkpoint, 14 November 2026

Research only. No competition submission or synthesis order has been made. No new candidate has measured SPR, competition or thermal results. The initial setup and benchmark were merged in PRs #2 and #3; this checkpoint adds the subsequent experiments.

## Results and locations

| Experiment | Observation and limit | Evidence |
| --- | --- | --- |
| FACS supervision | 524,346 labels; published test ROC AUC 0.9859 and AP 0.9737. This predicts the FACS class, not KD or Tm. | [Original benchmark](../results/experiment-20261114/benchmark_report.md) |
| Direct SPR supervision | Nested sequence-group validation on 158 eligible quantified variants: ridge Spearman 0.613 versus frozen FACS 0.347. Missing KD is not a negative label. | [SPR report](../results/spr-benchmark/report.md) |
| Distant SPR subset | 75 variants at edit distance >=3 from both relevant training sources: 0.425 versus 0.197. The paired difference interval includes zero. | [Distance diagnostic](../results/spr-distance-diagnostic/report.json) |
| Optimisation proposals | 10,000 proposals; 743 deliberately screened, 124 strict local passes and 10 provisional priorities. Absolute KD calibration around WT is inadequate. | [Ranking](../results/spr-expanded-ranking/report.md), [top 10](../results/spr-expanded-ranking/shortlist.csv), [all passes](../results/spr-expanded-ranking/candidates.csv) |
| RFD4 pilots | Four HER2 optimisation designs and twelve exploratory antibody designs. The k=1 pilot produces four variants below 70% paired identity to two local references, one with an extra VH cysteine. Adaptyv database novelty remains unknown. | [Identity audit](../results/experiment-20261114/denovo_identity.csv) |
| Structural diagnostics | OpenFold and Protenix controls, then ten target-only-template Protenix predictions. VH/VL assembly and antibody/target confidence must be separate. A FACS-low inversion is not an SPR nonbinding result. | [Refold scores](../results/denovo-refold/scores.csv), [geometry](../results/denovo-refold/geometry.json) |
| Thermal proxy | ThermoMPNN-D on two paired Fv structures fails some measured control directions and does not improve the small CDR panel ranking. No thermal score is used to promote candidates. | [Thermal pilot](../results/thermal-pilot/README.md) |
| Literature and tool audit | Sources distinguish design from cofolding, folding from binding ddG and measured Tm from proxies. AF2-IG means initial guess. | [Provider papers](research/provider-papers.md), [thermal assessment](research/thermal-assessment.md) |
| HER2 numbering | Nine structure cases now map author residues through mmCIF labels and SIFTS. For 1N8Z, author 557–605 becomes canonical P04626 579–627. Other structures can already use canonical numbering. | [Mapping provenance](../results/epitope-audit/provenance.json), [run output](../results/epitope-audit/run.txt) |
| Negative controls | Inspected supplied supplement does not establish the claimed 454 confirmed negative pairs. Missing fits, expression failures and absent antigen are distinct outcomes. | [Source audit](research/structural-negative-audit.md) |

The corrected numbering and negative-control audit supersede the corresponding older preparation notes in AGENTS.md and STRATEGY.md. Those historical notes have not been silently rewritten.

## Reproduce and continue

Start with [SETUP.md](SETUP.md), activate `.venv`, then follow [BENCHMARK.md](BENCHMARK.md). The organizer CSVs remain external inputs with hashes in the reports. On the team VM they are under `/home/hackathon/almaz/data/hackathon-data/`; obtain them before the VM shutdown. Do not commit credentials or the full organizer dataset.

The compact checked candidate inputs are versioned in [results/screen-inputs](../results/screen-inputs/). They allow both ranking reports to be replayed without the original local `/tmp` directory. For `rank_spr_candidates`, provide `--data-dir /path/to/csvs`, `--preflight results/screen-inputs/expanded-preflight.json`, `--facs-model results/experiment-20261114/selected_model.json`, `--old-top20 results/experiment-20261114/provisional_top20.csv` and a new `--output` directory. Use `initial-preflight.json` to replay the earlier 88-pass report. Input digests must match the recorded benchmark.

Raw working directories are excluded from Git: `results/campaign/` and `results/benchmark/`. A [versioned structural archive](../results/archives/README.md) preserves all 307 exported structural run files, including the completed context panel and RF3 outputs, so these results survive VM shutdown. Curated scores, model parameters, input specifications, controls, provenance and replay scripts are versioned separately. Provider job receipts are in [results/provider-jobs](../results/provider-jobs/).

## Current follow-ups

The completed Protenix panel contains all 21 predictions and compares seven antibodies with a split crop, continuous domain IV and continuous extracellular HER2. [Input specifications](../configs/her2-context/) and [protocol](research/her2-target-context.md) fix the controls, sequences, template coordinates and seed. Raw VM output is `/home/hackathon/almaz/out/campaign-20261114/context-panel/`. A completed seven-target JapanFold RF3 job uses the domain IV sequences without templates or MSA, so its conditioning differs from Protenix. Its job ID is `a6b1ac5596c635805e1f3312d28fd797` and recorded total charge is USD 0.022406. Completed job status alone does not validate its predicted interfaces.

The [measured-data thermal audit](research/thermal-training.md) examined 21 tables and 1,324 records. Confirmed scFv and purified Fab strata each came from only one study, while the target CDRH3 panel was measured as IgG1. No model was fitted and no candidate thermal score was produced. Synthetic mutant sequences can augment proposals, but predicted Tm labels do not become experimental ground truth. A measured stabilising mutation remains a hypothesis when combined with a new CDRH3 or transferred to another construct.

## Explanation and selection policy

[Russian visual primer](../output/explainer/her2-intuition.pdf), [whole-task decisions](EXPERIMENT_DECISIONS.md) and [stability combination hypotheses](research/stability-combination.md). Keep affinity, epitope, thermal evidence and eligibility separate. The local shortlist is an experiment queue, not an approved lab submission or a calibrated probability of S-tier success.
