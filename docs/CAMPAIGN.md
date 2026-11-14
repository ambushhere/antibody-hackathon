# HER2 design campaign

## Outcome and evidence boundary

Produce a ranked reserve list of up to 20 candidates across optimisation and an exploratory de novo arm. Four guaranteed experimental slots are shared, not four per track. Predicted scores are not experimental binding, epitope preservation or thermal stability.

## Rules audit, 14 November 2026

Primary local sources: `2026-11_Oslo_hackathon_Apheris_challenge.pdf`, slides 5–7, and `ai-antibodies-hackathon-apheris.pdf`, pages 2–3 and 7–8. The slides state optimisation CDRH3 identity <=70% and distance >=2; task text states <70% and strictly more than 2. Keep permissive and conservative gates separate. For ten residues conservative WT distance is >=4; conservative distance from provided sequences is >=3. Do not silently treat either as the confirmed validator implementation.

De novo means type-aware identity <70% to known antibody sequences in Adaptyv's database, not just CDRH3 novelty. VHH compares to VHH/heavy chains; paired formats compare as paired chains. Additional de novo requirement in task text: at least two CDRH3 changes from trastuzumab. Known scaffold CDR redesign can still be classified as optimisation. No local surrogate establishes database novelty. The official validator is required before lab submission.

The shared length limit is 250 residues, with task page 7 explicitly including the full scFv or combined VH/VL length. A linker can push a paired-domain design over the limit. No unpaired cysteines and immunum confidence >=0.1. Treat missing official validation as pending, never passed.

## Optimisation arms

1. Position-frequency baseline from training high binders, with pseudocounts.
2. Constrained WT mutations, including conservative WT-distance gate.
3. Two or more substitutions around diverse training high-binder parents. No held-out labels in generation.
4. RFD4 HER2 novel, followed by novelty and sequence filters.
5. HER2 exploratory as an additional arm only if the pilot supports throughput.

Compare accepted yield, duplicate rate, conservative novelty yield, predicted enrichment, liabilities and sequence diversity by arm. Screen cheaply before independent refolding. Do not select candidates solely by maximum classifier confidence, which is vulnerable to optimisation outside its training distribution.

## Benchmark contract

Use published remove_overlap clonotype train/validation/test splits. Audit sequence overlaps and label mapping against the provided FACS table. Compare dummy, position one-hot logistic and pairwise-interaction logistic models. Select using validation average precision, freeze selection, then report test ROC AUC, average precision, Brier score and top-fraction precision/enrichment. Report distance-to-WT distribution and candidate distance to training data where available. Single split estimates are not universal performance guarantees.

SPR is an external assay check, not a synonym for FACS classes. Aggregate repeated sequences, exclude WT reference and unrelated antibody backgrounds and remove training/validation overlaps. Report association with -log10(KD nM) for quantified values only. Missing, inconclusive and unreported KD are distinct, not automatic negatives. Report measured WT range rather than selecting a convenient reference row. No blind test claim for pretrained HER2 generation weights whose training examples are not disclosed.

The supplied Champloo supplement does not establish the previously claimed 454 confirmed negatives or 21 high-confidence negative predictions. Its missing fits, expression failures and unavailable antigens are distinct outcomes. Row-level AF3 confidence exists for only three cognate pairs. See [the source audit](research/structural-negative-audit.md). No false-positive rate or calibrated HER2 threshold can be inferred from these tables. Never transfer AF3 thresholds to Protenix/OpenFold without validation.

## Structural checks

Use the same HER2 crop, settings, seeds and format for reference and candidates. Record template use: template-conditioned refolding is not independent proof of the pose. Include WT, measured weak/strong controls and candidates if throughput permits. Include Apheris OpenFold for the sponsor prize requirement (participant guide).

Compute antibody-to-HER2 contacts explicitly, excluding VH–VL contacts. Organizer docs warn that the generic Fv interface score can reflect VH–VL rather than target binding. Report target-contact residues, overlap with WT footprint, heavy atom clashes and per-chain confidence where available. A predicted contact footprint cannot prove competition. Do not use CDRH3 RMSD alone as a rejection rule: organizer docs state Protenix fails to reproduce WT CDRH3.

## De novo arm

Start with the general antibody checkpoint and a supplied antibody scaffold, regenerate all CDRs with HER2 domain-IV hotspots. Keep candidate origins and scaffold identity. This is exploratory new-target design; track eligibility remains pending. If scaffold identity alone prevents <70%, do not relabel it de novo. A second scaffold or framework redesign requires checking antibody numbering and folding; global novelty cannot be solved by adding random mutations.

## Stability and selection

No validated thermal predictor has yet been identified in the installed tested Nebius catalog. `stab-ddg` estimates binding free energy change, not folding stability or melting temperature. Sequence liabilities, charge and hydrophobicity are risk descriptors only. Keep thermal evidence explicitly missing. Seek an appropriate predictor with WT-relative validation before asserting any thermal ranking.

Use separate columns for official eligibility, binding prediction, structural target evidence, stability evidence and liabilities. Preserve Pareto trade-offs and diversity across parent families; do not sum incomparable scores into a claimed probability of S tier. Reserve candidates with missing evidence must remain marked provisional.

## Delivery

Save exact commands, seed, image digest, input hashes, model files, score sources and failures. Export results off the temporary VM during work. Submit only after official validation. Check the fresh-repository rule before publishing the competition submission; the development PR is not itself the submission repository.

## Current experiment status

User requested experiments and explanation before any competition submission. No submission has been made. See BENCHMARK.md and results/experiment-20261114 for observed results and limitations. Current follow-ups test direct SPR supervision, sequence diversity in full-antibody generation, target-only-template refolding and published thermal controls. Each follow-up addresses a specific weakness before scaling generation.
