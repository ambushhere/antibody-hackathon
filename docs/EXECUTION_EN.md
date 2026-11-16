# Execution guide for judges

This repository contains custom code, recorded model outputs and a reproducible selection pipeline for both HER2 challenge tracks. This page is the English entry point for the Execution criterion. The handoff contains the final team-named candidate file and the final two-slide pitch. Submission and organizer acceptance are unverified; laboratory outcomes are not recorded.

## Presentation and deliverables

- [Final two-slide presentation, editable PowerPoint](../output/selection-slide/her2-selection-and-evidence-en.pptx), with [exact previews](../output/selection-slide/README.md) and an [English pitch summary](../PITCH.md).
- [English explanation and speaking notes](PIPELINE_EXPLAINED_EN.md), [English PDF](../output/pdf/her2-pipeline-explained-en.pdf) and [Russian PDF](../output/pdf/her2-pipeline-explained-ru.pdf).
- [Ranked CSV](../output/submission-draft/epitope-hunters.csv): 30 paired designs, 10 optimization and 20 denovo, with explicit `track`. The first four rows contain three optimization and one de novo design.
- [Exact-row validation receipt](../results/abguard-ranked/receipt.json) and [selection diagnostics](../output/submission-draft/team-3-ranked-draft30.selection.json).
- [Experiment index](EXPERIMENT_INDEX.md): model comparisons, data sizes, metrics, limitations and linked evidence.

## Optimization: generate, score and filter

Keep the trastuzumab heavy/light framework and change the ten-residue library window inside CDRH3. Generate local sequence proposals and score them with separate FACS classification and SPR regression signals. The FACS library contains 524,346 variants; the measured SPR evaluation contains 158 eligible variants with a nested distant subset of 75. The published FACS splits and grouped nested SPR validation are documented in the [experiment index](EXPERIMENT_INDEX.md).

The final optimization shortlist used the 200-feature positional ridge fit on 158 eligible variants with `alpha=10`, alongside the separate FACS logistic score. Later pairwise ridge was a retrospective comparison, not the model used to order that shortlist. Its paired difference intervals versus positional ridge span zero on both evaluated cohorts. On the reused SPR evaluation, position ridge has Spearman correlation 0.613 and position-plus-pair ridge 0.648. These are retrospective ranking results, not calibrated candidate KD or proof of improvement over WT. Transformer, MLP, AbLang2 and synthetic-label experiments did not justify changing the selected candidate order.

The statistical campaign generated 10,000 proposals. The expanded screen contained 739 statistical candidates and four earlier RFD4 optimization-pilot designs, with 124 strict local passes among those 743. The installed abguard independently accepted the same 124/743. Ten review picks were retained from the FACS/SPR Pareto-diversity shortlist. This preserves disagreement between the two signals instead of treating either score as a measured affinity. See the [ranking report](../results/spr-expanded-ranking/report.md), [ranking code](../scripts/rank_spr_candidates.py), [expanded screen](../scripts/expand_spr_screen.py) and [full-pool validation](research/abguard-validation.md).

## De novo: generate, score and filter

Apheris RFD4-Proteina 0.3.0 generated two arms of 32 new VH/VL pairs in the trastuzumab HER2 binding-site context. The arms used the antibody and HER2-novel checkpoints with `sequence_tsr_k=1`, seeds 20261118 and 20261119. Exact checkpoint names, image digests, input hashes and archived outputs are in the [generation audit](research/denovo-scale-results.md). Because seeds and checkpoints differ, arm yields are descriptive rather than a controlled checkpoint comparison.

Our local paired-sequence checks retained 43/64 designs. Installed abguard accepted 42 of those 43. Across the complete, unfiltered arms it accepted 21/32 and 26/32; these are different denominators. Eligibility checks include canonical amino acids, length, cysteines, antibody numbering and novelty. The installed tool uses its bundled reference of 610,109 CDRH3s. Its internal `novel` and `de_novo` classes are separate from the CSV's declared competition track.

The ranked draft excludes a declared-track/tool-class mismatch. It selects 20 de novo rows using fewer severe motif warnings, greedy paired-sequence diversity, nearest-reference CDRH3 distance and a deterministic ID tie break. These are selection heuristics. We do not transfer the trastuzumab SPR ranker to redesigned de novo frameworks. See the [selection policy](research/ranked-draft-policy.md), [generation audit code](../scripts/audit_denovo_scale.py) and [ranking/export code](../scripts/prepare_ranked_draft.py).

## Evaluation before laboratory testing

Abguard 0.1.0 with immunum 2.0.0 accepted all 30 exact ranked rows. The earlier local preflight uses immunum 2.1.0; these are separately recorded environments. Validator acceptance establishes eligibility under that installed build, not binding or final organizer acceptance.

We investigated structural predictions using OpenDDE, Protenix and RF3, including repeated WT and strong/weak SPR controls. The latest complete-template WT control passed its geometry criterion in only 1/3 seeds without MSA and 0/3 with target MSA, below the required 2/3. Candidate structural screening for these setups was withheld. ipTM is retained as structural confidence, not affinity. [Structural evidence](research/template-controls-results.md).

Thermal data audits and a ThermoMPNN-D folding proxy pilot did not establish a transferable candidate Tm ranking. No thermal score was used to claim improved stability. [Thermal evidence](research/thermal-assessment.md).

Laboratory SPR, competition against trastuzumab and direct thermal measurements remain necessary. No new candidate has demonstrated these improvements. Model agreement, sequence eligibility and laboratory success are distinct outcomes.

## Reproduce the selection from recorded evidence

Follow the lightweight replay in the [README](../README.md#reproduce-the-final-selection), or use the full development environment below in a fresh clone with `uv` installed and a local `.venv`. The bootstrap refuses a symlinked `.venv`; an existing shared environment should be activated directly instead. From the repository root:

```bash
bash scripts/bootstrap.sh
source .venv/bin/activate
PYTHONPATH=.:src pytest -q tests
python -m scripts.prepare_ranked_draft --output /tmp/reviewed-ranked30.csv
```

The last command reconstructs the ranked draft from archived local evidence without submitting anything. Compare it to the committed ranked CSV. [Export and validation instructions](../output/submission-draft/README.md) explain the six-column schema and the installed validator's handling of `track`.

Replaying a recorded selection is different from retraining models or regenerating structures. Full training requires the organizer datasets and the appropriate compute environment; new provider predictions require access credentials and may incur cost. Follow the per-experiment commands, input hashes, versions and seeds linked from the experiment index. [Setup](SETUP.md) describes access configuration; secrets are excluded from Git. The [final presentation README](../output/selection-slide/README.md) identifies the preserved manual edit. Earlier [presentation builders](../output/pipeline-draft/README.md) and [PDF builders](../output/pdf/README.md) produce historical explanatory artifacts. [Data retention](DATA_RETENTION.md) describes what was preserved and how to restore raw evidence.

The final team-facing filename is `epitope-hunters.csv`. The older `team-3` filenames are retained as frozen source evidence. Preparing this handoff does not submit designs; no submission receipt is recorded here.
