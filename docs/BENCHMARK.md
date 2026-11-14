# HER2 experiment benchmark

This is an experiment, not a competition submission. No new candidate has measured affinity, competition or thermal stability. The provisional list is an experiment queue, not a set of approved lab designs.

## Reproduce

Run `scripts/bootstrap.sh` on a fresh checkout and activate `.venv`. Obtain the organizer dataset from the team VM. The archive and CSV files are external inputs; hashes are recorded in the report.

```bash
python src/her2_benchmark.py benchmark \
  --archive /path/to/zenodo-10831512-trastuzumab-ngs/data.zip \
  --facs /path/to/tz_her2_facs_3classes_524346.csv \
  --spr /path/to/trastuzumab_her2_spr_all_768.csv \
  --output-dir results/my-benchmark --seed 0
```

The model is saved as JSON numeric parameters, not executable pickle. Prediction requires a CSV containing `cdrh3` and no existing `p_facs_high` column:

```bash
python src/her2_benchmark.py predict --model results/my-benchmark/selected_model.json \
  --candidates /path/to/candidates.csv --output results/my-benchmark/scores.csv
```

Extract the archive's published `data/affinity_data/her2/her2_aff_large/processed/remove_overlap/clonotype_split/0.7_0.15_0.15/train.csv` for generation. Only its high labels influence proposals; the full provided tables are used solely for novelty exclusion.

```bash
python src/her2_candidates.py generate --facs /path/to/facs.csv --spr /path/to/spr.csv \
  --train /path/to/train.csv --pool-size 10000 --seed 42 --export-vh-vl --output candidates.csv
PYTHONPATH=src python src/her2_preflight.py --candidates scored-shortlist.csv \
  --facs /path/to/facs.csv --spr /path/to/spr.csv --output preflight.json
```

The generator implements the slides' permissive >=2 edit rule. Preflight reports the conservative >=3 rule plus identity under every supported CDR scheme; do not mistake generated proposals for strict passes. Use `her2_candidates select` with a `cdrh3,predicted_score` CSV for explicit score-plus-diversity selection. The benchmark output column `p_facs_high` must be renamed to `predicted_score` for that generic selection command.

## Findings

The pairwise logistic model was selected using validation average precision. The held-out FACS test gives ROC AUC 0.9859 and average precision 0.9737. Position-only logistic reaches 0.9440 and 0.9028. All 524,346 rows match organizer labels, with no exact sequence or published clonotype overlap among splits. Similarity below exact identity is examined separately, not ruled out by this audit.

On 158 eligible quantified SPR variants, the selected FACS score correlates with -log10(KD nM) at Spearman 0.3473. This modest association is a limit on affinity-ranking claims. Missing KD was not converted into negative labels. FACS classification is not a surrogate for thermal stability.

The root reproduced the full benchmark and generated 10,000 candidates. A score/diversity subset with >=4 WT changes contains 500 statistical candidates; with four RFD4 proposals, 88 of 504 pass conservative local preflight. A provisional 20-entry experiment queue is recorded. This is a selected-subset yield, not the strict-pass rate of all 10,000.

Four RFD4 HER2-novel designs and eight exploratory antibody designs were generated, the latter split between fixed-framework CDR redesign and full-sequence redesign retaining canonical cysteines. All eight exploratory designs are numberable but database novelty is unverified. They remain close to known antibody sequences; do not claim a successful de novo entry.

OpenFold and Protenix evaluated WT, a prespecified FACS-low control and two RFD4 variants with matched target crop, template and seed. Protenix also evaluated two statistical candidates. Overall ipTM hides large VH–VL confidence contributions. One statistical candidate has VH–VL ipTM 0.182 compared with WT 0.948 in Protenix, despite its high FACS classifier score. It is flagged for structural review. This is a single prediction, not experimental rejection.

Target numbering is explicitly restored: predictions renumber target fragments from 1, while the source uses C489–580 and D591–607. Comparing unremapped residue IDs would produce a false epitope-overlap result. Contact reports use the checked mapping in configs/her2_target_residue_map.json.

## Limits and next experiment

Thermal improvement is not assessed. Binding ddG and sequence liability descriptors cannot fill that gap. Database-wide de novo novelty and official abguard behavior remain unknown. Candidate selection must preserve these missing fields.

Before scaling generation, test whether candidate rankings survive independent seeds, a target crop with a complete domain and appropriate control baselines. Investigate the affinity relationship on disjoint SPR groups and validate a thermal predictor against suitable measured examples. No model is entitled to a confidence threshold merely because another model uses it.

## Cross-assay follow-up

The direct SPR benchmark uses nested GroupKFold, with connected components at Levenshtein distance <=2 kept intact. On 158 primary variants without exact FACS train/validation overlap, position ridge gives out-of-fold Spearman 0.613 versus 0.347 for the frozen FACS score. The paired group bootstrap difference is 0.265, interval [0.110, 0.515]. These intervals resample fixed predictions and do not include all model-selection uncertainty. Ridge RMSE is 0.741 in -log10(KD nM), a substantial residual error. The FACS-plus-residual model does not improve on direct position ridge. See [report](../results/spr-benchmark/report.md), numeric parameters and input hashes.

This is a retrospective ranking result conditional on reportable KD, not a binding success classifier or proof of improvement over WT. The table supplies no assay batch identifier. Eight WT values span 0.0285–0.178 nM; the observed range is not an equivalence margin. Only one primary-cohort variant beats its lowest value. These limitations prevent precise claims about selecting improved antibodies.

```bash
python -m src.her2_spr_benchmark --data-dir /path/to/extracted-csv-directory \
  --facs-model results/experiment-20261114/selected_model.json \
  --output results/spr-benchmark
```

The input directory contains train.csv, val.csv, test.csv, facs.csv and spr.csv with the hashes in the report. The root independently reran this benchmark and its group isolation tests.

A frozen-prediction [distance diagnostic](../results/spr-distance-diagnostic/report.json) leaves 75 variants in 67 groups at edit distance >=3 from FACS train/validation and their respective outer-fold SPR training set. On this subset, SPR ridge rho is 0.425 versus FACS 0.197. The paired difference interval is [-0.064, 0.512], including zero. The direction is encouraging but does not establish superiority in the novelty regime. No model was tuned on this diagnostic.

## De novo diversity and independent-pose diagnostic

Four additional full-sequence designs use `sequence_tsr_k=1`, with the same input, seed and antibody checkpoint as the default k=4 pilot. The organizer parameter guide documents default 4 and says lowering to 1 removes sequence sharpening. All four k=1 outputs have IMGT-position paired identity below 70% to trastuzumab and adalimumab, versus zero of eight prior exploratory designs. The denominator is the union of numbered positions and the numerator is identical residues. This is a local two-reference diagnostic, not the official database alignment. One k=1 variant introduces a third VH cysteine and fails the conservative cysteine check. All four remain within the paired-domain length limit; adding a linker must be checked separately.

Protenix evaluated all eight full-sequence designs plus WT and a prespecified FACS-low control with the same seed and a target-only template. No antibody coordinates were supplied. [Scores and geometry](../results/denovo-refold/) show the tradeoff: VH–VL ipTM ranges 0.891–0.928 for k=4 and 0.098–0.860 for k=1. One locally novel candidate has no predicted contacts on the WT footprint; another retains 68% of the cropped WT contacts but has low target-pair confidence. The extra-cysteine candidate also has VH–VL ipTM 0.098. None is established as a successful de novo binder.

The FACS-low control has maximum antibody–target pair ipTM 0.556, higher than WT at 0.382. This ordering disagrees with the FACS screen. The control has no quantified SPR result, so it cannot establish a KD ranking failure. Maximum pair ipTM remains unsupported for affinity prioritisation here; the next comparison uses quantified SPR controls. The target is still the organizer's two-fragment crop, one seed was used and the short fragment affects confidence. A complete-domain control experiment has higher priority than generating a much larger pool.

The standard refold preparation rejected full-sequence designs because conditioned residues differ. A simple atom-only target CIF also failed Protenix parsing. The final target template preserves polymer metadata and removes both antibody coordinates and antibody entities. The first failed attempt and cancelled intermediate attempt are retained in raw local provenance. The completed run is `denovo-target-only-protenix-v3`.

## Thermal and source audits

[Thermal assessment](research/thermal-assessment.md) identifies measured trastuzumab scFv controls and a small HER2 CDRH3 stress panel. Folding ΔΔG is distinct from binding ΔΔG and from a measured Tm. A paired-Fv proxy needs control validation before candidate prioritisation. NanoMelt applies to VHH, not to our paired VH/VL designs.

The [completed CPU thermal pilot](../results/thermal-pilot/README.md) uses pinned ThermoMPNN-D single and double heads on two paired Fv backbones. The root reproduced the output tables exactly. R66G has a discordant direction on both structures; H91Y moves from -0.927 to a near-zero +0.016 kcal/mol. Ranking five measured HER2 CDRH3 variants is no better than a mutation-count baseline. Thermal ranking remains unsupported for our candidates.

The [expanded optimisation screen](../results/spr-expanded-ranking/report.md) adds SPR priorities from the existing pool, avoiding the earlier restriction to FACS-selected candidates. Of 239 newly checked proposals, 36 pass conservative local gates; the combined checked set has 124 passes from 743. This is a selected-subset yield. The SPR model predicts WT at -0.591 score units despite the measured WT range corresponding to much larger -log10(KD). Thus its absolute KD calibration near WT is inadequate. Use relative ranking with explicit extrapolation warnings, not claims of measured or calibrated KD superiority.

JapanFold's current catalog identifies `af2ig` as AlphaFold2 **initial guess**, not an immunoglobulin-specific model. It requires a designed complex and is pose-conditioned. It must not be counted as an independent antibody cofolding vote. `opendde-abag` is an antibody checkpoint on the Protenix-v2 stack, so agreement with Protenix also has shared-model dependence. RF3 refers to RoseTTAFold3 structure prediction and is distinct from RFdiffusion3/4 design.

[Structural negative audit](research/structural-negative-audit.md) withdraws the premise that the inspected supplement provides 454 confirmed negative pairs suitable for calibrating our model. Missing KD, non-expression and absent antigen cannot be recoded as nonbinding.

## Visual explanation

[Ten-slide Russian explainer](../output/explainer/her2-intuition.pdf) explains the target, epitope, paratope, formats, objectives and evidence boundary. Rebuild with `python -m pip install -r requirements-explainer.txt` then `python scripts/build_explainer.py`. It also requires Poppler (`pdftoppm`) and macOS supplemental fonts or Linux DejaVu fonts. Figures are schematic, not atomic structures.

## Novelty-distance diagnostic

A fixed random sample of 2,000 held-out test sequences was compared against every training CDRH3 using exact Levenshtein distance. This is a post-benchmark diagnostic and did not select or retrain the model. At distance >=2, 867 examples remain and average precision is 0.9526. At distance >=3, only 109 examples remain, including 9 high binders; average precision is 0.6532 versus prevalence 0.0826. This smaller cohort is closer to the novelty challenge but is too small for a precise performance claim. Its different prevalence also means average precision should not be compared without the baseline. The high overall test result alone does not establish prospective quality for novel designs.
