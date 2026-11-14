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

## Novelty-distance diagnostic

A fixed random sample of 2,000 held-out test sequences was compared against every training CDRH3 using exact Levenshtein distance. This is a post-benchmark diagnostic and did not select or retrain the model. At distance >=2, 867 examples remain and average precision is 0.9526. At distance >=3, only 109 examples remain, including 9 high binders; average precision is 0.6532 versus prevalence 0.0826. This smaller cohort is closer to the novelty challenge but is too small for a precise performance claim. Its different prevalence also means average precision should not be compared without the baseline. The high overall test result alone does not establish prospective quality for novel designs.
