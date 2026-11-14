# External HER2 SPR transfer results

The fixed models scored 36 sequence-compatible controls from the Absci human HER2 SPR panel. Ridge ranked the 25 quantified control binders at Spearman ρ 0.492 versus 0.342 for FACS. The paired group bootstrap interval for the difference was [−0.294, 0.647]. Binding detection gave ROC AUC 0.833 versus 0.756, with paired difference interval [−0.266, 0.317]. The estimates are too imprecise to establish an advantage for either model. This retrospective result is conditional on the study-level trastuzumab Fab scaffold. [Frozen protocol](external-spr-protocol.md) · [Metrics and input hashes](../../results/external-spr/metrics.json)

## Input audit

The published reference VH and VL start with the exact 120- and 107-residue benchmark variable domains. The columns give WT IMGT HCDR1 `GFNIKDTY` and HCDR2 `IYPTNGYT`, at linear VH positions 26–33 and 51–58. The WT IMGT HCDR3 `SRWGGDGFYAMDY` spans linear VH positions 97–109. The fitted models vary positions 99–108, `WGGDGFYAMD` at WT. Only 13-residue loops retaining `SR` and `Y` at the boundaries enter the ten-position projection. Individual alignments, row-level paired chains, assay construct IDs and run membership are unavailable. The selected binder CSV also lacks per-row HCDR1 and HCDR2. [Absci reference and CSV documentation](https://github.com/AbSciBio/unlocking-de-novo-antibody-design/blob/cd2bd4b6932b5301d435555dc22084d407e4b92b/README.md)

Absci describes cloning HCDR3 or HCDR1–HCDR3 inserts into trastuzumab Fab plasmids and measuring binding of captured Fab samples against human HER2 ECD on Carterra LSA. Chinery's FACS source used scFv yeast display. Its SPR source used Carterra LSA on constructs expressed in HEK293 or CHO and purified using MabSelect SuRe. The SPR methods do not specify a row-level constant region, linker or valency. Across-study format equivalence and full per-antibody scaffold identity remain unverified. The excluded WT illustrates the assay difference. The Absci control reads 1.94 nM, while eight internal reads span 0.0285–0.178 nM. [Absci methods, Cloning and SPR](https://absci.gcs-web.com/static-files/a53c2b92-0e77-4899-a52c-4485be691cf6) · [Chinery v2, Methods](https://storage.prod.researchhub.com/uploads/papers/2026/01/14/2024.03.26.586756v2.full.pdf) · [Internal SPR report](../../results/spr-benchmark/report.md)

| Sequential audit | Controls | Selected binders |
| --- | ---: | ---: |
| Published rows | 1,855 | 422 |
| Observed HCDR1 or HCDR2 differs from WT | 247 | Row annotations unavailable |
| Unsupported HCDR3 length among remaining rows | 589 | 221 |
| Changed HCDR3 anchors among remaining rows | 981 | 197 |
| Canonical 13-residue HCDR3 with WT anchors after the HCDR1/HCDR2 gate | 38 | 4 |
| WT excluded, also present in the full internal SPR file | 1 | 1 |
| Other exact internal match excluded | 1 | 0 |
| Other match to a compatible control excluded | 0 | 3 |
| Unique windows scored | 36 | 0 |

The other exact internal match was `WGGPGLYAMD`, in FACS validation. The three selected variants shared both HCDR3 and original KD with eligible controls. They count once in the control results. Of all 1,855 control rows, 758 have a positive binder call and a finite positive KD. All 1,097 reported negatives have blank KD. The scored controls had 25 binders and 11 negatives after gating. There were no duplicate or conflicting rows among these 36 windows. For repeated rows, the evaluator preserves source IDs, counts agreeing labels once, deduplicates identical KD values and takes the median of transformed quantitative values. It excludes conflicting binder calls from both endpoints. Empty or censored KD contributes no affinity target. [Published-row audit](../../results/external-spr/source_rows.csv) · [Scored windows](../../results/external-spr/window_scores.csv)

The follow-up file has 13 rows, including WT. All 13 HCDR3s link to the selected binder file. Original KD matches all 13 after rounding the selected file values to the follow-up file's displayed decimal precision, with zero exact numeric matches. The original and follow-up Fab and mAb values remain distinct. [Follow-up provenance audit](../../results/external-spr/followup_audit.csv)

Nearest Levenshtein distances were computed against 158 actual ridge fit windows, 367,044 FACS train windows and 78,651 FACS validation windows. Validation guided FACS model selection. Exact membership was also checked against 78,651 FACS test windows and all 754 unique internal SPR HCDR3s, including rows without numeric KD. The 36 controls had nearest distances of 1–6, 1–4 and 1–5 edits to ridge fit, FACS train and FACS validation respectively. The 13 controls in the distance sensitivity were at least three edits from each of those three reference sets. [Distances and membership by row](../../results/external-spr/source_rows.csv)

## Binding detection

The endpoint is Absci's published SPR `Binder` call. The detection panel has 36 controls in 19 groups, with at most seven rows per group and a positive fraction of 0.694. The distant panel has 13 controls in 10 groups, with at most two rows per group and a positive fraction of 0.538. These fractions give the corresponding average precision baselines. Both scores are used for ranking without a decision threshold.

| Panel | Score | ROC AUC [95% interval] | Average precision [95% interval] |
| --- | --- | ---: | ---: |
| All 36 | Ridge | 0.833 [0.578, 0.976] | 0.941 [0.770, 0.995] |
| All 36 | FACS | 0.756 [0.497, 0.996] | 0.860 [0.542, 0.999] |
| All 36 | Negative distance to WT | 0.689 [0.394, 0.933] | 0.832 [0.465, 0.987] |
| All 36 | Ridge minus FACS | 0.076 [−0.266, 0.317] | 0.081 [−0.081, 0.304] |
| Distant 13 | Ridge | 0.548 [0.133, 0.944] | 0.694 [0.367, 0.977] |
| Distant 13 | FACS | 0.738 [0.374, 1.000] | 0.852 [0.510, 1.000] |
| Distant 13 | Negative distance to WT | 0.345 [0.060, 0.813] | 0.501 [0.239, 0.926] |
| Distant 13 | Ridge minus FACS | −0.190 [−0.825, 0.500] | −0.158 [−0.584, 0.343] |

Across leave-one-group-out runs, the ROC AUC difference ranged from −0.074 to 0.138 in the full panel and −0.393 to 0.067 in the distant panel. The average precision differences ranged from −0.021 to 0.121 and −0.306 to 0.007 respectively.

## Affinity ranking

The 25 quantified controls span 1.07–1,052.08 nM across 13 groups, with the largest group containing seven. The distant subset contains seven quantified controls in six groups, with a largest group of two. Its KD spans 79.66–1,052.08 nM, so its estimates concern the weaker portion of the observed range.

| Panel | Score | Spearman ρ [95% interval] | Kendall τ-b [95% interval] |
| --- | --- | ---: | ---: |
| All 25 | Ridge | 0.492 [−0.263, 0.788] | 0.367 [−0.203, 0.641] |
| All 25 | FACS | 0.342 [−0.399, 0.563] | 0.233 [−0.292, 0.420] |
| All 25 | Negative distance to WT | 0.245 [−0.640, 0.622] | 0.185 [−0.490, 0.496] |
| All 25 | Ridge minus FACS | 0.150 [−0.294, 0.647] | 0.133 [−0.224, 0.483] |
| Distant 7 | Ridge | 0.643 [−0.800, 1.000] | 0.429 [−0.636, 1.000] |
| Distant 7 | FACS | −0.393 [−0.946, 0.800] | −0.238 [−0.862, 0.636] |
| Distant 7 | Negative distance to WT | 0.208 [−0.857, 0.962] | 0.218 [−0.745, 0.918] |
| Distant 7 | Ridge minus FACS | 1.036 [−0.800, 1.892] | 0.667 [−0.769, 1.724] |

The leave-one-group-out Spearman differences were 0.079–0.286 for all 25 and 0.400–1.314 for the distant seven. The corresponding Kendall differences were 0.095–0.239 and 0.200–0.933.

| Ridge error in −log10(KD nM) | All 25 [95% interval] | Distant 7 [95% interval] |
| --- | ---: | ---: |
| MAE | 0.765 [0.535, 1.057] | 0.588 [0.376, 0.823] |
| RMSE | 0.888 [0.661, 1.115] | 0.654 [0.430, 0.853] |
| Mean predicted minus measured | 0.686 [0.411, 1.032] | 0.555 [0.283, 0.818] |

The ridge predicts stronger binding than measured on average in this panel. The positive signed error is consistent with an assay or selection shift, with its cause unresolved. A numerical cross-assay correction was outside the protocol.

## Interpretation

Intervals use paired percentile resampling of intact Levenshtein ≤2 connected components. Group labels are preserved across endpoints and the distance sensitivity. One of 2,000 full detection draws and seven distant detection draws contained only one class. Eight distant KD draws had constant WT-distance scores. The other reported metrics had 2,000 valid draws. These intervals describe variation across the observed sequence groups with fixed scores. Shared study procedures, unavailable assay batches and possible higher-level dependencies remain unresolved.

The control panel was assembled for the source study and the quantitative endpoint is limited to observed binders. The window gate removes most of the source sequences. These are retrospective estimates within one antibody family and target. The results support no automatic promotion and provide no measured evidence about competition or thermal stability.

## Reproduction

Use the benchmark's original `train.csv`, `val.csv`, `test.csv` and `spr.csv`, as in the [internal SPR instructions](../BENCHMARK.md#cross-assay-follow-up). The evaluator verifies their hashes and rebuilds the 158 fitted IDs to check artifact provenance. The model parameters are read from repository JSON artifacts. The evaluator fits no model.

```sh
mkdir -p data/raw/external-spr
base=https://raw.githubusercontent.com/AbSciBio/unlocking-de-novo-antibody-design/cd2bd4b6932b5301d435555dc22084d407e4b92b
for file in README.md LICENSE spr-controls.csv zero-shot-binders.csv functionality-developability-cross_reactivity-data.csv; do
  curl -fsSL "$base/$file" -o "data/raw/external-spr/$file" || exit 1
done

source .venv/bin/activate
PYTHONDONTWRITEBYTECODE=1 python -m scripts.benchmark_external_spr \
  --data-dir /path/to/verified/internal/benchmark-inputs
PYTHONPATH=src:. PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider tests/test_external_spr.py
```

External measurement and sequence data are attributed to Absci Corporation (2023) under the source Clear BSD license. The observation audit records CSV row numbers counting the header as row 1. [Pinned source license](https://github.com/AbSciBio/unlocking-de-novo-antibody-design/blob/cd2bd4b6932b5301d435555dc22084d407e4b92b/LICENSE)
