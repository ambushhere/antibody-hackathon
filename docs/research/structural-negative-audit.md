# Structural negative audit

Date: 2026-11-14

## Finding

The supplied Champloo SPR supplement does not substantiate a set of 454 measured true negatives or a `21/454` AF3 threshold result. In the two named CSVs there are 29 unique VHH–target pairs in the pairwise table and 12 in the partial cross-reactivity table. Neither table has an explicitly confirmed nonbinder label. The pairwise table has 12 `no reportable fit` results, 7 `not expressed` rows and 5 `antigen not produced` rows. The source README explicitly warns that no reportable affinity is not confirmed nonbinding.

The pairwise CSV contains three numeric `af3_iptm` values, all on cognate pairs: CD4 0.12, leptin receptor 0.15 and CNPase 0.46. Thus the exact in-table threshold count is **0/3 numeric values at or above 0.7**, but these are positive controls, not negatives. The README separately reports that AF3 ipTM was 0.75–0.90 for **11/12** non-cognate pairs in Supplementary Figure S1B. Per-pair non-cognate values are absent from the CSV, so this aggregate cannot be joined to targets or independently counted from the table. It is not a false-positive rate because the pairs had no reportable affinity rather than confirmed negative binding calls.

## Assay fields and counts

In `spr_table_s1_pairwise.csv`, `pair_type` (`cognate` or `non-cognate`) describes pair identity, not the SPR result. The outcome fields are `spr_outcome`, `outcome_class` and, where fitted, `kd_nM`/`kd_method`. Across 29 rows, the measured or disposition classes are 4 reported KDs, 1 response without fitted KD, 12 no reportable fits, 7 not expressed and 5 antigens not produced. CSV `pair_type` counts are 6 cognate and 23 non-cognate; the README's original selected cohort summary is 5 cognate plus 12 non-cognate among 17 reported rows, with one provider-added pair also present in the 29-row CSV.

In `spr_crossreactivity_14target.csv`, the observation fields are `response_RU`, `surface_load_RU` and `call`. This extracted table is partial: 12 unique VHH–target pairs across 10 targets, not a full 14-target matrix. It records 1 confirmed binder, 2 off-target binders, 2 weak under-loaded responses, 5 weak responses with low-loading background, 1 silent response and 1 response not quantified in the text. The two off-target binder calls are positive binding observations. Weak, silent and unquantified calls do not provide a confirmed negative set under the source's stated low-loading caveat.

The two supplementary tables are separate experiments and are not pooled as one denominator. No valid per-target or pooled false-positive fraction can be computed, so Wilson 95% intervals are **not estimable** from these sources. There are zero matched confirmed-negative pairs with a row-level confidence score in the supplied tables.

## Match to current HER2 scores

`results/experiment-20261114/structure_comparison.csv` contains 10 rows for six HER2 Fv query IDs. It has OpenFold3 and Protenix v2 scores, including `max_antibody_target_pair_iptm`; it has no AF3 rows. None of its query IDs or model outputs are the Champloo VHH–antigen pair rows. The SPR tables have no matched OpenFold3 or Protenix scores. Therefore these records cannot benchmark or calibrate the current HER2 Fv refold scores, including a Protenix or OpenFold threshold.

The narrow supported claim is that high AF3 interface confidence can occur for non-cognate antibody–antigen pairs that have no reportable affinity in this selected dataset: the supplied README reports 11/12 at 0.75–0.90. This supports caution against treating structural confidence alone as evidence of binding specificity. It does not establish a false-positive rate, a universal threshold, or performance on HER2 Fv refolds.

## Data not available for this audit

- A row-level dataset for the claimed 454 measured negatives and the claimed 21 high-AF3 cases, including the definition of a confirmed negative and the exact pair identifiers.
- Per-pair AF3 values for the 12 non-cognate rows in `spr_table_s1_pairwise.csv` and a target-level mapping for the reported 11/12 aggregate.
- A complete 14-target cross-reactivity response matrix.
- Matched Protenix v2 or OpenFold3 confidence scores for the SPR pairs.
- A model-score join between either supplement table and the current HER2 Fv score CSV.

The VM dataset README lists a 1.3 GB AF3 real-pair archive and a 152 GB AF3 shuffled archive. The archives were not opened or extracted. The 152 GB archive was not present in the top-level directory listing. No archive was downloaded.

## Reproducibility and source hashes

The CSV counts were computed on the VM after activating `/home/hackathon/almaz/work/.venv`. Only aggregate counts and source hashes were returned; no raw records are included here. Machine-readable summaries and hashes are in [`results/negative-audit/`](../../results/negative-audit/).

| Source | SHA-256 |
|---|---|
| `spr_table_s1_pairwise.csv` | `6d28ccec9ed5d4f15cd88e4edd063bfcaaa1d80d48323447441ad03dd35e00f1` |
| `spr_crossreactivity_14target.csv` | `5ac9147ea62e9b84e641887e183a0d22bc592a97d6b13b718e1b8791a5ea33fb` |
| supplementary `README.md` | `6f26633b36da366f492e0d896d1853e7ee5d54aaa9a42e9fcf05c7a688487bc7` |
| VM dataset `README.txt` | `a7c1c811708c90d55c08b187e3d11b51abde24cea244dac4506dd993118f601d` |
| `results/experiment-20261114/structure_comparison.csv` | `a2819c7e652de3b196a9ddfc4db0586835910a1703082576d8babc0fe0a3a86d` |
