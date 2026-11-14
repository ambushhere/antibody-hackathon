# HER2 FACS supervised benchmark

Target is the published high FACS class versus mid and low. Models use CDRH3 only.
Seed 0. Selected pairwise by validation average precision.

| Split | Rows | CDRH3s | Clonotypes | High | Mid | Low |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| train | 367,044 | 367,044 | 13,829 | 120,505 | 143,278 | 103,261 |
| val | 78,651 | 78,651 | 251 | 25,822 | 18,629 | 34,200 |
| test | 78,651 | 78,651 | 706 | 25,822 | 28,104 | 24,725 |

All split pairs had zero CDRH3, clonotype centre and sequence-to-centre overlap.

The optional FACS mapping agreed on all 524,346 split rows.

| Split | Model | ROC AUC | PR AUC (AP) | Balanced accuracy | Brier | Precision top 1% | Enrichment top 1% | Precision top 5% | Enrichment top 5% |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| validation | dummy | 0.5000 | 0.3283 | 0.5000 | 0.2205 | 0.3283 | 1.0000 | 0.3283 | 1.0000 |
| validation | position | 0.9644 | 0.9365 | 0.9000 | 0.0759 | 1.0000 | 3.0459 | 0.9962 | 3.0343 |
| validation | pairwise | 0.9929 | 0.9868 | 0.9614 | 0.0299 | 1.0000 | 3.0459 | 0.9995 | 3.0443 |
| test | dummy | 0.5000 | 0.3283 | 0.5000 | 0.2205 | 0.3283 | 1.0000 | 0.3283 | 1.0000 |
| test | position | 0.9440 | 0.9028 | 0.8566 | 0.0904 | 0.9962 | 3.0343 | 0.9939 | 3.0273 |
| test | pairwise | 0.9859 | 0.9737 | 0.9404 | 0.0429 | 1.0000 | 3.0459 | 0.9997 | 3.0451 |

Balanced accuracy uses a probability threshold of 0.5. Precision at a tied cutoff is the expected precision over tied rows. Enrichment is relative to split prevalence.

## SPR evaluation

Correlations compare predicted high FACS probability to the median of −log10(KD in nM) for quantified replicates of each eligible CDRH3.
Eligibility requires the trastuzumab VH framework and VL, a ten-residue CDRH3 distinct from wild type and no exact train or validation CDRH3 match.

| Cohort | CDRH3s | Spearman ρ | Pearson r |
| --- | ---: | ---: | ---: |
| Eligible SPR | 158 | 0.3473 | 0.2297 |
| Also excludes FACS test CDRH3s | 152 | 0.3929 | 0.2560 |

Of 768 assay rows, 190 had quantified, finite positive KD.
Exclusions from those rows were 0 other backgrounds, 8 wild type, 14 train matches, 7 validation matches and 3 unsupported loop lengths.
Aggregation collapsed 0 replicate rows.
These correlations describe the cohort with a reportable KD.

## Input hashes

- train b639077ff116daa8b8c3a2122b6ed9291361892469f6237963967ca0772cb262 (/tmp/her2-benchmark.NjMJ71m1/train.csv)
- val ea5e2e0908f018f579e476cef4452f9de11ba9144513bffa783588d49123f5d1 (/tmp/her2-benchmark.NjMJ71m1/val.csv)
- test b350c10c7df4eb5e1b966f76909a805fc7d8011486f1ddb5dd1d7a6e19556a06 (/tmp/her2-benchmark.NjMJ71m1/test.csv)
- facs 7ff5a95bcc1942bb6e34b1e4ed6c02a9c177ee30ca160440f8abe82ed787735e (/tmp/her2-benchmark.NjMJ71m1/facs.csv)
- spr 4d6f53a6183f377670e9577b3233a7c4c14888965ac63824863219d13424fa9b (/tmp/her2-benchmark.NjMJ71m1/spr.csv)
