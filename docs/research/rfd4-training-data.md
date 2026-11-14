# RFD4 training data inspection

14 November 2026. Downloaded two metadata tables totalling **46,415,290 bytes (44.27 MiB)**, not the 2.13 TB structure corpus. The local cache is `data/external/rfd4-metadata/`. The committed [audit](../../results/rfd4-training-metadata/audit.json) pins upstream revision, file sizes, hashes, URLs, schemas and observed counts. The raw tables remain a reproducible local cache and are not added to Git.

Source: [NVIDIA RFD4-Proteina dataset at aa86df1](https://huggingface.co/datasets/nvidia/RFD4-Proteina/tree/aa86df1eb9e3529a0f39279f216dce3a95bb3dc2). The upstream card describes public release training recipes, not the undisclosed row selection of the Apheris HER2 fine-tune.

| Table | Download | Observed rows | Use |
| --- | ---: | ---: | --- |
| PDB immune interfaces | 20.93 MB | 81,895 antibody and 9,526 TCR entries | Inspect experimental structural contexts, cluster IDs and training overlap |
| ESMFold2 antibody-antigen distillation metadata | 25.49 MB | 75,576 | Inspect model confidence, interface contacts and possible teacher features |

These are sampling rows, not counts of unique independent complexes. We did not apply upstream recipe filters. The card reports 72,562 antibody rows and 6,189 distilled rows passing source filters; that is a different denominator. The distilled component's structure blobs occupy about 9.7 GB, and the PDB antibody view shares a roughly 182.3 GB structure store. Neither was downloaded.

## Findings that change interpretation

Both **1N8Z and 6OGE appear in the PDB immune sampling table**: two and four rows respectively. The six rows pass the two reported bond and ligand identity flags, but this audit does not establish every recipe filter or actual checkpoint exposure. We retain them as known-complex reproduction controls. They cannot substantiate an independence claim for a model that could have trained on these structures. The previous target-template proposal remains useful to diagnose inference geometry, not to establish novel-target generalization.

The distillation table contains only **13,652 paired Fv entries**, versus 38,631 single-heavy, 23,259 single-light, two scFv and 32 unannotated entries. These are provider format labels, not a reannotation. Of all rows, 71,091 have source `patents`. This is not an interchangeable set of 75,576 paired-antibody laboratory measurements. It does provide useful candidate features such as antibody-antigen ipTM, directional PAE, CDR contact fraction and clash counts. These can train a cheap surrogate of the structural teacher after sequence/target separation; they cannot supply new measured affinity or thermal labels.

## Decision

Keep the compact metadata for overlap checks and selection of a future targeted structural subset. Do not expand today's training with unfiltered rows. First qualify our HER2 structural setup on WT, and audit the external measured SPR cohort described in [benchmark evidence](benchmark-evidence.md). Any future teacher-distillation experiment should retain a separate real-measurement test and report format/target shift. No new model was trained on these tables and no candidate was promoted from them.

## Replay

```bash
source .venv/bin/activate
python -m pip install pyarrow
# Download the two exact URLs in results/rfd4-training-metadata/audit.json
# to data/external/rfd4-metadata/ using the manifest keys as filenames.
python -m scripts.audit_rfd4_training_metadata
```

The auditor checks exact byte counts and SHA-256 before reading selected columns. A second root run produced byte-identical audit JSON. The upstream annotations and packaging use CC BY 4.0; underlying structures retain their upstream terms. Attribute NVIDIA/RFD4-Proteina and any selected underlying source separately when reusing records.
