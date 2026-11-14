# Candidate evaluation board

Separate evidence axes. No weighted universal score, inferred binding probability or automatic submission.

Failure of a control removes support from that setup; it is not proof that every candidate is a nonbinder.

The frozen optimisation shortlist is unchanged. All candidates remain experimentally unmeasured.

| Candidate | Track | Local gate | Structural runs | Status |
| --- | --- | --- | ---: | --- |
| opt_WGAPGFYLMT | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_YTGDGFPAMN | optimization | pass | 7 | provisional research priority, not submission approval |
| opt_WEKGGFYRMD | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_WYGPGFPVMD | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_FTSDRFYAMD | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_WTAVGFFAMD | optimization | pass | 7 | provisional research priority, not submission approval |
| opt_YQADGFYHMD | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_WVGGWFYLMD | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_WSGRRFYHMD | optimization | pass | 0 | provisional research priority, not submission approval |
| opt_YTADGVYHLK | optimization | pass | 0 | provisional research priority, not submission approval |
| denovo_fullseq_k1_1 | de_novo | see design-specific validation | 7 | exploratory design |

## WT control diagnostics

Contact recovery is a diagnostic, not full pose qualification. Read target distortion and interface RMSD alongside it.

| Model | Target context | Seed | WT contacts out of 19 |
| --- | --- | ---: | ---: |
| opendde-abag | domain-iv | 42 | 2 |
| opendde-abag | domain-iv | 43 | 0 |
| opendde-abag | domain-iv | 44 | 5 |
| openfold3_of3ob | initial_split | 42 | not in this audit |
| protenix-v2 | domain-iv | 42 | 3 |
| protenix-v2 | ecd | 42 | 0 |
| protenix-v2 | initial_split | 42 | not in this audit |
| protenix-v2 | split | 42 | 7 |
| rf3 | domain-iv | 42 | 14 |

Full numeric geometry and confidence are in `structures.csv`, measured controls in `controls.csv` and separate candidate evidence fields in `candidates.csv`.
The file manifest records source hashes. Rebuild after a new audited model arm; raw provider outputs do not enter directly.
