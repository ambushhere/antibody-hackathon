# HER2 candidate screen with target MSA results

The target MSA control audit validated all 24 poses. Zero of three confidence-selected WT poses cleared the frozen reproduction gate, which requires at least two passing seeds. The candidate screen was not run. The derived report records `not_run_admission_failed`, with 57 unknown rows and no ranks. [Control audit](../../results/template-controls/target-msa/audit.json) · [Candidate report](../../results/candidate-msa-screen/report.md)

| Seed | Selected sample | Mean pair ipTM | WT contacts (of 19) | Target fit RMSD (Å) | WT interface RMSD (Å) | Crop clash pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 42 | 2 | 0.2105 | 3 | 2.019 | 19.605 | 0 |
| 43 | 2 | 0.2780 | 4 | 0.932 | 21.368 | 0 |
| 44 | 1 | 0.2900 | 7 | 2.210 | 16.541 | 0 |

All three selected WT poses missed the contact and interface RMSD limits. Seed 43 met the target fit limit and all three met the in-crop clash criterion. The [separate no MSA result](candidate-screen-results.md) remains at one passing seed. These controls evaluate reproduction for the two tested configurations. The [prepared fallback](candidate-msa-protocol.md) remains available for audit.
