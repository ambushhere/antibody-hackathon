# HER2 template ablation results

**Complete.** The audit verified 48/48 planned Protenix v2 no-MSA poses. It covered 24 arm, control and seed groups and retained both draws per completed group. [Per sample scores](../../results/template-controls/scores.csv), [Full audit](../../results/template-controls/audit.json), [Raw archive manifest](../../results/archives/template-controls-manifest.json)

The complete template fell below the two-seed threshold with 1/3 passing WT seeds (42).

The masked template fell below the two-seed threshold with 0/3 passing WT seeds.

The no-MSA candidate screen was not admitted. The conditional alternative requires the same WT gate on the separate target-MSA control. [Candidate fallback protocol](candidate-fallback-protocol.md)

## WT pose comparison

All WT draws are shown below. **Top 1** uses mean VH–HER2 and VL–HER2 pair ipTM within each arm and seed. Best reference geometry uses the lowest interface Cα RMSD in the 1N8Z target frame and is marked as a diagnostic. The operational rule uses thresholds of 13/19 contacts, 2 Å target RMSD, 4 Å interface RMSD and zero predicted antibody–target heavy atom pairs below 1.5 Å. The 6OGE ECD overlay is a separate steric diagnostic. [Protocol](template-controls-protocol.md), [Operational rule](structure-admission-rule.md)

| Template | Seed | Draw | Selection | Mean pair ipTM | Contacts of 19 | Target RMSD (Å) | Interface RMSD (Å) | Clashes | ECD overlay clashes | Gate |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| complete | 42 | 1 | **Top 1**, best reference geometry | 0.2985 | 14 | 1.30 | 3.22 | 0 | 0 | pass |
| complete | 42 | 2 |  | 0.2980 | 11 | 1.28 | 6.09 | 0 | 2 | fail |
| complete | 43 | 1 | **Top 1** | 0.3715 | 11 | 1.26 | 7.06 | 0 | 14 | fail |
| complete | 43 | 2 | best reference geometry | 0.3630 | 12 | 1.07 | 6.85 | 0 | 12 | fail |
| complete | 44 | 1 |  | 0.2495 | 9 | 6.42 | 17.18 | 0 | 66 | fail |
| complete | 44 | 2 | **Top 1**, best reference geometry | 0.3040 | 11 | 5.60 | 6.91 | 0 | 0 | fail |
| masked | 42 | 1 | **Top 1**, best reference geometry | 0.2995 | 17 | 1.32 | 2.53 | 1 | 0 | fail |
| masked | 42 | 2 |  | 0.2970 | 12 | 1.31 | 4.85 | 0 | 0 | fail |
| masked | 43 | 1 |  | 0.2335 | 0 | 0.99 | 25.75 | 0 | 0 | fail |
| masked | 43 | 2 | **Top 1**, best reference geometry | 0.2975 | 5 | 1.55 | 18.69 | 0 | 0 | fail |
| masked | 44 | 1 | best reference geometry | 0.2775 | 16 | 1.47 | 2.10 | 0 | 0 | pass |
| masked | 44 | 2 | **Top 1** | 0.2890 | 10 | 6.00 | 7.75 | 0 | 0 | fail |

A WT draw satisfying all four conditions ranked below top 1 for masked seed 44.

Paired top 1 deltas are complete minus masked at the same seed.

| Seed | Contacts | Target RMSD (Å) | Interface RMSD (Å) |
| ---: | ---: | ---: | ---: |
| 42 | -3 | -0.02 | +0.69 |
| 43 | +6 | -0.28 | -11.63 |
| 44 | +1 | -0.41 | -0.83 |

Supplying coordinates for residues 603–612 improved the selected interface fit at seeds 43 and 44, with the larger change at seed 43. At seed 42 the complete template passed all four criteria. The masked pose recovered more contacts and had lower interface RMSD, but its one clash failed the pose criterion. At seed 44 the masked arm sampled a close pose and the confidence ranking favored a lower quality pose. Both masks fell short of repeatable top 1 recovery.

## Other structural controls

The table reports the range of the top 1 geometry across observed seeds for each control and arm. SPRstrong and SPRweak are labels from the frozen panel. These structural scores provide no affinity estimate.

| Template | Control | Seeds | Contacts of 19 | Interface RMSD (Å) |
| --- | --- | ---: | ---: | ---: |
| complete | FACS-low | 3 | 4–10 | 17.82–25.36 |
| complete | SPRstrong | 3 | 5–16 | 1.86–15.42 |
| complete | SPRweak | 3 | 0–10 | 6.67–20.21 |
| masked | FACS-low | 3 | 9–13 | 8.12–17.99 |
| masked | SPRstrong | 3 | 11–12 | 3.52–15.50 |
| masked | SPRweak | 3 | 3–16 | 1.64–21.56 |

Both 1N8Z and 6OGE appear in upstream RFD4 sampling metadata. The results describe known-geometry reproduction and sensitivity to supplying coordinates for positions 603–612. The sampling metadata leave exact checkpoint exposure unresolved. The runtime validates the parsed template coverage and aggregate assignments. Per-query applied masks were unavailable. [Upstream sampling audit](../../results/rfd4-training-metadata/audit.json)

## Target-only MSA control

The audit verified 24/24 poses using the complete 6OGE template and the frozen target alignment. VH and VL used null MSA paths. The threshold and selection rule are unchanged. [MSA protocol](target-msa-protocol.md), [Per sample scores](../../results/template-controls/target-msa/scores.csv), [Full audit](../../results/template-controls/target-msa/audit.json)

The target MSA fell below the gate for the conditional candidate screen with 0/3 passing WT seeds. The complete no-MSA template had 1/3. [Conditional candidate screen](candidate-fallback-protocol.md)

All WT draws use the same 1N8Z fit and footprint as the baseline.

| Seed | Draw | Selection | Mean pair ipTM | Contacts of 19 | Target RMSD (Å) | Interface RMSD (Å) | Clashes | ECD overlay clashes | Gate |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 42 | 1 | best reference geometry | 0.2075 | 6 | 2.31 | 17.93 | 0 | 0 | fail |
| 42 | 2 | **Top 1** | 0.2105 | 3 | 2.02 | 19.61 | 0 | 0 | fail |
| 43 | 1 | best reference geometry | 0.2050 | 3 | 1.58 | 17.53 | 0 | 0 | fail |
| 43 | 2 | **Top 1** | 0.2780 | 4 | 0.93 | 21.37 | 0 | 13 | fail |
| 44 | 1 | **Top 1** | 0.2900 | 7 | 2.21 | 16.54 | 0 | 86 | fail |
| 44 | 2 | best reference geometry | 0.2645 | 8 | 1.87 | 16.41 | 0 | 0 | fail |

Paired WT top 1 deltas are target-only MSA minus complete no-MSA.

| Seed | Contacts | Target RMSD (Å) | Interface RMSD (Å) |
| ---: | ---: | ---: | ---: |
| 42 | -11 | +0.72 | +16.38 |
| 43 | -7 | -0.33 | +14.31 |
| 44 | -4 | -3.39 | +9.63 |

The target MSA reduced WT top 1 footprint recovery and increased interface RMSD at every seed. The best reference geometry among the six MSA WT draws had 16.41 Å interface RMSD and recovered 8/19 contacts. The conditional candidate screen stayed closed. These results describe sensitivity to one frozen target alignment in one known complex.
