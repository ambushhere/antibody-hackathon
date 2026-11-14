# HER2 candidate refold screen status

The complete 6OGE template with no MSA reproduced the WT pose at seed 42. Across the three seeds, 1 of 3 WT top poses passed all four criteria. The predeclared launch threshold was 2 of 3. The 57-query candidate run was therefore withheld. [Admission rule](structure-admission-rule.md) · [Baseline control audit](../../results/template-controls/audit.json)

| Seed | Top sample | Mean pair ipTM | WT contacts of 19 | Target fit Cα RMSD (Å) | Reference interface Cα RMSD (Å) | Crop clash pairs below 1.5 Å | Pass |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 42 | 1 | 0.2985 | 14 | 1.302 | 3.221 | 0 | Yes |
| 43 | 1 | 0.3715 | 11 | 1.265 | 7.062 | 0 | No |
| 44 | 2 | 0.3040 | 11 | 5.595 | 6.911 | 0 | No |

The frozen panel retains 43 de novo designs from the pre-existing 64-design VM batch, 10 optimization picks and four controls. The candidate audit has 57 unknown pose rows and no candidate Pareto fronts. Missing candidate measurements are null. Local de novo eligibility still requires official Adaptyv validation. [Frozen inputs](../../configs/her2-candidate-screen/manifest.json) · [Candidate audit](../../results/candidate-screen/audit.json) · [All expected rows](../../results/candidate-screen/scores.csv)
