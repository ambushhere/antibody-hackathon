# OpenDDE structural replicate results

14 November 2026. All 21 OpenDDE antibody–antigen predictions completed and passed the local artifact audit. WT recovered **2, 0 and 5 of 19** reference footprint positions at seeds 42, 43 and 44. All three fell below the prespecified 13/19 diagnostic. WT had 11.74–14.44 Å target Cα RMSD and 20.48–27.53 Å reference-interface Cα RMSD in the target frame. These control failures rule out using this setup as affirmative structural evidence for the candidate poses. [Protocol](structural-replicate-protocol.md) · [Audit JSON](../../results/structural-replicates/audit.json) · [All metrics CSV](../../results/structural-replicates/scores.csv)

## Controls

- **WT.** The assay manifest records eight WT SPR measurements from 0.02851802 to 0.178 nM, with 0.152 nM first. FACS-low recovered more contacts than WT at every matched seed (3, 6 and 9 of 19). The FACS-low sequence has an unknown SPR outcome.
- **Quantified SPR pair.** Strong (0.02343889 nM) recovered 4/19 at every seed. Weak (13,863.16 nM) recovered 5, 0 and 0 of 19. Contact recovery favored the weak control at seed 42 and the strong control at seeds 43 and 44. Neither control reliably recovered the WT interface in this comparison. The two quantified variants are controls for this panel and cannot calibrate general binding or affinity.
- **Confidence.** WT VH–VL ipTM stayed near 0.955 across seeds while VH–target and VL–target ranged from 0.219 to 0.365. WT mean PAE ranged from 2.41 to 2.63 Å between the antibody chains and from 17.17 to 27.92 Å for antibody–target directions. The fixed 1N8Z interface gave 18.70–22.01 Å with antibody rows and 22.94–24.17 Å with target rows. ipTM describes predicted structural confidence. Affinity, competition and thermal stability require assays.
- **Clashes.** Fourteen of 21 structures had at least one antibody–target heavy-atom pair below 1.5 Å. WT had one at each of seeds 43 and 44. The de novo draw at seed 44 also had one VH–VL heavy-atom pair below 1.5 Å. The OpenDDE result rows omit the provider clash flag.

## Results at each seed

RMSDs are in Å. Target RMSD uses the 109-position target fit. Antibody RMSD uses all mapped VH and VL Cα positions after that same fit. Interface RMSD uses the fixed 19 target and 20 antibody Cα positions. The contact column gives recovered WT positions out of 19, followed by Jaccard over the shared region. The three pair ipTM columns use inferred chain roles.

| Seed | Query | Contacts / Jaccard | Target RMSD | Antibody RMSD | Interface RMSD | Overall ipTM | VH–VL | VH–target | VL–target | Ab–target clashes |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 42 | WT | 2/19 / 0.047 | 12.98 | 60.32 | 27.53 | 0.665 | 0.955 | 0.356 | 0.365 | 0 |
| 42 | SPR strong | 4/19 / 0.093 | 12.77 | 55.41 | 24.97 | 0.630 | 0.952 | 0.312 | 0.305 | 1 |
| 42 | SPR weak | 5/19 / 0.139 | 13.08 | 48.86 | 24.40 | 0.654 | 0.945 | 0.361 | 0.354 | 1 |
| 42 | FACS-low | 3/19 / 0.075 | 14.11 | 62.21 | 29.29 | 0.603 | 0.943 | 0.273 | 0.273 | 0 |
| 42 | `denovo_fullseq_k1_1` | 7/19 / 0.200 | 13.95 | — | — | 0.719 | 0.939 | 0.624 | 0.590 | 0 |
| 42 | `opt_WTAVGFFAMD` | 4/19 / 0.091 | 15.08 | 56.71 | 27.04 | 0.666 | 0.950 | 0.375 | 0.371 | 1 |
| 42 | `opt_YTGDGFPAMN` | 3/19 / 0.086 | 14.19 | 60.47 | 27.66 | 0.637 | 0.949 | 0.314 | 0.323 | 1 |
| 43 | WT | 0/19 / 0.000 | 11.74 | 42.08 | 20.65 | 0.584 | 0.955 | 0.219 | 0.225 | 1 |
| 43 | SPR strong | 4/19 / 0.100 | 16.35 | 58.33 | 29.44 | 0.625 | 0.952 | 0.291 | 0.295 | 2 |
| 43 | SPR weak | 0/19 / 0.000 | 11.88 | 52.67 | 27.78 | 0.560 | 0.943 | 0.182 | 0.192 | 0 |
| 43 | FACS-low | 6/19 / 0.167 | 11.14 | 42.51 | 21.84 | 0.612 | 0.943 | 0.297 | 0.289 | 1 |
| 43 | `denovo_fullseq_k1_1` | 2/19 / 0.054 | 13.35 | — | — | 0.584 | 0.930 | 0.206 | 0.221 | 1 |
| 43 | `opt_WTAVGFFAMD` | 10/19 / 0.294 | 12.50 | 37.70 | 20.63 | 0.619 | 0.951 | 0.279 | 0.290 | 0 |
| 43 | `opt_YTGDGFPAMN` | 4/19 / 0.095 | 14.12 | 59.62 | 28.39 | 0.672 | 0.950 | 0.389 | 0.388 | 1 |
| 44 | WT | 5/19 / 0.143 | 14.44 | 34.32 | 20.48 | 0.630 | 0.956 | 0.297 | 0.302 | 1 |
| 44 | SPR strong | 4/19 / 0.105 | 15.44 | 56.72 | 26.85 | 0.620 | 0.951 | 0.294 | 0.290 | 1 |
| 44 | SPR weak | 0/19 / 0.000 | 14.31 | 54.26 | 29.00 | 0.620 | 0.945 | 0.288 | 0.291 | 0 |
| 44 | FACS-low | 9/19 / 0.243 | 11.82 | 33.87 | 18.59 | 0.595 | 0.942 | 0.264 | 0.259 | 3 |
| 44 | `denovo_fullseq_k1_1` | 8/19 / 0.258 | 14.54 | — | — | 0.614 | 0.935 | 0.254 | 0.268 | 0 |
| 44 | `opt_WTAVGFFAMD` | 5/19 / 0.116 | 14.40 | 55.26 | 24.99 | 0.661 | 0.951 | 0.362 | 0.363 | 1 |
| 44 | `opt_YTGDGFPAMN` | 5/19 / 0.128 | 19.72 | 30.11 | 20.35 | 0.603 | 0.947 | 0.259 | 0.262 | 2 |

## Range across seeds

The JSON records min, max and count for each metric, including every directional PAE value. De novo has no predefined WT antibody residue correspondence. Its target alignment, target-footprint contacts, chain confidence, chain-pair PAE and clash measurements remain defined.

| Query | Contacts out of 19 | Jaccard | Target RMSD | Antibody RMSD | Interface RMSD | VH–VL ipTM | VH–target ipTM | VL–target ipTM | Ab–target clashes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| WT | 0–5 | 0.000–0.143 | 11.74–14.44 | 34.32–60.32 | 20.48–27.53 | 0.955–0.956 | 0.219–0.356 | 0.225–0.365 | 0–1 |
| SPR strong | 4 | 0.093–0.105 | 12.77–16.35 | 55.41–58.33 | 24.97–29.44 | 0.951–0.952 | 0.291–0.312 | 0.290–0.305 | 1–2 |
| SPR weak | 0–5 | 0.000–0.139 | 11.88–14.31 | 48.86–54.26 | 24.40–29.00 | 0.943–0.945 | 0.182–0.361 | 0.192–0.354 | 0–1 |
| FACS-low | 3–9 | 0.075–0.243 | 11.14–14.11 | 33.87–62.21 | 18.59–29.29 | 0.942–0.943 | 0.264–0.297 | 0.259–0.289 | 0–3 |
| `denovo_fullseq_k1_1` | 2–8 | 0.054–0.258 | 13.35–14.54 | — | — | 0.930–0.939 | 0.206–0.624 | 0.221–0.590 | 0–1 |
| `opt_WTAVGFFAMD` | 4–10 | 0.091–0.294 | 12.50–15.08 | 37.70–56.71 | 20.63–27.04 | 0.950–0.951 | 0.279–0.375 | 0.290–0.371 | 0–1 |
| `opt_YTGDGFPAMN` | 3–5 | 0.086–0.128 | 14.12–19.72 | 30.11–60.47 | 20.35–28.39 | 0.947–0.950 | 0.259–0.389 | 0.262–0.388 | 1–2 |

The de novo draw at seed 42 had the largest inferred antibody–target pair ipTM in the panel (0.624 VH, 0.590 VL) and recovered 7/19 contacts. Across seeds its pair scores varied by 0.419 for VH and 0.370 for VL. These three draws and the failed WT controls leave candidate pose support unresolved.

## Provenance and measurement

- **Coverage.** Each finished job matches its submission receipt. Each request contains all seven exact T/H/L FASTAs, one sample per query, the declared seed, no MSA server and CIF plus PAE output. There are no supplied templates or restraints in the request, and the accepted settings disable TFG guidance. All results report one sample and MSA depth 1 with MSA disabled. Seeds are established at job level. The exported rows omit per-prediction seed and checkpoint revision.
- **Integrity.** The audit matched all 66 extracted files byte for byte to their provider archives and checked SHA-256 and sizes for those files and the three ZIPs against the local download receipts. Each provider result receipt matches the archived `results.json` and repeats identical metadata for `control_wt_pae.json`, with one verified file per archive. The three usage receipts report one attempt per target and $0.120819 USD in total historical charges. [Local receipts](../../results/provider-jobs/opendde-replicates/) · [Artifact hashes](../../results/structural-replicates/audit.json)
- **Confidence-axis inference.** OpenDDE exports index-based pTM and pair ipTM with no explicit chain labels or RF3-style error summaries. The PAE sidecar declares tokens in structure file order. The audit verified one Cα token per independently sequenced CIF residue, contiguous residue indices and matching PAE/PDE dimensions. It matched descending index pTM rank to ascending within-chain mean PAE and independently to ascending within-chain mean PDE. Both gave **T/H/L** for all 21 structures. The pair ipTM diagonal matches indexed chain pTM. This rank-based assignment assumes that lower mean error tracks higher pTM on these chains. The exporter provides no numerical reconstruction of pair ipTM from its released arrays. The raw indexed values, chain error means and inferred order remain in the JSON.
- **Geometry.** Contacts use a strict heavy-atom separation below 4.5 Å. Reference and prediction target coordinates use canonical P04626 numbering. The common target positions are 511–602 and 613–629. Jaccard compares contacts with the 19-position crystallographic footprint inside that common set. Matched antibody scaffolds use variable-domain residue indices after verifying equal-length CDRH3 substitutions and an identical VL. The same target-only rigid alignment is used for the antibody and interface RMSDs. The CSV keeps directional PAE values separate by matrix row and column labels.
- **Comparison boundary.** The previous RF3 domain IV WT draw recovered 14/19 with 10.15 Å target and 9.75 Å interface RMSD. OpenDDE WT fell below that contact count and had higher RMSD at all three seeds. The RF3 arm contains one draw, and the earlier Protenix domain IV arm used a target template. This panel leaves checkpoint-wide performance and laboratory outcomes unresolved. [Prior context results](context-results.md)

Raw predictions are preserved in `results/archives/opendde-replicates-20261114.tar.gz`; extract under `results/`. The adjacent `opendde-replicates-manifest.json` records archive and member SHA-256 values.

## Reproduction

```bash
source .venv/bin/activate
PYTHONDONTWRITEBYTECODE=1 python -m scripts.audit_structural_replicates \
  --reference /tmp/her2-ablang/assets/1N8Z.cif
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:src python -m pytest -q -p no:cacheprovider \
  tests/test_structural_replicates.py tests/test_context_panel_audit.py \
  tests/test_her2_structure_audit.py tests/test_her2_context.py tests/test_denovo_refold_audit.py
```
