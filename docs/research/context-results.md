# HER2 context and RF3 structural results

14 November 2026. The prespecified Protenix domain IV rescue failed at seed 42. WT recovered 3 of 19 crystallographic footprint residues in domain IV versus 7 in the matched split arm. All seven ECD predictions missed the WT footprint. RF3 recovered 14 of 19 for WT, but the aligned target remained distorted and measured strong and weak SPR controls reversed on footprint recovery. [Protocol](interaction-protocol.md) · [Context design](her2-target-context.md) · [Audited results](../../results/context-panel/audit.json)

## Controls

Each result cell gives **recovered WT positions out of 19 / Jaccard** within the same 109 resolved target positions.

| Control | Assay evidence | Protenix split | Protenix domain IV | Protenix ECD | RF3 domain IV |
| --- | --- | ---: | ---: | ---: | ---: |
| WT `WGGDGFYAMD` | 0.152 nM SPR, first of 8 WT reads spanning 0.02851802–0.178 nM | 7 / 0.226 | 3 / 0.111 | 0 / 0 | 14 / 0.452 |
| SPR strong `WGSQGFYAYD` | 0.02343889 nM | 2 / 0.062 | 4 / 0.143 | 0 / 0 | 10 / 0.333 |
| SPR weak `YLLNGVFVFD` | 13,863.16 nM, quantified binding | 0 / 0 | 3 / 0.103 | 0 / 0 | 12 / 0.375 |
| FACS-low `WSMFGHYMFA` | FACS class only, SPR status unknown | 3 / 0.103 | 9 / 0.300 | 0 / 0 | 13 / 0.371 |

- **Protenix split to domain IV.** WT had no new heavy-atom clashes below 1.5 Å and the strong SPR control recovered one more contact than the weak control in domain IV. WT lost four contacts, fell ten short of the 13/19 criterion and recovered six fewer than FACS-low. Its target Cα RMSD improved from 3.89 to 1.28 Å, but its antibody Cα RMSD in the target frame stayed above 50 Å (53.23 to 50.94 Å). WT contacted 13 positions outside the shared region in domain IV, including 12 in the extended tail. The better fit to the target did not recover the antibody pose.
- **Protenix ECD.** All seven poses made contacts outside the shared domain IV region and recovered 0/19 WT positions. The WT target fit was 0.63 Å over the same 109 Cα atoms. Its antibody centroid was 83.88 Å from the reference and its reference-interface Cα RMSD was 49.57 Å.
- **RF3 domain IV.** WT recovered the most reference contacts among its seven draws but had 10.15 Å target Cα RMSD, 15.52 Å antibody Cα RMSD and 9.75 Å reference-interface Cα RMSD after fitting the shared target. The strong control was closer than the weak control in the target frame by antibody Cα RMSD (23.97 vs 37.54 Å), while both footprint recall and target pair ipTM ranked the weak control above the strong control. FACS-low recovered 13/19 and had three sub-1.5 Å atom pairs.

### Confidence and clashes

The WT overall / VH–VL / maximum antibody–target pair ipTM was **0.673 / 0.947 / 0.467** in Protenix split, **0.583 / 0.936 / 0.246** in Protenix domain IV, **0.275 / 0.906 / 0.164** in Protenix ECD and **0.543 / 0.838 / 0.259** in RF3. For RF3 WT, published pair-average PAE was 5.98 Å for VH–VL and 21.20/21.42 Å for VH–target/VL–target. In split, the largest WT and FACS-low antibody–target scores were both VH–D (0.467 and 0.557). Split has four antibody–target pairs and the continuous arms have two. The maximum is a within-arm diagnostic. The CSV contains each exposed pair and the two directional PAE means for the four controls on fixed crystallographic interface positions. [Confidence and geometry CSV](../../results/context-panel/scores.csv)

All 28 provider clash flags were false. The separate strict heavy-atom check counted 1 below 1.5 Å for Protenix domain IV SPR weak, 4 for ECD WT and 1 for ECD SPR strong. RF3 counted 1 each for SPR strong and weak, 3 for FACS-low and 3 for `opt_YTGDGFPAMN`. The other 21 draws had zero at this cutoff. The closest distances and counts are in the CSV.

## Exploratory designs

| Query | Protenix split | Protenix domain IV | Protenix ECD | RF3 domain IV |
| --- | ---: | ---: | ---: | ---: |
| `denovo_fullseq_k1_1` | 13 / 0.351 | 13 / 0.520 | 0 / 0 | 11 / 0.314 |
| `opt_WTAVGFFAMD` | 4 / 0.095 | 4 / 0.138 | 0 / 0 | 11 / 0.289 |
| `opt_YTGDGFPAMN` | 4 / 0.093 | 16 / 0.696 | 0 / 0 | 12 / 0.273 |

The de novo design recovered the same number of footprint positions in both short Protenix contexts, with fewer off-footprint positions in the shared region in domain IV. Its different VH and VL sequences lack a predefined WT antibody residue correspondence, so the audit provides target fit and footprint without antibody or interface RMSD. `opt_YTGDGFPAMN` achieved 2.31 Å reference-interface Cα RMSD in the Protenix domain IV target frame (5.68 Å across the antibody Cα atoms). Its other three draws gave 19.03, 43.39 and 19.09 Å interface RMSD. ECD assigned maximum target pair ipTM of 0.689 and 0.712 to the de novo and `opt_YTGDGFPAMN` poses despite 0/19 footprint recall and VH–VL ipTM of 0.056 and 0.072. The prespecified controls failed to support candidate selection from these draws.

## Provenance and measurement

- **Coverage and hashes.** All 21 Protenix and seven RF3 predictions completed. The audit verifies the published archive digest, 106 in-scope raw member hashes and sizes, staged input hashes, complete query IDs and exact output chain sequences and numbering. Protenix result records identify seed 42 and sample 1 per prediction. The RF3 job and results receipts record seed 42 and one diffusion sample per target. The RF3 receipt lists `control_wt_pae.json` twice with identical metadata. One archive sidecar exists and its hash matches. RF3 output sequences match the exact domain IV queries. The provider receipts lack the submitted sequence payload and model weights revision. [Archive](../../results/archives/README.md) · [RF3 receipts](../../results/provider-jobs/)
- **Protenix template controls.** Apheris Protenix 1.3.0, `protenix-v2.pt`, image ID pinned in the audit, target-only templates, null MSA paths and no contact restraints. Logs show 14 assignments for split and 7 each for domain IV and ECD, with 109, 109 and 581 resolved template positions. The logs report aggregate counts and do not export each applied template mask. Canonical sequence and coordinate indices in the supplied template files match the manifest, including the gaps. All 836 split atoms occur at identical XYZ in domain IV and ECD. ECD adds 3,654 original 1N8Z coordinates upstream. Thus split versus domain IV changes topology and modeled sequence without changing their resolved source coordinates. ECD changes both sequence context and the extent of the template. The historical split template differs by 0.991 Å RMSD over the same 836 atoms after a rigid fit and yielded 11/19 WT contacts in that earlier run. The new split run is the primary comparator. [Input manifest](../../configs/her2-context/manifest.json)
- **Chain order and confidence.** Protenix query, CIF and pair score roles are C/D/H/L in split and T/H/L in domain IV and ECD. Detailed Protenix atom tokens and lDDT were checked against CIF residue and atom order. RF3 was run without a supplied template or MSA under the recorded protocol. Exact sequence matching identifies T/H/L in the RF3 CIF and token sidecar. Each RF3 summary has a unique H/L/T order when the exposed directional minima and bidirectional pair means for PAE and PDE are reproduced from the file-order arrays. The receipt exposes indexed pair ipTM and no index labels. The audit uses the inferred summary order for those values. RF3 and Protenix differ in model and template policy on the same domain IV sequence.
- **Geometry.** Contacts use protein heavy atoms at a strict distance below 4.5 Å. Target IDs are mapped to P04626 by the manifest, preserving the split gap at 603–612. Recall divides the recovered positions by the 19 deposited Fab footprint positions. Jaccard uses the intersection and union within 511–602 and 613–629. The audit fits the same 109 target Cα atoms to the original 1N8Z chain C with a proper rigid transform. Antibody Cα metrics use 1N8Z B for VH and A for VL, with the WT variable residues matched by index after exact sequence verification. For the matched variants, substitution-only changes of the same length in CDRH3 preserve this index correspondence. The fixed reference interface uses 19 target and 20 antibody Cα positions defined by crystallographic contacts. Its RMSD and the antibody-only RMSDs use the target transform throughout. There is no second fit on the antibody. Reference-interface PAE means report both matrix directions over these same 19 × 20 positions for the controls. [Audit code](../../scripts/audit_context_panel.py) · [1N8Z deposit](https://www.rcsb.org/structure/1N8Z)

Each arm contains one draw per query at one seed. Seed 42 does not estimate pose reproducibility. The runs use paired VH and VL and omit glycosylation, assay construct effects and expression conditions. The outputs provide structural hypotheses. Binding, competition and thermal stability still require measurements.

## Reproduction

Use the extracted raw files described by the archive README and a local copy of the original [1N8Z mmCIF](https://files.rcsb.org/download/1N8Z.cif). The script requires the deposit SHA-256 from the manifest and uses no network or accelerator.

```bash
source .venv/bin/activate
PYTHONDONTWRITEBYTECODE=1 python -m scripts.audit_context_panel --reference /path/to/1N8Z.cif
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:src python -m pytest -q -p no:cacheprovider \
  tests/test_context_panel_audit.py tests/test_her2_structure_audit.py \
  tests/test_her2_context.py tests/test_denovo_refold_audit.py
```
