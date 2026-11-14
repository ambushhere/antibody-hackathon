# Stability combination hypothesis

14 November 2026. This is an experimental hypothesis, not a submission or an admissibility determination.

## Evidence and choice

Okazaki et al. measured trastuzumab scFv VL R66G at Tm 66.7 ± 0.2 °C versus WT 61.8 ± 0.2 °C. Their Table 1 gives HER2 SPR KD 4.8 ± 0.7 nM versus WT 2.0 ± 0.4 nM. Although the abstract describes retained affinity, the numerical point estimates show a 2.4-fold weaker interaction. R66G is a framework substitution outside the CDRs. The work expressed scFv in E. coli and measured dye-based DSF. We use it as a measured thermal/binding tradeoff control, not a free stability gain. [Primary article and Table 1](https://pubs.acs.org/doi/10.1021/acsomega.3c01948).

Morrison et al. report a +4 °C thermal shift for VL H91Y. H91 lies in the light-chain binding loop. Their selection and affinity measurements involve the HER2-mimetic H98 peptide, so these binding results do not establish retention of binding to full HER2. H91Y remains a thermal control rather than our first combination choice. [Primary article](https://www.nature.com/articles/s41467-021-26279-8).

The full task text permits point mutations to an existing antibody (page 1) and applies identity and distance conditions to classify and admit optimisation designs (pages 2–3). We found no explicit sentence freezing every residue outside CDRH3 in that text. The slides call the track CDRH3 optimisation. This supports an exploratory combination experiment, but does not establish official eligibility. The stricter H3 novelty gate remains necessary in the combination. No known mutation alone is a qualifying design.

## Factorial experiment

For one chosen new H3, compare four constructs in the same exact scFv orientation, linker, tags and assay batch:

| Construct | H3 | VL66 | Purpose |
| --- | --- | --- | --- |
| WT | WT | R | Matched reference |
| Thermal control | WT | G | Test transfer of the published R66G effect |
| H3 parent | New | R | Binding/novelty proposal |
| Combined | New | G | Test the interaction between the two changes |

Measure SPR, competition and antigen-free thermal transitions for all four. The combination effect is not the sum of the two measured parent effects. A useful interaction diagnostic is `(Tm combined − Tm H3 parent) − (Tm R66G control − Tm WT)`, with replicate uncertainty and the same thermal transition. The hypothesis fails if the apparent thermal benefit disappears, binding deteriorates beyond the prespecified comparison margin or competition is lost. No such measurements exist for our candidates.

The sequence-only FACS and SPR models were trained on WT VL. Their H3-parent scores must not be copied onto the R66G constructs as validated predictions. Record them, if needed, explicitly as parent evidence. Cofolding can test for gross geometry failures but cannot verify ΔTm or the absence of an affinity penalty.

Our ThermoMPNN-D pilot predicts R66G in the wrong thermal direction on both evaluated backbones. We therefore prioritize this hypothesis from measured literature, not from that model. [Thermal pilot](../../results/thermal-pilot/README.md).
