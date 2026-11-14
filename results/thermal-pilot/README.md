## Paired Fv thermal proxy pilot

14 November 2026. **Decision: retain `thermal evidence missing` for campaign designs.** The CPU run demonstrates that the pinned model produces reproducible folding scores for paired VH/VL. R66G was directionally discordant on both structures, H91Y changed sign and the CDRH3 panel ranked as poorly as or worse than mutation count. The scores have an exploratory role as fixed-backbone folding hypotheses for follow-up. This pilot gives no basis for ranking, filtering or promoting designs by thermal stability.

### scFv controls

Positive experimental ΔTm means a higher measured apparent melting temperature. Negative model folding ΔΔG favours folding under the upstream convention. The columns use their original units. [ThermoMPNN-D paper](https://doi.org/10.1002/pro.70003), [pinned inference documentation](https://github.com/Kuhlman-Lab/ThermoMPNN-D/tree/df9a75aaddb674a7c4c193005031fc0536d325fb#inference).

| VL control | WT → mutant Tm (°C) | Measured ΔTm (°C) | 1N8Z folding ΔΔG (kcal/mol) | 4HKZ folding ΔΔG (kcal/mol) |
| --- | ---: | ---: | ---: | ---: |
| H91Y | 68.5 → 72.5 | +4.0 | −0.927 | +0.016 |
| A34D + Y49S, additive | 68.5 → 62.5 | −6.0 | +3.686 | +3.936 |
| A34D + Y49S, joint | 68.5 → 62.5 | −6.0 | +2.004 | +2.643 |
| R66G | 61.8 → 66.7 | +4.9 | +1.463 | +1.357 |

- All three contrasts scored on each template. The additive comparison agreed in direction for **2 of 3** on 1N8Z and **1 of 3** on 4HKZ. These are counts on three shared controls. On 4HKZ, H91Y is effectively near zero, with its raw positive sign reported as discordant. The joint double preserves the additive direction on both templates.
- A34D and Y49S contributed +1.904 and +1.783 on 1N8Z, and +1.342 and +2.594 on 4HKZ, respectively. Their Cα separation was 4.344 Å and 4.208 Å, within the 5 Å cutoff. The joint models used only this pair. The joint outputs also matched the native Siamese forward average to within 10⁻⁶ kcal/mol on both templates.
- Morrison variants 1.1 and 3.2 were trastuzumab scFvs measured with a fluorescent thermal shift dye, using two protein preparations and four technical replicates each. Table 1 implies −6.0 °C for 3.2, while the narrative gives −5 °C. Okazaki used a separate trastuzumab scFv and SYPRO Orange DSF protocol. Table 1 gives 61.8 ± 0.2 and 66.7 ± 0.2 °C for WT and R66G. Comparisons are paired within each study. [Morrison Table 1 and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8514426/#Tab1), [Okazaki Table 1 and Figure 3](https://doi.org/10.1021/acsomega.3c01948).

### HER2 CDRH3 stress test

The two pinned FLAb tables joined by exact VH/VL sequence pairs and matched the authors' sequence IDs and both inflections. All five equal-length variants among the 11 entries in Supplementary Table 8 are represented. They share the WT VL and VH outside the 13-residue study loop. Six entries changed length. FLAb also has ID 12, which is outside Supplementary Table 8. The study loop `SRWGGDGFYAMDY` spans source B97–109 and model B204–216. It includes two N-terminal flanking positions and one C-terminal position beyond the challenge's B99–108 loop. Mutation counts are positionwise comparisons to that 13-residue WT segment.

The assay measured full-length, CHO-produced IgG1 by nanoDSF on a Prometheus Panta in PBS at 0.5 mg/mL and 0.5 °C/min. The primary comparison uses the first inflection. IgG constant domains and multiple thermal transitions limit comparison with isolated Fv folding and fragment thermostability. The later inflection is recorded alongside the first without assigning either transition to a domain. The WT first and second inflection entries are both 79.52 °C. Candidate 4 is indeterminate for the direction check at +0.72 °C given the lack of repeat-level uncertainty. Its Table 8 value is above WT despite the narrative summary. [Absci methods and Supplementary Table 8](https://absci.gcs-web.com/static-files/a53c2b92-0e77-4899-a52c-4485be691cf6).

| Published Absci ID | Changes | First inflection (°C) | Second inflection (°C) | Δ first inflection (°C) | 1N8Z additive ΔΔG (kcal/mol) | 4HKZ additive ΔΔG (kcal/mol) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 9 | 70.79 | 79.42 | −8.73 | +3.941 | +1.399 |
| 4 | 12 | 80.24 | 80.24 | +0.72 | +6.338 | +8.729 |
| 5 | 7 | 71.86 | 71.86 | −7.66 | +2.967 | +2.744 |
| 7 | 6 | 70.66 | 78.36 | −8.86 | +2.582 | +2.164 |
| 11 | 8 | 75.66 | 75.66 | −3.86 | +3.442 | +2.643 |

Both templates gave positive additive ΔΔG for all five variants. Both agreed with all four clear negative shifts, as did the always-destabilising rule. Both ranked ID 4 as most destabilising. ID 1 moved from second most destabilising on 1N8Z to least on 4HKZ. The two templates share the same WT sequence but come from different crystal complexes.

| Compared set | Measured order, strongest Tm loss first | 1N8Z order | 4HKZ order | Mutation-count order | Spearman ρ vs measured order, 1N8Z / 4HKZ / count |
| --- | --- | --- | --- | --- | --- |
| All five | 7, 1, 5, 11, 4 | 4, 1, 11, 5, 7 | 4, 5, 11, 7, 1 | 4, 1, 11, 5, 7 | −0.70 / −0.80 / −0.70 |
| Four clear shifts | 7, 1, 5, 11 | 1, 11, 5, 7 | 5, 11, 7, 1 | 1, 11, 5, 7 | −0.40 / −0.60 / −0.40 |

The rank statistics describe this selected panel and use negative measured ΔTm as destabilisation severity. The 6–12 substitutions per loop are summed on a fixed WT backbone. Mutant relaxation and interactions among substitutions are outside that calculation. The scFv controls change VL, and the IgG panel has a single scaffold and was selected for binding and other properties.

### Coordinates and execution

- 1N8Z is HER2-bound trastuzumab Fab. 4HKZ is trastuzumab Fab bound to protein A and protein L. Both inputs retained only the original variable-domain N, Cα, C and O coordinates, without partners, constant domains or a scFv linker. The resulting antigen-free inputs still inherit their partner-bound conformations. [1N8Z deposition](https://www.rcsb.org/structure/1N8Z), [4HKZ deposition](https://www.rcsb.org/structure/4HKZ).
- Both ATOM-derived sequences matched the SEQRES prefixes and the pinned FLAb WT pair. Source VL A1–107 maps to model A1–107. Source VH B1–120 maps to model B108–227. Both source domains had contiguous residue numbers, blank insertion codes and 227/227 complete backbones per template. The chain boundary was retained with separate chain IDs and TER records. Within-chain C–N distances across the two templates were 1.306–1.350 Å. At the 11 alternate Cα sites in 4HKZ, the highest occupancy was selected, then alphabetical alternate ID for ties. The mapping records those choices. 1N8Z had zero alternate backbone sites.
- The model saw both chains in mapping order, and its returned WT letters matched the mapping at all 227 positions. Each scan yielded all 4,540 canonical scores, including 227 WT identity entries, at the upstream threshold of 100. Downstream tables retain the four prespecified control substitutions and the published panel's substitutions. No campaign proposal was scored.
- The single head used `ThermoMPNN-ens1.ckpt` and the targeted joint head used `ThermoMPNN-D-ens1.ckpt`, each with `v_48_020.pt` from the pinned tree. The scratch patch only changes execution devices and checkpoint map location. The wrapper sets the upstream checkpoint directory, uses evaluation and inference mode, seeds torch at 0 and caps thread counts at two intra-op and one inter-op.
- The single-only run took 3.57 s wall clock, and the single plus targeted joint run took 3.89 s, both including interpreter and model startup. A repeat of the joint run produced byte-identical CSV and JSON outputs. Four focused tests passed. Upstream emitted its Biopython `pairwise2` deprecation and TorchMetrics buffer warnings during import.

### Reproduction

From the repository root, run the block in Bash with failure propagation and a fresh private scratch directory. Inference verifies the pinned input and patched entry-point hashes before importing upstream code and refuses to overwrite existing results. The commands below assume the existing root `.venv` provides Python 3.13.18 on arm64, torch 2.15.0, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.1, Biopython 1.88 and pytest 9.2.0. Preflight confirmed these were available and that the five packages installed below were absent. The scratch venv exposes the parent site packages through a local `.pth` because `--system-site-packages` alone only exposes base interpreter packages.

```bash
set -euo pipefail
repo=$PWD
scratch=$(mktemp -d "${TMPDIR:-/tmp}/her2-thermal-replay.XXXXXX")
commit=df9a75aaddb674a7c4c193005031fc0536d325fb
flab=e05e3f0135783ffb7de412df6a1467e7a0031911
absci=cd2bd4b6932b5301d435555dc22084d407e4b92b
mkdir -p "$scratch"

curl -fsSL --max-time 60 "https://codeload.github.com/Kuhlman-Lab/ThermoMPNN-D/tar.gz/$commit" -o "$scratch/thermompnn-d.tar.gz"
(cd "$scratch" && head -n 1 "$repo/results/thermal-pilot/inputs.sha256" | shasum -a 256 -c -)
tar -xzf "$scratch/thermompnn-d.tar.gz" -C "$scratch"
for id in 1N8Z 4HKZ; do
    curl -fsSL --max-time 30 "https://files.rcsb.org/download/$id.pdb" -o "$scratch/$id.pdb"
done
for i in 1 2; do
    curl -fsSL --max-time 30 "https://raw.githubusercontent.com/Graylab/FLAb/$flab/data/thermostability/shanehsazzadeh2023unlocking_DSF$i.csv" -o "$scratch/dsf$i.csv"
done
curl -fsSL --max-time 30 "https://raw.githubusercontent.com/AbSciBio/unlocking-de-novo-antibody-design/$absci/functionality-developability-cross_reactivity-data.csv" -o "$scratch/absci.csv"
(cd "$scratch" && shasum -a 256 -c "$repo/results/thermal-pilot/inputs.sha256")
(cd "$scratch/ThermoMPNN-D-$commit" && patch -p1 < "$repo/results/thermal-pilot/upstream-cpu.patch")

source .venv/bin/activate
export PYTHONDONTWRITEBYTECODE=1
parent_site=$(python -c 'import sysconfig; print(sysconfig.get_path("purelib"))')
python -m venv --system-site-packages --without-pip "$scratch/venv"
source "$scratch/venv/bin/activate"
python -c 'import sysconfig,pathlib,sys; pathlib.Path(sysconfig.get_path("purelib"), "research-parent.pth").write_text(sys.argv[1]+"\n")' "$parent_site"
python -m pip --isolated install --index-url https://pypi.org/simple --no-cache-dir --disable-pip-version-check \
    omegaconf==2.3.0 pytorch-lightning==2.5.5 torchmetrics==1.8.2 tqdm==4.67.1 wandb==0.22.3

export WANDB_MODE=disabled
/usr/bin/time -p python scripts/thermal_proxy_pilot.py \
    --scratch "$scratch" --output "$scratch/reproduced" --joint
python -m pytest -q -p no:cacheprovider tests/test_thermal_proxy_pilot.py
(cd "$scratch" && shasum -a 256 -c "$repo/results/thermal-pilot/scratch-derived.sha256")
for file in controls.csv absci-stress.csv substitution-contributions.csv rank-stats.json 1N8Z-mapping.csv 4HKZ-mapping.csv; do
    cmp "$scratch/reproduced/$file" "$repo/results/thermal-pilot/$file"
done
shasum -a 256 -c results/thermal-pilot/results.sha256
```

The manifests contain SHA-256 for the source archive, checkpoints, downloaded structures and dataset inputs, patched upstream files, derived backbones, scripts and score tables. `substitution-contributions.csv` gives the source residue ID, model index and raw single score for every summed substitution. The downloaded archive, checkpoints, raw datasets and coordinate files are in scratch.
