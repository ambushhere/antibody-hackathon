# HER2 target context experiment

Test whether supplying HER2 as one continuous target chain improves recovery of the known trastuzumab footprint and separates the controls in Protenix v2. The historical target-only refold used seed 42 and one diffusion sample. WT had a maximum antibody–target pair ipTM of 0.382 and recovered 11/19 reference contacts. The prespecified FACS-low control scored 0.556 and recovered 4/19. Its largest pair score was H–D. The two statistics already disagree on which pose resembles the crystallographic interface. [Historical scores](../../results/denovo-refold/scores.csv)

## Numbering and continuity

RCSB and PDBe SIFTS independently map 1N8Z entity 3, author and label chain C, positions 1–607 to human HER2 UniProt P04626 positions 23–629. Both coordinate residue numbers and deposited polymer indices follow this offset. **P04626 position = 1N8Z author residue + 22.** The local 1N8Z PDB matches all 4,491 deposited HER2 protein atom coordinates in the RCSB mmCIF. UniProt annotates signal peptide 1–22, extracellular region 23–652 and transmembrane helix 653–675. This study uses the domain IV construct boundary 511–652, equivalent to mature 489–630. [RCSB chain and SIFTS alignment](https://data.rcsb.org/rest/v1/core/polymer_entity/1N8Z/3), [PDBe SIFTS mapping](https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/1n8z), [P04626 annotation](https://rest.uniprot.org/uniprotkb/P04626.json), [Domain IV construct definition](https://patents.google.com/patent/US20060275305A1/en)

| Arm and query order | P04626 sequence | Deposited coordinate coverage |
| --- | --- | --- |
| Split, C/D/H/L | C 511–602 (92 aa), D 613–629 (17 aa) | All 109 queried residues have coordinates |
| Domain IV, T/H/L | T 511–652 (142 aa) | 109 residues. P04626 603–612 are unobserved and 630–652 are beyond the deposited polymer |
| ECD, T/H/L | T 23–652 (630 aa) | 581 residues. P04626 124–132, 325–327, 383–386 and 603–612 are unobserved. P04626 630–652 are beyond the deposited polymer |

Positions 603–612 restore the sequence `GVKPDLSYMP` between `…CPS` and `IWK…`. The sequence past the deposited polymer is `CTHSCVDLDDKGCPAEQRASPLT`. The two continuous queries retain all P04626 cysteines. The deposited Cys600–Cys623 disulfide spans the artificial split between C and D. The continuous templates retain this deposited bond annotation. UniProt also annotates Cys626–Cys634 and Cys630–Cys642. The restored sequence contains those residues, with no template atoms for the three cysteines beyond P04626 629. [Deposited structure](https://files.rcsb.org/download/1N8Z.cif), [UniProt disulfide features](https://rest.uniprot.org/uniprotkb/P04626.json)

At a strict 4.5 Å heavy atom cutoff, the 1N8Z Fab and organizer Fv produce the same 19 target contact positions. The organizer C and D retain mature numbering even though D derives from deposited C. The 1N8Z footprint in P04626 coordinates is **579, 580, 582, 583, 591–595, 601, 613, 615, 618, 620, 622, 624–627**. `src/epitope_map.py` uses a structure's residue IDs directly in `main()` and sends them to precursor-numbered domain bins. Its 1N8Z values, and the matching AGENTS.md values 557–605, are therefore in mature numbering. The mapping is structure specific. For example 6OGE author chain A starts at 23 and already uses precursor numbering. [Epitope map source](../../src/epitope_map.py), [PDBe 6OGE mapping](https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/6oge)

The four RFD4 conditioning atoms are below. These describe the previous design input. The refold queries use ordinary sequence and supply no contact restraints. [Design specification](../../results/experiment-20261114/input-specs/denovo-full-sequence.json)

| Organizer site | 1N8Z chain C author position | P04626 position | Domain IV T position | ECD T position |
| --- | ---: | ---: | ---: | ---: |
| D Lys593 NZ | 593 | Lys615 | 105 | 593 |
| C Phe573 CD1 | 573 | Phe595 | 85 | 573 |
| D Ile591 CG2 | 591 | Ile613 | 103 | 591 |
| C Pro572 CD | 572 | Pro594 | 84 | 572 |

The target-only mmCIF files have query-length polymer sequences in `_entity_poly` and `_entity_poly_seq`, linked through `_struct_asym`. Author IDs preserve mature numbering at resolved residues. Label IDs index the complete polymer, leaving gaps in the atom table at unobserved positions. The sequence scheme lists all query positions. Deposited gaps are listed in `_pdbx_unobs_or_zero_occ_residues`. The extra tail appears in the polymer sequence metadata. All atom rows come from 1N8Z chain C. The terminal OXT at mature 607 is omitted uniformly because it is internal to the continuous arms. The split and domain IV templates contain exactly the same 836 atom keys and XYZ coordinates, and the ECD template adds 3,654 deposited atoms upstream. All three templates supply target chains only and encode the same polymer sequences as their respective query target chains. The open Protenix parser indexes atoms by label sequence ID and represents gaps as missing. [Protenix parser](https://github.com/bytedance/Protenix/blob/main/protenix/data/template/template_parser.py)

The organizer design crop's 837 atoms match the source deposit exactly, including OXT. The template actually passed to the historical refold was `her2_target_refold.cif` (SHA-256 `ddfebcd1406208fb58eccd4f7039c449259e092844c561f702be340ba7334f3d`). It shares 836 atom keys with the organizer crop but differs by 0.99 Å all-atom RMSD after optimal rigid superposition. Use the newly prepared `split-target.cif` for the primary matched comparison. Treat a repeat using the old split template as a historical sensitivity check. [Historical wrapper command](../../scripts/run_denovo_refold.sh), [Historical split template](../../results/experiment-20261114/input-specs/her2_target_refold.cif)

## Panel

Every context uses the same paired variable heavy and light sequences. `manifest.json` records sources, the full index map, canonical sequence SHA-256, file SHA-256 and per-chain hashes. `targets.fasta` contains all query target sequences. [Inputs and provenance](../../configs/her2-context/manifest.json), [Target sequences](../../configs/her2-context/targets.fasta)

| Query | Selection |
| --- | --- |
| `control_wt` | Trastuzumab `WGGDGFYAMD`, paired sequences from the historical control. The first reported SPR KD is 0.152 nM. The source contains eight WT entries, 0.02851802–0.178 nM |
| `control_facs_low` | `WSMFGHYMFA`, the first FACS-low CDRH3 without cysteine in source order. Exact historical pair. FACS label only |
| `control_spr_strong` | `WGSQGFYAYD`, 0.02343889 nM |
| `control_spr_weak` | `YLLNGVFVFD`, 13,863.16 nM |
| `denovo_fullseq_k1_1` | `denovo-fullseq-k1/her2_denovo_full_sequence_1` in the identity CSV. Historical footprint recall 13/19. 239 aa across VH/VL, canonical cysteines and 64.4%/66.1% paired IMGT identity to the two local references. Official novelty unverified |
| `opt_WTAVGFFAMD` | First pick from the conservative optimization shortlist, official validation pending |
| `opt_YTGDGFPAMN` | Second pick from the conservative optimization shortlist, official validation pending |

The SPR extremes are drawn from 182 quantified non-WT rows with the exact WT VL and exact WT VH outside the designated CDRH3 window. The weakest still has reportable binding. FACS-low and the measured weak control represent different evidence. [Identity audit](../../results/experiment-20261114/denovo_identity.csv), [Optimization shortlist](../../results/spr-candidate-ranking/shortlist.csv)

## Launch and readout

Start with the seven pairs in split and domain IV at seed 42, one draw each. The first four queries in each JSON form the prespecified control panel if compute is limited. Use ECD as a second target context test, starting with the four controls if needed. Preserve query order, model, template policy, seed, sample count and null MSA paths across compared arms. Before using scores, check logs for two target assignments per split query and one per continuous query, with none for H or L. The seven-pair runs should report 14 assignments in split and 7 in each continuous arm. Confirm T query length and coverage, particularly that the loop and extended tail have no template coordinates. Inspect score matrix and output chains for the intended query order. If the image fails to assign a continuous template or shifts mapping at the gap, the templated comparison is invalid.

Commands for the root operator after staging `configs/her2-context/` under the existing VM campaign directory. The image, seed, sample and detailed score flags match the historical runner. Record the local image ID and compare it to the historical `sha256:d0d1d4a8a911b8001a0868aa30b9f88f598ce12236cdd0b75a4b3465e7a3da58`. The commands were not executed for this preparation.

```bash
set -euo pipefail
campaign=/home/hackathon/almaz/out/campaign-20261114
inputs="$campaign/her2-context"
image=quay.io/apheris/oslo-hackathon:apheris-protenix-v1.3.0-v2
docker image inspect "$image" --format '{{.Id}}'
for arm in split domain-iv; do
  out="$campaign/target-context-${arm}-seed42"
  test ! -e "$out"
  docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v /home/hackathon/almaz:/home/hackathon/almaz \
    "$image" module-run predict --input "$inputs/$arm-queries.json" --output "$out" \
    --seeds 42 --num-diffusion-samples 1 --write-detailed-scores true \
    --template-paths "$inputs/$arm-target.cif" > "$out.log" 2>&1
done
```

Set the loop to `ecd` for the ECD arm. For a control-only run, subset each query JSON with `jq '{queries: .queries[:4]}'`, preserve the arm's template and use a distinct output path. Check the run outputs against the input hashes in the manifest.

- Convert target positions to P04626 using the manifest before calculating geometry. Use predicted atoms at the expected chain and sequence indices and validate the numbering against the predicted mmCIF. Score footprint recall over the same 19 positions and Jaccard within the shared target panel, 511–602 plus 613–629. Separately report all contacts outside that shared mask, minimum atom distance and H/L-to-target heavy atom clashes below 1.5 Å. Use the same strict 4.5 Å cutoff. Compare WT pose after superposing on the same resolved target Cα atoms in the shared mask.
- Report H–C, L–C, H–D and L–D for split, and H–T and L–T for continuous arms. Keep overall and VH–VL ipTM separate. The maximum over four pairs and maximum over two pairs have different opportunities for large values. Interpret pair ipTM as a within-arm confidence diagnostic. When detailed PAE is available, report both PAE directions over the same reference interface positions for the four WT-background controls.
- Prespecified context signal. Domain IV should increase WT footprint recall relative to the new split baseline, reach at least 13/19, avoid new sub-1.5 Å clashes and exceed FACS-low recall. The strongest measured SPR control should recover more of the footprint than the weakest. A flat or reversed WT result, controls moving together, or a strong/weak reversal falsifies this operational rescue at seed 42. Repeat paired arms with seeds 43 and 44 before treating a seed-42 rescue as repeatable. Record all draws and failed queries. De novo and optimization poses are exploratory and do not adjudicate the control comparison.

The experiment changes chain topology, restores missing sequence, removes the two artificial ends at the break and restores the disulfide annotation across it. Domain IV versus ECD further changes the templated target context and the number of residues in confidence summaries. A positive result cannot attribute the effect to any one of these changes. All arms represent one unglycosylated HER2 molecule using coordinates from a Fab-bound crystal. UniProt marks domain IV glycosylation at P04626 N530, N571 and N629. Glycan occupancy, HER2 homo- or hetero-oligomerization, membrane orientation and the assay antigen construct can change accessibility. Pairing VH and VL in silico omits scFv linker or Fab constant domains and the expression, immobilization and competitive binding conditions. Footprint recovery and confidence measure a predicted pose. Binding, competition and thermostability require their own measurements. [UniProt features](https://rest.uniprot.org/uniprotkb/P04626.json), [1N8Z structure and publication](https://www.rcsb.org/structure/1N8Z)

## Reproducing the inputs

From the repository root with the current environment. The script reads public 1N8Z and SIFTS APIs by default and accepts local copies of those inputs. Supply the original assay tables for regeneration. The candidate defaults are shortlist picks 1 and 2, with `--optimizer-cdrh3 FIRST SECOND` available for an explicit replacement from that shortlist.

```bash
source .venv/bin/activate
PYTHONDONTWRITEBYTECODE=1 python scripts/prepare_her2_context.py \
  --facs /tmp/her2-benchmark.NjMJ71m1/facs.csv \
  --spr /tmp/her2-benchmark.NjMJ71m1/spr.csv \
  --local-pdb /tmp/her2-thermal-pilot/1N8Z.pdb
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider tests/test_her2_context.py
```

Local tests parse the generated templates by author and label indices and check paired sequences, file hashes, atom identity, numbering gaps, cysteines and the restored cross-gap disulfide annotation. Regeneration reproduced the manifest and all artifact hashes. The seven new tests and five adjacent structure and de novo tests passed with `PYTHONPATH=src`. Apheris v1.3.0 template assignment and inference are pending.
