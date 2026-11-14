# HER2 target template coordinate ablation

The experiment tests whether supplying coordinates for HER2 P04626 residues 603–612 helps Protenix v2 reproduce the trastuzumab pose. It uses one continuous target chain, 511–629 (119 residues), and the four frozen Fv controls WT, FACS-low, SPRstrong and SPRweak. Both arms use the same query file, null MSA paths for all chains, no pocket or contact constraints and one chain A template from 6OGE. Antibody coordinates from 6OGE are excluded from the input.

| Arm | Template polymer positions | Positions with coordinates | Deposited atom rows |
| --- | ---: | ---: | ---: |
| Complete | 119 | 119 | 912 |
| Masked | 119 | 109 | 836 |

The 6OGE chain A author number equals the P04626 position and the original label index is 22 lower. The cropped label index is 1–119 and the author number stays 511–629. The manifest links all three indices to the query. The masked template removes only atom rows for 603–612. Both arms retain the eight disulfide records with both endpoints in the crop, including Cys600–Cys623. Cys626–Cys634 crosses the crop boundary, so the templates keep Cys626 and record the missing partner in the manifest. Glycans are outside both inference templates. [Original 6OGE mmCIF](https://files.rcsb.org/download/6OGE.cif)

## Run and fixed comparison

Run local `apheris-protenix` 1.3.0 with the Protenix v2 weights, seeds 42, 43 and 44 and two diffusion samples per query, seed and arm. The expected coverage is 48 poses. Seed 42 may be run first as an operational pilot (16 poses, eight per arm). Report all seed and sample rows as they arrive. The comparison remains incomplete until all 48 pass coverage checks.

For each query, arm and seed, select top 1 by the mean of VH–HER2 and VL–HER2 pair ipTM. Ties prefer the lower sample index. The score matrix order is target T, heavy H and light L. The audit verifies chain roles and order against the exact model sequences and the detailed score token order. Global ipTM and VH–VL pair ipTM are recorded as diagnostics. Best reference geometry uses the lower reference interface Cα RMSD of the two samples after the HER2 fit, with the same tie rule. Its rank by confidence shows whether the model preferred it. [Inference admission](template-launch-protocol.md)

Use 1N8Z chain C in P04626 coordinates (deposited residue +22) and Fv chains B (VH) and A (VL). Fit the predicted HER2 Cα coordinates to the 109 resolved positions 511–602 and 613–629. Use the same target frame for target Cα RMSD, Fv Cα RMSD, and interface Cα RMSD across reference target and antibody contact residues. Antibody Cα correspondence is by variable domain position, including the aligned substitution-only CDRH3. The interface RMSD uses the target fit without a second fit on the interface.

Define the contact footprint using heavy atom distances strictly below 4.5 Å. Report recall, precision and Jaccard against the 19 target positions from 1N8Z within the shared 109 positions, plus contacts with the 10 positions outside that common set. Record heavy atom clashes below 1.5 Å for the predicted antibody and target. As a steric check, fit the same 109 predicted HER2 Cα positions to 6OGE A and count antibody clashes against deposited HER2 polymer atoms outside the modeled interval. The overlay excludes glycans, other chains and the overlapping crop. Report directional PAE for the 1N8Z contact residues in both antibody/target directions. [1N8Z mmCIF](https://files.rcsb.org/download/1N8Z.cif), [Prior context audit](her2-target-context.md)

If a pose closer to 1N8Z ranks below top 1, report the geometry and the confidence rank. A similar outcome across masks leaves the missing coordinates unsupported as a sufficient explanation. FACS-low, SPRstrong and SPRweak are structural controls defined by the frozen panel. These runs measure recovery of known geometry, with experimental binding, affinity, epitope competition and thermostability requiring separate evidence.

The operational rule counts top 1 WT poses that simultaneously recover at least 13 of 19 target contacts, fit within 2 Å target Cα RMSD and 4 Å reference interface Cα RMSD and have zero antibody–target heavy atom clashes below 1.5 Å. The arm qualifies when at least two of three seeds pass. Every pose retains its four condition outcomes. The 6OGE overlay is a separate steric diagnostic. [Operational rule](structure-admission-rule.md)

Both 1N8Z and 6OGE occur in the published upstream RFD4 sampling metadata. This experiment is a reproduction control. The metadata do not establish exposure of a particular checkpoint, and the experiment provides no independent generalization estimate. [RFD4 sampling audit](../../results/rfd4-training-metadata/audit.json)

## Reproduce the inputs

```bash
source .venv/bin/activate
curl -fLsS https://files.rcsb.org/download/6OGE.cif -o /tmp/6OGE.cif
python -m scripts.prepare_her2_template_controls --source-cif /tmp/6OGE.cif
python -m pytest -q tests/test_template_controls.py
```

The manifest records the original RCSB file hash, upstream frozen input hashes, cropped file hashes, residue mapping, planned run settings and panel.

## Audit the raw outputs

```bash
source .venv/bin/activate
curl -fLsS https://files.rcsb.org/download/1N8Z.cif -o /tmp/1N8Z.cif
python -m scripts.audit_template_controls \
  --campaign results/campaign/template-controls \
  --reference /tmp/1N8Z.cif --template-source /tmp/6OGE.cif
```

The campaign directory holds `input.sha256`, `image-id.txt` and the six `{arm}/seed-{seed}/` directories with sibling logs. The audit verifies the staged input hashes, receipt keys, runtime template counts, model chain roles and detailed score token order. Results go to `results/template-controls/`. A valid partial campaign writes the verified rows and missing keys and exits with code 2. A malformed receipt or a missing artifact from a completed receipt fails validation.

For the [target-only MSA control](target-msa-protocol.md), pass `--msa-campaign` in place of `--campaign`. The audit expects the complete template, the same four queries with an alignment path on T only, the staged A3M and 24 samples. Its default output is `results/template-controls/target-msa/`.

For a combined seed run, Lightning logs the initialization seed. In the pinned Protenix image, `predict_step` reseeds from each batch's seed. The audit checks the first configured seed against the initialization log and verifies all configured diffusion seeds and samples against the result receipt and artifact paths.
