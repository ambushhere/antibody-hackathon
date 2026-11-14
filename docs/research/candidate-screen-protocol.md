# HER2 candidate refold screen protocol

The frozen panel comprises 43 unique locally eligible de novo VH/VL pairs from the 64 pre-existing team VM RFD4 designs, all 10 optimization pairs from the frozen SPR expanded shortlist and the four unchanged sequence controls. Local sequence eligibility is provisional. The panel uses all qualifying source records and was frozen before candidate refold scores. FACS and previous predicted structures did not affect membership. Candidate identifiers and generation metadata are in the manifest. [Inputs](../../configs/her2-candidate-screen/manifest.json)

Use one local Apheris Protenix v2 run, seed 42 with one diffusion sample per query and detailed scores. Every query uses canonical HER2 P04626 positions 511–629 as T and the submitted VH and VL sequences as H and L. All three MSA paths are null. The template is the byte-identical complete target-only 6OGE crop from the control experiment. The target-only MSA experiment is separate and has no effect on these inputs.

The root operator owns launch. Launch only after the complete, no MSA control arm qualifies under the frozen rule, with at least two of the three WT top poses passing all four reproduction criteria. Use the same image and resolved settings as the baseline, with the planned seed and sample count. Preserve the staged SHA-256 receipt, image ID, combined stdout and stderr, Foundry result and config, and every raw artifact. [Admission rule](structure-admission-rule.md)

## Audit

Validate the frozen source files and staged inputs, the image and runtime provenance, exact query indices and artifact paths, and the matching submitted chain sequences and confidence order. Report each expected row once, including any missing, failed or invalid artifact as unknown, with the reason. Unexpected or duplicate receipt entries, paths or raw artifacts invalidate attribution for the run. A missing expected row or artifact leaves unaffected verified rows auditable.

Report per-chain-pair and global ipTM, directional PAE, heavy atom contact and clash counts, overlap with the 19-position 1N8Z WT footprint and target fit on the 109 shared reference positions. The 6OGE overlay uses the same 109 positions as the control audit and counts antibody clashes against the resolved ECD outside 511–629. Report WT-mapped interface RMSD and directional reference-interface PAE for the substitution-only optimization and control scaffolds. De novo variable domains use their submitted sequence positions. WT antibody mapping and its dependent measures are inapplicable. Predicted-interface PAE summarizes the two directions between the contacting residue sets when the model places the two partners in contact.

Any descriptive Pareto front uses the same axes and units within each candidate track, with known support for every axis and a complete, admitted panel. The common axes are higher mean VH–target and VL–target pair ipTM, greater WT footprint recovery, lower target-fit RMSD and fewer in-crop and outside-crop clashes. Optimization additionally uses lower WT-mapped interface RMSD. Controls are displayed separately. The front gives the non-dominated rows under these axes without weights or a laboratory readiness decision.

This is a single-draw structural screen. Repeat draws and independent checks are required before advancing designs. Model confidence and geometry provide structural evidence. Affinity, competitive binding and melting temperature require experimental measurements. Both structure references occur in the upstream RFD4 sampling metadata, with checkpoint-level exposure unestablished.

## Local reproduction

```bash
source .venv/bin/activate
python -m scripts.prepare_candidate_screen
python -m scripts.audit_candidate_screen \
  --reference /tmp/her2-workers/template-controls-1n8z-source.cif \
  --template-source /tmp/her2-workers/template-controls-6oge-source.cif
```

The auditor defaults to `results/campaign/candidate-screen/` with `result.json`, `config.json`, `run.log`, `input.sha256` and `image-id.txt` at its root. The baseline control audit defaults to `results/template-controls/audit.json`. Derived JSON, CSV and Markdown are in `results/candidate-screen/`.
