# HER2 candidate screen with target MSA

The alternative panel follows the [predeclared fallback](candidate-fallback-protocol.md). It uses the original ordered 57 queries, comprising four controls, ten optimization candidates and 43 locally eligible de novo VH/VL pairs. The sole query change is the explicit `msa_path` for chain T. H and L retain null paths and the complete target-only 6OGE template is byte-identical to the original. The [variant manifest](../../configs/her2-candidate-msa/manifest.json) pins the parent input, alignment, control input and fallback protocol hashes.

## Admission

Require a complete local target MSA audit with all 24 control poses and no missing entries. It must cover the exact four queries at seeds 42, 43 and 44 with samples 1 and 2, using the complete 6OGE template, frozen target alignment and pinned Protenix v2 image. Select each WT top pose by mean VH–target and VL–target pair ipTM, with the lower sample index breaking ties. Require all four conditions from the [frozen reproduction rule](structure-admission-rule.md) in at least two seeds. Geometry-best poses are diagnostic. Root owns execution and uses seed 42 with one diffusion sample per query.

## Local reproduction

```bash
source .venv/bin/activate
PYTHONPATH=. python -m scripts.prepare_candidate_msa
PYTHONPATH=. python -m scripts.audit_candidate_screen --target-msa --check-admission
PYTHONPATH=. python -m scripts.audit_candidate_screen --target-msa
```

Preparation writes `configs/her2-candidate-msa/`. Admission reads `results/template-controls/target-msa/audit.json` and exits 0 only for verified qualification. The raw candidate audit expects `results/campaign/candidate-msa-screen/` with top-level `result.json`, `config.json`, `run.log`, `input.sha256` and `image-id.txt`. It verifies all four checksums at their exact paths under `/home/hackathon/almaz/work/candidate-msa-screen-20261114/input`. When supplied, `admission.json` ties the run to the control audit and candidate manifest hashes. Timing markers are checked against the Foundry receipt. The default 1N8Z and 6OGE original deposit paths are the same as for the [original candidate auditor](candidate-screen-protocol.md).

The audit always writes 57 expected rows to `results/candidate-msa-screen/`. A verified control failure with no candidate artifacts produces `not_run_admission_failed`, with 57 unknown poses and empty rank fields. Missing control evidence leaves admission pending. Missing individual pose files or receipts leave the affected rows unknown. An unqualified control or invalid run metadata prevents scoring. De novo poses use their submitted antibody residue positions and WT antibody measures remain inapplicable.

A complete admitted run can show the existing track-specific descriptive Pareto fronts. A single seed provides no final candidate promotion. Repeat draws and independent checks are required before advancing designs. Affinity, epitope competition and thermal stability require experimental measurements.
