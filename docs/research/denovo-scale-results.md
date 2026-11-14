# De novo scale audit

## Scope and status

This is a local sequence and provenance audit of the two completed RFD4 paired-Fv arms. It does not establish binding, trastuzumab epitope competition, thermal stability, or official Adaptyv novelty. The official novelty database was unavailable, so every design remains `unverified_database_unavailable`; no design is marked as an official pass. No SPR score was transferred to these designs.

The audit retained all 64 generated designs and their individual CIF SHA-256 values in [`designs.csv`](../../results/denovo-scale-audit/designs.csv). Run-level provenance and aggregate counts are in [`summary.json`](../../results/denovo-scale-audit/summary.json). Raw campaign directories were read only.

## Provider provenance

Both provider result manifests have start and end times, declare 32 designs and agree with exactly 32 copied generated CIFs. The provider-declared paths and local file paths match exactly in each arm. The local design settings match the provider resolved settings after omitting the `schema_version` metadata field. Input hashes match the local full-sequence design file. All 64 rows carry the structure hash, input hash, settings hash, image digest, model and seed.

| Arm | Checkpoint | Seed | Provider interval (UTC) | Coverage | Full-sequence input SHA-256 | Settings SHA-256 | Image SHA-256 | Result manifest SHA-256 |
|---|---|---:|---|---:|---|---|---|---|
| `denovo-scale-antibody-k1` | `antibody-2026-10-30` | 20261118 | 14:02:41–14:14:02 | 32/32 exact | `4acf71167985f6079fbc63d5d348d03b419520efe11ba07e6489f8b5e5b66a72` | `23830242ea2d1b90ff25b88f95e912266a0c84f5d450fe7d075276abc860240f` | `05690f7d570eaf29f5a0e344d0818e4e318e369ef7c8f2ef950dff791fa2036a` | `13f29955d8e902b9bb48826cb5a8166e58e012e54e02be1b90437d1ca0d3b6c4` |
| `denovo-scale-her2novel-k1` | `her2-novel-2026-11-09` | 20261119 | 14:14:04–14:25:27 | 32/32 exact | `4acf71167985f6079fbc63d5d348d03b419520efe11ba07e6489f8b5e5b66a72` | `a5d7bee7201f04839b3d31a9c03cc297b6e3ba740452d31fef9c44616b495d24` | `49d8ccd71fb480fb1edd83d6ec511032aefdbe02085ada398b6a6129f2757322` | `a5febd7f14b554c3d45fa4867595c6833bbdfd18d9c82cb0f70ac7ea6bb9672d` |

Both runs used `apheris-rfd4-proteina/0.3.0`, the same full-sequence input and `sequence_tsr_k=1`. The seeds differ and the checkpoint image digests differ. This is not a matched causal checkpoint experiment, so differences below are descriptive and cannot be attributed to checkpoint choice alone.

Each generation log also contains two RDKit atom-valence sanitization warnings and one skipped atom-name annotation warning. The provider manifests and 32 output files per arm confirm run completion, but this audit does not assess the effect of those warnings on the structures.

## Sequence audit

Each CIF was read on CPU. Antibody chains were selected by the explicit structure chain IDs `H` and `L`; each was numbered independently with `immunum` IMGT at minimum confidence 0.1. All 64 heavy chains mapped to H and all 64 light chains mapped to K. Numbering succeeded for all 64 designs. The HER2 target chain was excluded even where its sequence could receive an immunoglobulin-like classification.

Paired identity uses the existing `audit_denovo_identity.py` definition: matching residues at the same IMGT position divided by the union of numbered positions, summed over H and L. The local strict threshold is below 70% to both the trastuzumab and adalimumab references. The cysteine screen requires exactly two cysteines per chain, at IMGT positions 23 and 104. The length filter counts the complete generated VH plus VL sequences and allows at most 250 residues.

| Arm | Designs | Paired length range | Canonical cysteines | Paired identity range to trastuzumab | Paired identity range to adalimumab | Local sequence screen pass |
|---|---:|---:|---:|---:|---:|---:|
| `denovo-scale-antibody-k1` | 32 | 227–241 | 21/32 | 60.6–71.7% | 59.8–73.8% | 18/32 |
| `denovo-scale-her2novel-k1` | 32 | 224–243 | 25/32 | 55.6–69.4% | 57.0–71.0% | 25/32 |

Every design had an IMGT-segmented heavy-chain CDR3. Its edit distance from the trastuzumab IMGT CDR3 was at least 6 in the first arm and at least 9 in the second, so all passed the local minimum of two changes. CDR3 identity is populated only when candidate and reference CDR3 lengths are equal; otherwise it is left blank rather than treating unaligned residues as positional matches. The edit distance is retained for every numbered candidate.

There were 64 distinct full VH/VL pairs and no duplicate designs across the two arms. After cross-arm deduplication, 43 unique molecules pass the local sequence filters and have complete provider provenance. This is a local shortlist only. The official Adaptyv database and validator remain necessary before any eligibility claim or submission.

For simple sequence diversification, the audit reports the sum of VH and VL Levenshtein distances, keeping chain boundaries separate. Across the 43 locally eligible unique pairs, it calculated 903 pairwise distances: minimum 76, median 96, mean 96.38 and maximum 124. Within-arm summaries are in `summary.json`. These values describe sequence diversity only, not binding or ranking.

## Limits and handoff

Generated CIF availability and raw structure agreement are not proof of a correctly folded complex or binding. This audit did not calculate interface quality, epitope preservation, competition or stability; the structural panel remains separate. The report makes no affinity claim and no binding prediction.

Reproduce the audit from the repository root with the active project environment:

```bash
source .venv/bin/activate
python -m scripts.audit_denovo_scale \
  --campaign-root results/campaign/denovo-scale \
  --input results/campaign/campaign-20261114/denovo-input/full_sequence_design.json \
  --reference configs/adalimumab_reference.json \
  --output-dir results/denovo-scale-audit
```
