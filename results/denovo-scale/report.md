# De novo scale-up, 14 November 2026

Exploratory arm. These are predicted designs with local filters only. No official novelty
validation, no structural scoring, no measured binding.

## What was run

Full-sequence Fv redesign on the adalimumab framework against the HER2 domain IV crop,
`sequence_tsr_k=1`, 32 replicates per checkpoint, on the team VM H100.

| Run | Checkpoint | Seed | Designs |
| --- | --- | --- | --- |
| `denovo-scale-antibody-k1` | `antibody-2026-10-30` | 20261118 | 32 |
| `denovo-scale-her2novel-k1` | `her2-novel-2026-11-09` | 20261119 | 32 |

Condition spec `denovo-input/full_sequence_design.json`, sha256
`4acf71167985f6079fbc63d5d348d03b419520efe11ba07e6489f8b5e5b66a72`, identical to the
`denovo-fullseq-k1` pilot input. Hotspots are four domain IV side-chain atoms
(D/593/NZ, C/573/CD1, D/591/CG2, C/572/CD). Recycle strategy `high/front_loaded`.
Launcher `work/run_denovo_scale.sh` on the VM records image digest, input hash and
parameters per run.

## Local filter yield

| Checkpoint | n | <70% to both local refs | Canonical Cys only | Length <=250 | All local filters |
| --- | --- | --- | --- | --- | --- |
| `antibody-2026-10-30` | 32 | 28 | 21 | 32 | 18 |
| `her2-novel-2026-11-09` | 32 | 31 | 25 | 32 | 25 |

Paired lengths span 224–243 residues, so the 250-residue limit binds on none of them,
but it leaves no room for a scFv linker. Immunum confidence is at or above 0.70 for every
chain, well clear of the 0.1 floor.

The HER2 preference-trained checkpoint produced both lower identity and fewer
non-canonical cysteines than the general antibody checkpoint. This reverses the
expectation that HER2 conditioning would pull designs back toward known binders. It is a
single paired comparison at one seed each and one sampling setting; it is not a general
ranking of the two checkpoints.

CDRH3 differs from trastuzumab `WGGDGFYAMD` in all 64 designs, so the de novo requirement
of at least two CDRH3 changes is met trivially here. CDRH3 lengths run 11–22 residues.

Liability motif counts (`src.her2_candidates.liability_motifs`) are 2 for every passing
design, equal to the trastuzumab Fv reference count of 2. This is a sequence motif tally,
not a developability measurement.

## Cross-check against PR #6

PR #6 (`74a3c33`, merged 14 November 2026 at 14:15 UTC) lands two results that bear
directly on this arm.

**1N8Z is in the RFD4 published sampling table.**
[audit.json](../rfd4-training-metadata/audit.json) pins NVIDIA's `RFD4-Proteina` revision
`aa86df1eb9e3529a0f39279f216dce3a95bb3dc2` and finds both 1N8Z (two entries, chains A_1 and
B_1 against antigen C_1) and 6OGE (four entries) in `immune_interfaces.parquet`, both
passing the bond and ligand filters. The HER2 crop in our condition spec comes from 1N8Z,
and trastuzumab is one of the two references we measure identity against.

The consequence for this arm is specific: **low IMGT identity to trastuzumab is not
evidence that a design is independent of the model's training.** The model has seen the
trastuzumab–HER2 complex. Our identity numbers remain a rule-compliance diagnostic; they
are not a novelty-from-training argument, and the report above should not be read as one.
The audit's own wording is that presence in a sampling table does not establish exposure of
any particular checkpoint, so this is a caution, not a demonstrated leak.

The same audit shows `distilled_examples.parquet` holds 13,652 paired Fv entries, mostly
patent-derived. That is a far better local novelty reference than our two antibodies, and
it is now used — see the next section.

**The structural controls in PR #6 did not reproduce the known pose.** In
[scores.csv](../structural-replicates/scores.csv), across seeds 42/43/44 of the OpenDDE
replicates, `wt_contact_criterion_met` is False in all 21 rows. `control_wt` — trastuzumab,
whose pose is known from 1N8Z — recovered 2, 0 and 5 of 19 reference contacts
(footprint recall 0.11, 0.00, 0.26), with antibody CA RMSD in the target frame of 60.3,
42.1 and 34.3 Å. In seed 42 the de novo design `denovo_fullseq_k1_1` scored a higher
antibody–target pair ipTM (0.624) than WT (0.365), and the same design ranged 0.22–0.62
across the three seeds.

Two things follow for the scoring run below. A single-seed structural score cannot rank
these designs, so the Protenix batch uses three seeds. And the controls are not decoration:
if WT and the FACS-low control do not separate in the Protenix output either, then the
run yields no ranking, and no candidate may be promoted on it. That is the prespecified
reading, fixed before the scores were seen.

## Novelty against the RFD4 distillation table

The two pinned tables were downloaded to the ignored cache `data/external/rfd4-metadata/`.
Both sha256 digests match the values in `scripts/audit_rfd4_training_metadata.py`, so this
is the same revision `aa86df1e` that PR #6 audited.

All 13,652 paired Fv entries parse. Variable domains are cut at the `fw_start_seq_idx` and
`fw_end_seq_idx` boundaries in `cdr_annotations`; the raw `ig_seqs` strings are longer than
an Fv and comparing against them would understate identity. Reference VH run 90–143
residues, VL 72–121. Nearest neighbours of our designs come from patents (21), literature
(18) and REAP (4).

**Sixteen of the 43 designs that cleared the two-reference check do not clear this one.**

| Checkpoint | Local filters | Below 70% vs 13,652 paired Fv |
| --- | --- | --- |
| `antibody-2026-10-30` | 18 | 7 |
| `her2-novel-2026-11-09` | 25 | 20 |
| Total | 43 | 27 |

Two things follow. The two-antibody check was optimistic: Spearman between it and the
larger reference is 0.83, correlated but not a substitute, and the absolute threshold moves
a lot. And the gap between the checkpoints widens — 80% of the HER2 preference-trained
designs survive against 39% from the general antibody checkpoint. Combined with the earlier
yield, `her2-novel-2026-11-09` is the better de novo source in this setup, which is the
opposite of what the identity argument predicted before any of this was run.

Most novel designs by this measure are `hn_12` (0.624), `ab_6` (0.644), `hn_4` (0.647),
`hn_10` (0.654) and `hn_14` (0.657). Full table in
[training-novelty.csv](training-novelty.csv).

Limits. The search is a 4-mer Jaccard prefilter over both chains followed by exact global
BLOSUM62 alignment of the top 400 entries per design, because the full 43 × 13,652 × 2
alignment is about 78 hours. A reported maximum is therefore a lower bound if the prefilter
missed a closer entry.

That check was run: [training-novelty-shortlist2000.csv](training-novelty-shortlist2000.csv)
repeats everything at a shortlist of 2,000. Three designs of 43 found a closer entry, the
largest improvement is 0.0043 in paired identity, the nearest entry is unchanged for 40 of
43, and the set below 70% is identical — no design crosses the threshold either way. The
400-entry prefilter is adequate here. Join the two files on `(arm, design)`, not on
`design`: design names repeat across the two arms, and joining on the name alone silently
produces 67 rows and a fictitious maximum delta of 0.062.

The identity definition is alignment-column identity over the paired chains, which is not
the official type-aware definition.

Most importantly: 13,652 predicted-structure entries from one model's distillation set are
not Adaptyv's antibody database. Clearing this check is not clearing the official
validator, and `official_novelty` stays `unverified`.

## What this does not establish

Official novelty remains `unverified` in [identity.csv](identity.csv). The <70% figures are
IMGT-position identity against trastuzumab and adalimumab only, computed locally. The task
rule is type-aware identity below 70% against Adaptyv's antibody database, which no local
surrogate reproduces. Treat eligibility as pending, never passed.

Nothing here speaks to binding, to preservation of the domain IV epitope, or to thermal
stability.

## Structural scoring batch

45 queries: the 43 locally passing designs plus `control_wt` and `control_facs_low`, taken
verbatim from the pilot query file so the HER2 crop, chain layout and format are identical
across every row. Target-only template (`her2_target_refold.cif`); no antibody pose is
supplied. Protenix v2 with seeds 42, 43 and 44, one diffusion sample, detailed scores on;
OpenFold3 with seed 42 as the second voice and the sponsor requirement. Launcher
`work/run_scale_refold.sh` on the VM, outputs `scale-protenix` and `scale-openfold3`.

Scoring computes antibody-to-HER2 contacts explicitly and excludes VH–VL, per
[CAMPAIGN.md](../../docs/CAMPAIGN.md); the generic Fv interface score reflects VH–VL.
Nothing is rejected on CDRH3 RMSD alone — organizer docs state Protenix fails to reproduce
WT CDRH3.

### Result: the controls did not separate, so this batch yields no ranking

Both runs completed: Protenix 135 predictions (45 queries × 3 seeds), OpenFold3 225
(45 × 1 seed × 5 diffusion samples). The sponsor requirement to include Apheris OpenFold is
met. The prespecified reading was fixed before the scores were seen, and it fails.

| Metric | Model | `control_wt` (known binder) | `control_facs_low` (negative) |
| --- | --- | --- | --- |
| Antibody–target ipTM | Protenix v2 | 0.376 [0.363–0.384] | **0.438** [0.375–0.556] |
| WT footprint recall | Protenix v2 | 0.456 [0.211–0.579] | 0.351 [0.211–0.579] |
| Antibody–target ipTM | OpenFold3 | 0.121 [0.078–0.183] | **0.175** [0.118–0.278] |
| WT footprint recall | OpenFold3 | 0.211 [0.000–0.421] | **0.442** [0.211–0.684] |

The negative control scores higher than trastuzumab on three of the four measures. The one
measure where WT leads, Protenix footprint recall, has ranges that overlap completely: both
controls span 0.211 to 0.579.

Three further checks point the same way.

Replicate noise swamps the control gap. The three-seed range of antibody–target ipTM has a
median of 0.100 across the 45 Protenix queries, and **31 of 45 queries vary more between
their own seeds than the entire WT-to-negative-control mean gap of 0.062**. One design,
`ab_18`, ranges 0.506.

The two models do not agree. Spearman between mean Protenix and mean OpenFold3
antibody–target ipTM over the 43 designs is 0.243.

VH–VL dominates, as the organizer docs warn. Mean VH–VL ipTM is 0.478 against 0.344 for the
antibody–target interface, which is why the pair matrix is split by chain here rather than
read as a single Fv interface score.

**No candidate is promoted on these scores.** Eleven designs have a mean antibody–target
ipTM above WT and five above the negative control, but with controls inverted and seed
spread this large, those orderings carry no information. Per-prediction values are in
[protenix/scores.csv](protenix/scores.csv) and [openfold3/scores.csv](openfold3/scores.csv),
joined with novelty in [combined.csv](combined.csv), so the batch stays reproducible and
auditable; it is a negative result about the scoring setup, not about the designs.

One thing did improve on PR #6's OpenDDE replicates: there, WT recovered 0–5 of 19 contacts
and never met the contact criterion. With Protenix and a target-only template WT recovers a
mean 8.7 of 19. The pose is better; the discrimination still is not there.

What this does not license: concluding that the designs do not bind. The controls show the
measurement cannot tell binder from non-binder here, which is a statement about the assay,
not about any sequence.

## Reproduce

```bash
ssh team3 '/home/hackathon/almaz/work/run_denovo_scale.sh <tag> \
  /home/hackathon/almaz/out/campaign-20261114/denovo-input/full_sequence_design.json \
  <new output dir> 32 <seed> 1'
```

Then copy `*_generated.cif` under `results/campaign/<arm>/` and run:

```bash
PYTHONPATH=. .venv/bin/python scripts/audit_denovo_identity.py \
  --campaign results/campaign --arms <arm> ... --output results/denovo-scale/identity.csv
```

Raw VM output is `/home/hackathon/almaz/out/campaign-20261114/denovo-scale-*`. Export it
before VM shutdown; `results/campaign/` is excluded from Git.
