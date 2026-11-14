# Draft for both challenge tracks

**Current review file:** [30 ranked designs](team-3-ranked-draft30.csv), with [selection policy and row diagnostics](team-3-ranked-draft30.selection.json). It contains **10 optimization and 20 denovo**, ordered by priority, repeated as `priority=NN` in `design_description`. The first four contain three optimization and one de novo. **All 30 exact rows passed strict installed abguard validation again**; [receipt](../../results/abguard-ranked/receipt.json). One declaration/class mismatch was excluded. This is an unsent review draft, not a claim that these are the best binders.

The optimization ordering preserves the earlier FACS/SPR Pareto-diversity shortlist. De novo ordering uses motif warnings, paired-sequence diversity and database distance only; there is no validated affinity or Tm ranking for that arm. [Full policy](../../docs/research/ranked-draft-policy.md).

[team-3-draft.csv](team-3-draft.csv) contains the current research pool in the six-column format supplied by Davide Campagna on 14 November 2026. It has 10 intended optimization variants and 43 intended de novo designs. Every row uses `kind=paired`, with separate heavy and light variable domains. The organizer update at 16:33 requires an explicit sixth column `track`, with `optimization` or `denovo`. The current export follows it. `kind=paired` describes molecular format, not track. Installed abguard still produces its own novelty classes separately; these do not replace the declared track.

**This is not the final slot selection and has not been submitted.** The latest organizer update requires at least 28 designs, recommends about 30 and explicitly takes the first N rows. It supersedes the earlier PDF limit of 20. Both tracks belong in the same ranked file. The guaranteed four laboratory slots are shared across both tracks, with extra slots awarded through pitching and judge nominations. Do not submit the entire research pool. No candidate has measured affinity, epitope competition or thermal stability. The installed abguard check is complete; organizer-side acceptance remains unverified. The current structural control setup failed its reproduction criterion, so no new candidate scores from that setup are presented as evidence of readiness.

The local export checks the exact column schema, unique IDs and paired sequences, canonical amino acids, the combined length limit of 250, even cysteine counts in each domain and `immunum` 2.1.0 heavy/light assignment with confidence at least 0.1. [Per-design checks and provenance](team-3-draft.validation.json). Local sequence checks do not reproduce the complete official `abguard` tool or its Adaptyv database search. All ten optimization picks have nearest provided-window edit distance 3 and pass the stricter wording "more than 2" as well as the earlier "at least 2" wording. [PDF rule audit](../../docs/research/submission-rules-audit.md).

## Installed abguard result

[team-3-abguard-pass.csv](team-3-abguard-pass.csv) retains the **52 accepted rows** in the same six-column format. A strict second `team-validate` run accepted all 52. The six-column accepted pool declares **10 optimization and 42 denovo** rows. The installed tool assigns 11 `novel` and 41 `de_novo` classes; one RFD4-origin design is classified `novel`. These are tool classes, not proof that a redesign satisfies every PDF interpretation for optimization. The reference includes 610,109 CDRH3s. [Per-design results, software provenance and limitations](../../docs/research/abguard-validation.md).

The 53-row current export preserves the research pool, now with the required track field. The original five-column input remains under results/abguard-validation/raw for reproducibility. Its sidecar only describes that exporter's local checks; `not_run_by_this_exporter` does not mean the tool is unavailable. The accepted CSV preserves exact input descriptions, so "official novelty pending" there refers to final organizer acceptance. The six-column CSV also passed strict validation (52/52), but this installed version silently drops `track` from its output. Our exporter preserves the declared field. `denovo-scale-her2novel-k1_30` has declared denovo but internal novel class and remains flagged for review. Neither research pool is the final ranked submission. Sequence order preserves the frozen research panel, not a new final ranking.

The previously unavailable tool was found after the user's installation at `/home/hackathon/.local/bin/abguard` on the VM. It uses `immunum 2.0.0`; the earlier local preflight used 2.1.0. The actual package and bundled reference are archived with hashes. No Slack message or laboratory submission was made.

## Recreate or select a draft subset

```bash
source .venv/bin/activate
PYTHONPATH=.:src python -m scripts.export_submission_draft
PYTHONPATH=.:src python -m scripts.export_submission_draft \
  --ids opt_WGAPGFYLMT \
  --output /tmp/reviewed-subset.csv
```

Use actual IDs from the CSV for a subset; unknown or duplicate IDs fail before any file is written. `team-3` is the VM team label. Replace the filename with the confirmed team channel name for a final manual submission.

Davide's supplied instructions request one CSV in the team’s own Slack channel, tagging [Davide Campagna](https://aixantibodyha-v4i1692.slack.com/team/U0CG5RVQ4SW) in that message. Do not send it to another person or channel. Preparing these files does not perform that action.

To reproduce the ranked draft without overwriting it, use a new output filename:

```bash
source .venv/bin/activate
python -m scripts.prepare_ranked_draft --output /tmp/reviewed-ranked30.csv
```
