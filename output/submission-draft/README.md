# Draft for both challenge tracks

[team-3-draft.csv](team-3-draft.csv) contains the current research pool in the five-column format supplied by Davide Campagna on 14 November 2026. It has 10 intended optimization variants and 43 intended de novo designs. Every row uses `kind=paired`, with separate heavy and light variable domains. The intended track and design origin are in `design_description`, because the supplied schema has no track column. **Abguard assigns the actual track automatically**, using the installed tool’s sequence identity rule; a description does not override classification.

**This is not the final slot selection and has not been submitted.** These 53 rows exceed the PDF's limit of a ranked list of up to 20. The PDFs allow any proportion of optimization and de novo designs in that list. The guaranteed four laboratory slots are shared across both tracks, with extra slots awarded through pitching and judge nominations. Do not submit the entire research pool. No candidate has measured affinity, epitope competition or thermal stability. The installed abguard check is complete; organizer-side acceptance remains unverified. The current structural control setup failed its reproduction criterion, so no new candidate scores from that setup are presented as evidence of readiness.

The local export checks the exact column schema, unique IDs and paired sequences, canonical amino acids, the combined length limit of 250, even cysteine counts in each domain and `immunum` 2.1.0 heavy/light assignment with confidence at least 0.1. [Per-design checks and provenance](team-3-draft.validation.json). Local sequence checks do not reproduce the complete official `abguard` tool or its Adaptyv database search. All ten optimization picks have nearest provided-window edit distance 3 and pass the stricter wording "more than 2" as well as the earlier "at least 2" wording. [PDF rule audit](../../docs/research/submission-rules-audit.md).

## Installed abguard result

[team-3-abguard-pass.csv](team-3-abguard-pass.csv) retains the **52 accepted rows** in the same five-column format. A strict second `team-validate` run accepted all 52. The installed tool assigns 11 `novel` and 41 `de_novo` classes; one RFD4-origin design is classified `novel`. These are tool classes, not proof that a redesign satisfies every PDF interpretation for optimization. The reference includes 610,109 CDRH3s. [Per-design results, software provenance and limitations](../../docs/research/abguard-validation.md).

The 53-row original export and its descriptions preserve the pre-abguard research pool. Its sidecar only describes that exporter's local checks; `not_run_by_this_exporter` does not mean the tool is unavailable. The accepted CSV preserves exact input descriptions, so "official novelty pending" there refers to final organizer acceptance. Neither file is a final ranked ≤20 submission. Sequence order preserves the frozen research panel, not a new final ranking.

The previously unavailable tool was found after the user's installation at `/home/hackathon/.local/bin/abguard` on the VM. It uses `immunum 2.0.0`; the earlier local preflight used 2.1.0. The actual package and bundled reference are archived with hashes. No Slack message or laboratory submission was made.

## Recreate or select a draft subset

```bash
source .venv/bin/activate
PYTHONPATH=.:src python -m scripts.export_submission_draft
PYTHONPATH=.:src python -m scripts.export_submission_draft \
  --ids opt_WGAPGFYLMT \
  --output /tmp/reviewed-subset.csv
```

Use actual IDs from the CSV for a subset; unknown or duplicate IDs fail before any file is written. `team-3` is the VM team label. Replace the filename with the confirmed team name for a final manual submission.

Davide's supplied instructions request one CSV in the same Slack channel as his announcement, tagging [Davide Campagna](https://aixantibodyha-v4i1692.slack.com/team/U0CG5RVQ4SW) in that message. Do not send it to another person or channel. Preparing these files does not perform that action.
