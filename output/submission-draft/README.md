# Draft for both challenge tracks

[team-3-draft.csv](team-3-draft.csv) contains the current research pool in the five-column format supplied by Davide Campagna on 14 November 2026. It has 10 optimization variants and 43 de novo designs. Every row uses `kind=paired`, with separate heavy and light variable domains. The track and design origin are in `design_description`, because the supplied schema has no track column.

**This is not the final slot selection and has not been submitted.** These 53 rows exceed the currently planned laboratory panel. Their descriptions explicitly retain draft status. No candidate has measured affinity, epitope competition or thermal stability. Official novelty checks remain pending. The current structural control setup failed its reproduction criterion, so no new candidate scores from that setup are presented as evidence of readiness.

The local export checks the exact column schema, unique IDs and paired sequences, canonical amino acids, a conservative combined length limit of 250, even cysteine counts in each domain and `immunum` 2.1.0 heavy/light assignment with confidence at least 0.1. [Per-design checks and provenance](team-3-draft.validation.json). Local sequence checks do not reproduce the complete official `abguard` service or its Adaptyv database search.

## Official tool availability

On 14 November at approximately 16:04 Oslo, `abguard` was not found as a command or Python module in the active local or team VM environments. It was not present as a named local Docker image or in the inspected `/app`, `/opt` and `/usr/local/bin` paths of the supplied Apheris data and HER2-novel RFD4 images. The supplied VM guide mentions filtering with it but provides no command or installation link at that location. This is an availability finding within the inspected environments, not proof that no distributable tool exists. `immunum` 2.1.0 is installed locally and was used for every exported row.

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
