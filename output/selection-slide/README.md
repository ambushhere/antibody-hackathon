# Selection and evidence slides

[English presentation](her2-selection-and-evidence-en.pptx), three editable slides:

1. What we selected and why.
2. Tools and data behind our choices.
3. The architecture used for optimization, including position-one-hot SPR ridge and the separate FACS score.

The third slide distinguishes the actual selection model from the experimental Transformer. Candidate affinity and thermal improvement remain unmeasured. Sources are recorded in English speaker notes.

The builder is [selection-slide.mjs](../../scripts/presentation/selection-slide.mjs). Use the same runtime environment variables described in the [pipeline build instructions](../pipeline-draft/README.md), with a fresh `OUTPUT_STEM`, and run this builder instead. Omit `SOURCE_PPTX` to generate all three slides. Set it to a two-slide version to preserve those slides and append the architecture slide. All three final previews were reviewed; package and layout checks passed.
