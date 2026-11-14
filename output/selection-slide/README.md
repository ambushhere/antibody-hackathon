# Selection and evidence slides

[English presentation](her2-selection-and-evidence-en.pptx), four editable slides:

1. What we selected and why.
2. Tools and data behind our choices.
3. The architecture used for optimization, including position-one-hot SPR ridge and the separate FACS score.
4. The complete selection pipeline, with both tracks converging on eligibility checks and ranking.

The third slide distinguishes the actual selection model from the experimental Transformer. Candidate affinity and thermal improvement remain unmeasured. Sources are recorded in English speaker notes.

The builder is [selection-slide.mjs](../../scripts/presentation/selection-slide.mjs). Use the same runtime environment variables described in the [pipeline build instructions](../pipeline-draft/README.md), with a fresh `OUTPUT_STEM`, and run this builder instead. Omit `SOURCE_PPTX` to generate all four slides. Set it to a two- or three-slide version to preserve those slides and append the remaining diagrams. All four final previews were reviewed; package and layout checks passed.
