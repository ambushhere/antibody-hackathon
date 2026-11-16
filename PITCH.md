# Pitch summary

[Final two-slide presentation](output/selection-slide/her2-selection-and-evidence-en.pptx) · [Pipeline and evidence](README.md) · [Execution guide](docs/EXECUTION_EN.md)

## Start with measured evidence

We entered the challenge with more experience in software and machine learning than in antibody biology. Our literature review showed why a convincing structure or a high model score is not enough to establish binding. Laboratory testing is expensive, so we focused our optimization effort on evidence from experiments that had already been done.

The supplied FACS and SPR datasets studied variants of trastuzumab, an existing HER2 antibody. We followed their experimental design: keep the antibody framework and vary a ten-residue window within CDRH3. FACS gave us a large screening dataset; SPR gave us fewer, quantitative binding measurements. We trained separate scorers and selected candidates using both signals and sequence diversity.

We tried more complex models, embeddings and synthetic data. Those experiments did not establish a reason to replace the simple final selector. Our optimization shortlist therefore uses pairwise logistic regression for the FACS screening label and positional ridge regression for the SPR ranking signal. This was a choice to stay close to measured evidence, not proof that our new variants would bind better.

## Explore new designs with a relevant generator

For de novo, we used the supplied RFD4-Proteina workflow and antibody/HER2 checkpoints. It allowed us to generate new heavy/light chain pairs in the target binding-site context without training a generator from scratch. Public training metadata helped us understand the model's scope and possible overlap with known structures. It did not establish the exact fine-tuning data or superiority over other generators.

We selected this arm using sequence checks, motif warnings, novelty checks under the installed abguard build and diversity. We did not apply our trastuzumab-specific affinity scorer to new frameworks.

## Spend laboratory slots on a reasoned shortlist

Our final ranked file contains 30 candidates: ten optimization and twenty de novo. All passed the installed abguard build. Three of the first four positions are optimization designs, reflecting our stronger measured-data basis for that track.

The structural controls and thermal pilot did not support candidate ranking, and we report those negative results. The remaining test is laboratory measurement of HER2 binding, competition with trastuzumab and thermal stability. Our deliverable is a reproducible selection pipeline and a queue of hypotheses for those experiments.
