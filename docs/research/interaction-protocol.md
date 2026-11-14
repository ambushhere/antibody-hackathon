# Interaction experiment protocol

14 November 2026. Recorded before the new runs. This is retrospective method development using an already inspected dataset, not a new blind validation set. The user supplied an external-chat observation of pairwise ridge Spearman 0.648 overall and 0.472 on the existing distant subset. Those figures are unverified here and must not drive feature or hyperparameter selection.

## Fixed comparisons

Use the existing 158-variant primary SPR cohort, source hashes and outer sequence-group folds. Preserve all WT/background exclusions and reportable-KD conditioning. Keep FACS probability and position ridge as frozen comparators. All supervised tuning, centering, scaling and any feature reduction occur inside the corresponding training fold.

1. Position-plus-pair one-hot ridge using the same amino-acid feature definition as the FACS code. Tune alpha over 0.1, 1, 10 and 100 by inner group-fold pooled RMSE. Do not tune on the 75-variant distant diagnostic.
2. Frozen AbLang2 paired VH/VL representations with a small ridge head, subject to verified public code, weights, license and exact token-to-residue mapping. Prespecify CDRH3 mean pooling as the primary representation. Preserve the model revision, weight checksum and extraction inputs; no fine-tuning on SPR.
3. A structural-prior ablation using the same frozen residue representations and a prespecified WT contact mask. Determine contacts from the deposited WT structure without SPR labels. This is a WT-context hypothesis, not proof that generated loops preserve those contacts. Compare with the same representation without the mask; do not select the mask from observed test performance.

Report pooled OOF Spearman and log-KD RMSE, paired sequence-group bootstrap differences, fold-level results and the original 75-variant distance subset. A difference interval crossing zero is inconclusive. Report all planned variants, including failures. Do not pick a winner solely by the highest observed correlation among these correlated comparisons. Existing candidate ranks remain unchanged until reviewed.

## Separate structural experiment

Analyse the completed 21 Protenix context predictions against the same seven queries and the seven completed RF3 predictions. Validate exact sequences, chain identities, sample counts and source hashes first. Map target contacts to canonical P04626 and separate VH/VL assembly from antibody/target interfaces. Compare WT and quantified SPR strong/weak controls before interpreting candidate geometry. RF3 has no supplied template or MSA, while Protenix uses a target template. One seed cannot establish sampling reliability or measured binding.

No competition submission, synthesis order or new thermal ground-truth label is authorised by this protocol. FACS and SPR labels do not measure Tm. Root owns GPU and paid-provider admissions; workers may perform local CPU analysis and public source inspection within their assigned scopes.
