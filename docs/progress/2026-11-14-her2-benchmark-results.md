## 2026-11-14: benchmark and exploratory design results

Implemented and reproduced a leakage-audited FACS benchmark, candidate generator and conservative immunum preflight. Pairwise logistic regression reached test ROC AUC 0.9859 and average precision 0.9737. Its association with measured affinity was weaker, Spearman 0.3473 on 158 quantified SPR variants.

Generated 10,000 statistical proposals and 12 RFD4 designs across optimisation and two exploratory new-antibody arms. Of 504 shortlisted candidates, 88 passed conservative local checks. Official novelty and thermal improvement remain unverified. No designs were submitted to the competition.

Compared OpenFold and Protenix with WT and FACS-low controls. Fixed target residue mapping before contact comparison and recorded an apparent VH–VL assembly problem in a classifier-favoured candidate. Shared evidence retains separate predicted and measured status.

Replaced executable pickle model exchange with JSON parameters. Added bootstrap symlink protection, clearer secret-presence wording and MCP credential endpoint validation. The setup PR and manual GitHub VM check had already passed.
