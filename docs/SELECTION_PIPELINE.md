# Current selection pipeline

The pipeline supports local eligibility checks, retrospective affinity prioritisation in the trastuzumab background and structural failure diagnosis. It does not yet validate a reliable structural ranking of new binders. The latest WT cofolding controls failed pose recovery, and thermal evidence for new candidates remains missing.

```mermaid
flowchart TD
    A[WT sequence, 1N8Z and HER2 epitope]:::data --> B[Generate optimisation and de novo candidates]:::action
    B --> C{Local eligibility gates}:::gate
    C -->|Fail| X[Exclude and retain reason]:::blocked
    C -->|Pass| D[Optimisation: SPR and FACS priorities<br/>De novo: do not transfer the H3 predictor]:::action
    D --> E[Predict complexes alongside WT and measured controls<br/>Match target, template policy and sampling budget]:::action
    E --> F[Check chains, target fold, VH/VL assembly,<br/>clashes, epitope contacts and interface uncertainty]:::action
    F --> G{Does the method recover<br/>the reference pose?}:::gate
    G -->|Currently fails| H[Revise prediction setup<br/>No candidate promotion from these scores]:::blocked
    H --> E
    G -->|Supported by controls| I[Repeat seeds and compare another model<br/>Check pose and rank sensitivity]:::future
    I --> J[Diverse shortlist with separate evidence axes<br/>Official novelty and thermal evidence pending]:::future
    J --> K[Matched SPR, competition and thermal assays]:::future
    classDef data fill:#e7f0f8,stroke:#215987,color:#132639;
    classDef action fill:#e7f2ec,stroke:#276856,color:#132639;
    classDef gate fill:#f7edd2,stroke:#986c1c,color:#132639;
    classDef blocked fill:#faebdd,stroke:#b64c20,color:#132639;
    classDef future fill:#f1f3f6,stroke:#6b7785,color:#132639,stroke-dasharray:4 3;
```

Local gates cover format, length, numbering, cysteines and known-sequence distance. Passing them does not establish acceptance by the unavailable official database gate. Interface scores are not KD. Structural epitope overlap does not prove competition. There is no calibrated universal score combining these axes.

The matched domain IV Protenix WT recovered 3/19 footprint positions, while `YTGDGFPAMN` recovered 16/19 in one draw. That attractive candidate pose is a follow-up hypothesis, not sufficient evidence for promotion when the controls fail. Full ECD did not rescue the panel. RF3 recovered more WT contacts but distorted the target. [Structural audit](research/context-results.md).

Pairwise ridge and frozen AbLang2 provide additional retrospective affinity hypotheses, with no established improvement over the simple position baseline in the distant subset. [Comparison](research/ablang-benchmark.md). ThermoMPNN-D did not support thermal ranking. [Thermal controls](../results/thermal-pilot/README.md).
