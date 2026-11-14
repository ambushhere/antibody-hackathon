# Frozen candidate inputs

`initial-preflight.json` contains the 504 checked statistical/RFD4 proposals used by the first ranking, with 88 strict local passes. `expanded-preflight.json` contains those plus 239 SPR-prioritised proposals, with 124 total strict local passes. These are selected inputs, not estimates of whole-pool acceptance rate.

The JSON files are byte-identical copies of the run inputs whose hashes appear in the ranking models and prescreen audit. Paths in those original provenance records identify the run location. Use these versioned copies and a caller-supplied organizer data directory to replay the ranking, as described in [the checkpoint](../../docs/EXPERIMENT_CHECKPOINT.md).
