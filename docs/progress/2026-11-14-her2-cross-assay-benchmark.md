## 2026-11-14: cross-assay benchmarks and experiment handoff

Added nested sequence-group SPR benchmarking, candidate ranking and a novelty-distance diagnostic. The primary retrospective ridge correlation is 0.613 versus 0.347 for the frozen FACS score; the distant-subset difference remains uncertain. The expanded local screen has 124 strict passes and ten provisional priorities. No candidate has measured improvement or official eligibility approval.

Preserved de novo identity and target-only refolding results, all 21 completed Protenix context-panel predictions and seven completed RF3 predictions. The structural archive contains 307 files with checksums. RF3 provider usage records USD 0.022406 for this run. Its no-template conditioning differs from the Protenix template-based panel; biological interpretation is still pending.

Corrected the HER2 epitope map with per-structure mmCIF/SIFTS mapping and verified all nine cases. Audited the proposed structural negatives and withdrew unsupported calibration claims. ThermoMPNN-D failed measured control checks. A separate audit of 1,324 public thermal records found insufficient independent format-matched sources for the proposed transfer benchmark, so no thermal model or synthetic Tm ground truth was introduced.

Added the visual primer, source reviews, an experiment index and portable checked candidate inputs. Fixed missing/duplicate refold detection, sequence mismatch detection, nonfinite replay validation, FACS score provenance, split-hash verification and CSV formula escaping. Thermal replay now checks input hashes before importing its upstream entry point and writes to a new output directory.

Validation: 92 tests and 22 subtests passed locally. The root reproduced thermal tables and the thermal data audit, the expanded shortlist order, frozen distance statistics and all nine epitope mappings. Shell scripts passed syntax checks. CI verification and the merge are recorded in the pull request. No competition submission or synthesis order was made.
