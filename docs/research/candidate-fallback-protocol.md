# Conditional target-MSA candidate screen

14 November 2026, 14:53 UTC. This amendment is recorded after inspecting 40 of the 48 no-MSA control poses and before launching or viewing any target-MSA prediction. Complete-template WT top-1 passed the fixed reproduction gate in only one of three seeds. The no-MSA candidate screen is therefore not admitted. Masked-template WT passed neither of its first two seeds and cannot reach two of three.

The target-MSA control experiment remains unchanged: four controls, the complete 6OGE target template, seeds 42/43/44 and two samples per seed. After all 24 poses and their provenance have been audited, admit an alternative candidate screen only if WT top-1 passes the same pre-existing reproduction gate in at least two seeds. Do not relax the contact, RMSD or clash thresholds, and retain every control result.

If admitted, use the same frozen 57-query panel (four controls, ten optimization candidates and 43 locally eligible de novo candidates), with the sole biological input change being the already frozen target alignment. Preserve a separate input directory and manifest for this variant, including the original panel hashes. Use one sample at seed 42 as a first screen; this does not establish repeatability or laboratory success. The original no-MSA candidate inputs remain unchanged and are not launched.

If the target-MSA control also fails, preserve the prepared candidate pipeline but do not spend the batch on an unqualified configuration or promote candidate scores. This is explicitly control-based method development, not an independent estimate of performance on novel binders.
