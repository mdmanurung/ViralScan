# F-028 — Spike-in grid: flank index v2 retired

**Date:** 2026-10-23 · **Status:** confirmed (grid of 66 conditions, two backgrounds, one source library per virus)
**Run:** Slurm array 25864985 (66/66 COMPLETED), code tag `spikegrid-code-26c1ac5`, rescored with the viral-only
scorer (SPIKE-03, `scores_v2.tsv`). Receipts: `analysis/spikein_detection_limit/results_2026-10-23/` (with `SHA256SUMS`).

## Result
- **Verdict (locked rule, Q3/Q15): retire** in both backgrounds. The flank index's detection threshold equals the old
  index's everywhere (D03 10 = 10 spiked UMIs; R01 ≤3 = ≤3), but it adds up to **12.57** (D03) and **9** (R01) unsupported
  UMIs per unspiked null library, above the 5-UMI limit. Old index: mean 0.33, max 1.
- **False calls:** 0 in every null library for both kb arms.
- **Recovery:** flank recovers 1.13–1.38× (D03) and 1.12–1.23× (R01) the old index's net kb molecules per spiked read
  (`recovery_ratio.tsv`); no ratio below 0.80. The rule weighs the null floor, so the extra recovery does not keep it.
- **STAR arm:** threshold ≤3 (left-censored) in both backgrounds; nulls carry 0–8 non-spike STAR reads (D03 mean 6.7,
  R01 0). Input for STAR-GATE-01 (user decision).

## Limits (also in `report.md`)
One source library per virus (seeds resample reads, not libraries); NH policy for EBV (`--max-nh 2`) unvalidated;
thresholds are grid levels; STAR "absent" lumps too-short with host-filter loss.
