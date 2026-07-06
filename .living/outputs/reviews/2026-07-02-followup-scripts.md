# Verification review — EBV follow-up analysis scripts — 2026-07-02

**Scope**: correctness of the 5 load-bearing follow-up scripts (depth_confounder_check,
depth_matched_reanalysis, cpm_label_crosscheck, go_enrichment, threshold_hill_sweep).
Blind adversarial reviewer, focused on bugs that would change conclusions.

## Verdict: conclusions hold; 2 real issues fixed, minor caveats added

- **All CV loops are leakage-free** (StandardScaler fit on train folds only); depth
  matching is mechanically correct; depth-alone AUCs (0.803 / 0.967) verified.

## Findings and dispositions

- **[Major → FIXED, and conclusion CONFIRMED] %mito circularity for MT-ND4L.** The %mito
  covariate included MT-ND4L, so controlling for it when testing MT-ND4L was partly
  self-suppressive. Recomputed %mito with MT-ND4L excluded (leave-one-out): MT-ND4L still
  drops (FDR 4e-7 → 0.067). So the mito-QC verdict is real and now non-circular. Logged as
  a learning (composite-covariate circularity).
- **[Major → FIXED] cpm_label_crosscheck printed "agreeing near ~0.7"** but the computed
  CPM AUC is 0.636. The downstream report/findings already used the honest 0.64; fixed the
  script's hardcoded text to state the methods disagree on magnitude (0.64 vs 0.72).
- **[Minor → CAVEAT ADDED] depth-matched 0.72 is mildly optimistic** — the gene panel was
  selected on all cells (overlapping the matched subset); a fully nested selection would
  land slightly lower, still >> the 0.48 control. Noted in the report.
- **[Minor → FIXED] "depth-alone ~0.5 across thresholds" overstated** — it rises to ~0.59
  at stringent (30%) prevalence; corrected the threshold-sweep note.
- **[Minor] the 0.866 comparator** is the leakage-corrected headline (from
  hostresponse_summary.tsv), distinct from the older manuscript_draft.md's 0.845; same
  balanced+depth-filtered design, so the depth-vs-host comparison is apples-to-apples.
- **[Negligible] t-test df = n−2** ignoring covariates (Δ<0.01% at n=1906);
  **[Minor] double BH filter** is conservative (false-negatives only). No action.

## Net effect
No headline conclusion changed. The MT-ND4L-artifact claim is now cleanly established;
the CPM/threshold claims are stated more precisely. Report at 7pp, 0 drift.
