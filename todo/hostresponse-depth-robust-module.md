# Depth-robust, %mito-aware `hostresponse` module (v2.5 SH1.1–1.3)

- **Priority**: high
- **Status**: done (2026-07-06)
- **Category**: feature (package)
- **Date**: 2026-07-02
- **Author**: mdmanurung
- **Source**: Feature-completeness gap analysis 2026-07-02 (`.living/decisions.md`);
  findings F-001, F-003.

## Motivation

`hostresponse` is the package's flagship interpretive feature, and today it produces
**depth-confounded results**. It labels a cell virus-positive at `counts >= detection_threshold`
(raw UMI, default 10; `hostresponse.py:500`) and "controls" for depth only by keeping the
top-50%-depth cells and balancing class *counts* (`_balanced_split`). Findings F-001/F-003
showed this does **not** remove the confound — it *widens* it: sequencing depth alone
predicts the EBV label at AUC 0.967 vs the 0.866 model; the honest, depth-independent
effect is AUC ~0.64–0.72. The correct methodology exists but lives entirely in external
scripts under `analysis/hostresponse_ebv_matched/scripts/`, so a package user cannot
reproduce the honest result.

## Scope (three PLAN items, one module)

### SH1.1 — depth-robustness reporting
Add to the `hostresponse` output:
- a **depth-alone** AUC/MCC baseline (logistic on `log(depth)` only) — the null a real
  signal must beat;
- an optional **`log(depth)` covariate** model, reporting metrics before vs after adjustment;
- per-gene **E-values** (Ding & VanderWeele 2016): minimum confounder strength to explain
  away each gene's odds ratio; flag `E<1.5` (fragile) vs `E>3` (robust).
- Source: `depth_confounder_check.py`.
- **Risk (carry over from depth-confounder-check.md):** verify the depth covariate is
  **host-only** UMI, not total (host+viral) — the wrong denominator reintroduces the
  confound it is meant to block.

### SH1.2 — depth-independent label option
- Add a **CPM / fraction** label (`--label {raw,cpm,fraction}`) and/or a **depth-matched
  case-control** split option, replacing raw `>=10` + top-50% balancing as the default for
  publication-grade runs.
- Sources: `depth_matched_reanalysis.py`, `cpm_label_crosscheck.py`.

### SH1.3 — %mito control
- Add `pct_mito` as a covariate/filter, with **leave-one-out MT handling** (MT-ND4L was
  self-suppressing; FDR 4e-7 → 0.07 after control).
- Source: `go_enrichment.py` mito-control path.

## Definition of done
- `hostresponse` emits a depth-alone baseline and a depth-adjusted metric alongside the
  headline AUC/MCC, plus per-gene E-values, in its standard output.
- A `--label {raw,cpm,fraction}` (and/or depth-matched) option reproduces the honest
  AUC ~0.64–0.72 on the EBV showcase run without any external script.
- `%mito` covariate available; MT genes handled leave-one-out.
- Tests: a synthetic case where depth *is* the only signal → depth-adjusted AUC ≈ 0.5
  and E-values ≈ 1 (guards against silently re-confounding).
- Docs + CHANGELOG updated; PLAN SH1.1–1.3 flipped `[x]` in the same commit.

## Notes
- This is post-release (does not block v2.4.0) but is the **#1 correctness gap** in the
  package. Until it lands, any host-response result the package prints should be read as
  depth-inflated.
- Related: [depth-confounder-check.md](depth-confounder-check.md) (the analysis that
  established the confound), finding F-004 (host transcriptome is a continuum, not a
  discrete intermediate state).
