# covid_viralscan survey — no SARS-CoV-2, anellovirus dominant

**ID**: F-005
**Status**: preliminary (Stage-4 CellRanger overlap + low-alignment caveat pending)
**Date**: 2026-07-02

## Claim

The two 10x 5′ v3 GEX libraries `LUM-SJ-x213-g` and `LUM-SJ-x216-g`, scanned against a
combined human + SARS-CoV-2 + Serratus/anellovirus reference, show **no SARS-CoV-2
(NC_045512.2) signal** in either sample. The dominant viral signal is
**Alphatorquevirus (Torque teno virus / anellovirus)**.

## Evidence

`multimap_evidence.tsv` (host-conservative multimapping), both samples:

| Target | x213-g | x216-g |
|---|---|---|
| SARS-CoV-2 (NC_045512.2_gene1) | 0.0 UMI, 0 cells, `not_detected` | 0.0 UMI, 0 cells, `not_detected` |
| SARS-CoV-1 (sarsp1, neg. control) | 0.0 UMI, 0 cells, `not_detected` | 0.0 UMI, 0 cells, `not_detected` |

`viral_summary.tsv` — top detections:
- x213: Alphatorquevirus 57,138 UMI / 9,611 cells (5.9%); Beta/Gamma/Samektorque + Anelloviridae trace; MPXV_gp132 1 UMI.
- x216: Alphatorquevirus 85,471 UMI / 13,414 cells (6.6%); trace HHV-1/2/6 (≤3 UMI); MPXV_gp132 1 UMI.

The negative control being clean (SARS-CoV-1 = 0) alongside a strong, consistent
anellovirus signal indicates the specificity machinery is working, not that viral detection
is globally failing.

## Caveats (must resolve before firm conclusions)

- **Pseudoalignment is low: 6.4 % (x213) / 8.9 % (x216)** — consistent across samples;
  77–83 M reads aligned (usable) but the low rate is unexplained. Candidates: deep-library
  intronic/ambient fraction, chemistry/whitelist edge, D-list masking. No known-good baseline
  from the same index was locatable this session.
- **total_cells is high: 163,203 / 203,892 barcodes** at `min_counts=1000` on 1.2 B / 0.9 B
  read libraries — pulls in ambient barcodes, deflates `pct_infected`. Stage 4
  (`summarize_survey.py` CellRanger called-cell overlap) is needed to separate real cells from
  ambient before trusting per-cell fractions.
- Clinical expectation for these samples (were they expected COVID-positive?) is unconfirmed.

## Implications

If the libraries were expected COVID-positive, the low pseudoalignment + zero SARS-CoV-2
warrants a data/chemistry check before concluding true-negative. If they are general PBMC/
blood samples, anellovirus dominance with no SARS-CoV-2 is biologically unremarkable.

Related: [[raw-count-thresholds-confound-with-sequencing-depth]] (depth/threshold caveats),
pipeline bugs fixed en route logged in [[learnings.md]].
