# covid_viralscan survey — no SARS-CoV-2, anellovirus dominant

**ID**: F-005
**Status**: 🔧 RE-RUN VALIDATED — corrected whitelist (job 25140008); per-cell matrix now real (see update)
**Date**: 2026-07-02 (invalidated, then re-run + validated same day)

## ✅ RE-RUN VALIDATION (2026-07-02): corrected whitelist works end-to-end

Re-ran covid quant (job 25140008) with the correct barcode whitelist — CellRanger's raw
barcode universe (2,974,869 barcodes; raw R1 match 68.1% vs 0.4% for v3). The matrix is now
real, validated against CellRanger's called cells for x213 (same sample as the CellRanger run):

- **CellRanger's 28,922 real cells now overlap the ViralScan matrix at 100%** (was 0.5% broken).
- Per-barcode UMI recovered: x216 max 30,125 / 4,985 barcodes ≥1000 UMI (was max 5,311 / 195).
- **Viral signal concentrates in real cells**: among the 28,922 CellRanger cells, 99.7% are
  viral+ (≥1 UMI); among empty droplets only 13.3%. Real cells median 398 total UMI vs empty 1.
  → **The viruses ARE in non-empty droplets**, ~7.5× enriched vs empty. (≥1 UMI is permissive —
  much of the 99.7% is ambient anellovirus; meaningful per-virus rates need a ≥2–5 UMI threshold
  + proper cell-calling, the motivation for the report-both-denominators feature.)

**Net**: the F-005 SARS-CoV-2=0 result holds; the per-cell anellovirus story is now on a valid
matrix and answerable. Fix committed in `slurm_viralscan_quant.sh` (WHITELIST → CellRanger-derived).

## ⛔ INVALIDATION (2026-07-02): wrong 10x barcode whitelist  _(kept for the record — resolved above)_

The per-cell results below are **not trustworthy**. The covid quant used the 10x **v3**
whitelist (`10x_version3_whitelist.txt.gz`), but the data's barcodes do not match it:

- `bustools correct`: **96.5% of BUS records were "uncorrected"** (off-whitelist).
- ViralScan matrix: 163,203 barcodes but **median 1 UMI/barcode, max 5,311, only 195 > 1000 UMI**
  — essentially all empty droplets. CellRanger called **28,922 real cells** from the same FASTQs.
- Only **156/28,922 (0.5%)** CellRanger cells are in the v3 whitelist. Raw R1 barcodes match
  the v3 whitelist **0.4%** (RC 0%), and match **no** bundled ngs_tools whitelist > ~5%
  (best: `10x_version4`/GEM-X 4.7%), yet match CellRanger's called cells **56.5%** directly.

Conclusion: this library's chemistry (likely **GEM-X 5′** or another CellRanger-auto-detected
set) is not covered by the v3 whitelist, so ~96% of reads were dropped and the resulting
per-cell matrix is noise. **Re-run required with the correct whitelist** (source the GEM-X/5′
barcode list CellRanger used, or pass CellRanger's whitelist to `viralscan -w`) before any
per-cell or "% infected" claim. Bulk-level "SARS-CoV-2 = 0 UMI" is the only relatively robust
takeaway, and even that should be re-confirmed post-fix.

---
_Original (now-invalidated) preliminary write-up follows:_

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
