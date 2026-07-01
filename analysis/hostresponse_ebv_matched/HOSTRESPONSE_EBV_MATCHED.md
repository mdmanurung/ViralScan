# hostresponse_ebv_matched

EBV matched-cell host-response analysis for the manuscript: within a single LCL
sample (SRR12682296), test whether the host transcriptome discriminates
EBV-positive from EBV-negative cells, using paper-defined barcodes to match cells.

Registered into mycelium 2026-07-01 (analysis was already run; this documents it
under the living-repo layer). Conventions applied: `robust-analysis` (validation,
CV-based uncertainty) + `bioinformatics` (scRNA-seq).

## Inputs

- FASTQ / run: `SRR12682296` (see dataset `reference_strategy_fastqs`); ViralScan
  full-depth run dir on the showcase filesystem
  (`.../out_full_depth_wl/lcl_5lines/SRR12682296`).
- Paper barcodes: GEO `GSE158275` (`GSM4796271_LCL_777_B958_UMI_barcodes.tsv.gz`) —
  used to match cells; 10x `-1` suffix normalized before intersection.

## Method

`scripts/hostresponse_ebv_matched.py` (261 lines):
1. Load host + EBV AnnData; normalize barcodes; intersect tool cells with paper
   barcodes (matched-cell design).
2. Label cells EBV-positive at **≥10 EBV UMI** (detection_threshold_umi).
3. `run_hostresponse` (from `viralscan.scripts.hostresponse`) fits a host-gene
   classifier with cross-validation over `DEFAULT_SEEDS`, reporting mean±sd of
   sensitivity / specificity / balanced accuracy / AUC, plus a stability-selection
   gene ranking (min prob 0.6).

Reproduce: `python scripts/hostresponse_ebv_matched.py` (defaults point at the
showcase run dir; see argparse for overrides). Outputs land in
`results/hostresponse_ebv_matched/`.

## Headline results (registered in `outputs/numbers.json`)

| Metric | Value |
|--------|-------|
| Matched cells | 1906 (1179 EBV+ / 727 EBV−) |
| Detection threshold | ≥10 EBV UMI |
| AUC (mean) | 0.866 ± 0.036 |
| MCC (mean) | 0.570 ± 0.077 |
| Sensitivity (mean) | 0.822 ± 0.054 |
| Specificity (mean) | 0.744 ± 0.054 |
| Balanced accuracy (mean) | 0.783 ± 0.038 |
| Stable host genes | 15 (min prob 0.6) |

**Leakage fix (2026-07-01, post-review):** HVG feature selection was moved
*inside* the CV split (train cells only) after the mycelium review found it was
fit on all cells (feature-selection leakage). The corrected numbers above differ
from the pre-fix values (AUC 0.845, MCC 0.539) by **less than one seed-SD on every
metric** — unchanged within noise. HVG is unsupervised, so this leak class is
expected to be negligible; the fix is for methodological correctness, not because
it moves the answer. The 15-gene stable set was unchanged. Pre-fix table:

Matthews correlation coefficient (MCC) added 2026-07-01 by extending
`_run_l2_regression` in `src/viralscan/scripts/hostresponse.py` and re-running the
identical pipeline (fixed `DEFAULT_SEEDS`) on the saved matched h5ads. Sensitivity,
specificity, and balanced accuracy reproduced bit-for-bit; AUC reproduced to 3
decimals (0.8446→0.8449, cross-environment numerical noise); the 15-gene stable set
was identical. MCC 0.539 indicates moderate correlation between predicted and true
EBV status — a more conservative single-number summary than balanced accuracy under
this class balance.

Top stable genes (Ensembl IDs) in `results/.../Epstein-Barr_virus_stability.csv`
and `..._gene_weights.csv`.

## Outputs

- `hostresponse_summary.tsv`, `hostresponse_metrics.csv` — headline metrics (tracked)
- `Epstein-Barr_virus_stability.csv`, `..._gene_weights.csv` — gene rankings (tracked)
- `matched_barcodes.tsv` — matched cell barcodes (tracked)
- `ebv_burden_matched.h5ad`, `host_only_matched.h5ad` — AnnData (**gitignored**;
  regenerate via the script)

## Depth-confounder check (2026-07-01, follow-up — MAJOR caveat)

`scripts/depth_confounder_check.py` (results in `results/hostresponse_ebv_matched/depth_confounder.txt`)
tested whether the AUC 0.866 is real host biology or a sequencing-depth artifact.
`host_depth` here is host-only raw UMI (EBV burden is a separate matrix). Findings:

- **The ≥10-UMI EBV label is strongly depth-dependent**: EBV+ rate rises 27.5%→96.3%
  across host-depth quintiles; Spearman(depth, EBV UMI)=0.53.
- **Depth alone predicts EBV status** at AUC 0.803 (all cells) and **0.967 within the
  classifier's own balanced+depth-filtered design — higher than the 0.866 host-gene model.**
  The top-50%-depth filter widens the between-class depth gap (EBV+ median 34,981 vs
  EBV− 17,868; 3% overlap) rather than closing it.
- **Only 5 of the 15 stable genes** retain a depth-adjusted association (E-value ≥2,
  p<1e-7); the two top-stability genes are explained away by depth (p≈0.8).

**Conclusion**: the headline is **substantially depth-confounded**. A residual, genuine
host-response signal exists in ~5 genes, but it is much weaker than the raw AUC implies,
and the classifier underperforms depth alone. Root cause: the raw-UMI positivity label
(see finding F-003). Principled fix (not yet run): depth-normalized EBV label (CPM/
fraction), a depth-matched case-control design, or depth in the label definition. See
findings F-001 (revised to *contradicted*) and F-003.

## Robustness notes / open questions

- **Uncertainty**: metrics are mean±sd across `DEFAULT_SEEDS` CV splits — good.
- **scilintr**: clean after waiving 2 `unchecked-cache` false positives in
  `_h5ad_path` (input-location resolution, not output caching).
- **Single sample**: results are within one LCL sample; generalization across the
  other LCL lines / SRRs is untested here.
- **Threshold sensitivity**: EBV+ defined at exactly ≥10 UMI — a sensitivity sweep
  over the threshold would strengthen the claim (candidate follow-up, see todo).
- **Specificity (0.71)** is notably lower than sensitivity (0.82) — worth noting in
  any report.
