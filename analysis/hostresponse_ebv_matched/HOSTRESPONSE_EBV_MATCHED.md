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
| AUC (mean) | 0.845 ± 0.032 |
| MCC (mean) | 0.539 ± 0.059 |
| Sensitivity (mean) | 0.824 ± 0.044 |
| Specificity (mean) | 0.710 ± 0.052 |
| Balanced accuracy (mean) | 0.767 ± 0.030 |
| Stable host genes | 15 (min prob 0.6) |

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
