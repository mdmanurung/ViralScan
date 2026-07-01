# Todo List

Master list of future work items. Each item can have a detailed writeup
in a separate `.md` file in this directory.

## Items

<!-- Add todo items below. Link to detailed writeups as needed. -->

| Priority | Item | Category | Status | Date |
|----------|------|----------|--------|------|
| high | Depth-confounder check: re-run EBV classifier with `log(_raw_depth)`+`log(total_host_UMI)` covariates and compute per-gene E-values → quantify how much of AUC 0.866 survives depth adjustment. See [depth-confounder-check.md](depth-confounder-check.md) (idea 5a) — **DONE 2026-07-01: headline substantially depth-confounded; depth-alone AUC 0.80–0.97; 5/15 genes robust. Follow-up: depth-normalized-label re-analysis.** | validation | done | 2026-07-01 |
| high | Re-analysis with a depth-independent EBV label (CPM/fraction or depth-matched case-control), re-estimating host-response on the 5 depth-robust genes — the principled fix for the depth confound (from depth-confounder-check) | analysis | done | 2026-07-01 |
| medium | Cross-validate the depth-matched AUC ~0.72 with a CPM/fraction EBV label; characterize the 5 depth-robust genes (GO/pathway) — do they form a coherent EBV host-response program? | analysis | done | 2026-07-01 |
| low | Powered GO/pathway enrichment: re-select the depth-robust gene set on a depth-independent label (not just the top 5) then run enrichment; add a mitochondrial-fraction control for MT-ND4L | analysis | open | 2026-07-01 |
| high | Intermediate-attractor test: UMAP+HDBSCAN on host transcriptome; do the misclassified EBV cells form a distinct primed/intermediate state? Tests the convergent hypothesis for the specificity gap (0.744<0.822). See [intermediate-attractor-test.md](intermediate-attractor-test.md) (ideas 1a/4a/4b) | analysis | open | 2026-07-01 |
| medium | EBV host-response: sensitivity sweep over the ≥10 UMI detection threshold (currently a single fixed cutoff) to show metrics are stable — extend with Hill/EC50 dose-response fit per gene (idea 3a) to test whether 10 UMI is near the biological EC50 | validation | open | 2026-07-01 |
| medium | Reference-strategy Selectivity Index: SI = on-target EBV UMI ÷ off-target (HHV-6B+HSV-1) per (aligner×reference) → ranked, quantitative reference recommendation from the existing benchmark TSV (idea 3b) | method | open | 2026-07-01 |
| medium | Reference strategy as a noisy channel: per-pipeline I(true;measured) in bits; exact bit-cost of the HHV-6A/6B ambiguity (BSC) and HSV-1 denominator artifact (idea 7b) | method | open | 2026-07-01 |
| low | Minimal EBV gene panel via mutual-information bottleneck — can 3–5 genes match the 15-gene signature's I(Z;Y)? (idea 7a) | analysis | open | 2026-07-01 |
| low | Consider a pre-commit hook rejecting staged files > ~50 MB (see learnings: 64 GB references/ was not gitignored) | infrastructure | open | 2026-07-01 |
