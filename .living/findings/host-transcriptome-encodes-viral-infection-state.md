---
topic: host-transcriptome-encodes-viral-infection-state
description: Whether a cell's own (host) gene expression carries enough signal to predict its intracellular viral infection status at single-cell resolution.
created: 2026-07-01
last_updated: 2026-07-01
status: active
---

# Host transcriptome encodes viral infection state

## F-001: A host-gene classifier predicts per-cell EBV status well above chance in one LCL sample
**Status:** preliminary
**Claim:** Within a single Epstein–Barr-virus-transformed lymphoblastoid cell line sample (SRR12682296, GSE158275), a cross-validated L2 logistic classifier trained only on host (non-viral) highly-variable genes distinguishes EBV-positive from EBV-negative cells (EBV-positive defined as ≥10 viral UMI) at mean ROC AUC 0.866 and mean Matthews correlation coefficient 0.570 ± 0.077 across 6 seeds (n=1906 cells; 1179 positive / 727 negative), with the signal concentrated in 15 stably-selected host genes. Sensitivity (0.822) exceeds specificity (0.744).
**Implications:** EBV activity is legibly, though not perfectly, written into the host transcriptome at single-cell resolution — turning a prior qualitative observation (SoRelle et al. 2021) into a quantitative one. The compact 15-gene signature suggests specific host programs rather than a global shift. The sensitivity>specificity asymmetry raises the possibility that some threshold-negative cells carry real low-level viral activity (a candidate "primed" intermediate state). This is association, not causation — sequencing depth is a candidate confounder.
**Tags:** virology, scrna-seq, host-response, ebv, classifier, single-cell, latency

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-01 | mycelium-lifecycle-2026-07-01 | SRR12682296 (GSE158275, LCL) | ViralScan | AUC 0.866 / MCC 0.570 (leakage-corrected, per-fold HVG) | supports |

### Open Questions
- Is the specificity gap (0.744 < 0.822) a real primed/intermediate cell state or noise? (todo: intermediate-attractor-test)
- How much of the AUC survives adjustment for sequencing depth as a confounder? (todo: depth-confounder-check)
- Does the 15-gene signature generalize across the other LCL lines / donors in GSE158275? (single-sample result)
- Is the ≥10-UMI positive label near the biological EC50 of the host response? (todo: Hill/threshold sweep)
