---
topic: host-transcriptome-encodes-viral-infection-state
description: Whether a cell's own (host) gene expression carries enough signal to predict its intracellular viral infection status at single-cell resolution.
created: 2026-07-01
last_updated: 2026-07-01
status: active
---

# Host transcriptome encodes viral infection state

## F-001: A host-gene classifier predicts per-cell EBV status, but the headline is substantially depth-confounded
**Status:** contradicted (a major confounder overturns the naive interpretation; a weaker residual signal survives)
**Claim:** Within a single EBV-transformed LCL sample (SRR12682296, GSE158275), a cross-validated host-gene classifier separates EBV-positive from EBV-negative cells (EBV+ = ≥10 viral UMI) at AUC 0.866 / MCC 0.570 (n=1906; 1179+/727−; 15 stable genes). **However, a follow-up confounder analysis shows this is substantially driven by sequencing depth, not host biology:** the ≥10-raw-UMI label is strongly depth-dependent (EBV+ rate rises 27.5%→96.3% across host-depth quintiles); sequencing depth *alone* predicts EBV status at AUC 0.803 (all cells) and **0.967 within the classifier's own balanced+depth-filtered design — higher than the 0.866 host-gene model.** Because host features are depth-normalized, the classifier recovers a degraded proxy of the depth-driven label. Yet a residual genuine signal remains: **5 of the 15 stable genes retain a depth-adjusted association (E-value ≥2, p<1e-7)**, while the two top-stability genes are fully explained away by depth (p≈0.8 after adjustment).
**Implications:** The naive reading ("host transcriptome strongly predicts EBV status") is overstated — much of the apparent signal is a sequencing-depth artifact introduced by a raw-UMI positivity threshold. The pipeline's top-50%-depth + class-balancing step does NOT fix this (it matches class *counts*, not the between-class depth *distributions*; the filter actually widens the depth gap). A real but modest host-response signal exists in ~5 genes. Any downstream use (biomarker panel, mechanistic follow-up) should restrict to the depth-robust genes and re-define the label depth-independently.
**Tags:** virology, scrna-seq, host-response, ebv, classifier, single-cell, confounding, sequencing-depth

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-01 | mycelium-lifecycle-2026-07-01 | SRR12682296 (GSE158275, LCL) | ViralScan | AUC 0.866 / MCC 0.570 (leakage-corrected, per-fold HVG) | supports |
| 2026-07-01 | depth-confounder-check | SRR12682296 (GSE158275, LCL) | ViralScan | depth-alone AUC 0.803 (all)/0.967 (headline design); EBV+ rate 27%→96% by depth quintile; only 5/15 genes E≥2 | contradicts |

### Open Questions
- With a depth-normalized EBV label (CPM/fraction) or a depth-matched case-control design, how much host-response signal remains? (the principled re-analysis)
- Do the 5 depth-robust genes form a coherent biological program?
- Does any residual signal generalize across the other LCL lines / donors in GSE158275? (single-sample result)
