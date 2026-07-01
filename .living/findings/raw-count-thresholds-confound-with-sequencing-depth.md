---
topic: raw-count-thresholds-confound-with-sequencing-depth
description: Defining a positive/detected label by a raw (un-normalized) count threshold makes the label a proxy for sequencing depth, confounding any downstream predictor.
created: 2026-07-01
last_updated: 2026-07-01
status: active
---

# Raw-count thresholds confound with sequencing depth

## F-003: A raw-UMI positivity threshold makes the label depth-dependent, confounding predictors of it
**Status:** preliminary
**Claim:** Defining a binary label by a fixed raw-count threshold (here: EBV-positive = ≥10 raw viral UMI) makes the label strongly dependent on per-cell sequencing depth, because deeper cells reach the threshold more easily. In the EBV LCL data, EBV+ rate rose from 27.5% (shallowest host-depth quintile) to 96.3% (deepest), and sequencing depth alone predicted the label at AUC 0.803 (all cells) / 0.967 (in a top-50%-depth balanced design). Any feature that correlates with depth then appears predictive of the label. Standard mitigations that operate on class *counts* — class balancing, top-depth filtering — do NOT remove this: they do not equalize the between-class depth *distributions*, and depth-filtering can even widen the gap.
**Implications:** Positive/detected/expressed labels in count data (scRNA-seq, bulk RNA-seq, amplicon) should be defined depth-independently — as a normalized rate (CPM/fraction), via a depth-matched case-control design, or with depth entered into the label definition — before training or interpreting any predictor. When a raw-threshold label is unavoidable, always report depth-alone predictive performance as the confounding baseline, and E-values (or depth-adjusted effects) for individual features.
**Tags:** methodology, confounding, sequencing-depth, count-data, thresholding, label-definition, scrna-seq

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-01 | depth-confounder-check | SRR12682296 (GSE158275, LCL) | ViralScan | EBV+ rate 27%→96% across depth quintiles; depth-alone AUC 0.803/0.967; balancing+depth-filter did not fix it | supports |

### Open Questions
- Is a normalized-fraction (EBV CPM) label enough to decouple the label from depth, or does capture efficiency still bias it?
- General rule vs case-specific: how strong must the feature–depth correlation be before the confound dominates?
