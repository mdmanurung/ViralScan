---
topic: unsupervised-feature-selection-leakage-is-often-negligible
description: When feature selection performed before a train/test split is unsupervised (does not use labels), the resulting leakage often changes held-out metrics only within noise.
created: 2026-07-01
last_updated: 2026-07-01
status: active
---

# Unsupervised feature-selection leakage is often negligible

## F-002: Correcting unsupervised HVG-selection leakage moved classifier metrics within seed noise
**Status:** preliminary
**Claim:** In a single-cell classifier that selected highly-variable genes on all cells before an 80/20 cross-validation split (a feature-selection leak), moving HVG selection inside the fold (training cells only) changed every held-out metric by less than one seed-to-seed standard deviation (e.g. AUC 0.845→0.866 with SD≈0.036; specificity 0.710→0.744 with SD≈0.053). Because HVG selection is unsupervised — variance-based and never using the class labels — the leak carried little label information and its practical effect was negligible.
**Implications:** "Remove the leak ⇒ the number goes down" is not a safe assumption. The magnitude and even the direction of a leakage correction depend on whether the leaked step used the outcome labels. Unsupervised selection/normalization leaks are worth fixing for methodological correctness but rarely change conclusions; supervised leaks (target-aware selection, label-informed thresholds) are the ones that inflate. Always compare any post-fix metric shift to the estimator's own variability before narrating a cause.
**Tags:** methodology, data-leakage, cross-validation, feature-selection, hvg, scrna-seq, reproducibility

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-01 | mycelium-lifecycle-2026-07-01 | SRR12682296 (GSE158275, LCL) | ViralScan | 6 metrics all shifted < 1 seed-SD after per-fold HVG | supports |

### Open Questions
- Does the same near-null hold for other unsupervised steps (per-fold PCA, scaling) in this pipeline?
- Would a supervised feature-selection step (e.g. picking genes by EBV correlation before splitting) show the expected large inflation? (a useful positive control)
