# Findings Registry

> Per-project index of all scientific findings. See topic files for full evidence ledgers.

| ID | Claim | Status | Topic | Implications | Tags | Last Updated |
|----|-------|--------|-------|--------------|------|--------------|
| F-001 | Host-gene→EBV-status classifier (AUC 0.866) is substantially sequencing-depth-confounded; ~5/15 genes carry a residual real signal | contradicted | [host-transcriptome-encodes-viral-infection-state](host-transcriptome-encodes-viral-infection-state.md) | Headline overstated by a depth-driven raw-UMI label; depth alone beats the model; restrict to depth-robust genes + re-define label | virology, scrna-seq, host-response, ebv, confounding | 2026-07-01 |
| F-002 | Correcting unsupervised HVG-selection leakage moved metrics < 1 seed-SD (negligible) | preliminary | [unsupervised-feature-selection-leakage-is-often-negligible](unsupervised-feature-selection-leakage-is-often-negligible.md) | Unsupervised-selection leaks rarely change conclusions; compare shifts to estimator variance before narrating | methodology, data-leakage, cross-validation, feature-selection | 2026-07-01 |
| F-003 | A raw-UMI positivity threshold makes the label depth-dependent, confounding predictors of it | preliminary | [raw-count-thresholds-confound-with-sequencing-depth](raw-count-thresholds-confound-with-sequencing-depth.md) | Define positive/detected labels depth-independently (CPM/fraction, depth-matched, or depth covariate); report depth-alone baseline | methodology, confounding, sequencing-depth, count-data, thresholding | 2026-07-01 |
