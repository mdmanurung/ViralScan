# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

<!-- Document any project-specific convention overrides here. -->

## Cross-validation: fit feature selection inside the split

**Rule**: In any cross-validated classifier/regressor in this repo (e.g. the
host-response models in `src/viralscan/scripts/hostresponse.py`), highly-variable-gene
selection, feature filtering, and any fitted normalization/scaling must be computed on
the **training fold only** and applied to the held-out fold — never fit on the full
dataset before splitting. Prefer an explicit per-fold path (see `_hvg_mask` +
`use_hvg=True` in `_run_l2_regression`) and, where practical, a test asserting the
selection is recomputed per fold.

**Why**: Selecting features on all cells leaks held-out expression into the feature
space and biases reported metrics. The effect can be small for *unsupervised* selection,
so also check whether any metric shift exceeds the seed-to-seed SD before interpreting it.

**Source**: learnings.md — "Feature selection (HVG) before the CV split leaks into
held-out metrics" (2026-07-01, review finding F1). Promoted on severity + a shipped
structural mitigation, ahead of the usual 3-instance threshold. Complements the
`robust-analysis` and `bioinformatics` convention packs (general leakage guidance).
