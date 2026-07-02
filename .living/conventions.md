# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

<!-- Document any project-specific convention overrides here. -->

## scRNA-seq associations: control depth AND %mito, and define labels depth-independently

**Rule**: For any per-cell association / classifier / DE in this repo (host-response and
similar), (1) define positive/detected labels **depth-independently** (a normalized rate
like CPM/fraction, a depth-matched design, or depth in the label definition) — never a
raw-count threshold alone; (2) adjust associations for **both log sequencing depth and
percent-mitochondrial content** (add cell-cycle where relevant); (3) always report the
**depth-alone predictive baseline** for any count-threshold label; (4) flag any
mitochondrial-encoded hit for an explicit %mito control before calling it biology.

**Why**: A raw-count label is a depth proxy (EBV+ rate ran 27%→96% across depth
quintiles; depth alone hit AUC 0.80–0.97). And %mito is a *separate* confounder: it
removed 217/518 depth-robust genes here and turned MT-ND4L from FDR 4e-7 to 0.07.
Controlling one but not the other still yields artifacts.

**Source**: learnings.md — "Raw-count positivity thresholds silently confound with
sequencing depth" (2026-07-01) and "Mitochondrial genes need a %mito control" (2026-07-02);
findings F-001, F-003. Two related confounder instances → promoted to a convention.
Complements the `robust-analysis` + `bioinformatics` packs.

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
