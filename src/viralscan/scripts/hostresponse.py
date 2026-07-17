"""
Associate viral presence with host gene expression via logistic regression.

Implements the approach from Luebbert et al. 2026 (Nature Biotechnology):
  - Align virus and host h5ad on shared cell barcodes.
  - For each detected virus: train multi-seed L2 logistic regression on
    normalized host gene expression to predict virus presence/absence. When
    use_hvg is set, highly-variable-gene selection is performed inside each
    cross-validation fold (on training cells only) to avoid feature-selection
    leakage into the held-out metrics.
  - Run randomized Lasso stability selection (Meinshausen & Bühlmann 2010)
    to identify genes whose association is reproducible across sub-samples
    and random penalty perturbations.
  - Optionally run gget.enrichr pathway enrichment on stable genes.

Because the raw ">=N viral UMI" positive-call label tracks sequencing depth
(deeper cells carry more viral AND more host counts), a naive host-gene AUC is
partly a depth artifact (findings F-001/F-003). Every run therefore also reports
a depth-confound baseline: the AUC obtainable from sequencing depth ALONE under
the identical split, and per-gene depth-adjusted E-values (Ding & VanderWeele
2016) on the stably selected genes. If the depth-alone AUC approaches the model
AUC, the headline number should be read as depth-confounded.

Outputs (per virus, under <output>/hostresponse/):
  - <virus>_gene_weights.csv   — mean/SD of L2 coefficients (per-fold-HVG: genes
    selected in >=1 fold, with n_folds_selected)
  - <virus>_stability.csv      — per-gene stability probability + merged weights
  - <virus>_depth_diagnostics.csv — per stable gene: depth- (and, by default,
    %mito-) adjusted odds ratio + E-value (confounder strength needed to explain
    the association away)
  - <virus>_differential.csv — (when --differential) genome-wide depth-(and %mito-)
    adjusted partial correlation of every gene with virus status: partial_r, p_value,
    fdr, direction
  - hostresponse_metrics.csv   — sensitivity/specificity/balanced-acc/AUC/MCC +
    depth_alone_auc + n_genes_evalue_ge2 (depth-confound baseline) summary
  - <virus>_enrichment_<db>.csv (when --enrichment is set)
"""

import argparse
import contextlib
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import sklearn as _sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

from viralscan.kb_outputs import KbCountOutputs
from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging

# sklearn 1.8 deprecated the `penalty` kwarg; use l1_ratio=1 + saga instead.
# On older sklearn, l1_ratio without penalty='elasticnet' is silently ignored,
# so we must keep the explicit penalty kwarg there.
_SKLEARN_VER = tuple(int(x) for x in _sklearn.__version__.split(".")[:2])
_L1_LR_KWARGS: dict = (
    {"solver": "saga", "l1_ratio": 1.0}
    if _SKLEARN_VER >= (1, 8)
    else {"penalty": "l1", "solver": "liblinear"}
)

log = setup_script_logging()

DEFAULT_SEEDS = [0, 1, 10, 42, 100, 1234]
MIN_VIRUS_CELLS = 10
TOP_DEPTH_FRAC = 0.5
_TRAIN_FRAC = 0.8


def _safe_name(name: str) -> str:
    """Make a filesystem-safe version of a virus accession."""
    return name.replace("/", "_").replace(" ", "_").replace(".", "_").replace(":", "_")


def _load_viral_accessions(analysis_txt: str) -> set:
    """Read the viral gene-ID list produced by the analysis rule."""
    accessions: set = set()
    with open(analysis_txt) as fh:
        for line in fh:
            v = line.strip()
            if v:
                accessions.add(v)
    return accessions


def _detect_and_normalize(host_adata) -> None:
    """Normalize and log-transform the host matrix if it appears to be raw counts.

    Stores raw per-cell total counts in obs["_raw_depth"] BEFORE any normalization
    so that depth-based filtering is always computed from library size, not from
    the normalized feature subset.

    Detection heuristic:
      - If "log1p" is in adata.uns → already processed, skip.
      - Otherwise sample the top-left 100×100 block; if all values are integers
        AND the max exceeds 50, treat as raw counts.
    This avoids double-normalising while still handling the common case where the
    user passes a raw-count h5ad.  If the heuristic is ambiguous the user should
    pre-process and set adata.uns["log1p"] = {} themselves.
    """
    X = host_adata.X
    # np.asarray(...).flatten() handles both scipy sparse (via np.matrix) and ndarray.
    raw_depth = np.asarray(X.sum(axis=1)).flatten()
    host_adata.obs["_raw_depth"] = raw_depth

    if "log1p" in host_adata.uns:
        log.info(
            "Host h5ad: already log1p-normalised (found 'log1p' in uns). Skipping normalization."
        )
        return

    # Sample a small block to decide whether data is raw
    r = min(100, host_adata.n_obs)
    c = min(100, host_adata.n_vars)
    sample = X[:r, :c]
    if sp.issparse(sample):
        sample = sample.toarray()
    sample = np.asarray(sample, dtype=float)
    is_integer = np.all(sample == np.floor(sample))
    if is_integer and sample.max() > 10:
        log.info("Host h5ad: looks like raw counts (integer values, max > 10). Normalising.")
        sc.pp.normalize_total(host_adata, target_sum=1e4)
        sc.pp.log1p(host_adata)
    else:
        log.info("Host h5ad: appears already normalized. Skipping normalization.")


def _select_features(host_adata, use_hvg: bool):
    """Return (X, feature_names) after HVG selection or full-gene densification.

    When use_hvg=False the full matrix is densified — this may require tens of GB
    of RAM for large datasets; a warning is emitted.
    """
    if use_hvg:
        sc.pp.highly_variable_genes(host_adata, flavor="seurat")
        mask = host_adata.var["highly_variable"]
        n_hvg = int(mask.sum())
        log.info("Selecting %d highly variable genes as features.", n_hvg)
        sub = host_adata[:, mask]
        X = sub.X
        feature_names = sub.var_names.tolist()
    else:
        n_all = host_adata.n_vars
        log.warning(
            "use_hvg=False: densifying full gene matrix (%d genes × %d cells). "
            "This may require substantial RAM.",
            n_all,
            host_adata.n_obs,
        )
        X = host_adata.X
        feature_names = host_adata.var_names.tolist()

    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32), feature_names


def _balanced_split(pos_idx, neg_idx, depth, X, seed: int, top_depth_frac: float = TOP_DEPTH_FRAC):
    """Depth-filtered balanced train/test split.

    Filters each class to the top ``top_depth_frac`` by RAW sequencing depth
    (stored in obs["_raw_depth"] before normalization) to reduce false-negative
    viral-absence labels in low-coverage cells. Pass ``top_depth_frac=1.0`` to
    keep all cells (used by the depth-matched design, where the cohort has
    already been equalized on depth and further top-depth filtering would
    re-introduce a between-class depth gap).

    Returns (X_train, y_train, X_test_pos, X_test_neg) or None if too few cells.
    """
    rng = np.random.default_rng(seed)

    def top_depth(idx):
        if top_depth_frac >= 1.0:
            return idx
        d = depth[idx]
        cutoff = np.percentile(d, (1 - top_depth_frac) * 100)
        return idx[d >= cutoff]

    pos_top = top_depth(np.asarray(pos_idx))
    neg_top = top_depth(np.asarray(neg_idx))
    n = min(len(pos_top), len(neg_top))
    if n < 5:
        return None

    pos_s = rng.choice(pos_top, size=n, replace=False)
    neg_s = rng.choice(neg_top, size=n, replace=False)

    split = max(1, int(_TRAIN_FRAC * n))
    pos_train, pos_test = pos_s[:split], pos_s[split:]
    neg_train, neg_test = neg_s[:split], neg_s[split:]

    if len(pos_test) == 0 or len(neg_test) == 0:
        return None

    X_train = np.vstack([X[pos_train], X[neg_train]])
    y_train = np.array([1] * len(pos_train) + [0] * len(neg_train), dtype=np.int8)
    return X_train, y_train, X[pos_test], X[neg_test]


def _hvg_mask(x_train: np.ndarray) -> np.ndarray:
    """Seurat highly-variable-gene mask computed on TRAINING cells only.

    Selecting features from the full dataset before the train/test split leaks
    test-cell expression into the feature space (the held-out AUC/MCC become
    upward-biased). Computing the HVG mask inside the fold, on the training
    cells alone, removes that leak. Returns a boolean mask over x_train's columns.
    """
    import anndata as ad  # noqa: PLC0415

    tmp = ad.AnnData(np.asarray(x_train, dtype=np.float32))
    try:
        sc.pp.highly_variable_genes(tmp, flavor="seurat")
    except ValueError:
        # Degenerate training fold (e.g. an empty dispersion bin) — fall back to
        # using every gene for this fold rather than dropping the fold entirely.
        return np.ones(x_train.shape[1], dtype=bool)
    mask = tmp.var["highly_variable"].to_numpy()
    if not mask.any():
        return np.ones(x_train.shape[1], dtype=bool)
    return np.asarray(mask)


def _run_l2_regression(
    X, virus_presence, depth, seeds, feature_names, use_hvg=False, top_depth_frac=TOP_DEPTH_FRAC
):
    """Multi-seed L2 logistic regression on balanced, depth-filtered data.

    When ``use_hvg=True``, ``X`` is the full normalized gene matrix and highly
    variable genes are selected *inside each fold* from the training cells only
    (leakage-free). When ``use_hvg=False`` (default), ``X`` is used as-is and the
    behavior matches the original all-features path.

    Returns (weights_df, metrics_dict) or (None, None) when there are too few
    positive cells or every balanced split fails.
    """
    pos_idx = np.where(virus_presence)[0]
    neg_idx = np.where(~virus_presence)[0]

    if len(pos_idx) < MIN_VIRUS_CELLS:
        return None, None

    metric_lists: dict = {
        "sensitivity": [],
        "specificity": [],
        "balanced_acc": [],
        "auc": [],
        "mcc": [],
    }
    n_genes = X.shape[1]
    # Per-fold HVG sets differ, so weights are accumulated per global gene index.
    coef_sum = np.zeros(n_genes)
    coef_sq = np.zeros(n_genes)
    coef_cnt = np.zeros(n_genes, dtype=int)
    simple_weights: list = []  # used only on the use_hvg=False path (fixed feature set)
    n_fit = 0

    for seed in seeds:
        split = _balanced_split(pos_idx, neg_idx, depth, X, seed, top_depth_frac=top_depth_frac)
        if split is None:
            continue
        X_train, y_train, X_test_pos, X_test_neg = split

        # Leakage-free feature selection: fit the HVG mask on training cells only.
        mask = _hvg_mask(X_train) if use_hvg else None
        if mask is not None:
            X_train, X_test_pos, X_test_neg = (
                X_train[:, mask],
                X_test_pos[:, mask],
                X_test_neg[:, mask],
            )

        model = LogisticRegression(solver="lbfgs", max_iter=1000, C=1.0, random_state=seed)
        model.fit(X_train, y_train)

        sensitivity = float((model.predict_proba(X_test_pos)[:, 1] >= 0.5).mean())
        specificity = float((model.predict_proba(X_test_neg)[:, 1] < 0.5).mean())

        X_test = np.vstack([X_test_pos, X_test_neg])
        y_test = np.array([1] * len(X_test_pos) + [0] * len(X_test_neg))
        y_prob = model.predict_proba(X_test)[:, 1]
        y_pred = (y_prob >= 0.5).astype(int)

        metric_lists["sensitivity"].append(sensitivity)
        metric_lists["specificity"].append(specificity)
        metric_lists["balanced_acc"].append(float(balanced_accuracy_score(y_test, y_pred)))
        # Matthews correlation coefficient: single-number summary of the 2x2
        # confusion matrix, robust to the class balance produced by _balanced_split.
        metric_lists["mcc"].append(float(matthews_corrcoef(y_test, y_pred)))
        with contextlib.suppress(ValueError):
            metric_lists["auc"].append(float(roc_auc_score(y_test, y_prob)))

        coef = model.coef_[0]
        if mask is not None:
            gi = np.where(mask)[0]
            coef_sum[gi] += coef
            coef_sq[gi] += coef**2
            coef_cnt[gi] += 1
        else:
            simple_weights.append(coef)
        n_fit += 1

    if n_fit == 0:
        return None, None

    if use_hvg:
        safe = np.maximum(coef_cnt, 1)
        wmean = coef_sum / safe
        wsd = np.sqrt(np.maximum(coef_sq / safe - wmean**2, 0.0))
        selected = coef_cnt > 0
        weights_df = pd.DataFrame(
            {
                "gene": feature_names,
                "weight_mean": np.where(selected, wmean, 0.0),
                "weight_sd": np.where(selected, wsd, 0.0),
                "n_folds_selected": coef_cnt,
            }
        )
    else:
        weights_arr = np.stack(simple_weights)
        weights_df = pd.DataFrame(
            {
                "gene": feature_names,
                "weight_mean": weights_arr.mean(axis=0),
                "weight_sd": weights_arr.std(axis=0),
            }
        )
    metrics = {
        k: {"mean": float(np.mean(v)), "sd": float(np.std(v))} for k, v in metric_lists.items() if v
    }
    return weights_df, metrics


def _run_stability_selection(
    X, virus_presence, n_iter: int, seed: int = 42, alpha_min: float = 0.2
):
    """Randomized Lasso stability selection (Meinshausen & Bühlmann 2010).

    In each iteration:
      1. Draw a random half-subsample (balanced pos/neg).
      2. Standardize features (required for fair L1 penalty across genes).
      3. Apply random per-feature penalty scaling in [alpha_min, 1] — this
         perturbs the effective regularization per gene, improving selection
         stability over simple repeated sub-sampling.
      4. Fit L1 logistic regression; record which genes have non-zero coefficient.

    Returns a numpy array of per-gene selection probabilities in [0, 1].
    """
    pos_idx = np.where(virus_presence)[0]
    neg_idx = np.where(~virus_presence)[0]
    n_features = X.shape[1]
    selection_counts = np.zeros(n_features)
    rng = np.random.default_rng(seed)
    scaler = StandardScaler()
    valid_iters = 0

    for _ in range(n_iter):
        n_half = min(len(pos_idx), len(neg_idx)) // 2
        if n_half < 3:
            break
        pos_sub = rng.choice(pos_idx, size=n_half, replace=False)
        neg_sub = rng.choice(neg_idx, size=n_half, replace=False)
        idx = np.concatenate([pos_sub, neg_sub])
        y = np.array([1] * n_half + [0] * n_half)

        X_sub = scaler.fit_transform(X[idx])

        # Per-feature random penalty scaling: effective λ_j = λ / c_j
        c_per_feature = rng.uniform(alpha_min, 1.0, size=n_features)
        X_rand = X_sub * c_per_feature[np.newaxis, :]

        try:
            model = LogisticRegression(
                **_L1_LR_KWARGS,
                C=1.0,
                max_iter=300,
                random_state=int(rng.integers(0, 10000)),
            )
            model.fit(X_rand, y)
            selection_counts += (model.coef_[0] != 0).astype(float)
            valid_iters += 1
        except Exception:
            pass

    if valid_iters == 0:
        return np.zeros(n_features)
    return selection_counts / valid_iters


LABEL_CHOICES = ("raw", "cpm", "fraction")


def _virus_presence_label(counts, depth, detection_threshold: int, label: str) -> np.ndarray:
    """Compute the per-cell virus-positive label under the chosen definition.

    - ``raw``      : ``counts >= detection_threshold`` (raw UMI). The default;
      correlates with sequencing depth, so its downstream AUC is depth-confounded.
    - ``cpm`` /
      ``fraction`` : depth-NORMALIZED. Positive = the cells with the highest viral
      burden per host UMI, keeping prevalence identical to the raw label (same
      number of positives) so the metrics stay comparable. Viral-per-host-UMI
      (cpm) and viral-fraction rank cells identically, so the two names select the
      same cells; both decouple the label from depth (findings F-001/F-003).

    ``depth`` must be the host-only library size, so the CPM denominator does not
    include the viral counts in the numerator.
    """
    counts = np.asarray(counts, dtype=float)
    if label == "raw":
        return np.asarray(counts >= detection_threshold)
    if label not in ("cpm", "fraction"):
        raise ValueError(f"label must be one of {LABEL_CHOICES}, got {label!r}")
    n_pos = int((counts >= detection_threshold).sum())
    if n_pos == 0:
        return np.zeros(len(counts), dtype=bool)
    ratio = counts / np.maximum(np.asarray(depth, dtype=float), 1.0)
    # Positive = the n_pos cells with the highest depth-normalized viral burden.
    cutoff = np.sort(ratio)[::-1][n_pos - 1]
    return np.asarray(ratio >= cutoff)


def _depth_match_indices(virus_presence, depth, n_bins: int = 20, seed: int = 42) -> np.ndarray:
    """Coarsened exact matching on host-depth quantile bins.

    Keeps equal numbers of positive and negative cells within each depth bin, so
    the two classes end up with (near-)identical depth distributions and depth
    alone can no longer discriminate. The host transcriptome is then tested on a
    cohort where any surviving signal is depth-independent by construction
    (the principled fix for the confound; see depth_matched_reanalysis.py).
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(virus_presence).astype(int)
    depth = np.asarray(depth, dtype=float)
    edges = np.quantile(depth, np.linspace(0, 1, n_bins + 1))
    edges[-1] += 1.0
    binid = np.digitize(depth, edges[1:-1])
    keep: list = []
    for b in np.unique(binid):
        idx = np.where(binid == b)[0]
        pos, neg = idx[y[idx] == 1], idx[y[idx] == 0]
        k = min(len(pos), len(neg))
        if k == 0:
            continue
        keep.extend(rng.choice(pos, k, replace=False).tolist())
        keep.extend(rng.choice(neg, k, replace=False).tolist())
    return np.array(sorted(keep), dtype=int)


def _e_value(or_: float) -> float:
    """Ding & VanderWeele (2016) E-value from an odds ratio (approx risk ratio).

    The E-value is the minimum strength of association (on the risk-ratio scale)
    that an unmeasured confounder would need with BOTH the exposure and the
    outcome to fully explain away the observed odds ratio. E ~ 1 means the
    association is fragile (a weak confounder could explain it); E >= 2 means a
    confounder would have to be at least a 2-fold risk factor on both arms.
    """
    if not np.isfinite(or_) or or_ <= 0:
        return float("nan")
    rr = or_ if or_ >= 1 else 1.0 / or_
    return float(rr + np.sqrt(rr * (rr - 1.0)))


def _depth_alone_auc(pos_idx, neg_idx, depth, seeds, top_depth_frac: float = TOP_DEPTH_FRAC):
    """AUC using log(sequencing depth) as the ONLY predictor.

    Uses the same balanced, top-depth-filtered split as the host-gene model
    (:func:`_balanced_split`) so the number is directly comparable to the
    headline model AUC. This is the confound baseline: if depth alone predicts
    virus status about as well as the host transcriptome, the headline signal is
    largely a sequencing-depth artifact rather than biology (findings F-001/F-003).
    Under the depth-matched design (``top_depth_frac=1.0``) this should fall to
    ~0.5, confirming the match removed the depth signal.

    Returns {"mean", "sd"} over the seeds, or None if no split was valid.
    """
    logd = np.log1p(depth).reshape(-1, 1).astype(np.float32)
    aucs: list = []
    for seed in seeds:
        split = _balanced_split(pos_idx, neg_idx, depth, logd, seed, top_depth_frac=top_depth_frac)
        if split is None:
            continue
        x_train, y_train, x_test_pos, x_test_neg = split
        scaler = StandardScaler().fit(x_train)
        model = LogisticRegression(max_iter=1000, random_state=seed)
        model.fit(scaler.transform(x_train), y_train)
        x_test = np.vstack([x_test_pos, x_test_neg])
        y_test = np.array([1] * len(x_test_pos) + [0] * len(x_test_neg))
        prob = model.predict_proba(scaler.transform(x_test))[:, 1]
        with contextlib.suppress(ValueError):
            aucs.append(float(roc_auc_score(y_test, prob)))
    if not aucs:
        return None
    return {"mean": float(np.mean(aucs)), "sd": float(np.std(aucs))}


def _panel_depth_adjusted_auc(
    pos_idx, neg_idx, X_panel, depth, seeds, top_depth_frac: float = TOP_DEPTH_FRAC
):
    """AUC of the stable-gene panel logistic with log(depth) added as a covariate.

    Mirrors :func:`_depth_alone_auc` but augments the stable-gene feature matrix
    with ``log1p(host_depth)`` as an explicit last column, under the same balanced,
    depth-filtered split.  Comparing this to :func:`_depth_alone_auc` answers:
    once the depth channel is absorbed into the model, do the stable genes still add
    predictive value?  A value close to the depth-alone AUC means the panel is
    largely a depth proxy; a value substantially higher means genuine transcriptional
    biology survives depth adjustment.

    ``depth`` must be the host-only library size (``obs["_raw_depth"]``); using
    host+viral depth would reintroduce the confound being controlled.
    Returns ``{"mean", "sd"}`` over the seeds, or ``None`` if no split was valid.
    """
    logd = np.log1p(depth).reshape(-1, 1).astype(np.float32)
    X_aug = np.hstack([np.asarray(X_panel, dtype=np.float32), logd])
    aucs: list = []
    for seed in seeds:
        split = _balanced_split(pos_idx, neg_idx, depth, X_aug, seed, top_depth_frac=top_depth_frac)
        if split is None:
            continue
        x_train, y_train, x_test_pos, x_test_neg = split
        scaler = StandardScaler().fit(x_train)
        model = LogisticRegression(max_iter=1000, random_state=seed)
        model.fit(scaler.transform(x_train), y_train)
        x_test = np.vstack([x_test_pos, x_test_neg])
        y_test = np.array([1] * len(x_test_pos) + [0] * len(x_test_neg))
        prob = model.predict_proba(scaler.transform(x_test))[:, 1]
        with contextlib.suppress(ValueError):
            aucs.append(float(roc_auc_score(y_test, prob)))
    if not aucs:
        return None
    return {"mean": float(np.mean(aucs)), "sd": float(np.std(aucs))}


# 13 protein-coding mitochondrial genes (Ensembl gene IDs, version-stripped). Used
# to compute per-cell %mito for the QC covariate. Host h5ads that use gene *symbols*
# are handled separately by the "MT-"/"mt-" prefix rule in _mt_gene_mask.
_MT_ENSEMBL = frozenset(
    {
        "ENSG00000198899",  # MT-ATP6
        "ENSG00000198804",  # MT-CO1
        "ENSG00000198712",  # MT-CO2
        "ENSG00000198938",  # MT-CO3
        "ENSG00000198763",  # MT-ND2
        "ENSG00000198840",  # MT-ND3
        "ENSG00000198886",  # MT-ND4
        "ENSG00000212907",  # MT-ND4L
        "ENSG00000198786",  # MT-ND5
        "ENSG00000198695",  # MT-ND6
        "ENSG00000198727",  # MT-CYB
        "ENSG00000198888",  # MT-ND1
        "ENSG00000228253",  # MT-ATP8
    }
)


def _mt_gene_mask(var_names) -> np.ndarray:
    """Boolean mask of mitochondrial genes among ``var_names``.

    Matches either the known mitochondrial Ensembl gene IDs (version suffix
    stripped) or the conventional ``MT-`` / ``mt-`` gene-symbol prefix, so the
    same code works whether the host h5ad is annotated with Ensembl IDs or symbols.
    """
    out = np.zeros(len(var_names), dtype=bool)
    for i, v in enumerate(var_names):
        s = str(v)
        if s.split(".")[0] in _MT_ENSEMBL or s.upper().startswith("MT-"):
            out[i] = True
    return out


def _per_gene_evalues(
    x_genes: np.ndarray,
    gene_names: list,
    y: np.ndarray,
    depth,
    pct_mito=None,
    mt_total_counts=None,
    raw_total=None,
    mt_self_counts=None,
) -> pd.DataFrame:
    """Depth-adjusted odds ratio + E-value per gene, on all aligned cells.

    For each gene, fits ``y ~ std(gene) + std(log depth) [+ std(%mito)]``
    (lightly-penalized logistic) and reports ``adj_OR = exp(gene coef)`` and its
    E-value. Adjusting for ``log(depth)`` removes the shared-depth path that
    inflates the raw association; the E-value then quantifies how robust each
    gene's adjusted association is to any *remaining* unmeasured confounder.
    ``depth`` must be the host-only library size (obs["_raw_depth"]); using
    host+viral depth would reintroduce the very confound this adjustment removes.

    When ``pct_mito`` (per-cell mitochondrial fraction, %) is given, it is added
    as a further covariate so a mitochondrial-QC artifact (high-%mito stressed
    cells) cannot masquerade as a host-response gene. For a gene that is *itself*
    mitochondrial, including it in %mito would trivially suppress it (circularity,
    e.g. MT-ND4L); ``mt_self_counts[gene]`` (that gene's raw counts) plus
    ``mt_total_counts`` and ``raw_total`` are used to leave that gene out of its
    own %mito covariate. A zero-variance %mito covariate is skipped.
    """
    logd = StandardScaler().fit_transform(np.log1p(np.asarray(depth)).reshape(-1, 1))[:, 0]
    xs = StandardScaler().fit_transform(np.asarray(x_genes, dtype=float))
    pm = np.asarray(pct_mito, dtype=float) if pct_mito is not None else None
    mt_self_counts = mt_self_counts or {}
    rows: list = []
    for j, g in enumerate(gene_names):
        cols = [xs[:, j], logd]
        if pm is not None:
            pm_g = pm
            # Leave-one-out: drop this gene from its own %mito if it is mitochondrial.
            if g in mt_self_counts and mt_total_counts is not None and raw_total is not None:
                denom = np.maximum(np.asarray(raw_total, dtype=float), 1.0)
                pm_g = (np.asarray(mt_total_counts) - np.asarray(mt_self_counts[g])) / denom * 100.0
            if np.std(pm_g) > 0:
                cols.append(StandardScaler().fit_transform(pm_g.reshape(-1, 1))[:, 0])
        feat = np.column_stack(cols)
        try:
            # Default L2 (C=1.0). This mildly shrinks the odds ratio toward 1,
            # which makes the resulting E-value *conservative* (a robust gene may
            # look slightly less robust, but a depth-confounded gene never looks
            # spuriously robust) and — crucially — keeps the fit stable when a
            # gene is nearly collinear with depth, where an unpenalized fit would
            # diverge to a meaningless huge odds ratio via quasi-separation.
            model = LogisticRegression(max_iter=1000, C=1.0)
            model.fit(feat, y)
            or_g = float(np.exp(model.coef_[0][0]))
            rows.append((g, or_g, _e_value(or_g)))
        except Exception:  # noqa: BLE001 — a single separated gene is recorded as NaN, not dropped
            rows.append((g, float("nan"), float("nan")))
    df = pd.DataFrame(rows, columns=["gene", "adj_OR", "E_value"])
    df["evalue_flag"] = df["E_value"].apply(
        lambda e: "robust" if (np.isfinite(e) and e >= 3.0)
        else "moderate" if (np.isfinite(e) and e >= 1.5)
        else "fragile"
    )
    return df


def _bh_fdr(pvals) -> np.ndarray:
    """Benjamini–Hochberg FDR-adjusted p-values (no statsmodels dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    out = np.empty(n)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def _residualize(Y: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Residualize columns of ``Y`` on covariates ``C`` (intercept added)."""
    Cc = np.column_stack([np.ones(len(C)), C])
    beta, *_ = np.linalg.lstsq(Cc, Y, rcond=None)
    return np.asarray(Y - Cc @ beta)


def _genome_wide_differential(
    x_genes: np.ndarray, gene_names: list, y: np.ndarray, depth, pct_mito=None
) -> pd.DataFrame:
    """Covariate-adjusted partial correlation of every gene with virus status.

    Residualizes both the gene matrix and the label on ``log(depth)`` (and, when
    given, ``%mito``), then reports the per-gene partial correlation, its t-test
    p-value, and a BH-FDR. Unlike the stable-gene E-values this scans *all* genes
    in ``x_genes`` (genome-wide when use_hvg is off), so it recovers a
    depth-adjusted differential-expression table rather than only ranking the
    pre-selected panel. Folded from go_enrichment.py:partial_assoc.
    """
    from scipy import stats  # noqa: PLC0415

    cov = np.log1p(np.asarray(depth, dtype=float)).reshape(-1, 1)
    if pct_mito is not None and np.std(np.asarray(pct_mito, dtype=float)) > 0:
        cov = np.column_stack([cov, np.asarray(pct_mito, dtype=float)])
    xr = _residualize(np.asarray(x_genes, dtype=float), cov)
    yr = _residualize(np.asarray(y, dtype=float).reshape(-1, 1), cov).ravel()
    xr = xr - xr.mean(0)
    xr = xr / (xr.std(0) + 1e-9)
    yr = (yr - yr.mean()) / (yr.std() + 1e-9)
    r = (xr * yr[:, None]).mean(0)
    n = len(y)
    t = r * np.sqrt((n - 2) / np.maximum(1 - r**2, 1e-12))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    fdr = _bh_fdr(p)
    return (
        pd.DataFrame(
            {
                "gene": list(gene_names),
                "partial_r": r,
                "p_value": p,
                "fdr": fdr,
                "direction": np.where(r >= 0, "up", "down"),
            }
        )
        .sort_values("fdr", kind="stable")
        .reset_index(drop=True)
    )


def _run_enrichment(
    gene_names: list,
    background_names: list,
    virus_name: str,
    out_dir: str,
    database: str,
) -> None:
    """Run gget.enrichr pathway enrichment on a gene list.

    Requires gene symbols (not Ensembl IDs). If the host h5ad uses Ensembl IDs,
    enrichment will silently return empty results from the API.

    Install gget with: pip install 'ViralScan[enrichment]'
    """
    try:
        import gget  # noqa: PLC0415
    except ImportError:
        log.warning("gget is not installed. Install it with: pip install 'ViralScan[enrichment]'")
        return

    log.info(
        "[%s] Enrichment: %d genes vs %d background, db=%s.",
        virus_name,
        len(gene_names),
        len(background_names),
        database,
    )
    try:
        results = gget.enrichr(
            gene_names,
            database=database,
            background_list=background_names,
            plot=False,
            quiet=True,
        )
        if results is not None and len(results) > 0:
            out_path = Path(out_dir) / f"{_safe_name(virus_name)}_enrichment_{database}.csv"
            results.to_csv(out_path, index=False)
            log.info("[%s] Enrichment results written to %s", virus_name, out_path)
        else:
            log.info("[%s] Enrichment returned no results.", virus_name)
    except Exception as exc:
        log.warning("[%s] Enrichment failed: %s", virus_name, exc)


def _looks_like_ensembl(gene_names) -> bool:
    """True if the gene names look like Ensembl gene IDs (ENSG...)."""
    return any(str(g).upper().startswith("ENSG") for g in list(gene_names)[:50])


def _map_ensembl_to_symbols(ensembl_ids: list) -> dict:
    """Map Ensembl gene IDs → HGNC symbols via mygene.info (best-effort, network).

    Returns ``{version_stripped_ensembl_id: symbol}``. On any network/parse error
    returns ``{}`` so the caller silently falls back to Ensembl IDs. Folded from
    analysis/hostresponse_ebv_matched/scripts/go_enrichment.py:map_symbols.
    """
    import json  # noqa: PLC0415
    import urllib.parse  # noqa: PLC0415
    import urllib.request  # noqa: PLC0415

    ids = sorted({str(e).split(".")[0] for e in ensembl_ids})
    if not ids:
        return {}
    try:
        data = urllib.parse.urlencode(
            {"q": ",".join(ids), "scopes": "ensembl.gene", "fields": "symbol", "species": "human"}
        ).encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=data)
        res = json.load(urllib.request.urlopen(req, timeout=60))  # noqa: S310
    except Exception as exc:  # noqa: BLE001 — network is optional; fall back to IDs
        log.warning("Gene-symbol mapping failed (%s); keeping Ensembl IDs.", exc)
        return {}
    out: dict = {}
    for r in res:
        if isinstance(r, dict) and "symbol" in r and "query" in r:
            out[str(r["query"])] = r["symbol"]
    return out


def _add_symbol_column(df, symbol_map: dict):
    """Insert a ``symbol`` column next to ``gene`` from ``symbol_map`` (no-op if empty)."""
    if df is None or not symbol_map or "gene" not in df.columns or "symbol" in df.columns:
        return df
    pos = df.columns.get_loc("gene") + 1
    df.insert(pos, "symbol", df["gene"].map(lambda g: symbol_map.get(str(g).split(".")[0], "")))
    return df


def run_hostresponse(
    virus_h5ad: str,
    host_h5ad: str,
    viral_accessions_file: str,
    out_dir: str,
    use_hvg: bool = True,
    seeds: list | None = None,
    n_stab_iter: int = 100,
    stab_min_prob: float = 0.6,
    top_n_genes: int = 50,
    detection_threshold: int = 10,  # >=10 EBV UMI is the validated positive-call threshold
    do_enrichment: bool = False,
    enrichment_db: str = "GO_Biological_Process_2023",
    label: str = "raw",
    depth_match: bool = False,
    control_mito: bool = True,
    annotate_symbols: bool = False,
    differential: bool = False,
) -> None:
    """Main entry point: run per-virus logistic regression host-response analysis.

    ``label`` selects the positive-call definition: "raw" (default, depth-confounded
    ``counts >= detection_threshold``) or the depth-normalized "cpm"/"fraction"
    (prevalence-matched top viral-per-host-UMI). ``depth_match=True`` additionally
    restricts each virus's analysis to a coarsened-exact depth-matched cohort so any
    surviving host-gene signal is depth-independent by construction. Both are
    opt-in de-confounding controls (findings F-001/F-003); the defaults reproduce
    prior behaviour. ``control_mito`` (on by default) adds per-cell %mitochondrial
    content as a covariate to the per-gene E-values, so a mitochondrial-QC artifact
    cannot masquerade as a host-response gene; it is a no-op when the host h5ad has
    no mitochondrial genes.
    """
    import anndata as ad  # noqa: PLC0415

    if label not in LABEL_CHOICES:
        raise ValueError(f"label must be one of {LABEL_CHOICES}, got {label!r}")
    if seeds is None:
        seeds = DEFAULT_SEEDS

    os.makedirs(out_dir, exist_ok=True)
    log.info("hostresponse output directory: %s", out_dir)

    log.info("Loading virus h5ad: %s", virus_h5ad)
    virus_adata = ad.read_h5ad(virus_h5ad)
    log.info("Loading host h5ad: %s", host_h5ad)
    host_adata_full = ad.read_h5ad(host_h5ad)

    # Filter virus matrix to confirmed viral gene IDs (analysis.txt).
    # If the pipeline used a combined host+viral reference, the h5ad contains
    # host genes too; restricting here prevents training models that predict
    # host gene expression from other host gene expression.
    viral_accessions = _load_viral_accessions(viral_accessions_file)
    viral_vars = [v for v in virus_adata.var_names if v in viral_accessions]
    if not viral_vars:
        log.error(
            "No viral accessions from analysis.txt (%s) match var_names in the virus h5ad "
            "(%s). Verify that the h5ad was produced from the same viralscan run.",
            viral_accessions_file,
            virus_h5ad,
        )
        return
    log.info("Retained %d viral gene(s) from analysis.txt for hostresponse.", len(viral_vars))
    virus_adata = virus_adata[:, viral_vars].copy()

    # Align cells: intersect barcodes between virus and host matrices.
    shared_barcodes = sorted(set(virus_adata.obs_names) & set(host_adata_full.obs_names))
    if len(shared_barcodes) < MIN_VIRUS_CELLS * 2:
        log.error(
            "Only %d shared barcodes between virus h5ad and host h5ad (need >= %d). "
            "Ensure both h5ad files use identical barcode strings.",
            len(shared_barcodes),
            MIN_VIRUS_CELLS * 2,
        )
        return
    log.info("Aligned on %d shared barcodes.", len(shared_barcodes))
    virus_adata = virus_adata[shared_barcodes].copy()
    host_adata = host_adata_full[shared_barcodes].copy()

    # Compute %mito from RAW counts BEFORE normalization (so the QC covariate for
    # the E-values is a real fraction, not a normalized artifact). Keep the raw
    # per-gene MT counts so a mitochondrial gene can be left out of its own %mito
    # covariate (avoids the MT-ND4L self-suppression circularity; see go_enrichment.py).
    mt_total_counts = raw_total = pct_mito = None
    mt_self_counts: dict = {}
    if control_mito:
        mt_mask = _mt_gene_mask(host_adata.var_names)
        raw_total = np.asarray(host_adata.X.sum(axis=1)).ravel().astype(float)
        if mt_mask.any():
            mt_sub = host_adata[:, mt_mask].X
            mt_sub = mt_sub.toarray() if sp.issparse(mt_sub) else np.asarray(mt_sub)
            mt_total_counts = mt_sub.sum(axis=1).astype(float)
            pct_mito = mt_total_counts / np.maximum(raw_total, 1.0) * 100.0
            mt_self_counts = {
                str(g): mt_sub[:, k] for k, g in enumerate(host_adata.var_names[mt_mask].tolist())
            }
            log.info(
                "%%mito control: %d mitochondrial genes found (median %.1f%%).",
                int(mt_mask.sum()),
                float(np.median(pct_mito)),
            )
        else:
            log.info("%%mito control requested but no mitochondrial genes found; skipping.")

    # Normalise host data; stores raw depth in obs["_raw_depth"] first.
    _detect_and_normalize(host_adata)
    depth = host_adata.obs["_raw_depth"].values

    # Full normalized gene matrix — passed to the classifier so that HVG selection
    # happens INSIDE each cross-validation fold (leakage-free; see _run_l2_regression).
    X_full, all_gene_names = _select_features(host_adata, use_hvg=False)
    # HVG subset for stability selection / enrichment background (descriptive gene
    # ranking, not a held-out metric). When use_hvg is False these coincide.
    if use_hvg:
        X_stab, stab_feature_names = _select_features(host_adata, use_hvg=True)
    else:
        X_stab, stab_feature_names = X_full, all_gene_names
    background_names = stab_feature_names  # used as enrichment background

    # Optional Ensembl→symbol map (network, best-effort). Built once over the HVG
    # feature set so every output CSV can carry a human-readable `symbol` column.
    symbol_map: dict = {}
    if annotate_symbols:
        if _looks_like_ensembl(stab_feature_names):
            symbol_map = _map_ensembl_to_symbols(stab_feature_names)
            log.info("Annotated %d gene symbols from Ensembl IDs.", len(symbol_map))
        else:
            log.info("Gene symbols requested but genes are not Ensembl IDs; skipping.")

    all_metrics = []

    for virus in viral_vars:
        counts = virus_adata[:, virus].X
        counts = counts.toarray().flatten() if sp.issparse(counts) else np.asarray(counts).flatten()

        virus_presence_full = _virus_presence_label(counts, depth, detection_threshold, label)
        log.info(
            "[%s] %d / %d cells positive (label=%s, detection_threshold=%d).",
            virus,
            int(virus_presence_full.sum()),
            len(virus_presence_full),
            label,
            detection_threshold,
        )

        # Optionally restrict to a depth-matched cohort so any surviving signal is
        # depth-independent by construction; then the in-split top-depth filter is
        # disabled (top_depth_frac=1.0) because matching already equalized depth.
        def _sub(a, idx):
            return None if a is None else np.asarray(a)[idx]

        if depth_match:
            midx = _depth_match_indices(virus_presence_full, depth)
            X_full_v, X_stab_v = X_full[midx], X_stab[midx]
            vp_v, depth_v = virus_presence_full[midx], depth[midx]
            pct_mito_v, mt_total_v, raw_total_v = (
                _sub(pct_mito, midx),
                _sub(mt_total_counts, midx),
                _sub(raw_total, midx),
            )
            mt_self_v = {g: c[midx] for g, c in mt_self_counts.items()}
            top_depth_frac = 1.0
            log.info("[%s] Depth-matched cohort: %d cells (from %d).", virus, len(midx), len(depth))
        else:
            X_full_v, X_stab_v = X_full, X_stab
            vp_v, depth_v = virus_presence_full, depth
            pct_mito_v, mt_total_v, raw_total_v = pct_mito, mt_total_counts, raw_total
            mt_self_v = mt_self_counts
            top_depth_frac = TOP_DEPTH_FRAC

        n_pos = int(vp_v.sum())
        if n_pos < MIN_VIRUS_CELLS:
            log.info(
                "[%s] Skipping: fewer than %d positive cells in the analysis cohort.",
                virus,
                MIN_VIRUS_CELLS,
            )
            continue

        # ── L2 multi-seed regression (per-fold HVG when use_hvg) ───────────
        weights_df, metrics = _run_l2_regression(
            X_full_v,
            vp_v,
            depth_v,
            seeds,
            all_gene_names,
            use_hvg=use_hvg,
            top_depth_frac=top_depth_frac,
        )
        if weights_df is None:
            log.info("[%s] Skipping: balanced split produced no valid models.", virus)
            continue

        weights_df["virus"] = virus
        # With per-fold HVG, most genes are never selected (weight 0); drop them so
        # the CSV lists only genes that entered at least one fold's model.
        if "n_folds_selected" in weights_df.columns:
            weights_df = weights_df[weights_df["n_folds_selected"] > 0].reset_index(drop=True)
        weights_csv = Path(out_dir) / f"{_safe_name(virus)}_gene_weights.csv"
        _add_symbol_column(weights_df, symbol_map).to_csv(weights_csv, index=False)
        log.info("[%s] Gene weights written to %s", virus, weights_csv)

        # ── Randomized Lasso stability selection (descriptive gene ranking) ──
        stab_probs = _run_stability_selection(
            X_stab_v, vp_v, n_stab_iter, seed=seeds[0] if seeds else 42
        )
        stab_df = pd.DataFrame({"gene": stab_feature_names, "stab_prob": stab_probs})
        stab_df = stab_df.merge(
            weights_df[["gene", "weight_mean", "weight_sd"]], on="gene", how="left"
        )
        stab_df["stable"] = stab_df["stab_prob"] >= stab_min_prob
        stab_csv = Path(out_dir) / f"{_safe_name(virus)}_stability.csv"
        _add_symbol_column(stab_df, symbol_map).to_csv(stab_csv, index=False)
        log.info("[%s] Stability probabilities written to %s", virus, stab_csv)

        # ── Depth-confound diagnostics (always on; F-001/F-003) ────────────
        # The raw >=N-UMI label tracks sequencing depth, so a naive host-gene
        # AUC is partly a depth artifact. We always report (a) the AUC using
        # depth ALONE under the identical split, and (b) per-gene depth-adjusted
        # E-values on the stably selected genes — so a reader sees the confound
        # next to the headline number without having to run an external script.
        depth_alone = _depth_alone_auc(
            np.where(vp_v)[0], np.where(~vp_v)[0], depth_v, seeds, top_depth_frac=top_depth_frac
        )
        stable_genes_list = stab_df.loc[stab_df["stable"], "gene"].tolist()
        name_to_col = {g: i for i, g in enumerate(stab_feature_names)}
        ev_cols = [name_to_col[g] for g in stable_genes_list if g in name_to_col]
        n_evalue_ge2 = 0
        depth_adj_auc = None
        if ev_cols:
            ev_df = _per_gene_evalues(
                X_stab_v[:, ev_cols],
                [stab_feature_names[i] for i in ev_cols],
                vp_v.astype(int),
                depth_v,
                pct_mito=pct_mito_v,
                mt_total_counts=mt_total_v,
                raw_total=raw_total_v,
                mt_self_counts=mt_self_v,
            )
            ev_df.insert(0, "virus", virus)
            ev_csv = Path(out_dir) / f"{_safe_name(virus)}_depth_diagnostics.csv"
            _add_symbol_column(ev_df, symbol_map).to_csv(ev_csv, index=False)
            n_evalue_ge2 = int((ev_df["E_value"] >= 2.0).sum())
            depth_adj_auc = _panel_depth_adjusted_auc(
                np.where(vp_v)[0],
                np.where(~vp_v)[0],
                X_stab_v[:, ev_cols],
                depth_v,
                seeds,
                top_depth_frac=top_depth_frac,
            )
            log.info("[%s] Depth-adjusted E-values written to %s", virus, ev_csv)

        # ── Optional genome-wide depth-(and mito-)adjusted differential test ──
        n_diff_fdr05 = None
        if differential:
            diff_df = _genome_wide_differential(
                X_stab_v, stab_feature_names, vp_v.astype(int), depth_v, pct_mito=pct_mito_v
            )
            diff_df.insert(0, "virus", virus)
            diff_csv = Path(out_dir) / f"{_safe_name(virus)}_differential.csv"
            _add_symbol_column(diff_df, symbol_map).to_csv(diff_csv, index=False)
            n_diff_fdr05 = int((diff_df["fdr"] < 0.05).sum())
            log.info(
                "[%s] Genome-wide differential (%d genes, %d at FDR<0.05) written to %s",
                virus,
                len(diff_df),
                n_diff_fdr05,
                diff_csv,
            )
            if do_enrichment:
                sig_genes = diff_df.loc[diff_df["fdr"] < 0.05, "gene"].tolist()
                # Enrichr needs symbols; map when a symbol table is available.
                if symbol_map:
                    sig_genes = [symbol_map.get(str(g).split(".")[0], str(g)) for g in sig_genes]
                    bg = [symbol_map.get(str(g).split(".")[0], str(g)) for g in background_names]
                else:
                    bg = background_names
                if sig_genes:
                    _run_enrichment(sig_genes[:top_n_genes], bg, virus, out_dir, enrichment_db)

        # Collect summary metrics.
        row: dict = {
            "virus": virus,
            "n_positive": n_pos,
            "label": label,
            "depth_matched": depth_match,
            "mito_controlled": bool(control_mito and pct_mito_v is not None),
        }
        for metric, vals in (metrics or {}).items():
            row[f"{metric}_mean"] = vals["mean"]
            row[f"{metric}_sd"] = vals["sd"]
        if depth_alone is not None:
            row["depth_alone_auc_mean"] = depth_alone["mean"]
            row["depth_alone_auc_sd"] = depth_alone["sd"]
        row["n_stable_genes"] = len(ev_cols)
        row["n_genes_evalue_ge2"] = n_evalue_ge2
        if depth_adj_auc is not None:
            row["model_auc_depth_adjusted_mean"] = depth_adj_auc["mean"]
            row["model_auc_depth_adjusted_sd"] = depth_adj_auc["sd"]
        if n_diff_fdr05 is not None:
            row["n_differential_fdr05"] = n_diff_fdr05
        all_metrics.append(row)

        # Surface the confound verdict in the log so the honesty signal is
        # on-by-default, not buried in a CSV.
        if depth_alone is not None and metrics and "auc" in metrics:
            model_auc = metrics["auc"]["mean"]
            _depth_adj_str = (
                f"; panel+depth AUC {depth_adj_auc['mean']:.3f}"
                if depth_adj_auc is not None
                else ""
            )
            log.info(
                "[%s] Depth-confound check: model AUC %.3f vs depth-ALONE AUC %.3f "
                "(same balanced design)%s; %d/%d stable genes have E-value>=2.",
                virus,
                model_auc,
                depth_alone["mean"],
                _depth_adj_str,
                n_evalue_ge2,
                len(ev_cols),
            )
            if depth_alone["mean"] >= model_auc - 0.02:
                log.warning(
                    "[%s] Sequencing depth alone predicts virus status about as well as the "
                    "host-gene model (%.3f vs %.3f) — treat the headline AUC as depth-confounded "
                    "(F-001/F-003). Prefer the per-gene E-values, which are depth-adjusted.",
                    virus,
                    depth_alone["mean"],
                    model_auc,
                )

        # ── Optional pathway enrichment ───────────────────────────────────
        if do_enrichment:
            stable_genes = (
                stab_df[stab_df["stable"]]
                .assign(_abs_w=lambda df: df["weight_mean"].abs())
                .sort_values("_abs_w", ascending=False)["gene"]
                .tolist()
            )
            if stable_genes:
                _run_enrichment(
                    stable_genes[:top_n_genes],
                    background_names,
                    virus,
                    out_dir,
                    enrichment_db,
                )
            else:
                log.info("[%s] No stably selected genes; skipping enrichment.", virus)

    if all_metrics:
        metrics_csv = Path(out_dir) / "hostresponse_metrics.csv"
        pd.DataFrame(all_metrics).to_csv(metrics_csv, index=False)
        log.info("Summary metrics written to %s", metrics_csv)
    else:
        log.warning(
            "No virus met the minimum positive-cell threshold (%d). "
            "Consider lowering --detection-threshold or using a sample with higher viral load.",
            MIN_VIRUS_CELLS,
        )


# ── Snakemake entry point ─────────────────────────────────────────────────────

if "snakemake" in globals():
    cfg = RunConfig.from_yaml(snakemake.params.configfile)  # noqa: F821
    kb = KbCountOutputs.from_config_output(cfg.output)
    _virus_h5ad = str(kb.current_adata(multimapping=cfg.multimapping))
    _viral_acc_file = f"{cfg.output}log/analysis.txt"
    _out_dir = os.path.join(cfg.output, "hostresponse")
    _seeds = DEFAULT_SEEDS[: cfg.hostresponse_n_seeds]

    assert cfg.host_h5ad is not None, "hostresponse runs only when host_h5ad is set"
    run_hostresponse(
        virus_h5ad=_virus_h5ad,
        host_h5ad=cfg.host_h5ad,
        viral_accessions_file=_viral_acc_file,
        out_dir=_out_dir,
        use_hvg=cfg.hostresponse_use_hvg,
        seeds=_seeds,
        n_stab_iter=cfg.hostresponse_n_stab_iter,
        stab_min_prob=cfg.hostresponse_stab_min_prob,
        top_n_genes=cfg.hostresponse_top_n_genes,
        detection_threshold=cfg.detection_threshold,
        do_enrichment=cfg.hostresponse_enrichment,
        enrichment_db=cfg.hostresponse_enrichment_db,
        label=cfg.hostresponse_label,
        depth_match=cfg.hostresponse_depth_match,
        control_mito=cfg.hostresponse_control_mito,
        differential=cfg.hostresponse_differential,
    )

    Path(str(snakemake.output[0])).touch()  # noqa: F821


# ── Standalone CLI ────────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Associate viral presence with host gene expression (logistic regression)."
    )
    p.add_argument(
        "--virus-h5ad",
        required=True,
        metavar="PATH",
        help="Path to the virus count h5ad (from a viralscan run).",
    )
    p.add_argument(
        "--host-h5ad",
        required=True,
        metavar="PATH",
        help="Path to host gene-expression h5ad (cells × genes).",
    )
    p.add_argument(
        "--viral-accessions",
        required=True,
        metavar="PATH",
        help="Path to analysis.txt from the same viralscan run.",
    )
    p.add_argument(
        "--output", "-o", required=True, metavar="DIR", help="Output directory for results."
    )
    p.add_argument(
        "--use-hvg",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use highly variable genes (default: on).",
    )
    p.add_argument(
        "--n-seeds",
        type=int,
        default=6,
        help="Number of random seeds for multi-seed L2 regression.",
    )
    p.add_argument("--n-stab-iter", type=int, default=100, help="Stability selection iterations.")
    p.add_argument(
        "--stab-min-prob",
        type=float,
        default=0.6,
        help="Min selection probability to call a gene stably selected.",
    )
    p.add_argument(
        "--top-n-genes",
        type=int,
        default=50,
        help="Top N stable genes to pass to pathway enrichment.",
    )
    p.add_argument(
        "--detection-threshold",
        type=int,
        default=10,
        help="Min UMI count to call a cell virus-positive (validated default).",
    )
    p.add_argument(
        "--label",
        choices=LABEL_CHOICES,
        default="raw",
        help=(
            "Positive-call label. 'raw' (default) uses counts>=detection-threshold and is "
            "depth-confounded; 'cpm'/'fraction' use a depth-normalized, prevalence-matched "
            "label (viral burden per host UMI) that decouples the label from sequencing depth."
        ),
    )
    p.add_argument(
        "--depth-match",
        action="store_true",
        default=False,
        help=(
            "Restrict each virus's analysis to a coarsened-exact depth-matched cohort so any "
            "surviving host-gene signal is depth-independent by construction."
        ),
    )
    p.add_argument(
        "--mito-control",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Add per-cell %%mitochondrial content as a covariate to the per-gene E-values so a "
            "mito-QC artifact cannot pass as a host-response gene (default: on; --no-mito-control "
            "to disable)."
        ),
    )
    p.add_argument(
        "--gene-symbols",
        action="store_true",
        default=False,
        help=(
            "Annotate output CSVs with HGNC gene symbols mapped from Ensembl IDs via "
            "mygene.info (network; best-effort, falls back to IDs on failure)."
        ),
    )
    p.add_argument(
        "--differential",
        action="store_true",
        default=False,
        help=(
            "Write a genome-wide depth-(and %%mito-)adjusted differential table "
            "(<virus>_differential.csv: partial_r, p_value, fdr, direction) over ALL "
            "features, not just the stable panel."
        ),
    )
    p.add_argument(
        "--enrichment",
        action="store_true",
        default=False,
        help="Run pathway enrichment via gget (requires ViralScan[enrichment]).",
    )
    p.add_argument(
        "--enrichment-db",
        default="GO_Biological_Process_2023",
        help="gget.enrichr database (default: GO_Biological_Process_2023).",
    )
    p.add_argument("--verbose", action="store_true", default=False)
    return p


if __name__ == "__main__":
    import sys  # noqa: PLC0415

    args = _build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )
    run_hostresponse(
        virus_h5ad=args.virus_h5ad,
        host_h5ad=args.host_h5ad,
        viral_accessions_file=args.viral_accessions,
        out_dir=args.output,
        use_hvg=args.use_hvg,
        seeds=DEFAULT_SEEDS[: args.n_seeds],
        n_stab_iter=args.n_stab_iter,
        stab_min_prob=args.stab_min_prob,
        top_n_genes=args.top_n_genes,
        detection_threshold=args.detection_threshold,
        do_enrichment=args.enrichment,
        enrichment_db=args.enrichment_db,
        label=args.label,
        depth_match=args.depth_match,
        control_mito=args.mito_control,
        annotate_symbols=args.gene_symbols,
        differential=args.differential,
    )
