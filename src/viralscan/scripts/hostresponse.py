"""
Associate viral presence with host gene expression via logistic regression.

Implements the approach from Luebbert et al. 2025 (Nature Biotechnology):
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
import hashlib
import json
import logging
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import sklearn as _sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

from viralscan.kb_outputs import KbCountOutputs
from viralscan.run_safety import software_identity
from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging
from viralscan.virus_grouping import identity_path
from viralscan.virus_identity import TABLE_FILENAME, VirusIdentityTable

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


#: Per-virus CSV suffixes hostresponse writes. Cleared before a rerun because a
#: virus that falls below MIN_VIRUS_CELLS under a new multimap method is skipped
#: with `continue`, leaving the previous method's file in place.
_OWNED_OUTPUT_SUFFIXES = (
    "_gene_weights.csv",
    "_stability.csv",
    "_depth_diagnostics.csv",
    "_differential.csv",
    "_cv_folds.tsv",
    "_cv_metrics.tsv",
    "_cv_summary.tsv",
    "_cv_group_auc.tsv",
    "_permutation_metrics.tsv",
)

#: Run-level status table (one row per virus x stratum); always rewritten.
STATUS_FILENAME = "hostresponse_status.tsv"
MANIFEST_FILENAME = "hostresponse_manifest.json"
_STATUS_COLUMNS = (
    "stratum",
    "virus",
    "status",
    "reason",
    "n_cells",
    "n_positive",
    "n_negative",
    "n_groups",
)


def _status_row(
    virus: str, status: str, reason: str = "", *, stratum: str = "ALL", **counts
) -> dict:
    """One ``hostresponse_status.tsv`` row; unspecified counts stay empty."""
    row = dict.fromkeys(_STATUS_COLUMNS, "")
    row.update(stratum=stratum, virus=virus, status=status, reason=reason, **counts)
    return row


def clear_stale_virus_outputs(out_dir) -> list[str]:
    """Remove the previous run's per-virus CSVs before regenerating.

    Which viruses clear MIN_VIRUS_CELLS depends on the multimap method, so
    `rerun-multimap` can legitimately drop one from the set. The directory was
    only ever created with ``exist_ok=True`` and the skip path is a bare
    ``continue``, so the demoted virus's CSVs survived into a tree labelled with
    the new method (SW-04).

    Enrichment CSVs carry a ``_enrichment_<database>.csv`` suffix and are matched
    by prefix rather than by a fixed name, since the database is configurable.
    """
    directory = Path(out_dir)
    if not directory.is_dir():
        return []
    removed = []
    # Only owned artifacts inside our stratum namespace are removed; retain user files.
    for subdir in sorted(directory.glob("cell_type_*")):
        if subdir.is_dir() and not subdir.is_symlink():
            removed.extend(clear_stale_virus_outputs(subdir))
            (subdir / "hostresponse_metrics.csv").unlink(missing_ok=True)
            if not any(subdir.iterdir()):
                subdir.rmdir()
    (directory / MANIFEST_FILENAME).unlink(missing_ok=True)
    for path in sorted([*directory.glob("*.csv"), *directory.glob("*.tsv")]):
        name = path.name
        if (
            name.endswith(_OWNED_OUTPUT_SUFFIXES)
            or "_enrichment_" in name
            or name == STATUS_FILENAME
        ):
            path.unlink()
            removed.append(str(path))
    if removed:
        log.info("cleared %d stale per-virus hostresponse file(s)", len(removed))
    return removed


def _safe_name(name: str) -> str:
    """Make a filesystem-safe version of a virus accession."""
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in name)


def _load_viral_accessions(analysis_txt: str) -> set:
    """Read the Run's viral gene IDs.

    ``analysis_txt`` is either a Virus Identity table (``virus_identity.tsv``;
    the ``viral`` column decides) or the legacy gene-ID list of the analysis
    rule (one ID per line).
    """
    if os.path.basename(analysis_txt) == TABLE_FILENAME:
        return set(VirusIdentityTable.read_tsv(analysis_txt).viral_gene_ids())
    accessions: set = set()
    with open(analysis_txt) as fh:
        for line in fh:
            v = line.strip()
            if v:
                accessions.add(v)
    return accessions


def _sha256_file(path) -> str:
    """Hex SHA-256 of a file, streamed (input-hash provenance)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _resolve_raw_depth(adata):
    """Return ``(per-cell raw library size, source label)`` for the host matrix.

    Source preference: ``layers["counts"]`` (an explicit raw-count layer), else
    ``X``. When ``X`` carries ``uns["log1p"]`` and there is no counts layer the
    row sums are NOT library sizes; the label says so (``X_log1p_not_raw``) and
    the caller records it in provenance instead of presenting it as raw depth.
    Raises ``ValueError`` if the chosen source has negative values (scaled data
    cannot define a depth).
    """
    if "counts" in adata.layers:
        src, label = adata.layers["counts"], "layers['counts']"
    else:
        src = adata.X
        label = "X_log1p_not_raw" if "log1p" in adata.uns else "X"
    values = src.data if sp.issparse(src) else np.asarray(src)
    if not np.isfinite(values).all():
        raise ValueError(f"raw-depth source {label} has nonfinite values")
    if src.shape[0] and float(src.min()) < 0:
        raise ValueError(f"raw-depth source {label} has negative values (scaled data?)")
    # np.asarray(...).flatten() handles both scipy sparse (via np.matrix) and ndarray.
    return np.asarray(src.sum(axis=1)).flatten(), label


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
    raw_depth, depth_source = _resolve_raw_depth(host_adata)
    host_adata.obs["_raw_depth"] = raw_depth
    host_adata.uns["_raw_depth_source"] = depth_source

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
        lambda e: (
            "robust"
            if (np.isfinite(e) and e >= 3.0)
            else "moderate"
            if (np.isfinite(e) and e >= 1.5)
            else "fragile"
        )
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


# ── HR-01: group-disjoint evaluation ──────────────────────────────────────────
#
# Cell-split mode (``_run_l2_regression``) is untouched. The functions below add
# an opt-in group mode in which whole groups (donor/sample) are held out, so the
# reported metrics estimate generalisation to NEW groups, not to new cells of
# groups already seen in training.

HOSTRESPONSE_SCHEMA_VERSION = 1
CV_MODES = ("cell", "group")
DEFAULT_CV_FOLDS = 5
MIN_FOLD_CLASS_CELLS = 2

#: Frozen aggregation statistic (HR-01/HR-02). Headline AUC = mean over seeds of
#: the mean over folds of the per-fold held-out AUC, i.e. every fold has equal
#: weight regardless of how many cells it holds (a pooled-cell AUC would let the
#: largest group dominate). Equal-group AUC is reported alongside it.
AUC_STATISTIC = "mean_over_seeds_of_fold_mean_test_auc"

_METRIC_NAMES = ("sensitivity", "specificity", "balanced_accuracy", "auc", "mcc")
_FOLD_COLUMNS = (
    "model",
    "seed",
    "fold",
    "status",
    "reason",
    "n_train_cells",
    "n_test_cells",
    "n_train_groups",
    "n_test_groups",
    *_METRIC_NAMES,
)


class FoldDesignError(ValueError):
    """A valid group-disjoint fold design cannot be built.

    ``reason`` is a stable machine-readable code (``insufficient_groups``,
    ``single_class``, ``duplicate_cell_ids``, ``invalid_folds``,
    ``invalid_fold_classes``, ``length_mismatch``) so callers can store it as a
    failure status instead of parsing the message.
    """

    def __init__(self, reason: str, message: str):
        super().__init__(f"{reason}: {message}")
        self.reason = reason


def validate_obs_column(obs: pd.DataFrame, column: str, *, what: str = "group") -> np.ndarray:
    """Return ``obs[column]`` as a string array, failing on missing metadata.

    Raises ``ValueError`` when the column is absent or any value is missing/blank.
    Group CV never falls back to cell CV: the caller must let this propagate.
    """
    if column not in obs.columns:
        raise ValueError(f"{what} column {column!r} not found in obs (have: {sorted(obs.columns)})")
    series = obs[column]
    if series.isna().any():
        raise ValueError(
            f"{what} column {column!r} has {int(series.isna().sum())} missing value(s)"
        )
    values = np.asarray(series.astype(str).to_numpy())
    if (pd.Series(values).str.strip() == "").any():
        raise ValueError(f"{what} column {column!r} has blank value(s)")
    return values


def _canonical_order(cell_ids) -> np.ndarray:
    """Indices that sort cells by id: the row-order-invariant processing order."""
    ids = np.asarray(cell_ids).astype(str)
    if len(set(ids.tolist())) != len(ids):
        raise FoldDesignError("duplicate_cell_ids", "cell ids must be unique")
    return np.asarray(np.argsort(ids, kind="stable"))


@dataclass(frozen=True)
class FoldPlan:
    """Group-disjoint fold assignment for one seed.

    ``fold[i]`` is the held-out fold of input cell ``i`` (input order).
    ``method`` is the splitter that produced it; ``stratified`` is False when the
    GroupKFold fallback had to be used and ``note`` then says why stratification
    was not possible.
    """

    fold: np.ndarray
    n_folds: int
    method: str
    stratified: bool
    note: str
    seed: int


def _plan_class_problem(y, fold, n_folds: int, min_class_cells: int):
    for k in range(n_folds):
        test = fold == k
        for side, mask in (("test", test), ("train", ~test)):
            for cls in (0, 1):
                n = int((y[mask] == cls).sum())
                if n < min_class_cells:
                    return (
                        f"fold {k} {side} has {n} class-{cls} cell(s) (need >= {min_class_cells})"
                    )
    return None


def make_group_folds(
    y,
    groups,
    cell_ids,
    n_folds: int = DEFAULT_CV_FOLDS,
    seed: int = 0,
    min_class_cells: int = 1,
) -> FoldPlan:
    """Assign every cell to a held-out fold so that no group spans two folds.

    Cells and groups are processed in canonical (sorted-id) order, so the
    assignment depends only on *which* cell has which label/group, never on the
    row order of the inputs. ``StratifiedGroupKFold(shuffle=True,
    random_state=seed)`` is tried first with
    ``n_splits=min(n_folds, n_groups)``; if it cannot be built, or leaves a train
    or test side with fewer than ``min_class_cells`` cells of either class, the
    deterministic ``GroupKFold`` fallback is tried and the reason recorded in
    ``FoldPlan.note``. If neither is valid a :class:`FoldDesignError` is raised;
    invalid folds are never dropped.
    """
    from sklearn import model_selection  # noqa: PLC0415

    StratifiedGroupKFold = getattr(model_selection, "StratifiedGroupKFold", None)
    if not callable(StratifiedGroupKFold):
        raise FoldDesignError(
            "splitter_unavailable", "group CV requires sklearn StratifiedGroupKFold"
        )
    GroupKFold = model_selection.GroupKFold

    y_arr = np.asarray(y).astype(int)
    g_arr = validate_obs_column(pd.DataFrame({"group": groups}), "group")
    if not (len(y_arr) == len(g_arr) == len(cell_ids)):
        raise FoldDesignError("length_mismatch", "y, groups and cell_ids must have equal length")
    order = _canonical_order(cell_ids)
    if set(np.unique(y_arr).tolist()) != {0, 1}:
        raise FoldDesignError("single_class", "both classes are required to build folds")
    n_groups = len(set(g_arr.tolist()))
    n_splits = min(int(n_folds), n_groups)
    if n_splits < 2:
        raise FoldDesignError(
            "insufficient_groups",
            f"{n_groups} group(s) and n_folds={n_folds}: at least 2 folds/groups are required",
        )
    y_c, g_c = y_arr[order], g_arr[order]
    dummy = np.zeros((len(y_c), 1))
    splitters = (
        ("StratifiedGroupKFold", StratifiedGroupKFold(n_splits, shuffle=True, random_state=seed)),
        ("GroupKFold", GroupKFold(n_splits)),
    )
    notes: list = []
    for name, splitter in splitters:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                splits = list(splitter.split(dummy, y_c, g_c))
        except ValueError as exc:
            notes.append(f"{name} unavailable ({exc})")
            continue
        fold_c = np.full(len(y_c), -1, dtype=int)
        for k, (_, test) in enumerate(splits):
            fold_c[test] = k
        problem = _plan_class_problem(y_c, fold_c, n_splits, min_class_cells)
        if problem:
            notes.append(f"{name} invalid ({problem})")
            continue
        fold = np.empty(len(y_c), dtype=int)
        fold[order] = fold_c
        stratified = name == "StratifiedGroupKFold"
        if not stratified:
            log.info(
                "group CV seed %d: stratification not possible; using GroupKFold (%s)",
                seed,
                "; ".join(notes),
            )
        return FoldPlan(fold, n_splits, name, stratified, "; ".join(notes), int(seed))
    raise FoldDesignError("invalid_folds", "; ".join(notes))


def _sub_seed(seed: int, fold: int) -> int:
    """Deterministic per-(seed, fold) integer seed for fold-local randomness."""
    return int(np.random.SeedSequence([int(seed), int(fold)]).generate_state(1)[0])


def _side_cohort(
    idx: np.ndarray,
    y: np.ndarray,
    depth: np.ndarray,
    top_depth_frac: float,
    depth_match: bool,
    seed: int,
) -> np.ndarray:
    """Fold-local depth filter/matching for ONE side (train or test) of a fold.

    Mirrors the cell-mode cohort rules (per-class top-depth filter, or coarsened
    depth matching which replaces it) but is computed from that side's own cells
    only, so no held-out depth/label information shapes the training cohort and
    the test cohort is defined without reference to training cells. ``idx`` must
    already be in canonical order; the result keeps that order.
    """
    if len(idx) == 0:
        return idx
    if depth_match:
        return np.asarray(idx[_depth_match_indices(y[idx], depth[idx], seed=seed)])
    if top_depth_frac >= 1.0:
        return idx
    keep = []
    for cls in (0, 1):
        ii = idx[y[idx] == cls]
        if len(ii):
            d = depth[ii]
            keep.append(ii[d >= np.percentile(d, (1 - top_depth_frac) * 100)])
    kept = np.concatenate(keep) if keep else idx[:0]
    return np.asarray(idx[np.isin(idx, kept)])


def _class_counts_problem(y, idx, minimum: int):
    for cls in (0, 1):
        n = int((y[idx] == cls).sum())
        if n < minimum:
            return f"{n} class-{cls} cell(s) (need >= {minimum})"
    return None


def _fold_metrics(y_true, prob) -> dict:
    y_true = np.asarray(y_true)
    pred = np.asarray(prob) >= 0.5
    return {
        "sensitivity": float(pred[y_true == 1].mean()),
        "specificity": float((~pred[y_true == 0]).mean()),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)),
        "auc": float(roc_auc_score(y_true, prob)),
        "mcc": float(matthews_corrcoef(y_true, pred)),
    }


def summarize_fold_metrics(fold_metrics: pd.DataFrame) -> pd.DataFrame:
    """Per model and metric: mean, sd (ddof=0), median, min, max, n_valid_folds.

    Only ``status == "ok"`` rows with a finite value count. Folds are repeated
    over seeds, so the SD describes fold-to-fold variation, not independent
    biological replicates.
    """
    cols = ["model", "metric", "mean", "sd", "median", "min", "max", "n_valid_folds"]
    rows = []
    for model in dict.fromkeys(fold_metrics["model"]):
        sub = fold_metrics[(fold_metrics["model"] == model) & (fold_metrics["status"] == "ok")]
        for metric in _METRIC_NAMES:
            v = pd.to_numeric(sub[metric], errors="coerce").dropna().to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            nan = float("nan")
            rows.append(
                {
                    "model": model,
                    "metric": metric,
                    "mean": float(v.mean()) if len(v) else nan,
                    "sd": float(v.std()) if len(v) else nan,
                    "median": float(np.median(v)) if len(v) else nan,
                    "min": float(v.min()) if len(v) else nan,
                    "max": float(v.max()) if len(v) else nan,
                    "n_valid_folds": int(len(v)),
                }
            )
    return pd.DataFrame(rows, columns=cols)


@dataclass
class CVResult:
    """Outcome of :func:`evaluate_group_cv`.

    ``status`` is ``ok`` (every requested seed produced a valid fold design),
    ``partial`` (some seeds failed; their reasons are in ``failures``) or
    ``failed``. ``statistic`` holds the frozen headline numbers
    (:data:`AUC_STATISTIC`) and is ``None`` unless ``status == "ok"`` — inference
    must not use a partially valid result.
    """

    status: str
    fold_metrics: pd.DataFrame
    summary: pd.DataFrame
    folds: pd.DataFrame
    group_auc: pd.DataFrame
    failures: list
    statistic: Optional[dict]
    weights: Optional[pd.DataFrame]
    fold_features: dict
    n_groups: int
    n_cells: int
    n_seeds_requested: int
    n_seeds_valid: int
    plan_methods: dict = field(default_factory=dict)
    feature_names: list = field(default_factory=list)

    def provenance(self) -> dict:
        """JSON-serialisable description of what was evaluated."""
        digests = {}
        for seed, sub in self.folds.groupby("seed"):
            text = "\n".join(f"{c}\t{f}" for c, f in zip(sub["cell_id"], sub["fold"]))
            digests[int(seed)] = hashlib.sha256(text.encode()).hexdigest()
        return {
            "status": self.status,
            "auc_statistic": AUC_STATISTIC,
            "n_groups": self.n_groups,
            "n_cells": self.n_cells,
            "n_seeds_requested": self.n_seeds_requested,
            "n_seeds_valid": self.n_seeds_valid,
            "fold_methods": {int(k): v for k, v in self.plan_methods.items()},
            "fold_digests": digests,
            "failures": list(self.failures),
            "feature_names": list(self.feature_names),
            "fold_features": {
                f"{seed}:{fold}": values for (seed, fold), values in self.fold_features.items()
            },
        }


# ── HR-04: same-fold baselines ────────────────────────────────────────────────

#: Low-dimensional nuisance baselines (never a high-dimensional nuisance model).
BASELINE_MODELS = ("depth_only", "cell_type_only", "depth_plus_cell_type")
_BASELINE_KEYS = {
    "depth_only": "depth_auc",
    "cell_type_only": "cell_type_auc",
    "depth_plus_cell_type": "depth_cell_type_auc",
}


def _baseline_design(model: str, logd_tr, logd_te, ct_tr=None, ct_te=None):
    """Design matrices ``(train, test)`` for a nuisance baseline.

    Everything learned is fitted on the TRAINING rows only: the log-depth
    ``StandardScaler`` and the cell-type level set (sorted, so deterministic).
    A test cell whose type is absent from training gets an all-zero indicator
    row (the defined handling of unseen levels).
    """
    tr_blocks, te_blocks = [], []
    if model in ("depth_only", "depth_plus_cell_type"):
        sc = StandardScaler().fit(np.reshape(logd_tr, (-1, 1)))
        tr_blocks.append(sc.transform(np.reshape(logd_tr, (-1, 1))))
        te_blocks.append(sc.transform(np.reshape(logd_te, (-1, 1))))
    if model in ("cell_type_only", "depth_plus_cell_type"):
        levels = sorted(set(np.asarray(ct_tr).astype(str).tolist()))
        for ct, blocks in ((ct_tr, tr_blocks), (ct_te, te_blocks)):
            arr = np.asarray(ct).astype(str)
            blocks.append(np.stack([(arr == lv) for lv in levels], axis=1).astype(float))
    if not tr_blocks:
        raise ValueError(f"unknown baseline model {model!r}")
    return np.hstack(tr_blocks), np.hstack(te_blocks)


def _fit_prob(x_tr, y_tr, x_te, seed: int, class_weight="balanced") -> np.ndarray:
    """Fit the shared L2 logistic regression; return held-out P(y=1)."""
    model = LogisticRegression(
        max_iter=1000, class_weight=class_weight, random_state=int(seed)
    ).fit(x_tr, y_tr)
    return np.asarray(model.predict_proba(x_te)[:, 1])


def _panel_depth_fold_prob(x_tr, y_tr, x_te, logd_tr, logd_te, iters, min_prob, seed):
    """Stable-panel + depth probabilities with the panel selected in-fold.

    Randomized-lasso stability selection runs on the training cohort alone; the
    whole-cohort stability table stays descriptive. Returns ``None`` when no gene
    reaches ``min_prob`` (the fold is then reported ``not_estimable``).
    """
    stable = _run_stability_selection(x_tr, y_tr.astype(bool), iters, seed=seed) >= min_prob
    if not stable.any():
        return None
    f_tr = np.column_stack([x_tr[:, stable], logd_tr])
    f_te = np.column_stack([x_te[:, stable], logd_te])
    sc = StandardScaler().fit(f_tr)
    return _fit_prob(sc.transform(f_tr), y_tr, sc.transform(f_te), seed)


def _cell_mode_baselines(
    pos_idx, neg_idx, depth, cell_type, seeds, models=BASELINE_MODELS, top_depth_frac=TOP_DEPTH_FRAC
):
    """Cell-split nuisance baselines on the SAME balanced splits as the expression model.

    The splits come from :func:`_balanced_split` with identical arguments, so the
    train/test cells match :func:`_run_l2_regression` exactly (the split depends
    only on indices, depth and seed). Scaler and cell-type levels are fitted on
    each split's training cells. Returns ``{model: {"mean","sd"} | None}``;
    cell-type models are ``None`` when ``cell_type`` is not given.
    """
    depth = np.asarray(depth, dtype=float)
    logd = np.log1p(depth)
    ct = None if cell_type is None else np.asarray(cell_type).astype(str)
    index_col = np.arange(len(depth)).reshape(-1, 1)  # lets _balanced_split return indices
    aucs: dict = {m: [] for m in models}
    for seed in seeds:
        split = _balanced_split(pos_idx, neg_idx, depth, index_col, seed, top_depth_frac)
        if split is None:
            continue
        tr_col, y_tr, te_pos, te_neg = split
        tr = tr_col[:, 0]
        te = np.concatenate([te_pos[:, 0], te_neg[:, 0]])
        y_te = np.array([1] * len(te_pos) + [0] * len(te_neg))
        for m in models:
            if m != "depth_only" and ct is None:
                continue
            x_tr, x_te = _baseline_design(
                m,
                logd[tr],
                logd[te],
                None if ct is None else ct[tr],
                None if ct is None else ct[te],
            )
            prob = _fit_prob(x_tr, y_tr, x_te, seed, class_weight=None)
            with contextlib.suppress(ValueError):
                aucs[m].append(float(roc_auc_score(y_te, prob)))
    return {
        m: ({"mean": float(np.mean(v)), "sd": float(np.std(v))} if v else None)
        for m, v in aucs.items()
    }


def baseline_comparison(expression_auc, baseline_aucs: dict) -> dict:
    """Headline expression-vs-baseline AUC columns.

    ``baseline_aucs`` maps baseline model name to its mean AUC (or ``None``/NaN
    when not assessed). Returns ``expression_auc``, ``depth_auc``,
    ``cell_type_auc``, ``depth_cell_type_auc`` and
    ``expression_minus_best_baseline_auc``; unassessed entries are ``None`` and the
    difference is ``None`` when no baseline was assessed (never zero-filled).
    """

    def _num(v) -> Optional[float]:
        return float(v) if v is not None and np.isfinite(v) else None

    expr = _num(expression_auc)
    out: dict = {"expression_auc": expr}
    present = []
    for model, key in _BASELINE_KEYS.items():
        v = _num(baseline_aucs.get(model))
        out[key] = v
        if v is not None:
            present.append(v)
    out["expression_minus_best_baseline_auc"] = (
        expr - max(present) if present and expr is not None else None
    )
    return out


def _failed_fold_row(model, seed, reason) -> dict:
    row: dict = dict.fromkeys(_FOLD_COLUMNS, float("nan"))
    row.update(model=model, seed=seed, fold=-1, status="failed", reason=reason)
    return row


def evaluate_group_cv(
    X,
    y,
    groups,
    depth,
    cell_ids,
    *,
    seeds,
    n_folds: int = DEFAULT_CV_FOLDS,
    feature_names=None,
    use_hvg: bool = True,
    top_depth_frac: float = TOP_DEPTH_FRAC,
    depth_match: bool = False,
    min_class_cells: int = MIN_FOLD_CLASS_CELLS,
    collect_weights: bool = False,
    cell_type=None,
    baselines=(),
    panel_iters: int = 0,
    stab_min_prob: float = 0.6,
) -> CVResult:
    """Group-disjoint cross-validation of the host-expression L2 model.

    ``baselines`` is a subset of :data:`BASELINE_MODELS` (``depth_only``,
    ``cell_type_only``, ``depth_plus_cell_type``) plus ``panel_depth`` (stable-gene
    panel + depth, panel chosen by randomized-lasso stability selection with
    ``panel_iters`` iterations INSIDE each training fold; needs ``panel_iters > 0``).
    Each baseline is scored on exactly the expression model's folds/cohorts/seeds.
    ``cell_type`` (one label per cell) is required for the cell-type baselines;
    its categories are learned from the training cohort and an unseen test level
    maps to the all-zero indicator row.

    For every seed: build a group-disjoint :class:`FoldPlan`, derive each fold's
    train/test cohorts (fold-local depth filtering/matching, see
    :func:`_side_cohort`), require both classes on both sides of every fold, then
    fit and score. Anything that learns from expression (HVG selection) is fitted
    on the training cohort of that fold only. A seed whose design is invalid is
    recorded in ``failures`` and a ``failed`` row in ``fold_metrics``; its folds
    are never silently dropped or repaired.

    ``X`` is the normalised cell x gene matrix, ``y`` the 0/1 label, ``groups``
    the donor/sample per cell, ``depth`` the raw host depth and ``cell_ids`` the
    unique barcodes. Processing follows sorted barcode order so results do not
    depend on row order. Deterministic for fixed inputs and seeds.
    """
    seeds = list(seeds)
    y = np.asarray(y).astype(int)
    groups = np.asarray(groups).astype(str)
    depth = np.asarray(depth, dtype=float)
    ids = np.asarray(cell_ids).astype(str)
    n = len(y)
    if not (X.shape[0] == len(groups) == len(depth) == len(ids) == n):
        raise FoldDesignError("length_mismatch", "X, y, groups, depth, cell_ids differ in length")
    baselines = tuple(baselines)
    unknown = set(baselines) - set(BASELINE_MODELS) - {"panel_depth"}
    if unknown:
        raise ValueError(f"unknown baseline model(s) {sorted(unknown)}")
    ct = None if cell_type is None else np.asarray(cell_type).astype(str)
    if ct is not None and len(ct) != n:
        raise FoldDesignError("length_mismatch", "cell_type differs in length from y")
    if ct is None and {"cell_type_only", "depth_plus_cell_type"} & set(baselines):
        raise ValueError("cell-type baselines require cell_type")
    if "panel_depth" in baselines and panel_iters <= 0:
        raise ValueError("panel_depth baseline requires panel_iters > 0")
    if not np.isfinite(depth).all() or (depth < 0).any():
        raise ValueError("depth must be finite and nonnegative")
    if not 0 < top_depth_frac <= 1 or min_class_cells < 1:
        raise ValueError("top_depth_frac must be in (0, 1] and min_class_cells >= 1")
    if len(set(seeds)) != len(seeds) or any(int(seed) < 0 for seed in seeds):
        raise ValueError("seeds must be unique and nonnegative")
    logd = np.log1p(depth)
    canon = _canonical_order(ids)
    n_features = X.shape[1]
    fnames = (
        list(feature_names) if feature_names is not None else [f"f{j}" for j in range(n_features)]
    )

    metric_rows: list = []
    member_rows: list = []
    group_auc_rows: list = []
    failures: list = []
    fold_features: dict = {}
    plan_methods: dict = {}
    seed_stats: list = []
    coef_sum = np.zeros(n_features)
    coef_sq = np.zeros(n_features)
    coef_cnt = np.zeros(n_features, dtype=int)
    simple_weights: list = []

    for seed in seeds:
        try:
            plan = make_group_folds(y, groups, ids, n_folds, seed, min_class_cells=1)
            sides = []
            for k in range(plan.n_folds):
                in_test = plan.fold[canon] == k
                sub = _sub_seed(seed, k)
                tr = _side_cohort(canon[~in_test], y, depth, top_depth_frac, depth_match, sub)
                te = _side_cohort(canon[in_test], y, depth, top_depth_frac, depth_match, sub)
                for side, idx in (("train", tr), ("test", te)):
                    problem = _class_counts_problem(y, idx, min_class_cells)
                    if problem:
                        raise FoldDesignError(
                            "invalid_fold_classes",
                            f"seed {seed} fold {k} {side} cohort has {problem}",
                        )
                sides.append((tr, te))
        except FoldDesignError as exc:
            failures.append({"seed": int(seed), "reason": exc.reason, "message": str(exc)})
            metric_rows.append(_failed_fold_row("expression", int(seed), str(exc)))
            continue

        plan_methods[int(seed)] = plan.method
        fold_stats = []
        oof_prob = np.full(n, np.nan)
        for k, (tr, te) in enumerate(sides):
            x_tr, x_te = X[tr], X[te]
            mask = _hvg_mask(x_tr) if use_hvg else None
            if mask is not None:
                x_tr, x_te = x_tr[:, mask], x_te[:, mask]
            fold_features[(int(seed), k)] = (
                np.where(mask)[0].tolist() if mask is not None else list(range(n_features))
            )
            model = LogisticRegression(
                solver="lbfgs",
                max_iter=1000,
                C=1.0,
                class_weight="balanced",
                random_state=int(seed),
            ).fit(x_tr, y[tr])
            prob = model.predict_proba(x_te)[:, 1]
            oof_prob[te] = prob
            m = _fold_metrics(y[te], prob)
            fold_stats.append(m)
            common = {
                "seed": int(seed),
                "fold": k,
                "n_train_cells": len(tr),
                "n_test_cells": len(te),
                "n_train_groups": len(set(groups[tr].tolist())),
                "n_test_groups": len(set(groups[te].tolist())),
            }
            metric_rows.append({"model": "expression", "status": "ok", "reason": "", **common, **m})
            # Baselines: identical cohorts/folds/seed; every learned object (scaler,
            # category encoder, stable panel) is fitted on this fold's training cohort only.
            for name in baselines:
                if name == "panel_depth":
                    probs = _panel_depth_fold_prob(
                        x_tr, y[tr], x_te, logd[tr], logd[te], panel_iters, stab_min_prob, int(seed)
                    )
                    if probs is None:
                        metric_rows.append(
                            {
                                "model": name,
                                "status": "not_estimable",
                                "reason": "no_stable_panel_in_training_fold",
                                **common,
                            }
                        )
                        continue
                else:
                    b_tr, b_te = _baseline_design(
                        name,
                        logd[tr],
                        logd[te],
                        None if ct is None else ct[tr],
                        None if ct is None else ct[te],
                    )
                    probs = _fit_prob(b_tr, y[tr], b_te, int(seed))
                metric_rows.append(
                    {
                        "model": name,
                        "status": "ok",
                        "reason": "",
                        **common,
                        **_fold_metrics(y[te], probs),
                    }
                )
            coef = model.coef_[0]
            if collect_weights:
                if mask is not None:
                    gi = np.where(mask)[0]
                    coef_sum[gi] += coef
                    coef_sq[gi] += coef**2
                    coef_cnt[gi] += 1
                else:
                    simple_weights.append(coef)
        # Fold of every cell; test_eligible marks survivors of the test-cohort rule.
        eligible = np.zeros(n, dtype=bool)
        for _, te in sides:
            eligible[te] = True
        for i in canon:
            member_rows.append(
                (int(seed), ids[i], groups[i], int(plan.fold[i]), int(y[i]), bool(eligible[i]))
            )
        # Equal-group AUC from out-of-fold predictions (single-class group: not estimable).
        per_group = []
        for g in sorted(set(groups.tolist())):
            sel = np.where((groups == g) & eligible)[0]
            ok = len(sel) > 0 and len(np.unique(y[sel])) == 2
            auc = float(roc_auc_score(y[sel], oof_prob[sel])) if ok else float("nan")
            group_auc_rows.append(
                {
                    "seed": int(seed),
                    "group": g,
                    "n_test_cells": len(sel),
                    "auc": auc,
                    "status": "ok" if ok else "not_estimable",
                }
            )
            if ok:
                per_group.append(auc)
        seed_stats.append(
            {
                "auc": float(np.mean([f["auc"] for f in fold_stats])),
                "balanced_accuracy": float(np.mean([f["balanced_accuracy"] for f in fold_stats])),
                "mcc": float(np.mean([f["mcc"] for f in fold_stats])),
                "group_equal_auc": float(np.mean(per_group)) if per_group else float("nan"),
                "n_groups_auc_estimable": len(per_group),
            }
        )

    n_valid = len(seed_stats)
    status = "ok" if n_valid == len(seeds) and n_valid > 0 else ("partial" if n_valid else "failed")
    statistic = None
    if status == "ok":
        statistic = {
            k: float(np.mean([s[k] for s in seed_stats]))
            for k in ("auc", "balanced_accuracy", "mcc", "group_equal_auc")
        }
        statistic["n_groups_auc_estimable"] = int(
            round(float(np.mean([s["n_groups_auc_estimable"] for s in seed_stats])))
        )
    weights = None
    if collect_weights and n_valid:
        if use_hvg:
            safe = np.maximum(coef_cnt, 1)
            wmean = coef_sum / safe
            wsd = np.sqrt(np.maximum(coef_sq / safe - wmean**2, 0.0))
            weights = pd.DataFrame(
                {
                    "gene": fnames,
                    "weight_mean": wmean,
                    "weight_sd": wsd,
                    "n_folds_selected": coef_cnt,
                }
            )
        else:
            arr = np.stack(simple_weights)
            weights = pd.DataFrame(
                {"gene": fnames, "weight_mean": arr.mean(axis=0), "weight_sd": arr.std(axis=0)}
            )
    fold_metrics = pd.DataFrame(metric_rows, columns=list(_FOLD_COLUMNS))
    return CVResult(
        status=status,
        fold_metrics=fold_metrics,
        summary=summarize_fold_metrics(fold_metrics),
        folds=pd.DataFrame(
            member_rows, columns=["seed", "cell_id", "group", "fold", "y", "test_eligible"]
        ),
        group_auc=pd.DataFrame(
            group_auc_rows, columns=["seed", "group", "n_test_cells", "auc", "status"]
        ),
        failures=failures,
        statistic=statistic,
        weights=weights,
        fold_features=fold_features,
        n_groups=len(set(groups.tolist())),
        n_cells=n,
        n_seeds_requested=len(seeds),
        n_seeds_valid=n_valid,
        plan_methods=plan_methods,
        feature_names=fnames,
    )


# ── HR-02: structured permutation null ───────────────────────────────────────
#
# The null answers one question: how large is the group-CV AUC when the viral
# label carries no host-transcriptome information, *given the declared
# exchangeability structure*? The structure is never inferred (``--groups``
# alone does not imply sample-level exchangeability): the caller names the unit.

#: ``cell_within_block``: cell labels are shuffled inside each block (e.g. within
#: a sample or cell type), preserving each block's prevalence. ``group``: whole
#: per-group label vectors are exchanged between groups with the same number of
#: cells (sample-level permutation; within-group label structure is preserved).
NULL_UNITS = ("cell_within_block", "group")
NULL_ASSUMPTIONS = {
    "cell_within_block": (
        "viral labels are exchangeable among cells within the same block under the null; "
        "block prevalence is preserved"
    ),
    "group": (
        "per-group label vectors are exchangeable between groups with the same number of "
        "cells under the null; within-group label structure is preserved"
    ),
}
_NULL_COLUMNS = (
    "virus",
    "permutation_id",
    "seed",
    "permutation_unit",
    "auc",
    "balanced_accuracy",
    "mcc",
    "n_groups",
    "n_cells",
    "status",
    "reason",
)


def empirical_p_value(null_aucs, observed_auc: float, n_requested: int) -> float:
    """``(1 + #(null >= observed)) / (1 + n_requested)`` (ties count as >=).

    Valid only for a complete null: ``len(null_aucs) == n_requested``; the caller
    must not drop failed replicates (see :func:`run_structured_null`).
    """
    null = np.asarray(null_aucs, dtype=float)
    if (
        n_requested < 1
        or len(null) != n_requested
        or not np.isfinite(null).all()
        or not np.isfinite(observed_auc)
    ):
        raise ValueError("empirical p-value requires a complete finite null and observed AUC")
    return float((1 + int((null >= observed_auc).sum())) / (1 + int(n_requested)))


def _null_seed(seed: int, replicate: int) -> int:
    """Deterministic seed of null replicate ``replicate`` (1-based) of base ``seed``."""
    return int(np.random.SeedSequence([int(seed), int(replicate), 7919]).generate_state(1)[0])


def validate_null_spec(unit, has_block: bool, n_permutations: int) -> None:
    """Fail fast on an undeclared or inconsistent null specification.

    ``has_block`` says whether a permutation block (column/values) was supplied.
    ``n_permutations == 0`` (the default, null disabled) accepts anything.
    """
    if int(n_permutations) < 0:
        raise ValueError("permutations must be >= 0")
    if int(n_permutations) == 0:
        return
    if unit not in NULL_UNITS:
        raise ValueError(
            f"permutation_unit must be declared explicitly, one of {NULL_UNITS} (got {unit!r}); "
            "exchangeability is never inferred from the group column"
        )
    if unit == "cell_within_block" and not has_block:
        raise ValueError("permutation_unit='cell_within_block' requires a permutation block")
    if unit == "group" and has_block:
        raise ValueError("a block column only applies to permutation_unit='cell_within_block'")


def permute_labels(y, cell_ids, groups, unit: str, block=None, seed: int = 0):
    """One structured null draw of the label vector.

    Returns ``(y_perm, None)`` or ``(None, reason)`` when the declared structure
    admits no non-trivial exchange (``all_blocks_constant``,
    ``no_exchangeable_groups``). Cells/groups are visited in canonical (sorted
    id) order, so the draw depends on ``seed`` and cell identities, not row order.
    ``block`` is required (one value per cell, no missing) for
    ``cell_within_block``.
    """
    y = np.asarray(y).astype(int)
    if not (len(y) == len(cell_ids) == len(groups)):
        raise ValueError("labels, cell_ids and groups must have equal length")
    if not set(np.unique(y).tolist()) <= {0, 1}:
        raise ValueError("permutation labels must be binary")
    order = _canonical_order(cell_ids)
    rng = np.random.default_rng(seed)
    out = y.copy()
    if unit == "cell_within_block":
        if block is None:
            raise ValueError("cell_within_block requires block values")
        blk = np.asarray(block)
        if len(blk) != len(y) or pd.isna(pd.Series(blk)).any():
            raise ValueError("block values must cover every cell with no missing value")
        blk = blk.astype(str)[order]
        y_c = y[order]
        perm = y_c.copy()
        exchangeable = False
        for b in sorted(set(blk.tolist())):
            pos = np.where(blk == b)[0]
            if len(pos) < 2 or len(set(y_c[pos].tolist())) < 2:
                continue
            exchangeable = True
            perm[pos] = y_c[pos][rng.permutation(len(pos))]
        if not exchangeable:
            return None, "all_blocks_constant"
        out[order] = perm
        return out, None
    if unit == "group":
        g = validate_obs_column(pd.DataFrame({"group": groups}), "group")[order]
        y_c = y[order]
        members = {k: np.where(g == k)[0] for k in sorted(set(g.tolist()))}
        by_size: dict = {}
        for k, pos in members.items():
            by_size.setdefault(len(pos), []).append(k)
        perm = y_c.copy()
        exchangeable = False
        for _, keys in sorted(by_size.items()):
            vectors = {k: tuple(y_c[members[k]].tolist()) for k in keys}
            if len(keys) < 2 or len(set(vectors.values())) < 2:
                continue
            exchangeable = True
            donors = [keys[i] for i in rng.permutation(len(keys))]
            for k, d in zip(keys, donors):
                perm[members[k]] = y_c[members[d]]
        if not exchangeable:
            return None, "no_exchangeable_groups"
        out[order] = perm
        return out, None
    raise ValueError(f"unknown permutation unit {unit!r}")


@dataclass
class NullResult:
    """Outcome of :func:`run_structured_null`.

    ``replicates`` keeps EVERY requested replicate (failed ones with
    ``status="failed"`` and a reason). ``summary['status']`` is ``ok`` only when
    the observed statistic and all requested replicates are valid; otherwise it is
    ``not_estimable`` (or ``disabled`` for ``permutations=0``) and every
    inference field is ``None``.
    """

    replicates: pd.DataFrame
    summary: dict


def run_structured_null(
    procedure,
    y,
    cell_ids,
    groups,
    *,
    unit,
    block=None,
    n_permutations: int,
    seed: int = 0,
    virus: str = "",
    observed: Optional[dict] = None,
) -> NullResult:
    """Empirical null of a full evaluation procedure under a declared exchangeability.

    ``procedure(y) -> {"auc", "balanced_accuracy", "mcc", "n_groups", "n_cells"}``
    is the COMPLETE pipeline (cohort selection, folds, learned features, fit,
    scoring); it is re-run from scratch on every permuted label vector and may
    raise :class:`FoldDesignError` for an invalid design, which is recorded as a
    failed replicate rather than discarded. The p-value is
    :func:`empirical_p_value` over all ``n_permutations`` replicates, so any
    failed replicate makes the null not estimable (dropping invalid arrangements
    would condition on the labels and bias the null). ``observed`` may supply the
    already-computed observed statistic dict to avoid a recomputation.
    """
    validate_null_spec(unit, block is not None, n_permutations)
    n_req = int(n_permutations)
    base = {
        "permutation_unit": unit,
        "auc_statistic": AUC_STATISTIC,
        "n_requested": n_req,
        "n_valid": 0,
        "n_failed": 0,
        "observed_group_auc": None,
        "null_auc_mean": None,
        "null_auc_sd": None,
        "empirical_p_auc": None,
        "observed_minus_null_auc": None,
        "assumption": NULL_ASSUMPTIONS.get(unit, ""),
        "seed": int(seed),
    }
    empty = pd.DataFrame(columns=list(_NULL_COLUMNS))
    if n_req == 0:
        return NullResult(empty, {**base, "status": "disabled", "reason": "permutations=0"})

    if observed is None:
        try:
            observed = procedure(np.asarray(y).astype(int))
        except FoldDesignError as exc:
            rows = []
            for r in range(1, n_req + 1):
                failed_row: dict = dict.fromkeys(_NULL_COLUMNS, float("nan"))
                failed_row.update(
                    virus=virus,
                    permutation_id=r,
                    seed=_null_seed(seed, r),
                    permutation_unit=unit,
                    status="failed",
                    reason=f"observed_invalid: {exc}",
                )
                rows.append(failed_row)
            return NullResult(
                pd.DataFrame(rows, columns=list(_NULL_COLUMNS)),
                {
                    **base,
                    "n_failed": n_req,
                    "status": "not_estimable",
                    "reason": f"observed_invalid: {exc}",
                },
            )
    base["observed_group_auc"] = float(observed["auc"])
    if not np.isfinite(base["observed_group_auc"]):
        observed = None

    rows = []
    for r in range(1, n_req + 1):
        s = _null_seed(seed, r)
        row: dict = dict.fromkeys(_NULL_COLUMNS, float("nan"))
        row.update(virus=virus, permutation_id=r, seed=s, permutation_unit=unit, reason="")
        y_perm, why = permute_labels(y, cell_ids, groups, unit, block, seed=s)
        if observed is None:
            row.update(status="failed", reason="observed_statistic_not_finite")
        elif y_perm is None:
            row.update(status="failed", reason=why)
        else:
            try:
                stat = procedure(y_perm)
                if not all(np.isfinite(stat[k]) for k in ("auc", "balanced_accuracy", "mcc")):
                    raise FoldDesignError("nonfinite_null_statistic", "null metrics must be finite")
                row.update(
                    status="ok",
                    auc=float(stat["auc"]),
                    balanced_accuracy=float(stat["balanced_accuracy"]),
                    mcc=float(stat["mcc"]),
                    n_groups=int(stat["n_groups"]),
                    n_cells=int(stat["n_cells"]),
                )
            except FoldDesignError as exc:
                row.update(status="failed", reason=str(exc))
        rows.append(row)
    reps = pd.DataFrame(rows, columns=list(_NULL_COLUMNS))
    ok = reps[reps["status"] == "ok"]
    base["n_valid"] = len(ok)
    base["n_failed"] = n_req - len(ok)
    if base["n_failed"]:
        reasons = sorted(set(reps.loc[reps["status"] == "failed", "reason"].astype(str)))
        return NullResult(
            reps,
            {
                **base,
                "status": "not_estimable",
                "reason": "failed replicates: " + " | ".join(reasons),
            },
        )
    null = ok["auc"].to_numpy(dtype=float)
    base.update(
        null_auc_mean=float(null.mean()),
        null_auc_sd=float(null.std()),
        empirical_p_auc=empirical_p_value(null, base["observed_group_auc"], n_req),
        observed_minus_null_auc=float(base["observed_group_auc"] - null.mean()),
    )
    return NullResult(reps, {**base, "status": "ok", "reason": ""})


def group_cv_permutation_null(
    X,
    y,
    groups,
    depth,
    cell_ids,
    *,
    unit,
    n_permutations: int,
    seed: int = 0,
    block=None,
    virus: str = "",
    observed: Optional[dict] = None,
    **cv_kwargs,
) -> NullResult:
    """Structured null of the group-CV AUC with the full procedure refit per replicate.

    Every replicate calls :func:`evaluate_group_cv` on the permuted labels with the
    same ``cv_kwargs`` (seeds, folds, HVG, depth filtering/matching, ...), so
    label-dependent cohort selection, fold construction/validation and learned
    features (in-fold HVG) are recomputed from the permuted labels; nothing
    label-dependent is carried over from the observed run. A replicate whose design
    is invalid (any seed) is a failed replicate. ``cv_kwargs`` must not request
    baselines/weights (the null statistic is the expression AUC only).
    """
    for forbidden in ("baselines", "collect_weights", "panel_iters", "cell_type"):
        cv_kwargs.pop(forbidden, None)

    def procedure(y_in):
        res = evaluate_group_cv(X, y_in, groups, depth, cell_ids, **cv_kwargs)
        if res.status != "ok" or res.statistic is None:
            reasons = "; ".join(sorted({f["reason"] for f in res.failures})) or res.status
            raise FoldDesignError("invalid_null_replicate", reasons)
        return {**res.statistic, "n_groups": res.n_groups, "n_cells": res.n_cells}

    return run_structured_null(
        procedure,
        y,
        cell_ids,
        groups,
        unit=unit,
        block=block,
        n_permutations=n_permutations,
        seed=seed,
        virus=virus,
        observed=observed,
    )


def _summary_auc(summary: pd.DataFrame, model: str) -> Optional[dict]:
    """``{"mean","sd"}`` of a model's fold AUC from a CV summary, or ``None``."""
    sub = summary[(summary["model"] == model) & (summary["metric"] == "auc")]
    if sub.empty or not int(sub.iloc[0]["n_valid_folds"]):
        return None
    return {"mean": float(sub.iloc[0]["mean"]), "sd": float(sub.iloc[0]["sd"])}


def _write_cv_tables(out_dir, virus: str, cv: CVResult) -> None:
    """Write the per-virus fold membership, per-fold metric and summary TSVs."""
    base = Path(out_dir) / _safe_name(virus)
    for suffix, df in (
        ("_cv_folds.tsv", cv.folds),
        ("_cv_metrics.tsv", cv.fold_metrics),
        ("_cv_summary.tsv", cv.summary),
        ("_cv_group_auc.tsv", cv.group_auc),
    ):
        df.assign(virus=virus)[["virus", *df.columns]].to_csv(
            f"{base}{suffix}", sep="\t", index=False
        )


def _legacy_metrics(summary: pd.DataFrame, model: str) -> dict:
    """``{metric: {"mean","sd"}}`` in the cell-mode naming for one model's summary."""
    out = {}
    for _, r in summary[summary["model"] == model].iterrows():
        if r["n_valid_folds"]:
            key = "balanced_acc" if r["metric"] == "balanced_accuracy" else r["metric"]
            out[key] = {"mean": float(r["mean"]), "sd": float(r["sd"])}
    return out


@dataclass
class _RunOpts:
    """Analysis options shared by every cohort (pooled or per cell-type stratum)."""

    use_hvg: bool
    seeds: list
    n_stab_iter: int
    stab_min_prob: float
    top_n_genes: int
    detection_threshold: int
    do_enrichment: bool
    enrichment_db: str
    label: str
    depth_match: bool
    control_mito: bool
    annotate_symbols: bool
    differential: bool
    cv_mode: str
    groups_column: Optional[str]
    cv_folds: int
    panel_in_fold: bool
    permutations: int
    permutation_unit: Optional[str]
    min_negative_cells: int
    min_groups: int


def _run_cohort(
    virus_adata,
    host_adata,
    viral_vars,
    groups_arr,
    ct_full,
    block_arr,
    out_dir,
    stratum,
    o,
    prov_viruses,
    baseline_models,
):
    """Run the per-virus analysis on ONE aligned cohort (all cells, or one cell-type stratum).

    ``host_adata`` must hold raw counts (it is normalised in place). ``o`` is the
    :class:`_RunOpts` bundle. Returns ``(metric_rows, status_rows)``; per-virus
    outputs are written into ``out_dir``. ``stratum`` is ``"ALL"`` for the pooled
    analysis, otherwise the cell-type name (added to every metrics row).
    """
    (
        use_hvg, seeds, n_stab_iter, stab_min_prob, top_n_genes, detection_threshold,
        do_enrichment, enrichment_db, label, depth_match, control_mito, annotate_symbols,
        differential, cv_mode, groups_column, cv_folds, panel_in_fold, permutations,
        permutation_unit,
    ) = (
        o.use_hvg, o.seeds, o.n_stab_iter, o.stab_min_prob, o.top_n_genes, o.detection_threshold,
        o.do_enrichment, o.enrichment_db, o.label, o.depth_match, o.control_mito, o.annotate_symbols,
        o.differential, o.cv_mode, o.groups_column, o.cv_folds, o.panel_in_fold, o.permutations,
        o.permutation_unit,
    )  # fmt: skip
    status_rows: list = []
    cnt: dict = {}

    def _st(virus, status, reason=""):
        return _status_row(virus, status, reason, stratum=stratum, **cnt)

    # Compute %mito from RAW counts BEFORE normalization (so the QC covariate for
    # the E-values is a real fraction, not a normalized artifact). Keep the raw
    # per-gene MT counts so a mitochondrial gene can be left out of its own %mito
    # covariate (avoids the MT-ND4L self-suppression circularity; see go_enrichment.py).
    mt_total_counts = raw_total = pct_mito = None
    mt_self_counts: dict = {}
    if control_mito:
        mt_mask = _mt_gene_mask(host_adata.var_names)
        raw_counts = host_adata.layers.get("counts", host_adata.X)
        raw_total = np.asarray(raw_counts.sum(axis=1)).ravel().astype(float)
        if mt_mask.any():
            mt_sub = raw_counts[:, mt_mask]
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
            ct_v = _sub(ct_full, midx)
            top_depth_frac = 1.0
            log.info("[%s] Depth-matched cohort: %d cells (from %d).", virus, len(midx), len(depth))
        else:
            X_full_v, X_stab_v = X_full, X_stab
            vp_v, depth_v = virus_presence_full, depth
            pct_mito_v, mt_total_v, raw_total_v = pct_mito, mt_total_counts, raw_total
            mt_self_v = mt_self_counts
            ct_v = ct_full
            top_depth_frac = TOP_DEPTH_FRAC

        # Group mode gates on the full label (its folds build their own cohorts);
        # cell mode gates on the (optionally depth-matched) analysis cohort.
        n_pos = int(virus_presence_full.sum()) if cv_mode == "group" else int(vp_v.sum())
        n_neg = (
            len(virus_presence_full) - int(virus_presence_full.sum())
            if cv_mode == "group"
            else len(vp_v) - int(vp_v.sum())
        )
        cnt.update(
            n_cells=n_pos + n_neg,
            n_positive=n_pos,
            n_negative=n_neg,
            n_groups=len(set(groups_arr.tolist())) if groups_arr is not None else "",
        )
        if n_pos < MIN_VIRUS_CELLS:
            log.info(
                "[%s] Skipping: fewer than %d positive cells in the analysis cohort.",
                virus,
                MIN_VIRUS_CELLS,
            )
            status_rows.append(_st(virus, "skipped_min_positive", f"{n_pos} < {MIN_VIRUS_CELLS}"))
            continue

        if n_neg < o.min_negative_cells:
            status_rows.append(
                _st(virus, "skipped_min_negative", f"{n_neg} < {o.min_negative_cells}")
            )
            continue
        if groups_arr is not None and len(set(groups_arr.tolist())) < o.min_groups:
            status_rows.append(
                _st(
                    virus,
                    "group_cv_failed"
                    if len(set(groups_arr.tolist())) < 2
                    else "skipped_min_groups",
                    f"insufficient_groups: {len(set(groups_arr.tolist()))} < {o.min_groups}",
                )
            )
            continue

        cv_extra: dict = {"cv_mode": cv_mode}
        if cv_mode == "group":
            assert groups_arr is not None
            cv = evaluate_group_cv(
                X_full,
                virus_presence_full,
                groups_arr,
                depth,
                host_adata.obs_names.to_numpy(),
                seeds=seeds,
                n_folds=cv_folds,
                feature_names=all_gene_names,
                use_hvg=use_hvg,
                depth_match=depth_match,
                collect_weights=True,
                cell_type=ct_full,
                baselines=baseline_models,
                panel_iters=n_stab_iter if panel_in_fold else 0,
                stab_min_prob=stab_min_prob,
            )
            _write_cv_tables(out_dir, virus, cv)
            prov_viruses[virus] = cv.provenance()
            if cv.status == "failed":
                reasons = "; ".join(f["message"] for f in cv.failures)
                log.warning("[%s] Group CV invalid for every seed: %s", virus, reasons)
                status_rows.append(_st(virus, "group_cv_failed", reasons))
                continue
            weights_df = cv.weights
            metrics = _legacy_metrics(cv.summary, "expression")
            cv_extra.update(
                group_column=groups_column,
                cv_status=cv.status,
                cv_n_groups=cv.n_groups,
                cv_n_valid_seeds=cv.n_seeds_valid,
                cv_n_folds=int(
                    cv.fold_metrics.loc[cv.fold_metrics["status"] == "ok", "fold"].nunique()
                ),
                cv_group_equal_auc=(cv.statistic or {}).get("group_equal_auc"),
            )
            if permutations > 0:
                null = group_cv_permutation_null(
                    X_full,
                    virus_presence_full,
                    groups_arr,
                    depth,
                    host_adata.obs_names.to_numpy(),
                    unit=permutation_unit,
                    n_permutations=permutations,
                    seed=seeds[0] if seeds else 0,
                    block=block_arr,
                    virus=virus,
                    observed=cv.statistic,
                    seeds=seeds,
                    n_folds=cv_folds,
                    feature_names=all_gene_names,
                    use_hvg=use_hvg,
                    depth_match=depth_match,
                )
                null.replicates.to_csv(
                    Path(out_dir) / f"{_safe_name(virus)}_permutation_metrics.tsv",
                    sep="\t",
                    index=False,
                )
                prov_viruses[virus]["null"] = null.summary
                cv_extra.update(
                    {
                        k: null.summary[k]
                        for k in (
                            "observed_group_auc",
                            "null_auc_mean",
                            "null_auc_sd",
                            "empirical_p_auc",
                            "observed_minus_null_auc",
                        )
                        if null.summary[k] is not None
                    },
                    permutation_status=null.summary["status"],
                    permutation_unit=permutation_unit,
                    permutation_n_requested=null.summary["n_requested"],
                    permutation_n_valid=null.summary["n_valid"],
                    permutation_n_failed=null.summary["n_failed"],
                )
                if null.summary["status"] != "ok":
                    log.warning(
                        "[%s] Permutation null not estimable: %s", virus, null.summary["reason"]
                    )
        else:
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
            status_rows.append(_st(virus, "no_valid_models", "balanced split failed"))
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
        baseline_means: dict = {}
        cv_panel_auc = None
        if cv_mode == "group":
            # Baselines were scored on the group folds themselves (fold-local fits).
            depth_alone = _summary_auc(cv.summary, "depth_only")
            for m in (*BASELINE_MODELS, "panel_depth"):
                s = _summary_auc(cv.summary, m)
                if m == "panel_depth":
                    cv_panel_auc = s
                elif s is not None:
                    baseline_means[m] = s["mean"]
        else:
            depth_alone = _depth_alone_auc(
                np.where(vp_v)[0],
                np.where(~vp_v)[0],
                depth_v,
                seeds,
                top_depth_frac=top_depth_frac,
            )
            if ct_v is not None:
                cell_base = _cell_mode_baselines(
                    np.where(vp_v)[0],
                    np.where(~vp_v)[0],
                    depth_v,
                    ct_v,
                    seeds,
                    models=("cell_type_only", "depth_plus_cell_type"),
                    top_depth_frac=top_depth_frac,
                )
                baseline_means.update({m: v["mean"] for m, v in cell_base.items() if v})
        if depth_alone is not None:
            baseline_means["depth_only"] = depth_alone["mean"]
        stable_genes_list = stab_df.loc[stab_df["stable"], "gene"].tolist()
        name_to_col = {g: i for i, g in enumerate(stab_feature_names)}
        ev_cols = [name_to_col[g] for g in stable_genes_list if g in name_to_col]
        n_evalue_ge2 = 0
        depth_adj_auc = cv_panel_auc  # None in cell mode (whole-cohort panel computed below)
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
            # The panel is chosen from the WHOLE cohort, so its AUC is descriptive and
            # not a held-out estimate; it is therefore only reported in cell mode.
            if cv_mode != "group":
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
            "stratum": stratum,
            "n_positive": n_pos,
            "label": label,
            "depth_matched": depth_match,
            "mito_controlled": bool(control_mito and pct_mito_v is not None),
            **cv_extra,
        }
        status_rows.append(_st(virus, "ok", cv_extra.get("cv_status", "")))
        if baseline_means:
            expr_auc = ((metrics or {}).get("auc") or {}).get("mean")
            comp = baseline_comparison(expr_auc, baseline_means)
            row.update({k: v for k, v in comp.items() if v is not None})
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

    return all_metrics, status_rows


def run_hostresponse(
    virus_h5ad: str,
    host_h5ad: str,
    viral_accessions_file: str,
    out_dir: str,
    use_hvg: bool = True,
    seeds: Optional[list] = None,
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
    cv_mode: str = "cell",
    groups_column: Optional[str] = None,
    cv_folds: int = DEFAULT_CV_FOLDS,
    cell_type_column: Optional[str] = None,
    panel_in_fold: bool = False,
    permutations: int = 0,
    permutation_unit: Optional[str] = None,
    permutation_block: Optional[str] = None,
    cell_types: Optional[list] = None,
    min_negative_cells: int = MIN_VIRUS_CELLS,
    min_groups: int = 2,
) -> dict:
    """Main entry point: run per-virus logistic regression host-response analysis.

    Returns a JSON-serialisable provenance dict (see :func:`_sha256_file` and
    ``HOSTRESPONSE_SCHEMA_VERSION``): input hashes, raw-depth source, cv mode,
    group column, seeds, fold digests and per-virus status. It is written to
    disk as ``hostresponse_manifest.json`` with hashes of every owned artifact.

    ``cell_type_column`` (host or virus obs) adds the ``cell_type_only`` and
    ``depth_plus_cell_type`` baselines next to the always-on depth baseline, on the
    same splits/folds as the expression model (training-only encoding, see
    :func:`_baseline_design`); the metrics row then carries ``expression_auc``,
    ``depth_auc``, ``cell_type_auc``, ``depth_cell_type_auc`` and
    ``expression_minus_best_baseline_auc``. In group mode ``panel_in_fold=True``
    additionally reports stable-panel + depth AUC with the panel selected inside
    each training fold (costly: ``n_stab_iter`` fits per fold and seed); the
    whole-cohort stable panel is never scored as a held-out predictor there.

    ``permutations > 0`` (group mode only; default 0 = off) runs a structured null of
    the group-CV AUC (:func:`group_cv_permutation_null`): the exchangeability unit
    ``permutation_unit`` (``cell_within_block`` with ``permutation_block``, or
    ``group``) must be declared explicitly. The full procedure is refit for every
    replicate; ``<virus>_permutation_metrics.tsv`` keeps failed replicates and an
    invalid null is reported ``not_estimable`` (no p-value), never repaired.

    ``cv_mode="cell"`` (default) keeps the repeated balanced cell-holdout
    evaluation unchanged. ``cv_mode="group"`` holds out whole groups
    (``groups_column`` in the host or virus obs) via :func:`evaluate_group_cv`; a
    missing/blank group column raises ``ValueError`` and never falls back to cell
    CV. Group-mode and cell-mode numbers are separate estimands and are labelled
    by ``cv_mode`` in every output row. Whole-cohort stability selection, E-values
    and the differential table remain descriptive in both modes.

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
    Selected ``cell_types`` run separately in collision-checked ``cell_type_*``
    directories; no pooled fit is run when strata are requested. Each stratum
    validates its own class/group support. Differential p-values stay cell-wise.
    """
    import anndata as ad  # noqa: PLC0415

    if label not in LABEL_CHOICES:
        raise ValueError(f"label must be one of {LABEL_CHOICES}, got {label!r}")
    if cv_mode not in CV_MODES:
        raise ValueError(f"cv_mode must be one of {CV_MODES}, got {cv_mode!r}")
    if cv_mode == "group" and not groups_column:
        raise ValueError("cv_mode='group' requires groups_column (no silent fallback to cell CV)")
    if seeds is None:
        seeds = DEFAULT_SEEDS
    seeds = list(seeds)
    if not seeds or len(set(seeds)) != len(seeds) or any(int(seed) < 0 for seed in seeds):
        raise ValueError("seeds must be nonempty, unique and nonnegative")
    if cv_folds < 2 or min_groups < 2 or min_negative_cells < 1:
        raise ValueError("cv_folds/min_groups must be >= 2 and min_negative_cells >= 1")
    if panel_in_fold and (cv_mode != "group" or n_stab_iter < 1):
        raise ValueError("panel_in_fold requires group CV and n_stab_iter >= 1")
    if cell_types is not None and (not cell_types or not cell_type_column):
        raise ValueError("cell_types requires a nonempty list and cell_type_column")
    validate_null_spec(permutation_unit, bool(permutation_block), permutations)
    if permutations > 0 and cv_mode != "group":
        raise ValueError("permutations > 0 requires cv_mode='group' (the null is for group CV)")

    provenance: dict = {
        "hostresponse_schema_version": HOSTRESPONSE_SCHEMA_VERSION,
        "sklearn_version": _sklearn.__version__,
        "software": software_identity(),
        "cv_mode": cv_mode,
        "group_column": groups_column if cv_mode == "group" else None,
        "n_folds": int(cv_folds) if cv_mode == "group" else None,
        "seeds": [int(s) for s in seeds],
        "auc_statistic": AUC_STATISTIC,
        "label": label,
        "detection_threshold": int(detection_threshold),
        "depth_matching": bool(depth_match),
        "feature_mode": "hvg_in_training_fold" if use_hvg else "all_genes",
        "inputs": {
            "virus_h5ad": {"path": str(virus_h5ad), "sha256": _sha256_file(virus_h5ad)},
            "host_h5ad": {"path": str(host_h5ad), "sha256": _sha256_file(host_h5ad)},
            "viral_accessions": {
                "path": str(viral_accessions_file),
                "sha256": _sha256_file(viral_accessions_file),
            },
        },
        "raw_depth_source": None,
        "n_cells_aligned": 0,
        "n_groups": None,
        "permutation_count": int(permutations),
        "permutation_unit": permutation_unit if permutations > 0 else None,
        "permutation_block": permutation_block if permutations > 0 else None,
        "viruses": {},
        "strata": {},
        "cell_type_filter": list(cell_types) if cell_types is not None else None,
        "min_negative_cells": int(min_negative_cells),
        "min_groups": int(min_groups),
        "descriptive_inference_unit": "cell",
        "legacy_cell_panel_auc": "descriptive_whole_cohort_selection"
        if cv_mode == "cell"
        else None,
        "status": "started",
    }
    status_rows: list = []

    log.info("Loading virus h5ad: %s", virus_h5ad)
    virus_adata = ad.read_h5ad(virus_h5ad)
    log.info("Loading host h5ad: %s", host_h5ad)
    host_adata_full = ad.read_h5ad(host_h5ad)

    # Hard input validation happens BEFORE any existing output is cleared, so an
    # invalid request never destroys the previous run's results.
    for name, obj in (("virus", virus_adata), ("host", host_adata_full)):
        if not obj.obs_names.is_unique:
            raise ValueError(f"{name} h5ad has duplicate barcodes; cell alignment is ambiguous")
    _, provenance["raw_depth_source"] = _resolve_raw_depth(host_adata_full)
    if cv_mode == "group" and provenance["raw_depth_source"] == "X_log1p_not_raw":
        raise ValueError("group CV requires raw depth: supply layers['counts'] with log1p X")

    # Filter virus matrix to confirmed viral gene IDs (identity table or analysis.txt).
    # If the pipeline used a combined host+viral reference, the h5ad contains
    # host genes too; restricting here prevents training models that predict
    # host gene expression from other host gene expression.
    viral_accessions = _load_viral_accessions(viral_accessions_file)
    viral_vars = [v for v in virus_adata.var_names if v in viral_accessions]

    # Align cells: intersect barcodes between virus and host matrices.
    shared_barcodes = sorted(set(virus_adata.obs_names) & set(host_adata_full.obs_names))

    groups_arr = None
    if cv_mode == "group":
        assert groups_column is not None
        obs_src = host_adata_full.obs if groups_column in host_adata_full.obs else virus_adata.obs
        groups_arr = validate_obs_column(obs_src.loc[shared_barcodes], groups_column)
    ct_full = None
    if cell_type_column:
        ct_src = host_adata_full.obs if cell_type_column in host_adata_full.obs else virus_adata.obs
        ct_full = validate_obs_column(
            ct_src.loc[shared_barcodes], cell_type_column, what="cell-type"
        )
    block_arr = None
    if permutations > 0 and permutation_block:
        blk_src = (
            host_adata_full.obs if permutation_block in host_adata_full.obs else virus_adata.obs
        )
        block_arr = validate_obs_column(
            blk_src.loc[shared_barcodes], permutation_block, what="permutation block"
        )
    provenance["n_cells_aligned"] = len(shared_barcodes)
    provenance["cell_type_column"] = cell_type_column
    provenance["baseline_models"] = [
        "depth_only",
        *(["cell_type_only", "depth_plus_cell_type"] if ct_full is not None else []),
        *(["panel_depth"] if panel_in_fold and cv_mode == "group" else []),
    ]

    strata_names: dict = {}
    if cell_types is not None:
        assert ct_full is not None
        available = set(ct_full.tolist())
        for ct in cell_types:
            if not isinstance(ct, str) or not ct.strip() or ct not in available:
                raise ValueError(f"unknown/blank cell type {ct!r}; available: {sorted(available)}")
            safe = _safe_name(ct)
            if safe in strata_names.values():
                raise ValueError(f"cell-type output name collision for {ct!r}")
            strata_names[ct] = safe
    virus_names = [_safe_name(str(v)) for v in viral_vars]
    if len(set(virus_names)) != len(virus_names):
        raise ValueError("viral output names collide after path sanitization")
    opts = _RunOpts(
        use_hvg,
        seeds,
        n_stab_iter,
        stab_min_prob,
        top_n_genes,
        detection_threshold,
        do_enrichment,
        enrichment_db,
        label,
        depth_match,
        control_mito,
        annotate_symbols,
        differential,
        cv_mode,
        groups_column,
        cv_folds,
        panel_in_fold,
        permutations,
        permutation_unit,
        min_negative_cells,
        min_groups,
    )

    os.makedirs(out_dir, exist_ok=True)
    clear_stale_virus_outputs(out_dir)
    (Path(out_dir) / "hostresponse_metrics.csv").unlink(missing_ok=True)
    log.info("hostresponse output directory: %s", out_dir)

    if not viral_vars:
        log.error(
            "No viral accessions from analysis.txt (%s) match var_names in the virus h5ad "
            "(%s). Verify that the h5ad was produced from the same viralscan run.",
            viral_accessions_file,
            virus_h5ad,
        )
        provenance["status"] = "no_viral_accessions_matched"
        status_rows.append(_status_row("", provenance["status"]))
        _write_manifest(out_dir, provenance, status_rows)
        return provenance
    log.info("Retained %d viral gene(s) from analysis.txt for hostresponse.", len(viral_vars))
    virus_adata = virus_adata[:, viral_vars].copy()

    if len(shared_barcodes) < MIN_VIRUS_CELLS * 2:
        log.error(
            "Only %d shared barcodes between virus h5ad and host h5ad (need >= %d). "
            "Ensure both h5ad files use identical barcode strings.",
            len(shared_barcodes),
            MIN_VIRUS_CELLS * 2,
        )
        provenance["status"] = "too_few_shared_barcodes"
        status_rows.append(_status_row("", provenance["status"], n_cells=len(shared_barcodes)))
        _write_manifest(out_dir, provenance, status_rows)
        return provenance
    log.info("Aligned on %d shared barcodes.", len(shared_barcodes))
    virus_adata = virus_adata[shared_barcodes].copy()
    host_adata = host_adata_full[shared_barcodes].copy()
    provenance["n_cells_aligned"] = len(shared_barcodes)
    if groups_arr is not None:
        provenance["n_groups"] = len(set(groups_arr.tolist()))

    all_metrics = []
    cohorts = list(strata_names) if cell_types is not None else ["ALL"]
    for stratum in cohorts:
        idx = (
            np.where(ct_full == stratum)[0]
            if cell_types is not None
            else np.arange(len(shared_barcodes))
        )
        cohort_out = (
            Path(out_dir) / f"cell_type_{strata_names[stratum]}"
            if cell_types is not None
            else Path(out_dir)
        )
        cohort_out.mkdir(exist_ok=True)
        virus_prov: dict = {}
        cohort_metrics, cohort_status = _run_cohort(
            virus_adata[idx].copy(),
            host_adata[idx].copy(),
            viral_vars,
            None if groups_arr is None else groups_arr[idx],
            None if ct_full is None else ct_full[idx],
            None if block_arr is None else block_arr[idx],
            cohort_out,
            stratum,
            opts,
            virus_prov,
            provenance["baseline_models"],
        )
        all_metrics.extend(cohort_metrics)
        status_rows.extend(cohort_status)
        provenance["strata"][stratum] = {
            "output_directory": str(cohort_out.relative_to(out_dir)),
            "n_cells": len(idx),
            "viruses": virus_prov,
            "status": "ok" if cohort_metrics else "no_eligible_virus",
        }
        if cell_types is None:
            provenance["viruses"] = virus_prov

    metrics_csv = Path(out_dir) / "hostresponse_metrics.csv"
    if all_metrics:
        pd.DataFrame(all_metrics).to_csv(metrics_csv, index=False)
        log.info("Summary metrics written to %s", metrics_csv)
    else:
        # A no-eligible run must not leave the previous run's summary behind.
        metrics_csv.unlink(missing_ok=True)
        log.warning(
            "No virus met the minimum positive-cell threshold (%d). "
            "Consider lowering --detection-threshold or using a sample with higher viral load.",
            MIN_VIRUS_CELLS,
        )
    provenance["status"] = "ok" if all_metrics else "no_eligible_virus"
    provenance["viruses_status"] = {r["virus"]: r["status"] for r in status_rows}
    _write_manifest(out_dir, provenance, status_rows)
    return provenance


def _write_manifest(out_dir, provenance: dict, status_rows: list) -> None:
    """Write artifact-local provenance, including all owned output hashes/statuses."""
    directory = Path(out_dir)
    pd.DataFrame(status_rows, columns=list(_STATUS_COLUMNS)).to_csv(
        directory / STATUS_FILENAME, sep="\t", index=False
    )
    provenance["statuses"] = status_rows
    provenance["artifacts"] = {
        str(path.relative_to(directory)): _sha256_file(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
        and not path.is_symlink()
        and (
            path.name.endswith(_OWNED_OUTPUT_SUFFIXES)
            or "_enrichment_" in path.name
            or path.name in (STATUS_FILENAME, "hostresponse_metrics.csv")
        )
        and (path.parent == directory or path.parent.name.startswith("cell_type_"))
    }

    # Strict JSON records unavailable numbers as null, including invalid null draws.
    def finite(value):
        if isinstance(value, dict):
            return {str(k): finite(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [finite(v) for v in value]
        if isinstance(value, float) and not np.isfinite(value):
            return None
        return value

    cleaned = finite(provenance)
    provenance.clear()
    provenance.update(cleaned)
    target = directory / MANIFEST_FILENAME
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(provenance, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(target)


# ── Snakemake entry point ─────────────────────────────────────────────────────

if "snakemake" in globals():
    cfg = RunConfig.from_yaml(snakemake.params.configfile)  # noqa: F821
    kb = KbCountOutputs.from_config_output(cfg.output)
    _virus_h5ad = str(kb.current_adata(multimapping=cfg.multimapping))
    _identity_file = identity_path(cfg.output)
    _viral_acc_file = (
        str(_identity_file) if _identity_file.is_file() else f"{cfg.output}log/analysis.txt"
    )
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
        cv_mode=cfg.hostresponse_cv,
        groups_column=cfg.hostresponse_groups,
        cv_folds=cfg.hostresponse_cv_folds,
        cell_type_column=cfg.hostresponse_cell_type_column,
        cell_types=cfg.hostresponse_cell_types,
        min_negative_cells=cfg.hostresponse_min_negative_cells,
        min_groups=cfg.hostresponse_min_groups,
        panel_in_fold=cfg.hostresponse_panel_in_fold,
        permutations=cfg.hostresponse_permutations,
        permutation_unit=cfg.hostresponse_permutation_unit,
        permutation_block=cfg.hostresponse_permutation_block,
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
    p.add_argument("--cv", choices=CV_MODES, default="cell")
    p.add_argument("--groups", default=None)
    p.add_argument("--cv-folds", type=int, default=DEFAULT_CV_FOLDS)
    p.add_argument("--cell-type-column", default=None)
    p.add_argument("--cell-types", type=lambda value: value.split(","), default=None)
    p.add_argument("--min-negative-cells", type=int, default=MIN_VIRUS_CELLS)
    p.add_argument("--min-groups", type=int, default=2)
    p.add_argument("--panel-in-fold", action="store_true")
    p.add_argument("--permutations", type=int, default=0)
    p.add_argument("--permutation-unit", choices=NULL_UNITS, default=None)
    p.add_argument("--permutation-block", default=None)
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
        cv_mode=args.cv,
        groups_column=args.groups,
        cv_folds=args.cv_folds,
        cell_type_column=args.cell_type_column,
        cell_types=args.cell_types,
        min_negative_cells=args.min_negative_cells,
        min_groups=args.min_groups,
        panel_in_fold=args.panel_in_fold,
        permutations=args.permutations,
        permutation_unit=args.permutation_unit,
        permutation_block=args.permutation_block,
    )
