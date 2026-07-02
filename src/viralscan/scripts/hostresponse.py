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

Outputs (per virus, under <output>/hostresponse/):
  - <virus>_gene_weights.csv   — mean/SD of L2 coefficients (per-fold-HVG: genes
    selected in >=1 fold, with n_folds_selected)
  - <virus>_stability.csv      — per-gene stability probability + merged weights
  - hostresponse_metrics.csv   — sensitivity/specificity/balanced-acc/AUC/MCC summary
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
    if sp.issparse(X):
        raw_depth = np.asarray(X.sum(axis=1)).flatten()
    else:
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


def _balanced_split(pos_idx, neg_idx, depth, X, seed: int):
    """Depth-filtered balanced train/test split.

    Filters each class to the top TOP_DEPTH_FRAC by RAW sequencing depth
    (stored in obs["_raw_depth"] before normalization) to reduce false-negative
    viral-absence labels in low-coverage cells.

    Returns (X_train, y_train, X_test_pos, X_test_neg) or None if too few cells.
    """
    rng = np.random.default_rng(seed)

    def top_depth(idx):
        d = depth[idx]
        cutoff = np.percentile(d, (1 - TOP_DEPTH_FRAC) * 100)
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


def _run_l2_regression(X, virus_presence, depth, seeds, feature_names, use_hvg=False):
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
        split = _balanced_split(pos_idx, neg_idx, depth, X, seed)
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
) -> None:
    """Main entry point: run per-virus logistic regression host-response analysis."""
    import anndata as ad  # noqa: PLC0415

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

    all_metrics = []

    for virus in viral_vars:
        counts = virus_adata[:, virus].X
        counts = counts.toarray().flatten() if sp.issparse(counts) else np.asarray(counts).flatten()

        virus_presence = counts >= detection_threshold
        n_pos = int(virus_presence.sum())
        log.info(
            "[%s] %d / %d cells positive (detection_threshold=%d).",
            virus,
            n_pos,
            len(virus_presence),
            detection_threshold,
        )

        if n_pos < MIN_VIRUS_CELLS:
            log.info("[%s] Skipping: fewer than %d positive cells.", virus, MIN_VIRUS_CELLS)
            continue

        # ── L2 multi-seed regression (per-fold HVG when use_hvg) ───────────
        weights_df, metrics = _run_l2_regression(
            X_full, virus_presence, depth, seeds, all_gene_names, use_hvg=use_hvg
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
        weights_df.to_csv(weights_csv, index=False)
        log.info("[%s] Gene weights written to %s", virus, weights_csv)

        # ── Randomized Lasso stability selection (descriptive gene ranking) ──
        stab_probs = _run_stability_selection(
            X_stab, virus_presence, n_stab_iter, seed=seeds[0] if seeds else 42
        )
        stab_df = pd.DataFrame({"gene": stab_feature_names, "stab_prob": stab_probs})
        stab_df = stab_df.merge(
            weights_df[["gene", "weight_mean", "weight_sd"]], on="gene", how="left"
        )
        stab_df["stable"] = stab_df["stab_prob"] >= stab_min_prob
        stab_csv = Path(out_dir) / f"{_safe_name(virus)}_stability.csv"
        stab_df.to_csv(stab_csv, index=False)
        log.info("[%s] Stability probabilities written to %s", virus, stab_csv)

        # Collect summary metrics.
        row: dict = {"virus": virus, "n_positive": n_pos}
        for metric, vals in (metrics or {}).items():
            row[f"{metric}_mean"] = vals["mean"]
            row[f"{metric}_sd"] = vals["sd"]
        all_metrics.append(row)

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
    )
