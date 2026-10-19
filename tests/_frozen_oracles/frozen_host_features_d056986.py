"""Verbatim _select_features from pre-B21 commit d056986; never regenerate."""

# ruff: noqa
import logging
import numpy as np
import scipy.sparse as sp
import scanpy as sc

log = logging.getLogger("viralscan.scripts.hostresponse")


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
