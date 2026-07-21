"""D1: duplicate var_names are handled, not crashed on cryptically."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from viralscan.utils import matrix_for_genes


def _adata(var_names):
    a = ad.AnnData(sp.csr_matrix(np.ones((3, len(var_names)), dtype=float)))
    a.var_names = var_names
    return a


def test_matrix_for_genes_rejects_non_unique_var_names_with_clear_error():
    adata = _adata(["EBV_1", "EBV_1", "HOST"])  # duplicate accession
    with pytest.raises(ValueError, match="duplicate gene IDs"):
        matrix_for_genes(adata, adata.X, ["EBV_1"])


def test_matrix_for_genes_works_after_make_unique():
    adata = _adata(["EBV_1", "EBV_1", "HOST"])
    adata.var_names_make_unique()  # what Detection.preprocessing does up front
    out = matrix_for_genes(adata, adata.X, list(adata.var_names[:2]))
    assert out.shape == (3, 2)
