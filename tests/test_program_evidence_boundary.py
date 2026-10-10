"""A missing unique layer is unavailable evidence, not a programme-negative matrix."""

import anndata as ad
import pytest
from scipy import sparse

from viralscan.gene_programs import EVIDENCE_LAYER
from viralscan.scripts.gene_programs import _marker_matrix


def test_missing_unique_layer_rejected_even_with_no_resolved_markers():
    data = ad.AnnData(sparse.csr_matrix([[3.0]]))
    with pytest.raises(ValueError, match="counts_unique_viral.*Rebuild"):
        _marker_matrix(data, [])


def test_present_unique_layer_preserves_empty_marker_matrix():
    data = ad.AnnData(sparse.csr_matrix([[3.0]]))
    data.layers[EVIDENCE_LAYER] = sparse.csr_matrix([[2.0]])
    unique, selected = _marker_matrix(data, [])
    assert unique.shape == selected.shape == (1, 0)
