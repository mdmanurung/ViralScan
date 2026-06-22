"""Tests for _swap_multimap_layer and the rerun-multimap CLI subcommand."""

from __future__ import annotations

import argparse
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def _make_multimap_h5ad(path: Path) -> None:
    """Write a minimal multimap h5ad with all three pre-stored non-EM layers."""
    n_cells, n_genes = 4, 6
    base = sparse.csr_matrix(np.zeros((n_cells, n_genes), dtype=np.float32))

    # Distinct sentinel data so we can tell layers apart after the swap.
    equal = sparse.csr_matrix(np.ones((n_cells, n_genes), dtype=np.float32) * 1.0)
    hc = sparse.csr_matrix(np.ones((n_cells, n_genes), dtype=np.float32) * 2.0)
    uw = sparse.csr_matrix(np.ones((n_cells, n_genes), dtype=np.float32) * 3.0)

    adata = ad.AnnData(
        X=base,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(n_genes)]),
    )
    adata.layers["counts_corrected"] = base.copy()
    adata.layers["counts_multimap_equal"] = equal
    adata.layers["counts_multimap_host_conservative"] = hc
    adata.layers["counts_multimap_unique_weighted"] = uw
    adata.uns["multimap_method"] = "equal"
    adata.write_h5ad(str(path))


class TestSwapMultimapLayer:
    def test_swap_to_host_conservative(self, tmp_path: Path) -> None:
        from viralscan.menu import _swap_multimap_layer

        h5ad = tmp_path / "multimap.h5ad"
        _make_multimap_h5ad(h5ad)

        result = _swap_multimap_layer(h5ad, "host-conservative")

        assert result is True
        adata = ad.read_h5ad(str(h5ad))
        # counts_corrected must now contain the host-conservative data (sentinel = 2.0)
        assert adata.layers["counts_corrected"].toarray().max() == pytest.approx(2.0)
        assert adata.uns["multimap_method"] == "host-conservative"

    def test_swap_to_equal(self, tmp_path: Path) -> None:
        from viralscan.menu import _swap_multimap_layer

        h5ad = tmp_path / "multimap.h5ad"
        _make_multimap_h5ad(h5ad)

        result = _swap_multimap_layer(h5ad, "equal")

        assert result is True
        adata = ad.read_h5ad(str(h5ad))
        assert adata.layers["counts_corrected"].toarray().max() == pytest.approx(1.0)
        assert adata.uns["multimap_method"] == "equal"

    def test_swap_to_unique_weighted(self, tmp_path: Path) -> None:
        from viralscan.menu import _swap_multimap_layer

        h5ad = tmp_path / "multimap.h5ad"
        _make_multimap_h5ad(h5ad)

        result = _swap_multimap_layer(h5ad, "unique-weighted")

        assert result is True
        adata = ad.read_h5ad(str(h5ad))
        assert adata.layers["counts_corrected"].toarray().max() == pytest.approx(3.0)
        assert adata.uns["multimap_method"] == "unique-weighted"

    def test_returns_false_when_layer_missing(self, tmp_path: Path) -> None:
        from viralscan.menu import _swap_multimap_layer

        # Write h5ad without the pre-stored layers (simulates an older run).
        n_cells, n_genes = 2, 3
        adata = ad.AnnData(
            X=sparse.csr_matrix((n_cells, n_genes)),
            obs=pd.DataFrame(index=["c0", "c1"]),
            var=pd.DataFrame(index=["g0", "g1", "g2"]),
        )
        h5ad = tmp_path / "old.h5ad"
        adata.write_h5ad(str(h5ad))

        result = _swap_multimap_layer(h5ad, "host-conservative")

        assert result is False
        # File must be unchanged (no counts_corrected written).
        adata2 = ad.read_h5ad(str(h5ad))
        assert "counts_corrected" not in adata2.layers

    def test_other_layers_preserved_after_swap(self, tmp_path: Path) -> None:
        from viralscan.menu import _swap_multimap_layer

        h5ad = tmp_path / "multimap.h5ad"
        _make_multimap_h5ad(h5ad)

        _swap_multimap_layer(h5ad, "host-conservative")

        adata = ad.read_h5ad(str(h5ad))
        assert "counts_multimap_equal" in adata.layers
        assert "counts_multimap_unique_weighted" in adata.layers


class TestRerunMultimapParser:
    def test_rerun_multimap_help(self) -> None:
        with patch("sys.argv", ["viralscan", "rerun-multimap", "--help"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit) as exc:
                create_help()
        assert exc.value.code == 0

    def test_rerun_multimap_parses_method(self) -> None:
        with patch("sys.argv", ["viralscan", "rerun-multimap", "-o", "out/", "--multimap-method", "host-conservative"]):
            from viralscan.menu import create_help

            args = create_help()
        assert args._subcommand == "rerun-multimap"
        assert args.multimap_method == "host-conservative"
        assert args.output == "out/"

    def test_rerun_multimap_rejects_unknown_method(self) -> None:
        with patch("sys.argv", ["viralscan", "rerun-multimap", "-o", "out/", "--multimap-method", "bogus"]):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()
