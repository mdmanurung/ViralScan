"""Tests for _swap_multimap_layer and the rerun-multimap CLI subcommand."""

from __future__ import annotations

import os
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
    adata.layers["counts_unique"] = base.copy()
    adata.layers["counts_multimap_equal"] = equal
    adata.layers["counts_multimap_host_conservative"] = hc
    adata.layers["counts_multimap_unique_weighted"] = uw
    adata.uns["multimap_method"] = "equal"
    adata.uns["count_schema_version"] = "3.0.0"
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
        assert adata.layers["counts_ambiguous_allocated"].toarray().max() == pytest.approx(2.0)
        assert adata.X.toarray().max() == pytest.approx(2.0)
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

    def test_pre_v3_h5ad_is_not_numerically_migrated(self, tmp_path: Path) -> None:
        from viralscan.menu import _swap_multimap_layer

        h5ad = tmp_path / "legacy.h5ad"
        _make_multimap_h5ad(h5ad)
        adata = ad.read_h5ad(h5ad)
        del adata.uns["count_schema_version"]
        before = adata.X.copy()
        adata.write_h5ad(h5ad)

        assert _swap_multimap_layer(h5ad, "host-conservative") is False
        np.testing.assert_array_equal(ad.read_h5ad(h5ad).X.toarray(), before.toarray())

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
        with patch(
            "sys.argv",
            [
                "viralscan",
                "rerun-multimap",
                "--run-dir",
                "source/",
                "-o",
                "out/",
                "--multimap-method",
                "host-conservative",
            ],
        ):
            from viralscan.menu import create_help

            args = create_help()
        assert args._subcommand == "rerun-multimap"
        assert args.multimap_method == "host-conservative"
        assert args.run_dir == "source/"
        assert args.output == "out/"

    def test_rerun_multimap_rejects_unknown_method(self) -> None:
        with patch(
            "sys.argv",
            [
                "viralscan",
                "rerun-multimap",
                "--run-dir",
                "source/",
                "-o",
                "out/",
                "--multimap-method",
                "bogus",
            ],
        ):
            from viralscan.menu import create_help

            with pytest.raises(SystemExit):
                create_help()


class TestRerunLeavesSourceUntouched:
    """SW-19: drive the real ``_run_rerun_multimap``, with snakemake stubbed out."""

    @staticmethod
    def _source_run(root: Path) -> Path:
        from viralscan.runconfig import RunConfig

        sample = root / "source" / "SAMPLE"
        (sample / "log").mkdir(parents=True)
        (sample / "log" / "multimap.done").touch()
        h5ad = sample / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
        h5ad.parent.mkdir(parents=True)
        _make_multimap_h5ad(h5ad)
        RunConfig(output=str(sample) + "/", cell_calling="none").to_yaml(sample / "config.yaml")
        return root / "source"

    @staticmethod
    def _sha(path: Path) -> str:
        import hashlib

        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_fast_layer_swap_writes_the_copy_not_the_source(self, tmp_path: Path) -> None:
        import argparse

        from viralscan import menu

        source = self._source_run(tmp_path)
        rel = Path("SAMPLE/kb-python/counts_unfiltered/adata_multimap.h5ad")
        before = self._sha(source / rel)
        args = argparse.Namespace(
            run_dir=str(source),
            output=str(tmp_path / "rerun"),
            multimap_method="host-conservative",
            multimap_em_max_iter=None,
            multimap_em_tol=None,
            cores=1,
            verbose=False,
            quiet=True,
        )
        with (
            patch.object(menu, "_check_required_tools"),
            patch.object(menu, "_check_cell_caller_tools"),
            patch.object(menu.subprocess, "run") as run,
        ):
            menu._run_rerun_multimap(args)

        assert self._sha(source / rel) == before, "rerun-multimap modified the source run"
        swapped = ad.read_h5ad(tmp_path / "rerun" / rel)
        assert swapped.uns["multimap_method"] == "host-conservative"

        # snakemake gets the copy's sample dir, with the trailing separator the
        # Snakefile's f"{config['output']}log/..." paths depend on.
        snakemake_argv = run.call_args_list[0].args[0]
        expected = f"output={tmp_path / 'rerun' / 'SAMPLE'}{os.sep}"
        assert expected in snakemake_argv

    def test_rewritten_config_keeps_its_mtime_so_kb_count_is_not_rerun(
        self, tmp_path: Path
    ) -> None:
        """SW-22: config.yaml is an input of kb_count, host_filter and analysis.

        Rewriting it with a fresh mtime made snakemake's mtime trigger re-run
        kallisto and STAR on the FASTQs for what should be a layer swap.
        """
        import argparse

        from viralscan import menu

        source = self._source_run(tmp_path)
        old = 1_000_000_000
        os.utime(source / "SAMPLE" / "config.yaml", (old, old))
        args = argparse.Namespace(
            run_dir=str(source),
            output=str(tmp_path / "rerun"),
            multimap_method="host-conservative",
            multimap_em_max_iter=None,
            multimap_em_tol=None,
            cores=1,
            verbose=False,
            quiet=True,
        )
        with (
            patch.object(menu, "_check_required_tools"),
            patch.object(menu, "_check_cell_caller_tools"),
            patch.object(menu.subprocess, "run"),
        ):
            menu._run_rerun_multimap(args)

        copied = tmp_path / "rerun" / "SAMPLE" / "config.yaml"
        assert "host-conservative" in copied.read_text()  # it was rewritten
        assert copied.stat().st_mtime == old, "rewrite made config.yaml newer than outputs"
