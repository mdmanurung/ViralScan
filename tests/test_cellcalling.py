"""Tests for cell-calling and the report-both-denominators summary change.

Covers the dependency-free paths (external list, knee) and that compute_stats
reports viral rates over BOTH called cells and all barcodes. emptyDrops (R) is
exercised only via the dispatch contract, not a live R call.
"""

from __future__ import annotations

import gzip

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from viralscan.enrichment import cell_type_enrichment
from viralscan.runconfig import RunConfig
from viralscan.scripts import cellcalling
from viralscan.scripts.cellcalling import call_cells, external_cells, knee_cells
from viralscan.scripts.detection import compute_stats


def test_emptydrops_r_script_is_next_to_cellcalling_module():
    from pathlib import Path

    import viralscan.scripts.cellcalling as cellcalling

    script = Path(cellcalling.__file__).with_name("emptydrops.R")
    assert script.exists()
    assert "emptyDrops" in script.read_text()


# ---------------------------------------------------------------------------
# external_cells
# ---------------------------------------------------------------------------
class TestExternalCells:
    def test_matches_plain_barcodes(self, tmp_path):
        f = tmp_path / "cells.txt"
        f.write_text("AAAA\nCCCC\n")
        mask = external_cells(["AAAA", "GGGG", "CCCC"], f)
        assert mask.tolist() == [True, False, True]

    def test_strips_dash_one_suffix_on_both_sides(self, tmp_path):
        # CellRanger writes 'AAAA-1'; ViralScan obs_names are bare 'AAAA'
        f = tmp_path / "cells.txt"
        f.write_text("AAAA-1\nCCCC-1\n")
        mask = external_cells(["AAAA", "CCCC", "TTTT"], f)
        assert mask.tolist() == [True, True, False]

    def test_reads_gzip(self, tmp_path):
        f = tmp_path / "cells.txt.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("AAAA-1\n")
        mask = external_cells(["AAAA", "GGGG"], f)
        assert mask.tolist() == [True, False]

    def test_no_match_returns_all_false(self, tmp_path):
        f = tmp_path / "cells.txt"
        f.write_text("ZZZZ\n")
        mask = external_cells(["AAAA", "CCCC"], f)
        assert mask.sum() == 0


# ---------------------------------------------------------------------------
# knee_cells
# ---------------------------------------------------------------------------
class TestKneeCells:
    def test_separates_clear_bimodal_population(self):
        # 20 "cells" at ~5000 UMI, 2000 empties at ~2 UMI
        rng_cells = np.full(20, 5000.0)
        rng_empty = np.full(2000, 2.0)
        total = np.concatenate([rng_cells, rng_empty])
        mask = knee_cells(total, min_umi=10)
        # all real cells kept, no empties
        assert mask[:20].all()
        assert not mask[20:].any()

    def test_too_few_barcodes_falls_back_to_min_umi(self):
        total = np.array([500.0, 3.0])
        mask = knee_cells(total, min_umi=10)
        assert mask.tolist() == [True, False]

    def test_returns_bool_mask_aligned_to_input(self):
        total = np.array([9000.0, 8000.0, 7000.0, 1.0, 1.0, 1.0, 1.0])
        mask = knee_cells(total, min_umi=5)
        assert mask.dtype == bool and len(mask) == len(total)


class TestAutoCellCalling:
    def test_auto_prefers_external_list(self, tmp_path, monkeypatch):
        called = tmp_path / "called.txt"
        called.write_text("bc0\n")
        adata = _make_adata()
        monkeypatch.setattr(
            cellcalling,
            "emptydrops_cells",
            lambda *_args, **_kwargs: pytest.fail("EmptyDrops should not run"),
        )
        mask = call_cells(
            adata,
            RunConfig(cell_calling="auto", called_cells_file=str(called)),
            matrix_dir=tmp_path,
        )
        assert mask.tolist() == [True, False, False, False, False]

    def test_auto_uses_emptydrops_without_external_list(self, tmp_path, monkeypatch):
        adata = _make_adata()
        expected = np.array([True, True, False, False, False])
        monkeypatch.setattr(cellcalling, "emptydrops_cells", lambda *_args, **_kwargs: expected)
        result = call_cells(adata, RunConfig(cell_calling="auto"), matrix_dir=tmp_path)
        np.testing.assert_array_equal(result, expected)

    def test_auto_does_not_silently_fall_back_to_knee(self, monkeypatch):
        adata = _make_adata()
        monkeypatch.setattr(
            cellcalling,
            "emptydrops_cells",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("R unavailable")),
        )
        with pytest.raises(RuntimeError, match="R unavailable"):
            call_cells(adata, RunConfig(cell_calling="auto"), matrix_dir="matrix")


# ---------------------------------------------------------------------------
# compute_stats — report both denominators
# ---------------------------------------------------------------------------
def _make_adata():
    # 5 barcodes x 3 genes: v1,v2 viral; h1 host. barcodes 0,1 are "called cells".
    X = np.array(
        [
            [10, 0, 100],  # bc0 called, viral+
            [0, 5, 200],  # bc1 called, viral+
            [1, 0, 3],  # bc2 empty, viral+ (ambient)
            [0, 0, 2],  # bc3 empty, viral-
            [2, 1, 1],
        ],  # bc4 empty, viral+
        dtype=float,
    )
    a = ad.AnnData(sp.csr_matrix(X))
    a.obs_names = ["bc0", "bc1", "bc2", "bc3", "bc4"]
    a.var_names = ["v1", "v2", "h1"]
    return a


class TestReportBothDenominators:
    def setup_method(self):
        self.adata = _make_adata()
        self.group = {"virusA": ["v1", "v2"]}
        self.called = np.array([True, True, False, False, False])

    def test_all_barcode_stats_unchanged_when_no_mask(self):
        stats, _ = compute_stats(self.adata, {}, self.group, [])
        s = stats["virusA"]
        # 4 of 5 barcodes have viral UMI (bc0,1,2,4)
        assert s["infected_cells"] == 4
        assert s["total_cells"] == 5
        assert s["pct_infected"] == pytest.approx(80.0)
        # with no mask, called == all
        assert s["n_called_cells"] == 5
        assert s["pct_infected_called"] == pytest.approx(80.0)

    def test_called_cell_denominator_restricts_correctly(self):
        stats, _ = compute_stats(self.adata, {}, self.group, [], called_mask=self.called)
        s = stats["virusA"]
        # all-barcode view stays the same
        assert s["infected_cells"] == 4 and s["total_cells"] == 5
        # called view: 2 called cells, both viral+
        assert s["n_called_cells"] == 2
        assert s["infected_called"] == 2
        assert s["pct_infected_called"] == pytest.approx(100.0)

    def test_per_cell_df_flags_called_cells(self):
        _, per_cell = compute_stats(self.adata, {}, self.group, [], called_mask=self.called)
        assert "is_called_cell" in per_cell.columns
        called_bc = set(per_cell.loc[per_cell["is_called_cell"], "barcode"])
        assert called_bc == {"bc0", "bc1"}

    def test_primary_call_matrix_controls_viral_numerators(self):
        primary = sp.csr_matrix(
            np.array(
                [
                    [10, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                ],
                dtype=float,
            )
        )

        stats, per_cell = compute_stats(self.adata, {}, self.group, [], viral_count_matrix=primary)
        s = stats["virusA"]

        assert s["viral_molecules_total_est"] == 10
        assert s["infected_cells"] == 1
        assert s["pct_infected"] == pytest.approx(20.0)
        assert s["viral_molecules_per_10k_est"] == pytest.approx(10 / self.adata.X.sum() * 10_000)
        assert per_cell["barcode"].tolist() == ["bc0"]


def test_cell_type_enrichment_uses_primary_call_matrix(tmp_path):
    adata = _make_adata()
    primary = sp.csr_matrix(
        np.array(
            [
                [10, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
            dtype=float,
        )
    )
    labels = tmp_path / "cell_types.csv"
    labels.write_text("barcode,cell_type\nbc0,T\nbc1,B\nbc2,B\nbc3,B\nbc4,B\n")

    result = cell_type_enrichment(
        adata,
        {"virusA": ["v1", "v2"]},
        RunConfig(cell_types=str(labels)),
        viral_count_matrix=primary,
    )

    by_type = result.set_index("cell_type")
    assert by_type.loc["T", "n_infected"] == 1
    assert by_type.loc["B", "n_infected"] == 0
