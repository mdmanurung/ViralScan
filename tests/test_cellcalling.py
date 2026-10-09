"""Tests for cell-calling and the report-both-denominators summary change.

Covers the dependency-free paths (external list, knee) and that compute_stats
reports viral rates over BOTH called cells and all barcodes. emptyDrops (R) is
exercised only via the dispatch contract, not a live R call.
"""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from viralscan.enrichment import cell_type_enrichment
from viralscan.runconfig import RunConfig
from viralscan.scripts import cellcalling, detection
from viralscan.scripts.cellcalling import CellCallingError, call_cells, external_cells, knee_cells
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

    def test_no_match_is_fatal_rather_than_an_empty_mask(self, tmp_path):
        """SW-11: zero matches means the barcode spaces differ, not that there are no cells.

        Returning an all-false mask let the caller fall back to every barcode,
        which reports whole-droplet rates under a called-cell label.
        """
        f = tmp_path / "cells.txt"
        f.write_text("ZZZZ\n")
        with pytest.raises(CellCallingError, match="none of 2 barcodes matched"):
            external_cells(["AAAA", "CCCC"], f)

    def test_an_empty_called_cell_list_is_fatal(self, tmp_path):
        f = tmp_path / "cells.txt"
        f.write_text("\n  \n")
        with pytest.raises(CellCallingError, match="empty"):
            external_cells(["AAAA"], f)

    def test_canonical_barcode_collisions_are_fatal(self, tmp_path):
        """Two raw barcodes collapsing to one canonical form makes membership ambiguous."""
        f = tmp_path / "cells.txt"
        f.write_text("AAAA\n")
        with pytest.raises(CellCallingError, match="collide after suffix"):
            external_cells(["AAAA-1", "AAAA-2"], f)


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

    def test_configured_seed_and_niters_reach_emptydrops(self, tmp_path, monkeypatch):
        """The protocol freezes seeds.cell_calling; it must reach the R call.

        emptyDrops is a Monte-Carlo test, so a seed that stops at the config
        boundary makes the shared cell anchor irreproducible while every run
        still succeeds.
        """
        adata = _make_adata()
        seen: dict[str, object] = {}

        def _capture(*_args, **kwargs):
            seen.update(kwargs)
            return np.ones(adata.n_obs, dtype=bool)

        monkeypatch.setattr(cellcalling, "emptydrops_cells", _capture)
        call_cells(
            adata,
            RunConfig(cell_calling="emptydrops", emptydrops_seed=20260727002, emptydrops_niters=7),
            matrix_dir=tmp_path,
        )
        assert seen["seed"] == 20260727002
        assert seen["niters"] == 7

    def test_emptydrops_cells_requires_an_explicit_seed(self):
        """No signature default: an unsupplied seed must be an error, not 100."""
        with pytest.raises(TypeError, match="seed"):
            cellcalling.emptydrops_cells(
                ["bc0"],
                "matrix",
                rscript="Rscript",
                fdr=0.01,
                lower=100,
                niters=10,
            )

    def test_emptydrops_seed_is_forwarded_to_the_r_command(self, tmp_path, monkeypatch):
        """Guards the argv position the R script reads the seed from."""
        recorded: dict[str, list[str]] = {}

        def _fake_run(cmd, **_kwargs):
            recorded["cmd"] = cmd
            out_tsv = Path(cmd[3])
            out_tsv.write_text("barcode\tis_cell\nbc0\tTRUE\n", encoding="utf-8")
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr(cellcalling.subprocess, "run", _fake_run)
        cellcalling.emptydrops_cells(
            ["bc0"],
            tmp_path,
            rscript="Rscript",
            fdr=0.01,
            lower=100,
            niters=10,
            seed=20260727002,
        )
        assert recorded["cmd"][-1] == "20260727002"

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


class TestFailClosedContract:
    """SW-11: no failure path may silently turn every barcode into a cell."""

    def test_unknown_method_is_a_cell_calling_error(self):
        cfg = SimpleNamespace(cell_calling="magic")
        adata = ad.AnnData(np.ones((2, 2)))
        with pytest.raises(CellCallingError, match="Unknown cell_calling method"):
            call_cells(adata, cfg)

    def test_external_without_a_list_is_a_cell_calling_error(self):
        cfg = SimpleNamespace(cell_calling="external", called_cells_file=None)
        adata = ad.AnnData(np.ones((2, 2)))
        with pytest.raises(CellCallingError, match="requires config.called_cells_file"):
            call_cells(adata, cfg)

    def test_emptydrops_without_a_matrix_dir_is_a_cell_calling_error(self):
        cfg = SimpleNamespace(cell_calling="emptydrops")
        adata = ad.AnnData(np.ones((2, 2)))
        with pytest.raises(CellCallingError, match="counts_unfiltered"):
            call_cells(adata, cfg)

    def test_zero_barcodes_is_a_cell_calling_error(self):
        cfg = SimpleNamespace(cell_calling="knee")
        adata = ad.AnnData(np.zeros((0, 3)))
        with pytest.raises(CellCallingError, match="at least one barcode"):
            call_cells(adata, cfg)

    def test_knee_warns_that_it_is_sensitivity_only(self, caplog):
        # SW-23: knee lands at knee_min_umi on real libraries; decided
        # 2026-10-01 that reported numbers come from emptyDrops.
        cfg = SimpleNamespace(cell_calling="knee")
        adata = ad.AnnData(np.vstack([np.full((20, 2), 2500.0), np.ones((200, 2))]))
        with caplog.at_level(logging.WARNING, logger="viralscan"):
            call_cells(adata, cfg)
        assert "sensitivity-only" in caplog.text and "SW-23" in caplog.text

    def test_none_remains_the_explicit_way_to_use_every_barcode(self):
        """The escape hatch stays, but it must be asked for."""
        cfg = SimpleNamespace(cell_calling="none")
        adata = ad.AnnData(np.ones((4, 2)))

        mask = call_cells(adata, cfg)

        assert mask.all()
        assert mask.shape == (4,)


def test_emptydrops_calling_zero_cells_is_fatal(tmp_path, monkeypatch):
    """SW-11 claimed all four callers fail closed; emptyDrops did not.

    A zero-cell mask turns every called-cell rate into 0/0, which compute_stats
    reports as 0.0 — a number, not a failure.
    """

    def _fake_run(cmd, **_kwargs):
        Path(cmd[3]).write_text("barcode\tis_cell\nbc0\tFALSE\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(cellcalling.subprocess, "run", _fake_run)

    with pytest.raises(CellCallingError, match="called zero cells"):
        cellcalling.emptydrops_cells(
            ["bc0"], tmp_path, rscript="Rscript", fdr=0.01, lower=100, niters=10, seed=1
        )


def test_emptydrops_calling_at_least_one_cell_is_allowed(tmp_path, monkeypatch):
    def _fake_run(cmd, **_kwargs):
        Path(cmd[3]).write_text("barcode\tis_cell\nbc0\tTRUE\nbc1\tFALSE\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(cellcalling.subprocess, "run", _fake_run)
    mask = cellcalling.emptydrops_cells(
        ["bc0", "bc1"], tmp_path, rscript="Rscript", fdr=0.01, lower=100, niters=10, seed=1
    )

    assert mask.tolist() == [True, False]


# ---------------------------------------------------------------------------
# results/called_cells.tsv: the set layer 2 scores (PLAN PROG-17)
# ---------------------------------------------------------------------------
class TestCalledCellsFile:
    def _named(self):
        adata = ad.AnnData(X=sp.csr_matrix(np.ones((4, 1), dtype=np.float32)))
        adata.obs_names = ["AAA", "CCC", "GGG", "TTT"]
        return adata

    def test_written_set_round_trips_to_the_same_mask(self, tmp_path):
        adata = self._named()
        mask = np.array([True, False, True, False])
        cellcalling.write_called_cells(adata.obs_names, mask, tmp_path)
        got = cellcalling.load_called_mask(adata, SimpleNamespace(), tmp_path)
        assert got.tolist() == mask.tolist()

    def test_a_list_from_another_matrix_is_fatal(self, tmp_path):
        adata = self._named()
        (tmp_path / "results").mkdir()
        (tmp_path / "results" / "called_cells.tsv").write_text("barcode\nAAA\nNOT_HERE\n")
        with pytest.raises(CellCallingError, match="another run"):
            cellcalling.load_called_mask(adata, SimpleNamespace(), tmp_path)

    def test_pre_prog17_emptydrops_run_reuses_its_own_output(self, tmp_path, monkeypatch):
        adata = self._named()
        counts = tmp_path / "kb-python" / "counts_unfiltered"
        counts.mkdir(parents=True)
        (counts / "emptydrops_cells.tsv").write_text(
            "barcode\tis_cell\nAAA\tTRUE\nCCC\tFALSE\nTTT\tTRUE\n"
        )
        monkeypatch.setattr(cellcalling, "call_cells", lambda *a, **k: pytest.fail("re-called"))
        config = SimpleNamespace(cell_calling="auto", called_cells_file=None)
        got = cellcalling.load_called_mask(adata, config, tmp_path)
        assert got.tolist() == [True, False, False, True]

    def test_stale_emptydrops_output_is_ignored_for_an_external_run(self, tmp_path):
        adata = self._named()
        counts = tmp_path / "kb-python" / "counts_unfiltered"
        counts.mkdir(parents=True)
        (counts / "emptydrops_cells.tsv").write_text("barcode\tis_cell\nAAA\tTRUE\n")
        listed = tmp_path / "cells.txt"
        listed.write_text("GGG-1\nTTT-1\n")
        config = SimpleNamespace(cell_calling="auto", called_cells_file=str(listed))
        got = cellcalling.load_called_mask(adata, config, tmp_path)
        assert got.tolist() == [False, False, True, True]


class TestAccessionBreadth:
    """ANDET-01: breadth over the virus's index genes, not its detected genes."""

    def test_breadth_counts_undetected_index_genes(self):
        a = _make_adata()
        # Detection kept only v1, but the index holds four genes for virusA.
        detected = {"virusA": ["v1"]}
        index = {"virusA": ["v1", "v2", "v3_absent_from_matrix", "v4_absent_from_matrix"]}
        stats, _ = compute_stats(a, {}, detected, [], index_genes_by_virus=index)
        # v1 and v2 carry counts; four index genes -> 0.5, never the old 1.0.
        assert stats["virusA"]["accession_breadth"] == 0.5

    def test_without_the_index_map_the_old_definition_is_kept(self):
        stats, _ = compute_stats(_make_adata(), {}, {"virusA": ["v1"]}, [])
        assert stats["virusA"]["accession_breadth"] == 1.0


def _host_raw_fixture(tmp_path, layer="GeneFull"):
    raw = tmp_path / "host_filtered" / "star_tmp" / "Solo.out" / layer / "raw"
    raw.mkdir(parents=True)
    (raw / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n2 3 3\n1 1 900\n2 2 150\n1 3 3\n"
    )
    (raw / "barcodes.tsv").write_text("AAA\nCCC\nGGG\n")
    (raw / "features.tsv").write_text("host1\tHost1\nhost2\tHost2\n")
    return raw


def test_solo_raw_dir_only_for_twostep_runs(tmp_path):
    # DEF-05 replaces the old test that enshrined the audited Gene/raw reader.
    raw = _host_raw_fixture(tmp_path)
    assert cellcalling.solo_raw_dir(SimpleNamespace(output=str(tmp_path), host_index=None)) is None
    assert cellcalling.solo_raw_dir(SimpleNamespace(output=str(tmp_path), host_index="idx")) == raw


@pytest.mark.parametrize("layer", [None, "Gene"])
def test_twostep_missing_genefull_is_fatal_even_with_viral_counts(tmp_path, layer):
    if layer:
        _host_raw_fixture(tmp_path, layer)
    viral = tmp_path / "kb-python" / "counts_unfiltered"
    viral.mkdir(parents=True)
    (viral / "matrix.mtx").write_text("viral matrix must never be used")
    config = SimpleNamespace(output=str(tmp_path), host_index="idx", cell_calling="auto")
    with pytest.raises(CellCallingError, match="GeneFull/raw"):
        cellcalling.solo_raw_dir(config)


@pytest.mark.parametrize("missing", ["matrix.mtx", "barcodes.tsv", "features.tsv"])
def test_twostep_incomplete_host_matrix_is_fatal(tmp_path, missing):
    raw = _host_raw_fixture(tmp_path)
    (raw / missing).unlink()
    config = SimpleNamespace(output=str(tmp_path), host_index="idx", cell_calling="auto")
    with pytest.raises(CellCallingError, match=missing):
        cellcalling.solo_raw_dir(config)


@pytest.mark.parametrize("method", ["auto", "external", "none", "knee"])
def test_twostep_override_does_not_require_host_matrix(tmp_path, method):
    config = SimpleNamespace(
        output=str(tmp_path),
        host_index="idx",
        cell_calling=method,
        called_cells_file="cells.txt" if method in {"auto", "external"} else None,
    )
    assert cellcalling.solo_raw_dir(config) is None


def test_host_emptydrops_uses_genefull_fixture_and_records_input(tmp_path, monkeypatch):
    raw = _host_raw_fixture(tmp_path)
    config = SimpleNamespace(output=str(tmp_path), host_index="idx", cell_calling="auto")
    seen = []

    def fake_run(cmd, check):
        assert check is True
        seen.append(cmd)
        assert Path(cmd[2]) == raw
        Path(cmd[3]).write_text(
            "barcode\ttotal\tis_cell\nAAA\t900\tTRUE\nCCC\t150\tTRUE\nGGG\t3\tFALSE\n"
        )

    monkeypatch.setattr(cellcalling.subprocess, "run", fake_run)
    cells = cellcalling.host_called_cells(config, cellcalling.solo_raw_dir(config), 200)
    assert cells == {"called": {"AAA", "CCC"}, "comparable": {"AAA"}}
    assert len(seen) == 1
    receipt = json.loads((tmp_path / "results" / "cell_calling_input.json").read_text())
    assert receipt["matrix_dir"] == str(raw)
    assert receipt["count_layer"] == "GeneFull.raw"


@pytest.mark.parametrize("with_host_matrix", [False, True])
def test_external_twostep_keeps_cells_absent_from_viral_matrix(tmp_path, with_host_matrix):
    if with_host_matrix:
        _host_raw_fixture(tmp_path)
    listed = tmp_path / "cells.txt"
    listed.write_text("AAA-1\nCCC-1\nTTT-1\n")
    config = SimpleNamespace(output=str(tmp_path), called_cells_file=str(listed))
    host_cells = cellcalling.external_host_cells(config, 200)
    assert host_cells["called"] == {"AAA", "CCC", "TTT"}
    assert host_cells["comparable"] == ({"AAA"} if with_host_matrix else None)
    a = ad.AnnData(
        X=sp.csr_matrix([[5]]),
        obs=__import__("pandas").DataFrame(index=["AAA"]),
        var=__import__("pandas").DataFrame(index=["v1"]),
    )
    stats, _ = compute_stats(a, ["v1"], {"virusA": ["v1"]}, ["v1"], host_cells=host_cells)
    s = stats["virusA"]
    assert s["n_called_cells"] == 3
    assert s["pct_infected_called"] == pytest.approx(100 / 3, abs=0.0001)
    if with_host_matrix:
        assert s["n_comparable_cells"] == 1
        assert s["pct_infected_comparable"] == 100
    else:
        assert s["n_comparable_cells"] is None
        assert s["infected_comparable"] is None
        assert s["pct_infected_comparable"] is None


@pytest.mark.parametrize("partition", ["gene_role", "is_viral", "viral"])
def test_comparable_cells_use_host_counts_only(partition):
    a = ad.AnnData(X=sp.csr_matrix([[5000, 199], [1, 200], [5000, 500]]))
    a.var_names = ["v1", "host1"]
    a.var[partition] = ["target", "host"] if partition == "gene_role" else [True, False]
    called = np.array([True, True, False])
    stats, _ = compute_stats(a, ["v1"], {"virusA": ["v1"]}, ["v1"], called_mask=called)
    s = stats["virusA"]
    assert s["n_comparable_cells"] == 1
    assert s["infected_comparable"] == 1
    assert s["pct_infected_comparable"] == 100


@pytest.mark.parametrize("reason", ["virus_only", "missing_partition", "unavailable_identity"])
def test_comparable_cells_are_unavailable_without_host_depth(reason):
    a = ad.AnnData(X=sp.csr_matrix([[5000]]))
    a.var_names = ["v1"]
    if reason == "virus_only":
        a.uns["index_kind"] = "virus_only"
        a.var["is_viral"] = True
    elif reason == "unavailable_identity":
        a.uns["gene_identity_source"] = "unavailable"
        a.var["is_viral"] = False
    stats, _ = compute_stats(a, ["v1"], {"virusA": ["v1"]}, ["v1"])
    s = stats["virusA"]
    assert s["n_comparable_cells"] is None
    assert s["infected_comparable"] is None
    assert s["pct_infected_comparable"] is None


def _detection_fixture(tmp_path, monkeypatch, config):
    a = ad.AnnData(X=sp.csr_matrix([[5]]))
    a.obs_names = ["AAA"]
    a.var_names = ["v1"]
    monkeypatch.setattr(detection, "config", config)
    monkeypatch.setattr(
        detection, "preprocessing", lambda: (a, ["v1"], str(tmp_path), ["v1"], None, {})
    )
    monkeypatch.setattr(detection, "clear_stale_virus_plots", lambda _: None)
    monkeypatch.setattr(
        detection, "histogram", lambda *args, **kwargs: ({"virusA": ["v1"]}, ["v1"])
    )
    return a


def test_detection_missing_host_matrix_never_calls_viral_emptydrops(tmp_path, monkeypatch):
    config = SimpleNamespace(
        output=str(tmp_path), host_index="idx", visual=False, cell_calling="auto"
    )
    _detection_fixture(tmp_path, monkeypatch, config)

    def forbidden(*args, **kwargs):
        pytest.fail("two-step cell calling must not use the viral matrix")

    monkeypatch.setattr(cellcalling, "call_cells", forbidden)
    with pytest.raises(CellCallingError, match="GeneFull/raw"):
        detection.main()


def test_legacy_twostep_scoring_rejects_viral_emptydrops_fallback(tmp_path):
    a = ad.AnnData(X=sp.csr_matrix([[5]]))
    a.obs_names = ["AAA"]
    counts = tmp_path / "kb-python" / "counts_unfiltered"
    counts.mkdir(parents=True)
    (counts / "emptydrops_cells.tsv").write_text("barcode\tis_cell\nAAA\tTRUE\n")
    config = SimpleNamespace(output=str(tmp_path), host_index="idx", cell_calling="auto")
    with pytest.raises(CellCallingError, match="GeneFull/raw"):
        cellcalling.load_called_mask(a, config, tmp_path)


def test_legacy_twostep_scoring_uses_host_emptydrops_output(tmp_path):
    raw = _host_raw_fixture(tmp_path)
    (raw / "emptydrops_cells.tsv").write_text(
        "barcode\tis_cell\nAAA\tTRUE\nCCC\tTRUE\nGGG\tFALSE\n"
    )
    a = ad.AnnData(X=sp.csr_matrix([[5], [1]]))
    a.obs_names = ["AAA", "GGG"]
    config = SimpleNamespace(output=str(tmp_path), host_index="idx", cell_calling="auto")
    assert cellcalling.load_called_mask(a, config, tmp_path).tolist() == [True, False]


def test_detection_external_override_preserves_full_universe(tmp_path, monkeypatch):
    cells = tmp_path / "cells.txt"
    cells.write_text("AAA-1\nCCC-1\nTTT-1\n")
    config = SimpleNamespace(
        output=str(tmp_path),
        host_index="idx",
        visual=False,
        cell_calling="auto",
        called_cells_file=str(cells),
    )
    a = _detection_fixture(tmp_path, monkeypatch, config)

    class StopAfterCalling(Exception):
        pass

    def check_stats(*args, **kwargs):
        assert args[0] is a
        assert kwargs["called_mask"].tolist() == [True]
        assert kwargs["host_cells"] == {"called": {"AAA", "CCC", "TTT"}, "comparable": None}
        raise StopAfterCalling

    monkeypatch.setattr(detection, "compute_stats", check_stats)
    with pytest.raises(StopAfterCalling):
        detection.main()
    assert (
        tmp_path / "results" / "host_called_cells.tsv"
    ).read_text() == "barcode\nAAA\nCCC\nTTT\n"
    assert (tmp_path / "results" / "called_cells.tsv").read_text() == "barcode\nAAA\n"


def test_host_cells_from_tsv_splits_called_and_comparable(tmp_path):
    tsv = tmp_path / "ed.tsv"
    tsv.write_text(
        "barcode\ttotal\tFDR\tis_cell\tknee\tinflection\n"
        "AAA\t900\t0\tTRUE\t5\t5\n"
        "CCC\t150\t0\tTRUE\t5\t5\n"
        "GGG\t3\t1\tFALSE\t5\t5\n"
    )
    cells = cellcalling.host_cells_from_tsv(tsv, 200.0)
    assert cells == {"called": {"AAA", "CCC"}, "comparable": {"AAA"}}
    tsv.write_text("barcode\ttotal\tFDR\tis_cell\tknee\tinflection\nGGG\t3\t1\tFALSE\t5\t5\n")
    with pytest.raises(CellCallingError):
        cellcalling.host_cells_from_tsv(tsv, 200.0)


def test_compute_stats_twostep_denominator_is_the_host_cell_set():
    # 2 barcodes survive host filtering (both infected) out of 5 host-called cells:
    # the viral matrix alone would say 2/2 = 100 %.
    adata = ad.AnnData(
        X=sp.csr_matrix(np.array([[5.0], [3.0]])),
        obs=__import__("pandas").DataFrame(index=["AAA", "CCC"]),
        var=__import__("pandas").DataFrame(index=["V_gene1"]),
    )
    host_cells = {
        "called": {"AAA", "CCC", "GGG", "TTT", "ACG"},
        "comparable": {"AAA", "GGG", "TTT", "ACG"},
    }
    called = np.array([True, True])
    stats, _ = compute_stats(
        adata,
        ["V_gene1"],
        {"V": ["V_gene1"]},
        ["V_gene1"],
        called_mask=called,
        host_cells=host_cells,
    )
    v = stats["V"]
    assert v["n_called_cells"] == 5 and v["pct_infected_called"] == 40.0
    assert v["n_comparable_cells"] == 4 and v["infected_comparable"] == 1
    assert v["pct_infected_comparable"] == 25.0
