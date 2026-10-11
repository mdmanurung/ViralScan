"""Regression tests for the 2026-10 audit, bundle S2 (B-05..B-09)."""

from __future__ import annotations

import re
from types import SimpleNamespace

import anndata as ad
import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from viralscan.runconfig import RunConfig  # noqa: E402
from viralscan.scripts import detection as D  # noqa: E402


def _adata(X, genes):
    X = np.asarray(X, dtype=float)
    a = ad.AnnData(X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    return a


class TestSummaryHeadline:
    STATS = {
        "viral_molecules_total_est": 7,
        "infected_cells": 2,
        "total_cells": 3,
        "pct_infected": 66.67,
        "viral_molecules_per_10k_est": 1.5,
        "n_called_cells": 2,
        "infected_called": 2,
        "pct_infected_called": 100.0,
    }

    def test_called_cells_lead_and_all_barcodes_are_labelled(self):
        line = D._virus_summary_line("V", self.STATS)
        assert line.index("2/2 called cells") < line.index("2/3 all barcodes")
        assert "(100.00%)" in line and "(66.67%)" in line

    def test_no_semicolon(self):
        """compare_legacy_v2_v3._GENE_TOTAL_RE reads any ``x;number`` line as a gene total."""
        assert ";" not in D._virus_summary_line("V", self.STATS)


class TestBelowGateCountsContract:
    def test_default_tuple_index_viruses_does_not_crash(self):
        a = _adata([[5, 0, 500]], ["v1", "v2", "h"])
        assert D.below_gate_counts(a, None, (), {}) == {}

    def test_mapping_still_counts(self):
        a = _adata([[5, 0, 500]], ["v1", "v2", "h"])
        assert D.below_gate_counts(a, None, {"V": ["v1"], "W": ["v2"]}, {}) == {"V": 5.0}


class TestPositiveControlExpected:
    @pytest.mark.parametrize("expected", [0, -5])
    def test_non_positive_expected_is_a_clear_error(self, expected):
        a = _adata([[5, 0, 500]], ["v1", "v2", "h"])
        cfg = SimpleNamespace(
            positive_control_gene="v1",
            positive_control_expected_molecules=expected,
            positive_control_scope=None,
        )
        with pytest.raises(ValueError, match="must be > 0"):
            D.measure_positive_control(a, cfg)


class TestDetectGenesInputs:
    def test_nan_counts_raise_not_silently_dropped(self):
        with pytest.raises(ValueError, match="non-finite"):
            D.detect_genes(["h1", "v1"], np.array([[0, np.nan], [0, 3]]), ["v1"])

    def test_duplicate_var_names_raise(self):
        with pytest.raises(ValueError, match="duplicate var_names"):
            D.detect_genes(["v1", "v1"], np.array([[5, 0], [0, 0]]), ["v1"])


class TestSuperExpressorInputs:
    @pytest.fixture(autouse=True)
    def _cfg(self, monkeypatch):
        monkeypatch.setattr(D, "config", RunConfig(), raising=False)

    def test_zero_cell_adata_is_a_clear_error(self, tmp_path):
        a = _adata(np.zeros((0, 3)), ["v1", "v2", "h"])
        (tmp_path / "plots").mkdir()
        with pytest.raises(ValueError, match="no cells"):
            D.super_expressor(a, "V", ["v1"], str(tmp_path))

    def test_caller_adata_is_not_mutated(self, tmp_path):
        a = _adata([[3, 0, 5], [0, 0, 9], [1, 1, 1]], ["v1", "v2", "h"])
        (tmp_path / "plots").mkdir()
        D.super_expressor(a, "V", ["v1"], str(tmp_path))
        assert list(a.obs.columns) == []
        assert (tmp_path / "plots" / "SuperExpressor_V.png").exists()

    def test_docstring_states_the_real_rule(self):
        assert re.search(r">= ``config.se_threshold``", D.super_expressor.__doc__)
