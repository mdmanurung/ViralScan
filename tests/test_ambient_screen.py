"""M-03: descriptive ambient (soup) screen of `compute_stats`, with hand-computed expectations.

Design: docs/plans/ambient-and-barcode-swapping-controls.md section 3.1. The screen is the
calibration-free rho = 1 bound: with `f_U` the 95 % Poisson upper bound on the viral fraction of the
ambient barcodes, a called cell of N molecules receives at most Poisson(f_U * N) soup molecules.
Poisson tails below are closed forms (1 - e^-l, ...), not scipy calls.
"""

from __future__ import annotations

import csv
import math

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from viralscan.evidence import interpretation_flags  # noqa: E402
from viralscan.scripts import detection as D  # noqa: E402

LN20 = -math.log(0.05)  # 2.995732..., the 95 % upper bound on a Poisson mean given 0 events


def _adata(X, genes, cells, roles=None):
    a = ad.AnnData(np.asarray(X, dtype=float))
    a.var_names = genes
    a.obs_names = cells
    if roles:
        a.var["gene_role"] = roles
    return a


def _stats(a, called, **kw):
    stats, _ = D.compute_stats(a, {}, {"V": ["v1"]}, ["v1"], called_mask=called, **kw)
    return stats["V"]


def test_upper_bounds_are_the_exact_poisson_limits():
    assert D._poisson_ub(0) == pytest.approx(LN20, abs=1e-9)
    assert D._poisson_ub(1) == pytest.approx(4.743865, abs=1e-6)
    assert D._poisson_ub(2) == pytest.approx(6.295794, abs=1e-6)


def _single_run():
    # c0, c1 called (500 molecules each, c1 holds the one viral molecule); e0..e2 are empty
    # droplets of 40, 60 and 100 molecules (100 = the ambient ceiling, inclusive).
    X = [[0, 500], [1, 499], [0, 40], [0, 60], [0, 100]]
    a = _adata(X, ["v1", "h"], ["c0", "c1", "e0", "e1", "e2"], ["viral", "host"])
    return a, np.array([True, True, False, False, False])


def test_zero_ambient_viral_molecules_give_a_bound_that_cannot_reject_soup():
    """V_A = 0, T_A = 200: f_U = ln(20)/200. Each 500-molecule cell: l = 7.489, P(>=1) = 0.99944."""
    a, called = _single_run()
    s = _stats(a, called, cell_calling_method="emptydrops")
    lam = LN20 / 200 * 500
    e1 = 2 * (1 - math.exp(-lam))
    e2 = 2 * (1 - math.exp(-lam) * (1 + lam))
    assert s["ambient_viral_fraction"] == 0
    assert s["ambient_viral_fraction_ub95"] == pytest.approx(LN20 / 200, rel=1e-5)
    assert s["n_ambient_barcodes_with_virus"] == 0
    assert s["soup_expected_cells_max"] == pytest.approx(e1, rel=1e-5)
    assert s["soup_expected_cells_ge2_max"] == pytest.approx(e2, rel=1e-5)
    assert s["soup_tail_p_cells"] == pytest.approx(1 - math.exp(-e1), rel=1e-5)  # P(X >= 1)
    assert s["ambient_assessment"] == "soup_sufficient"
    assert s["viral_depth_spearman"] is None  # needs >= 3 infected called cells
    # Descriptive only: the gate is untouched by the screen.
    assert s["detected"] is True and s["infected_called"] == 1


def test_many_deep_infected_cells_over_a_deep_ambient_set_reject_soup():
    """1,000 ambient barcodes of 100 molecules, no viral molecule: T_A = 1e5, f_U = ln(20)/1e5.

    Three called cells of 1,000 / 2,000 / 3,000 molecules hold 1 / 2 / 3 viral molecules.
    E1 = sum(1 - e^-(f_U N)); the tail of 3 observed cells given E1 is 1 - e^-E1 (1 + E1 + E1^2/2).
    """
    n_amb = 1000
    X = [[1, 999], [2, 1998], [3, 2997]] + [[0, 100]] * n_amb
    cells = ["c0", "c1", "c2"] + [f"e{i}" for i in range(n_amb)]
    a = _adata(X, ["v1", "h"], cells, ["viral", "host"])
    called = np.array([True] * 3 + [False] * n_amb)
    s = _stats(a, called, cell_calling_method="emptydrops")
    f_u = LN20 / 1e5
    e1 = sum(1 - math.exp(-f_u * n) for n in (1000, 2000, 3000))
    tail = 1 - math.exp(-e1) * (1 + e1 + e1**2 / 2)
    assert s["soup_expected_cells_max"] == pytest.approx(e1, rel=1e-4)
    assert s["soup_tail_p_cells"] == pytest.approx(tail, rel=1e-3)
    assert s["soup_tail_p_cells"] < D.AMBIENT_TAIL_ALPHA
    assert s["ambient_assessment"] == "soup_rejected_rho1"
    assert s["viral_depth_spearman"] == 1.0  # viral and depth rise together over the called cells
    assert "not a clean result" in s["ambient_basis"] and "not assessed" in s["ambient_basis"]
    assert s["detected"] is True


def test_ambient_viral_molecules_enter_the_fraction_and_the_count_of_barcodes():
    """One viral molecule in an ambient barcode: V_A = 1, T_A = 41 (40 + 1), f_hat = 1/41."""
    X = [[0, 500], [1, 499], [1, 40]]
    a = _adata(X, ["v1", "h"], ["c0", "c1", "e0"], ["viral", "host"])
    s = _stats(a, np.array([True, True, False]), cell_calling_method="emptydrops")
    assert s["ambient_viral_fraction"] == pytest.approx(1 / 41, rel=1e-5)
    assert s["ambient_viral_fraction_ub95"] == pytest.approx(4.743865 / 41, rel=1e-5)
    assert s["n_ambient_barcodes_with_virus"] == 1


def _assert_blank_not_assessed(s, reason_part):
    assert s["ambient_assessment"] == "not_assessed"
    assert reason_part in s["ambient_basis"]
    assert "Swapping between libraries is not assessed" in s["ambient_basis"]
    for column in D.AMBIENT_COLUMNS:
        if column not in ("ambient_assessment", "ambient_basis"):
            assert s[column] is None, column  # blank, never 0


def test_virus_only_index_is_not_assessed_and_never_a_number():
    a = _adata([[1], [0], [1]], ["v1"], ["c0", "c1", "e0"])
    a.uns["index_kind"] = "virus_only"
    s = _stats(a, np.array([True, True, False]), cell_calling_method="emptydrops")
    _assert_blank_not_assessed(s, "no host genes")
    assert s["detected"] is True


def test_run_without_a_host_gene_column_is_not_assessed():
    a = _adata([[1, 500], [0, 40]], ["v1", "h"], ["c0", "e0"])  # no gene_role / is_viral
    s = _stats(a, np.array([True, False]), cell_calling_method="emptydrops")
    _assert_blank_not_assessed(s, "no host genes")


@pytest.mark.parametrize(
    "method, part", [("none", "cell calling is off"), ("knee", "sensitivity-only")]
)
def test_cell_calling_none_and_knee_are_not_assessed(method, part):
    a, called = _single_run()
    _assert_blank_not_assessed(_stats(a, called, cell_calling_method=method), part)


def test_no_empty_droplet_is_not_assessed():
    a, _ = _single_run()
    s = _stats(a, np.ones(5, dtype=bool), cell_calling_method="emptydrops")
    _assert_blank_not_assessed(s, "no non-called barcode")


def test_barcodes_deeper_than_the_ceiling_are_not_ambient():
    a, called = _single_run()
    s = _stats(a, called, cell_calling_method="emptydrops", ambient_max_molecules=50)
    # only e0 (40 molecules) is now ambient: T_A = 40
    assert s["ambient_viral_fraction_ub95"] == pytest.approx(LN20 / 40, rel=1e-5)


def _two_step(host_totals):
    # Viral matrix of a two-step run: only barcodes with a non-host read. Suffixes are raw, host sets
    # are stripped (B-10). cB is a called cell and amb2 an ambient droplet, neither in this matrix.
    a = _adata([[1], [1]], ["v1"], ["cA-1", "amb1-1"], ["viral"])
    a.uns["index_kind"] = "virus_only"  # host depth comes from STARsolo, not the kb matrix
    host_cells = {"called": {"cA", "cB"}, "comparable": {"cA", "cB"}}
    if host_totals is not None:
        host_cells["host_totals"] = host_totals
    return a, host_cells, np.array([True, False])


def test_two_step_run_counts_barcodes_the_viral_matrix_lacks():
    """Depths: cA 500 + 1, amb1 30 + 1, so V_A = 1 and T_A = 31 + amb2 50 = 81 (big, 5000, is not ambient).

    f_U = 4.743865 / 81. Called cells cA (501) and cB (700, absent from the viral matrix) both reach
    P(>=1) ~ 1, so E1 = 2 (a barcode-less cB would give 1). One infected called cell: tail = 1 - e^-2.
    """
    a, host_cells, called = _two_step(
        {"cA": 500.0, "cB": 700.0, "amb1": 30.0, "amb2": 50.0, "big": 5000.0}
    )
    s = _stats(a, called, host_cells=host_cells, cell_calling_method="emptydrops")
    assert s["ambient_viral_fraction"] == pytest.approx(1 / 81, rel=1e-5)
    assert s["n_ambient_barcodes_with_virus"] == 1
    assert s["soup_expected_cells_max"] == pytest.approx(2.0, abs=1e-6)
    assert s["soup_tail_p_cells"] == pytest.approx(1 - math.exp(-2.0), rel=1e-5)
    assert s["ambient_assessment"] == "soup_sufficient"
    assert s["n_called_cells"] == 2


def test_two_step_run_without_host_totals_is_not_assessed():
    a, host_cells, called = _two_step(None)
    s = _stats(a, called, host_cells=host_cells, cell_calling_method="external")
    _assert_blank_not_assessed(s, "GeneFull/raw")


def test_detected_does_not_depend_on_the_screen():
    """The same matrix with the screen on and off (`none` blanks it) detects the same."""
    a, called = _single_run()
    on = _stats(a, called, cell_calling_method="emptydrops")
    off = _stats(a, called, cell_calling_method="none")
    assert on["ambient_assessment"] == "soup_sufficient"
    assert off["ambient_assessment"] == "not_assessed"
    for key in ("detected", "infected_called", "background_umis", "pct_infected_called"):
        assert on[key] == off[key]


def test_columns_reach_viral_summary_tsv_and_blank_where_undefined(tmp_path):
    a, called = _single_run()
    stats = {"V": _stats(a, called, cell_calling_method="emptydrops")}
    stats["W"] = _stats(a, called, cell_calling_method="none")
    D.write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path))
    with open(tmp_path / "results" / "viral_summary.tsv", newline="") as fh:
        rows = {r["virus_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["V"]["ambient_assessment"] == "soup_sufficient"
    assert float(rows["V"]["soup_tail_p_cells"]) > 0
    assert rows["W"]["ambient_assessment"] == "not_assessed"
    assert rows["W"]["soup_tail_p_cells"] == "" and rows["W"]["ambient_viral_fraction"] == ""


def test_alignment_only_rows_are_not_assessed():
    row = {c: "" for c in ("virus_name", "n_called_cells", "n_comparable_cells", "total_cells")}
    template = D._alignment_only_template([row])
    assert template["ambient_assessment"] == "not_assessed"
    assert template["soup_tail_p_cells"] is None


def test_summary_line_names_the_screen_without_a_semicolon():
    a, called = _single_run()
    line = D._virus_summary_line("V", _stats(a, called, cell_calling_method="emptydrops"))
    assert "Ambient screen" in line and "soup_sufficient" in line and ";" not in line
    off = D._virus_summary_line("V", _stats(a, called, cell_calling_method="none"))
    assert "Ambient screen" not in off


def test_contamination_flag_carries_the_screen_but_never_reads_clean():
    default = {r["flag"]: r for r in interpretation_flags([], [], [])}["contamination"]
    assert default["status"] == "not_assessed" and default["metric"] == ""
    ambient = {
        "ambient_assessment": "soup_rejected_rho1",
        "ambient_basis": "rho = 1 soup bound",
        "soup_tail_p_cells": "0.0001",
    }
    flags = {r["flag"]: r for r in interpretation_flags([], [], [], ambient=ambient)}
    row = flags["contamination"]
    assert row["status"] == "not_assessed"  # swapping between libraries is never assessed in a run
    assert row["metric"] == "0.0001" and row["basis"].startswith("soup_rejected_rho1:")
    unassessed = {"ambient_assessment": "not_assessed"}
    flags = {r["flag"]: r for r in interpretation_flags([], [], [], ambient=unassessed)}
    assert flags["contamination"] == default
