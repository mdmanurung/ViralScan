"""CELLCALL-PCT-01: what ``pct_infected`` divides by under ``--cell-calling none`` and ``knee``.

These tests pin today's behaviour so a later policy decision shows up as a deliberate edit.
Nothing here is a fix. See the PLAN row for the open choice.
"""

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from viralscan.scripts.cellcalling import call_cells
from viralscan.scripts.detection import compute_stats

KW = dict(
    infected_min_molecules=1.0,
    detection_threshold=1,
    detection_min_cells=1,
    ambient_max_molecules=100,
)


def _stats(adata, method, host_cells=None):
    config = SimpleNamespace(cell_calling=method, knee_min_umi=10.0)
    called = call_cells(adata, config)
    stats, _ = compute_stats(
        adata,
        {"v1"},
        {"VirusX": ["v1"]},
        {"v1"},
        called_mask=called,
        host_cells=host_cells,
        cell_calling_method=method,
        **KW,
    )
    return stats["VirusX"], called


def _combined():
    """50 deep cells (2 infected) plus 950 empty droplets of 3 molecules each."""
    rows = [[5000, 100, 5 if i < 2 else 0] for i in range(50)] + [[3, 0, 0]] * 950
    return ad.AnnData(
        X=sp.csr_matrix(np.array(rows, float)),
        obs=pd.DataFrame(index=[f"b{i}" for i in range(1000)]),
        var=pd.DataFrame({"gene_role": ["host", "host", "viral"]}, index=["h1", "h2", "v1"]),
    )


def _two_step_matrix():
    """A kb matrix of a two-step run: only barcodes with a non-host read, viral genes only."""
    return ad.AnnData(
        X=sp.csr_matrix(np.array([[5.0], [3.0], [1.0], [1.0]])),
        obs=pd.DataFrame(index=["b0", "b1", "b2", "b3"]),
        var=pd.DataFrame(index=["v1"]),
    )


def test_none_divides_by_every_barcode_of_the_matrix():
    s, called = _stats(_combined(), "none")
    assert called.all()
    assert s["n_called_cells"] == s["total_cells"] == 1000
    assert s["pct_infected_called"] == s["pct_infected"] == 0.2  # 2 / 1000, empty droplets in


def test_knee_keeps_a_called_subset_on_a_combined_matrix():
    s, called = _stats(_combined(), "knee")
    assert 0 < called.sum() < 1000
    assert s["n_called_cells"] == int(called.sum())
    assert s["pct_infected"] == 0.2  # the all-barcode rate never uses the called set


def test_two_step_none_has_no_host_denominator_note_today():
    """OPEN (CELLCALL-PCT-01): with a host index and ``none``, ``host_cells`` stays None.

    The viral matrix holds only barcodes with a non-host read, so every rate is 100 % and the
    AUD3-03 NA guard (``if host_cells``) does not fire. Update this test when the choice is made.
    """
    s, called = _stats(_two_step_matrix(), "none")
    assert called.all()
    assert s["n_called_cells"] == 4
    assert s["pct_infected_called"] == s["pct_infected"] == 100.0
    assert s["denominator_note"] == ""
