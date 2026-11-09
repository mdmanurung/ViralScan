"""CELLCALL-PCT-01: what ``pct_infected`` divides by under ``--cell-calling none`` and ``knee``.

The compute_stats tests pin the all-barcode behaviour. With a host index, ``none`` and ``knee``
are refused in ``menu.errorhandler`` (decision (a), 2026-11-09), so the two-step matrix never
reaches that path from ``run``.
"""

import argparse
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan import menu
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


def _run_args(tmp_path, **kw):
    (tmp_path / "R1.fastq.gz").touch()
    (tmp_path / "R2.fastq.gz").touch()
    base = {
        "ncbi_accession": None, "reference": False, "gtf": None, "fasta": None,
        "index": str(tmp_path / "index.idx"), "transcripts": "/fake/t2g.txt", "f1": None,
        "host_filter": "starsolo", "host_index": "/fake/host",
        "sample1": str(tmp_path / "R1.fastq.gz"), "sample2": str(tmp_path / "R2.fastq.gz"),
        "star_confirm": False, "star_hostvirus_index": None,
        "discover": False, "kraken2_db": None, "cell_calling": "auto",
    }  # fmt: skip
    base.update(kw)
    return argparse.Namespace(**base)


@pytest.mark.parametrize("method", ["none", "knee"])
def test_run_refuses_none_and_knee_with_a_host_index(tmp_path, method, caplog):
    """CELLCALL-PCT-01 (a): a two-step matrix holds only non-host barcodes, so every rate is 100 %."""
    with patch("os.path.exists", return_value=True):
        with pytest.raises(SystemExit) as exc:
            menu.errorhandler(_run_args(tmp_path, cell_calling=method))
    assert exc.value.code == 2
    assert f"--cell-calling {method} cannot be used with --host-index" in caplog.text


@pytest.mark.parametrize("method", ["auto", "emptydrops", "external"])
def test_run_accepts_other_cell_calling_with_a_host_index(tmp_path, method):
    called = tmp_path / "called.tsv"
    called.touch()
    args = _run_args(tmp_path, cell_calling=method, called_cells_file=str(called))
    with patch("os.path.exists", return_value=True):
        menu.errorhandler(args)


def test_none_without_a_host_index_is_still_allowed(tmp_path):
    args = _run_args(tmp_path, cell_calling="none", host_filter=None, host_index=None)
    with patch("os.path.exists", return_value=True):
        menu.errorhandler(args)
