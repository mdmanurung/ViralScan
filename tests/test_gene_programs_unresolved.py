"""PROG-14: unresolved catalogue markers reach the summary, not just an info log."""

from __future__ import annotations

import anndata as ad
import pandas as pd
from scipy import sparse

from viralscan.gene_programs import EVIDENCE_LAYER, load_catalogue, summarise_programs
from viralscan.scripts.gene_programs import run_one

EBV = "Epstein-Barr virus"


def _adata(var_names):
    x = sparse.csr_matrix((2, len(var_names)))
    return ad.AnnData(x, var=pd.DataFrame(index=var_names), layers={EVIDENCE_LAYER: x})


def test_run_one_records_resolution_counts_and_warns(caplog) -> None:
    cat = load_catalogue()
    part = run_one(_adata(["EPSTEIN_HHV4_EBNA-1.1", "host"]), [EBV], cat, 2, "bundled")
    res = part.attrs["marker_resolution"][EBV]
    assert res["resolved"] >= 1 and res["unresolved"] > 0
    with caplog.at_level("WARNING"):
        none = run_one(_adata(["host"]), [EBV], cat, 2, "bundled")
    assert none.attrs["marker_resolution"][EBV]["resolved"] == 0
    assert any(r.levelname == "WARNING" and EBV in r.getMessage() for r in caplog.records)


def test_summary_carries_counts_and_the_right_caveat() -> None:
    cat = load_catalogue()
    empty = pd.DataFrame(columns=["barcode", "virus_name", "state"])
    none = summarise_programs(
        empty, cat, viruses=[EBV], marker_resolution={EBV: {"resolved": 0, "unresolved": 30}}
    ).iloc[0]
    assert (none["n_markers_resolved"], none["n_markers_unresolved"]) == (0, 30)
    assert "none of 30" in none["caveat"]
    assert "uniquely-placing" not in none["caveat"]

    part = summarise_programs(
        empty, cat, viruses=[EBV], marker_resolution={EBV: {"resolved": 20, "unresolved": 10}}
    ).iloc[0]
    assert (part["n_markers_resolved"], part["n_markers_unresolved"]) == (20, 10)
    assert "10 of 30" in part["caveat"] and "reduced" in part["caveat"]
