"""AUD3-03: report headline denominators (P3-03) and two-step all-barcode rates (P4-02).

``pct_infected_called`` is the primary rate everywhere it is shown. In two-step mode the
kb matrix holds only barcodes with a non-host read, so ``pct_infected`` and
``viral_molecules_per_10k_est`` have no meaningful denominator and are written as NA.
"""

import csv
import re

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection
from viralscan.scripts.detection import _in_barcode_set, compute_stats, write_tsv_outputs
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

KW = dict(
    infected_min_molecules=1.0,
    detection_threshold=1,
    detection_min_cells=1,
    ambient_max_molecules=100,
)


def _text(html):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))


def test_report_headlines_called_cell_rate(tmp_path):
    """100 barcodes, 10 called; Virus A in 1 called cell + 8 empty droplets (P3-03 repro)."""
    run = tmp_path / "run"
    (run / "results").mkdir(parents=True)
    (run / "log").mkdir()
    d = run / "kb-python" / "counts_unfiltered"
    d.mkdir(parents=True)

    def gi(g, viral, key="", name=""):
        return GeneIdentity(
            gene_id=g,
            genome_accession=g,
            status="catalogued" if viral else "host",
            viral=viral,
            virus_key=key,
            virus_name=name,
        )

    VirusIdentityTable((gi("h1", False), gi("vA_g1", True, "A", "Virus A"))).write_tsv(
        run / "results" / "virus_identity.tsv"
    )
    (run / "log" / "analysis.txt").write_text("vA_g1\n")
    X = np.zeros((100, 2), dtype=np.float32)
    X[:10, 0] = 500
    X[0, 1] = 1
    X[10:18, 1] = 1
    a = ad.AnnData(X=sp.csr_matrix(X))
    a.obs_names = [f"BC{i:03d}" for i in range(100)]
    a.var_names = ["h1", "vA_g1"]
    a.write_h5ad(d / "adata.h5ad")
    (run / "cells.txt").write_text("\n".join(f"BC{i:03d}" for i in range(10)) + "\n")
    cfg = RunConfig(
        output=str(run) + "/",
        multimapping=False,
        visual=False,
        cell_calling="external",
        called_cells_file=str(run / "cells.txt"),
        detection_threshold=1,
        detection_min_cells=1,
    )
    detection.run(
        RunContext.from_config(cfg), str(run / "log" / "analysis.txt"), str(run / "log" / "d")
    )

    row = next(csv.DictReader((run / "results" / "viral_summary.tsv").open(), delimiter="\t"))
    assert float(row["pct_infected_called"]) == 10.0 and float(row["pct_infected"]) == 9.0
    txt = _text((run / "report.html").read_text())
    summary = txt[txt.index("Run Summary") : txt.index("Per-Virus Statistics")]
    # headline: called-cell rate first, all-barcode rate present but labelled secondary
    assert "1 / 10" in summary and "10.00%" in summary
    assert re.search(r"all barcodes[^%]*9\.00%", summary, re.I)
    assert summary.index("10.00%") < summary.index("9.00%")
    table = txt[txt.index("Per-Virus Statistics") :]
    i = table.index("Virus A")
    cells = table[i : i + 200]
    assert cells.index("10.00%") < cells.index("9.00%")
    assert re.search(r"% Infected \(called cells\)", table)
    assert re.search(r"% Infected \(all barcodes", table)


def _row(pct_called, pct_all=0.5, **extra):
    return {
        "viral_molecules_total_est": 3,
        "infected_cells": 5,
        "total_cells": 1000,
        "pct_infected": pct_all,
        "viral_molecules_per_10k_est": 1.0,
        "infected_called": 5,
        "n_called_cells": 100,
        "pct_infected_called": pct_called,
        "detected": True,
        **extra,
    }


def _badge_html(tmp_path, stats):
    detection.generate_html_report(stats, pd.DataFrame(), None, None, {}, [], str(tmp_path))
    return (tmp_path / "report.html").read_text()


def test_badge_uses_called_rate_not_all_barcode_rate(tmp_path):
    # all-barcode 0.5% (no badge) but called 5% -> badge on the called figure
    html = _badge_html(tmp_path, {"V": _row(5.0, 0.5)})
    assert re.search(r'badge-high">\s*5\.00%', html)
    # all-barcode 50% but called 0.5% -> no badge
    html = _badge_html(tmp_path, {"V": _row(0.5, 50.0)})
    assert "badge-high" not in html.split("Per-Virus Statistics")[1]


def test_comparable_rate_is_headlined_when_no_called_cells(tmp_path):
    s = _row(
        0.0,
        0.5,
        n_called_cells=0,
        pct_infected_comparable=2.5,
        n_comparable_cells=40,
        infected_comparable=1,
    )
    html = _badge_html(tmp_path, {"V": s})
    assert re.search(r'badge-high">\s*2\.50%', html)


# --- P4-02 -------------------------------------------------------------------------------
def _combined_and_two_step():
    cells = [f"cell{i}" for i in range(1000)]
    empties = [f"empty{i}" for i in range(4000)]
    rows = [[5000, 100, 5 if i == 0 else 3 if i == 1 else 0] for i in range(1000)]
    rows += [[5, 1, 0] for _ in empties]
    A = ad.AnnData(
        X=sp.csr_matrix(np.array(rows, float)),
        obs=pd.DataFrame(index=cells + empties),
        var=pd.DataFrame({"gene_role": ["host", "host", "viral"]}, index=["h1", "h2", "v1"]),
    )
    callA = np.r_[np.ones(1000, bool), np.zeros(4000, bool)]
    sA, _ = compute_stats(
        A,
        {"v1"},
        {"VirusX": ["v1"]},
        {"v1"},
        called_mask=callA,
        cell_calling_method="external",
        **KW,
    )
    B = ad.AnnData(
        X=sp.csr_matrix(np.array([[5.0], [3.0], [1.0]])),
        obs=pd.DataFrame(index=["cell0", "cell1", "empty0"]),
        var=pd.DataFrame(index=["v1"]),
    )
    B.uns["index_kind"] = "virus_only"
    hc = {"called": set(cells), "comparable": set(cells)}
    sB, _ = compute_stats(
        B,
        {"v1"},
        {"VirusX": ["v1"]},
        {"v1"},
        called_mask=_in_barcode_set(B.obs_names, hc["called"]),
        host_cells=hc,
        cell_calling_method="emptydrops",
        **KW,
    )
    return sA["VirusX"], sB["VirusX"]


def test_two_step_all_barcode_rates_are_na():
    a, b = _combined_and_two_step()
    assert b["pct_infected"] is None and b["viral_molecules_per_10k_est"] is None
    assert "host-filtered" in b["denominator_note"]
    # the defended fields still agree across strategies
    assert a["pct_infected_called"] == b["pct_infected_called"] == 0.2
    assert a["pct_infected_comparable"] == b["pct_infected_comparable"]


def test_combined_mode_numbers_unchanged():
    a, _ = _combined_and_two_step()
    assert a["pct_infected"] == pytest.approx(2 / 5000 * 100)  # 0.04
    assert a["viral_molecules_per_10k_est"] == pytest.approx(8 / 5_124_008 * 10_000, abs=1e-3)
    assert a.get("denominator_note", "") == ""
    line = detection._virus_summary_line("VirusX", a)
    assert "all barcodes" in line and "host-filtered" not in line


def test_two_step_tsv_summary_and_report_render_na(tmp_path):
    _, b = _combined_and_two_step()
    write_tsv_outputs({"VirusX": b}, pd.DataFrame(), str(tmp_path))
    row = next(csv.DictReader((tmp_path / "results" / "viral_summary.tsv").open(), delimiter="\t"))
    assert row["pct_infected"] == "" and row["viral_molecules_per_10k_est"] == ""
    assert "host-filtered" in row["denominator_note"]
    assert float(row["pct_infected_called"]) == 0.2
    line = detection._virus_summary_line("VirusX", b)
    assert "NA" in line and "100.00" not in line and "10000" not in line
    # the barcode count of a two-step run is the host-filtered set, and the line says so
    assert "host-filtered barcodes" in line and "all barcodes" not in line
    html = _badge_html(tmp_path, {"VirusX": b})
    assert "100.00%" not in html and "10000.00" not in html
    assert re.search(r"0\.20%", html)
