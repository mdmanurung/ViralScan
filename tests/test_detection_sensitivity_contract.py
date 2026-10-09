"""Output concordance for the canonical SensitivityRecord (SENS-CORR-03).

``sensitivity.tsv``, ``positive_control.json``, the summary text and ``report.html``
must agree, for every target, on: whether capture is measured, whether the row is an
informative negative, and why. These tests build one fixture per scenario, render all
four outputs through the same functions ``detection.main`` uses, and compare them.

The fixtures are hand-computable. Every cell holds exactly 1000 molecules, so depth is
``n_cells * 1000`` exactly; a control planted at 100 molecules and recovered at 50 gives
capture 0.5; at detection threshold 1 the LOD95 needs ``2*ln(20) = 5.9915`` true
molecules, of which ``ln(20) = 2.9957`` are expected to be observed.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from tests._exact_control import exact_control_kwargs
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection as D

LN20 = math.log(20.0)
SPIKE = "SPIKEIN_planted"
TARGET, OTHER, BELOW, CALLED = "Target virus", "Other virus", "Below-gate virus", "Called virus"
GENES = {  # virus -> its index genes
    TARGET: ["TGT_gp1"],
    OTHER: ["OTH_gp1"],
    BELOW: ["BEL_gp1"],
    CALLED: ["CAL_gp1"],
}
VAR = ["H0", SPIKE, *[g for genes in GENES.values() for g in genes]]


def _adata(n_cells: int, spike: float = 50.0) -> sc.AnnData:
    """``n_cells`` cells of exactly 1000 molecules each; depth = n_cells * 1000.

    Cell 0 carries ``spike`` molecules of the control, cell 1 twenty molecules of the
    called virus, cell 2 half a molecule of the below-gate virus (a fractional
    multimapper share, nonzero but under the gate of 1).
    """
    x = np.zeros((n_cells, len(VAR)), dtype=np.float64)
    x[:, 0] = 1000.0
    for cell, gene, n in ((0, SPIKE, spike), (1, "CAL_gp1", 20.0), (2, "BEL_gp1", 0.5)):
        x[cell, VAR.index(gene)] = n
        x[cell, 0] -= n
    assert np.all(x.sum(axis=1) == 1000.0)
    return sc.AnnData(
        sparse.csr_matrix(x),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=VAR),
    )


CALLED_STATS = {
    CALLED: {
        "viral_molecules_total_est": 20.0,
        "infected_cells": 1,
        "total_cells": 10,
        "pct_infected": 10.0,
        "viral_molecules_per_10k_est": 1.0,
    }
}


def _scoped(tmp_path, monkeypatch, virus=TARGET) -> dict:
    """Verified exact_sequence control on ``virus`` (receipt, index, identity row)."""
    return exact_control_kwargs(tmp_path, monkeypatch, SPIKE, virus)


def _cfg(**kw) -> RunConfig:
    base = {"detection_threshold": 1}
    if kw.pop("control", True):
        base |= {"positive_control_gene": SPIKE, "positive_control_expected_molecules": 100.0}
    return RunConfig(**(base | kw))


def _render(adata, cfg, tmp_path, monkeypatch, stats=CALLED_STATS):
    """Render all four outputs from one fixture and return their parsed forms."""
    monkeypatch.setattr(D, "config", cfg)
    depth = float(D._sum_axis(adata.X, 1).sum())
    _, detail = D.measure_positive_control(adata, cfg, depth=depth)
    sens = D.build_sensitivity_outputs(
        adata, stats, cfg, depth, detail, index_viruses=GENES, count_matrix=None
    )
    D.write_sensitivity_table(sens.table, str(tmp_path))
    D.write_control_report(
        detail,
        str(tmp_path),
        control=sens.control,
        certified_targets=sens.certified_targets,
        informative_negative_targets=sens.informative_targets,
    )
    D.generate_html_report(
        stats,
        pd.DataFrame(),
        None,
        None,
        {},
        [],
        str(tmp_path),
        sensitivity_df=sens.table,
        sensitivity_statement=sens.interpretation,
        control_line=sens.control_line,
    )
    tsv = pd.read_csv(tmp_path / "results" / "sensitivity.tsv", sep="\t").set_index("virus_name")
    js = json.loads((tmp_path / "results" / "positive_control.json").read_text())
    html = (tmp_path / "report.html").read_text()
    return sens, tsv, js, html


def _html_cells(html: str, virus: str) -> list[str]:
    row = re.search(
        rf'<span class="badge [^"]*">{re.escape(virus)}</span></td>(.*?)</tr>', html, re.S
    )
    assert row, f"{virus} has no row in the HTML sensitivity table"
    return [c.strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", row.group(1), re.S)]


def assert_concordant(sens, tsv, js, html) -> None:
    """The invariant: TSV, JSON, statement and HTML say the same thing per target."""
    informative_tsv = sorted(tsv.index[tsv["informative_negative"]])
    measured_tsv = sorted(tsv.index[tsv["capture_measured"]])
    # JSON lists exactly the TSV's measured and informative targets.
    assert js["certified_targets"] == measured_tsv
    assert js["informative_negative_targets"] == informative_tsv
    assert js["certifies_negatives"] is bool(informative_tsv)
    # The statement names an informative negative iff the TSV has one, and only those.
    named = re.findall(r"Informative negative for (.+?) only", sens.interpretation)
    assert sorted(named) == informative_tsv
    # The HTML carries the same statement and the same per-row verdicts.
    plain = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    assert " ".join(sens.interpretation.split())[:200] in plain
    for virus, row in tsv.iterrows():
        cells = _html_cells(html, str(virus))
        assert (cells[1] == "✓") is bool(row["detected"])
        assert cells[2].startswith("not measured") is (not bool(row["capture_measured"]))
        assert (cells[7] == "✓") is bool(row["informative_negative"])
        assert cells[4] == f"{row['lod95_molecules']:.4g}"
        assert cells[5] == f"{row['lod95_expected_observed_molecules']:.4g}"
        blockers = "" if pd.isna(row["negative_blockers"]) else row["negative_blockers"]
        assert cells[8] == blockers
        if row["detected"]:  # a call is never an informative negative
            assert not row["informative_negative"]


DEEP = 10_000  # 1e7 molecules: capture 0.5 gives LOD95 0.00599/10k, inside the 0.01 band
SHALLOW = 1_000  # 1e6 molecules: capture 0.5 gives 0.0599/10k, outside it


def test_measured_in_scope_deep_is_the_one_informative_negative(tmp_path, monkeypatch) -> None:
    sens, tsv, js, html = _render(
        _adata(DEEP), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    assert_concordant(sens, tsv, js, html)
    t = tsv.loc[TARGET]
    assert bool(t["capture_measured"]) and t["capture"] == pytest.approx(0.5)
    assert t["capture_status"] == "measured"
    # The audited defect: the row used to report 2.996 "molecules" at capture 0.5.
    assert t["lod95_molecules"] == pytest.approx(2 * LN20)  # required TRUE molecules, 5.9915
    assert t["lod95_expected_observed_molecules"] == pytest.approx(LN20)  # 2.9957 observed
    assert t["lod95_per_10k"] == pytest.approx(2 * LN20 / 1e7 * 1e4)  # 0.0059915
    assert 1e7 * t["lod95_per_10k"] / 1e4 == pytest.approx(t["lod95_molecules"])
    assert bool(t["sensitivity_eligible"]) and bool(t["informative_negative"])
    assert js["informative_negative_targets"] == [TARGET] and js["certifies_negatives"] is True
    assert "5.991 true molecules (2.996 expected observed)" in sens.interpretation
    assert f"Informative negative for {TARGET} only" in sens.interpretation
    # Mixed scope: every other undetected row is the depth-only floor, labelled so.
    for other in (OTHER, BELOW):
        o = tsv.loc[other]
        assert not bool(o["capture_measured"]) and pd.isna(o["capture"])
        assert o["capture_status"] == "control-out-of-scope"
        assert not bool(o["informative_negative"])
        assert o["negative_blockers"] == "capture-not-measured"
        assert o["lod95_molecules"] == pytest.approx(LN20)  # capture 1: required == observed
    assert "CANNOT be read as absence" in sens.interpretation  # the unqualified rows


def test_below_gate_nonzero_count_is_reported_not_inferred_absent(tmp_path, monkeypatch) -> None:
    sens, tsv, js, html = _render(
        _adata(DEEP), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    b = tsv.loc[BELOW]
    assert b["observed_molecules"] == 0.5 and not bool(b["detected"])
    assert "below the detection gate" in b["notes"]
    assert f"{BELOW} (0.5)" in sens.interpretation
    assert tsv.loc[OTHER, "observed_molecules"] == 0  # a true zero stays zero


def test_detected_positive_is_never_an_informative_negative(tmp_path, monkeypatch) -> None:
    # Control scoped to the CALLED virus: capture measured, depth ample, still a call.
    cfg = _cfg(**_scoped(tmp_path, monkeypatch, CALLED))
    sens, tsv, js, html = _render(_adata(DEEP), cfg, tmp_path, monkeypatch)
    assert_concordant(sens, tsv, js, html)
    c = tsv.loc[CALLED]
    assert bool(c["detected"]) and bool(c["capture_measured"]) and bool(c["sensitivity_eligible"])
    assert not bool(c["informative_negative"])
    assert c["negative_blockers"] == "detected"
    assert js["informative_negative_targets"] == [] and js["certifies_negatives"] is False
    assert js["certified_targets"] == [CALLED]  # capture applies, but it certifies no negative
    assert f"Called by the detection gate: {CALLED}" in sens.interpretation


def test_shallow_run_does_not_certify_even_with_measured_capture(tmp_path, monkeypatch) -> None:
    sens, tsv, js, html = _render(
        _adata(SHALLOW), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    assert_concordant(sens, tsv, js, html)
    t = tsv.loc[TARGET]
    assert bool(t["capture_measured"]) and t["lod95_per_10k"] == pytest.approx(2 * LN20 / 1e6 * 1e4)
    assert not bool(t["depth_sufficient"]) and not bool(t["sensitivity_eligible"])
    assert not bool(t["informative_negative"]) and t["negative_blockers"] == "depth-insufficient"
    assert js["certified_targets"] == [TARGET] and js["certifies_negatives"] is False
    assert "Depth is NOT sufficient" in sens.interpretation


def test_capture_halves_the_depth_that_depth_alone_would_accept(tmp_path, monkeypatch) -> None:
    """4e6 molecules: depth-only floor passes the 0.01/10k band, capture 0.5 does not."""
    sens, tsv, js, html = _render(
        _adata(4_000), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    t = tsv.loc[TARGET]
    assert bool(t["depth_only_sufficient"]) and not bool(t["depth_sufficient"])
    assert not bool(t["informative_negative"])
    assert bool(tsv.loc[OTHER, "depth_sufficient"])  # rows with no capture: floor == depth-only


def test_unconfigured_control_certifies_nothing(tmp_path, monkeypatch) -> None:
    cfg = _cfg(control=False)
    sens, tsv, js, html = _render(_adata(DEEP), cfg, tmp_path, monkeypatch)
    assert_concordant(sens, tsv, js, html)
    assert js["status"] == "not-configured" and js["certified_targets"] == []
    assert not tsv["capture_measured"].any() and not tsv["informative_negative"].any()
    assert tsv["capture"].isna().all()
    assert "No positive control was supplied" in sens.control_line
    assert "Informative negative" not in sens.interpretation
    assert "No positive control was supplied" in re.sub(r"\s+", " ", html)


def test_failed_control_certifies_nothing(tmp_path, monkeypatch) -> None:
    sens, tsv, js, html = _render(
        _adata(DEEP, spike=0.0), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    assert_concordant(sens, tsv, js, html)
    assert js["status"] == "failed" and js["certifies_negatives"] is False
    assert not tsv["capture_measured"].any()
    assert set(tsv.loc[[TARGET, OTHER, BELOW], "capture_status"]) == {"control-failed"}
    assert "status: failed" in sens.control_line


def test_over_recovered_control_certifies_nothing(tmp_path, monkeypatch) -> None:
    sens, tsv, js, html = _render(
        _adata(DEEP, spike=118.0), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    assert_concordant(sens, tsv, js, html)
    assert js["status"] == "over-recovered" and js["certifies_negatives"] is False
    assert not tsv["capture_measured"].any() and tsv["capture"].isna().all()
    assert set(tsv.loc[[TARGET, OTHER], "capture_status"]) == {"control-over-recovered"}
    # Not tightened below the depth floor either.
    assert tsv.loc[TARGET, "lod95_molecules"] == pytest.approx(LN20)


@pytest.mark.parametrize(
    "extra",
    [{}, {"positive_control_scope": "panel_mechanics"}],
    ids=["legacy-no-scope", "panel-mechanics"],
)
def test_unrelated_panel_control_certifies_nothing(tmp_path, monkeypatch, extra) -> None:
    sens, tsv, js, html = _render(_adata(DEEP), _cfg(**extra), tmp_path, monkeypatch)
    assert_concordant(sens, tsv, js, html)
    assert js["status"] == "measured" and js["capture"] == pytest.approx(0.5)  # recovery recorded
    assert js["certified_targets"] == [] and js["certifies_negatives"] is False
    assert not tsv["capture_measured"].any() and not tsv["informative_negative"].any()
    assert set(tsv["capture_status"]) == {"control-out-of-scope"}
    assert "It supplied a capture term to no row" in sens.control_line


def test_zero_depth_is_not_estimable_and_certifies_nothing(tmp_path, monkeypatch) -> None:
    adata = _adata(10)
    adata.X = sparse.csr_matrix(np.zeros(adata.shape))
    sens, tsv, js, html = _render(
        adata, _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch, stats={}
    )
    assert sens.interpretation.startswith("Detection limit: NOT ESTIMABLE")
    assert not tsv["informative_negative"].any() and not tsv["depth_sufficient"].any()
    assert (tsv["negative_blockers"].str.contains("depth-insufficient")).all()
    assert np.isinf(tsv["lod95_per_10k"]).all()
    assert js["certifies_negatives"] is False
    assert "NOT ESTIMABLE" in re.sub(r"\s+", " ", html)


def test_html_no_longer_prints_the_unmeasured_sentence_for_a_measured_run(
    tmp_path, monkeypatch
) -> None:
    """The old HTML always called negative_result_statement(depth) with no capture."""
    sens, tsv, js, html = _render(
        _adata(DEEP), _cfg(**_scoped(tmp_path, monkeypatch)), tmp_path, monkeypatch
    )
    plain = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    assert f"Informative negative for {TARGET} only" in plain
    assert "depth-only floor" in plain  # the headline is the floor, not the measured row
    assert "Positive control SPIKEIN_planted: 50.0 of 100.0 planted molecules recovered" in plain


@pytest.mark.parametrize("breakage", ["no-receipt", "index-digest-mismatch"])
def test_unverified_exact_control_certifies_no_negative(tmp_path, monkeypatch, breakage) -> None:
    """Flags are not evidence: without a verified receipt the recovery is recorded, nothing certifies."""
    kwargs = _scoped(tmp_path, monkeypatch)
    if breakage == "no-receipt":
        kwargs.pop("positive_control_receipt")
    else:  # the index was rebuilt after the control was measured
        Path(kwargs["index"]).write_bytes(b"a different binary index")
    sens, tsv, js, html = _render(_adata(DEEP), _cfg(**kwargs), tmp_path, monkeypatch)
    assert_concordant(sens, tsv, js, html)
    assert js["status"] == "measured" and js["capture"] == pytest.approx(0.5)
    assert js["verification_status"] in {"unavailable", "failed"}
    assert js["certified_targets"] == [] and js["certifies_negatives"] is False
    assert not tsv["capture_measured"].any() and not tsv["informative_negative"].any()
    assert set(tsv["capture_status"]) == {"control-out-of-scope"}
