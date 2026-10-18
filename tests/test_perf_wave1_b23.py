"""B23 strict oracle: blocked ``run_one`` must equal the frozen unblocked baseline exactly.

The oracle is a verbatim copy of ``scripts/gene_programs.py`` at commit 5e551c9 (tests/_frozen_oracles), loaded from its own path.
Frames are compared with ``check_exact=True`` plus dtypes, attrs, log records
(level + message, so value-count tie order) and warnings.
"""

import importlib.util
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan import gene_programs as gp
from viralscan.scripts import gene_programs as cand

BASELINE = Path(__file__).resolve().parent / "_frozen_oracles" / "frozen_gene_programs_5e551c9.py"
_spec = importlib.util.spec_from_file_location("b23_baseline_sgp", BASELINE)
base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(base)

CAT = gp.load_catalogue()
EBV = "Epstein-Barr virus"
VARS = [r["gene_id_bundled"] for r in CAT if r["virus"] == EBV] + ["host"]


def _adata(n, density, seed, selected_density=None):
    rng = np.random.default_rng(seed)
    m = len(VARS)

    def mat(d):
        dense = np.where(rng.random((n, m)) < d, rng.choice([0.0, 0.5, 1.0, 2.0], (n, m)), 0.0)
        return sparse.csr_matrix(dense.astype(np.float32))

    unique = mat(density)
    sel = mat(density if selected_density is None else selected_density)
    return ad.AnnData(
        sel,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=VARS),
        layers={gp.EVIDENCE_LAYER: unique},
    )


def _run(mod, adata, caplog, block=None, **kw):
    caplog.clear()
    old = getattr(mod, "_CALL_ROW_BLOCK", None)
    if block is not None and old is not None:
        mod._CALL_ROW_BLOCK = block
    try:
        with caplog.at_level(logging.DEBUG, logger="viralscan"):
            df = mod.run_one(
                adata, [EBV, "no model"], CAT, kw.pop("min_breadth", 2), "bundled", **kw
            )
    finally:
        if old is not None:
            mod._CALL_ROW_BLOCK = old
    return df, [(r.levelname, r.getMessage()) for r in caplog.records]


@pytest.mark.parametrize("block", [1, 2, 5, 8, 9, 10_000])
@pytest.mark.parametrize(
    "n,density,sel_density",
    [(0, 0.1, None), (1, 0.3, None), (17, 0.02, 0.2), (23, 0.4, 0.0), (40, 0.6, None)],
)
@pytest.mark.parametrize(
    "kw", [dict(), dict(min_breadth=1, min_umi=1.0), dict(latent_min_breadth=2, min_umi=0.5)]
)
def test_blocked_equals_baseline(block, n, density, sel_density, kw, caplog):
    adata = _adata(n, density, 7 + n, sel_density)
    want, want_logs = _run(base, adata, caplog, **dict(kw))
    got, got_logs = _run(cand, adata, caplog, block=block, **dict(kw))
    pd.testing.assert_frame_equal(got, want, check_exact=True, check_dtype=True)
    assert [str(t) for t in got.dtypes] == [str(t) for t in want.dtypes]
    assert got.attrs == want.attrs
    assert got_logs == want_logs


@pytest.mark.parametrize(
    "order", ["PLIR", "IPLR", "RLIP", "IIPPLL", "LLIIPP", "IPLRM", "MRLIP", "IIIPPL", "PPLLLI"]
)
@pytest.mark.parametrize("block", [1, 2, 3, 100])
def test_log_tie_order_matches_baseline(order, block, caplog):
    pos = {n: i for i, n in enumerate(VARS)}
    kinds = {
        "P": [pos["EPSTEIN_HHV4_BALF5"], pos["EPSTEIN_HHV4_BBLF4"]],
        "L": [pos["EPSTEIN_HHV4_LMP-1"]],
        "R": [pos["EPSTEIN_HHV4_BZLF1"]],
        "I": [],
        "M": [pos["EPSTEIN_HHV4_BALF5"], pos["EPSTEIN_HHV4_BBLF4"], pos["EPSTEIN_HHV4_LMP-1"]],
    }
    rows = np.zeros((len(order), len(VARS)), dtype=np.float32)
    for i, ch in enumerate(order):
        rows[i, kinds[ch]] = 1.0
    adata = ad.AnnData(
        sparse.csr_matrix(rows.copy()),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(len(order))]),
        var=pd.DataFrame(index=VARS),
        layers={gp.EVIDENCE_LAYER: sparse.csr_matrix(rows)},
    )
    _, want = _run(base, adata, caplog)
    _, got = _run(cand, adata, caplog, block=block)
    assert got == want


def test_calls_are_bounded_and_cover_every_cell(monkeypatch, caplog):
    adata = _adata(25, 0.3, 3)
    sizes = []
    real = cand.call_cell_programme

    def spy(unique, markers, **kw):
        sizes.append(unique.shape[0])
        return real(unique, markers, **kw)

    monkeypatch.setattr(cand, "call_cell_programme", spy)
    monkeypatch.setattr(cand, "_CALL_ROW_BLOCK", 4)
    _run(cand, adata, caplog)
    assert sizes == [4, 4, 4, 4, 4, 4, 1] and sum(sizes) == 25


def test_zero_cells_still_validates_parameters(caplog):
    adata = _adata(0, 0.1, 1)
    with pytest.raises(ValueError, match="min_breadth must be >= 1"):
        _run(cand, adata, caplog, min_breadth=0)
    with pytest.raises(ValueError, match="min_umi must be >= 0"):
        _run(cand, adata, caplog, min_umi=-1.0)


def test_pruned_rows_still_counted_in_log(caplog):
    adata = _adata(30, 0.0, 1)  # every call is zero-evidence indeterminate
    df, logs = _run(cand, adata, caplog)
    assert [r for r in df["barcode"] if r] == []  # only the no-model marker row remains
    assert (("INFO", f"{EBV}: {{'indeterminate': 30}}")) in logs
