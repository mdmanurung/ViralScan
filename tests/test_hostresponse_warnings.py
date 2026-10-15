"""AUD3-06a: host-response warnings and labels (no computed value changes)."""

import logging

import numpy as np
import pandas as pd
import pytest

from viralscan.scripts import hostresponse as hr
from viralscan.scripts.hostresponse import (
    _detect_and_normalize,
    _run_l2_regression,
    run_hostresponse,
)

from .test_hostresponse import _make_host_adata, _make_virus_adata


def _files(tmp_path, donor: bool = True):
    host = _make_host_adata(n_obs=80, n_vars=50, raw=True)
    if donor:
        host.obs["donor"] = ["d1", "d2", "d3", "d4"] * 20
    virus = _make_virus_adata(n_obs=80)
    (tmp_path / "host.h5ad").parent.mkdir(exist_ok=True)
    host.write_h5ad(tmp_path / "host.h5ad")
    virus.write_h5ad(tmp_path / "virus.h5ad")
    (tmp_path / "analysis.txt").write_text("VIRUS_A\nVIRUS_B\n")
    return (
        str(tmp_path / "virus.h5ad"),
        str(tmp_path / "host.h5ad"),
        str(tmp_path / "analysis.txt"),
        str(tmp_path / "out"),
    )


def _run(tmp_path, donor=True, **kw):
    kw.setdefault("cv_mode", "cell")  # these tests pin the explicit legacy mode (AUD3-06b)
    v, h, a, o = _files(tmp_path, donor)
    return run_hostresponse(
        v, h, a, o, use_hvg=False, seeds=[0, 1], n_stab_iter=3, detection_threshold=1, **kw
    )


def _msgs(caplog, level=logging.WARNING):
    return [r.getMessage() for r in caplog.records if r.levelno >= level]


def test_groups_ignored_under_cell_cv_warns(tmp_path, caplog):
    caplog.set_level(logging.WARNING, logger="viralscan")
    _run(tmp_path, groups_column="donor")
    assert any("--groups" in m and "IGNORED" in m for m in _msgs(caplog))


def test_donor_like_column_under_cell_cv_warns_once(tmp_path, caplog):
    caplog.set_level(logging.WARNING, logger="viralscan")
    _run(tmp_path)
    hits = [m for m in _msgs(caplog) if "donor leakage possible" in m]
    assert len(hits) == 1 and "--cv group" in hits[0]


def test_no_donor_warning_without_donor_like_column(tmp_path, caplog):
    caplog.set_level(logging.WARNING, logger="viralscan")
    _run(tmp_path, donor=False)
    assert not [m for m in _msgs(caplog) if "donor leakage" in m or "--groups" in m]


def test_group_cv_does_not_warn_about_groups(tmp_path, caplog):
    caplog.set_level(logging.WARNING, logger="viralscan")
    _run(tmp_path, cv_mode="group", groups_column="donor", cv_folds=2)
    assert not [m for m in _msgs(caplog) if "donor leakage" in m or "IGNORED" in m]


def _l2(npos, nneg=400, seeds=(0, 1, 10, 42)):
    rng = np.random.default_rng(0)
    n = npos + nneg
    X = rng.normal(size=(n, 8)).astype(np.float32)
    y = np.zeros(n, bool)
    y[:npos] = True
    return _run_l2_regression(
        X, y, rng.integers(500, 5000, n).astype(float), list(seeds), list("abcdefgh")
    )


def test_low_support_test_folds_warn_once(caplog):
    caplog.set_level(logging.WARNING, logger="viralscan")
    weights, _ = _l2(10)
    assert weights is not None
    hits = [m for m in _msgs(caplog) if "positive" in m and "test fold" in m]
    assert len(hits) == 1  # aggregated over seeds


def test_adequate_support_does_not_warn(caplog):
    caplog.set_level(logging.WARNING, logger="viralscan")
    _l2(200)
    assert not [m for m in _msgs(caplog) if "test fold" in m]


def test_normalization_skip_is_a_warning_and_docstring_matches_code(caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    adata = _make_host_adata(raw=False)  # small, non-integer values: skipped
    _detect_and_normalize(adata)
    skipped = [r for r in caplog.records if "SKIPPING normalisation" in r.getMessage()]
    assert skipped and all(r.levelno == logging.WARNING for r in skipped)
    assert "exceeds 10" in _detect_and_normalize.__doc__
    assert "exceeds 50" not in _detect_and_normalize.__doc__


def test_stability_selection_many_positives_warns_uncalibrated(tmp_path, caplog, monkeypatch):
    caplog.set_level(logging.WARNING, logger="viralscan")
    _run(tmp_path)
    assert not [m for m in _msgs(caplog) if "uncalibrated" in m]
    caplog.clear()
    monkeypatch.setattr(hr, "STAB_SELECTION_MAX_POSITIVES", 5)
    _run(tmp_path)
    assert [m for m in _msgs(caplog) if "uncalibrated" in m and "positive cells" in m]
    assert hr.STAB_SELECTION_MAX_POSITIVES == 5


def test_stability_default_threshold_is_300():
    assert hr.STAB_SELECTION_MAX_POSITIVES == 300


# Recorded on the unmodified code (87c9e35-era main, c5b00a7) for this synthetic design.
GOLDEN = {"auc_mean": 0.5, "stab_sum": 10.666667, "n_stable": 8, "w_abs_sum": 9.158124}


def test_outputs_unchanged_for_fixed_seed(tmp_path):
    _run(tmp_path)
    out = tmp_path / "out"
    m = pd.read_csv(out / "hostresponse_metrics.csv")
    stab = pd.read_csv(out / "VIRUS_A_stability.csv")
    w = pd.read_csv(out / "VIRUS_A_gene_weights.csv")
    got = {
        "auc_mean": round(float(m["auc_mean"].iloc[0]), 6),
        "stab_sum": round(float(stab["stab_prob"].sum()), 6),
        "n_stable": int(stab["stable"].sum()),
        "w_abs_sum": round(float(w["weight_mean"].abs().sum()), 6),
    }
    # The regression weights move by ~1e-4 (relative) across BLAS builds and platforms (CI: Linux and macOS,
    # Python 3.10-3.12), so they get a relative tolerance; the rank-based and discrete values stay tight.
    assert got["w_abs_sum"] == pytest.approx(GOLDEN["w_abs_sum"], rel=1e-3)
    assert {k: v for k, v in got.items() if k != "w_abs_sum"} == pytest.approx(
        {k: v for k, v in GOLDEN.items() if k != "w_abs_sum"}, abs=1e-5
    )
    assert list(stab.columns) == ["gene", "stab_prob", "weight_mean", "weight_sd", "stable"]
