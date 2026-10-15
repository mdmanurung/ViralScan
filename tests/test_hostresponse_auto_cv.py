"""AUD3-06b (audit P5-02): ``--cv auto`` is the default; donor leakage is the exception, not the rule."""

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.hostresponse_contract import CV_MODES, validate_evaluation
from viralscan.scripts import hostresponse as hr
from viralscan.scripts.hostresponse import run_hostresponse

N_GROUPS, PER_GROUP = 6, 10
N = N_GROUPS * PER_GROUP


def _inputs(tmp_path, y=None, **obs_cols):
    """6 groups x 10 cells; 4 positives per group; obs columns given as name -> values."""
    ids = [f"c{i:03d}" for i in range(N)]
    y = np.tile([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], N_GROUPS) if y is None else y
    rng = np.random.default_rng(5)
    raw = rng.negative_binomial(5, 0.3, size=(N, 20)).astype(np.float32)
    raw[:, 0] += (y * 30).astype(np.float32)
    obs = pd.DataFrame(index=ids)
    for name, values in obs_cols.items():
        obs[name] = values
    host = ad.AnnData(
        X=sp.csr_matrix(raw), obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(20)])
    )
    counts = np.zeros((N, 1), dtype=np.float32)
    counts[y == 1, 0] = 5
    virus = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=ids),
        var=pd.DataFrame(index=["VIRUS_A"]),
    )
    (tmp_path / "acc.txt").write_text("VIRUS_A\n")
    host.write_h5ad(tmp_path / "host.h5ad")
    virus.write_h5ad(tmp_path / "virus.h5ad")
    return dict(
        virus_h5ad=str(tmp_path / "virus.h5ad"),
        host_h5ad=str(tmp_path / "host.h5ad"),
        viral_accessions_file=str(tmp_path / "acc.txt"),
        out_dir=str(tmp_path / "out"),
        use_hvg=False,
        seeds=[0, 1],
        n_stab_iter=5,
        detection_threshold=1,
        control_mito=False,
    )


GROUPS = np.repeat([f"g{k}" for k in range(N_GROUPS)], PER_GROUP)
ALT = np.tile([f"x{k}" for k in range(N_GROUPS)], PER_GROUP)  # a second, different split


def _msgs(caplog, level=logging.WARNING):
    return [r.getMessage() for r in caplog.records if r.levelno >= level]


def _metrics(kw):
    return pd.read_csv(Path(kw["out_dir"]) / "hostresponse_metrics.csv")


def test_auto_is_the_documented_default():
    import inspect

    assert "auto" in CV_MODES
    assert inspect.signature(run_hostresponse).parameters["cv_mode"].default == "auto"
    flags = ["--virus-h5ad", "v", "--host-h5ad", "h", "--viral-accessions", "a", "--output", "o"]
    assert hr._build_parser().parse_args(flags).cv == "auto"


def test_auto_with_groups_behaves_as_explicit_group(tmp_path):
    kw = _inputs(tmp_path, sample_id=GROUPS)
    prov = run_hostresponse(**kw, cv_mode="auto", groups_column="sample_id", cv_folds=3)
    auto = _metrics(kw)
    assert prov["cv_mode"] == "group" and prov["group_column"] == "sample_id"
    assert prov["cv_mode_requested"] == "auto"
    run_hostresponse(**kw, cv_mode="group", groups_column="sample_id", cv_folds=3)
    pd.testing.assert_frame_equal(auto, _metrics(kw))


def test_auto_with_named_column_is_fail_closed(tmp_path):
    kw = _inputs(tmp_path, sample_id=GROUPS)
    with pytest.raises(ValueError, match="nope"):
        run_hostresponse(**kw, cv_mode="auto", groups_column="nope")
    # Named column but group CV infeasible: no fallback to cell (explicit group semantics).
    run_hostresponse(**kw, cv_mode="auto", groups_column="sample_id", min_groups=N_GROUPS + 1)
    assert not (Path(kw["out_dir"]) / "hostresponse_metrics.csv").exists()
    status = pd.read_csv(Path(kw["out_dir"]) / "hostresponse_status.tsv", sep="\t")
    assert status["status"].iloc[0] in {"skipped_min_groups", "group_cv_failed"}


def test_auto_single_donor_like_column_is_used_as_groups(tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    kw = _inputs(tmp_path, sample_id=GROUPS)
    prov = run_hostresponse(**kw, cv_folds=3)  # default cv_mode
    assert prov["cv_mode"] == "group" and prov["group_column"] == "sample_id"
    assert prov["cv_mode_requested"] == "auto"
    assert _metrics(kw)["cv_mode"].eq("group").all()
    info = [m for m in _msgs(caplog, logging.INFO) if "sample_id" in m and "auto" in m]
    assert info, "the auto-selected column must be named in the log"
    assert not [m for m in _msgs(caplog) if "leakage" in m]
    # Same numbers as the explicit request on that column.
    auto = _metrics(kw)
    run_hostresponse(**kw, cv_mode="group", groups_column="sample_id", cv_folds=3)
    pd.testing.assert_frame_equal(auto, _metrics(kw))


def test_auto_two_candidate_columns_does_not_guess(tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    kw = _inputs(tmp_path, sample_id=GROUPS, donor=ALT)
    prov = run_hostresponse(**kw)
    assert prov["cv_mode"] == "cell" and prov["group_column"] is None
    assert _metrics(kw)["cv_mode"].eq("cell").all()
    warned = [m for m in _msgs(caplog) if "leakage" in m and "--cv group" in m]
    assert warned and "sample_id" in warned[0] and "donor" in warned[0]


def test_auto_nothing_found_is_cell_with_loud_warning(tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    kw = _inputs(tmp_path, batch=GROUPS)  # not donor-like
    prov = run_hostresponse(**kw)
    assert prov["cv_mode"] == "cell" and prov["cv_mode_requested"] == "auto"
    warned = [m for m in _msgs(caplog) if "leakage" in m and "--cv group --groups" in m]
    assert len(warned) == 1


def test_auto_infeasible_selected_column_falls_back_with_warning(tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    kw = _inputs(tmp_path, sample_id=GROUPS)
    prov = run_hostresponse(**kw, min_groups=N_GROUPS + 1)  # too few groups for group CV
    assert prov["cv_mode"] == "cell" and prov["cv_mode_requested"] == "auto"
    assert "sample_id" in prov["cv_auto_fallback"]
    assert _metrics(kw)["cv_mode"].eq("cell").all()
    assert not list(Path(kw["out_dir"]).glob("*_cv_*.tsv")), "no group-CV leftovers"
    warned = [m for m in _msgs(caplog) if "falling back to cell" in m.lower()]
    assert warned and "sample_id" in warned[0] and "group" in warned[0]
    # Explicit --cv group on the same data stays fail-closed.
    run_hostresponse(**kw, cv_mode="group", groups_column="sample_id", min_groups=N_GROUPS + 1)
    assert not (Path(kw["out_dir"]) / "hostresponse_metrics.csv").exists()


def test_auto_falls_back_when_no_virus_has_usable_positive_groups(tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    y = (GROUPS == "g0").astype(int)  # every positive sits in ONE group: no valid group folds
    kw = _inputs(tmp_path, y=y, sample_id=GROUPS)
    prov = run_hostresponse(**kw)
    assert prov["cv_mode"] == "cell" and prov["cv_mode_requested"] == "auto"
    assert "sample_id" in prov["cv_auto_fallback"]
    assert _metrics(kw)["cv_mode"].eq("cell").all()
    assert not list(Path(kw["out_dir"]).glob("*_cv_*.tsv"))
    assert [m for m in _msgs(caplog) if "falling back to cell" in m.lower()]
    run_hostresponse(**kw, cv_mode="group", groups_column="sample_id")  # explicit: fail-closed
    assert not (Path(kw["out_dir"]) / "hostresponse_metrics.csv").exists()


def test_auto_selected_column_with_missing_values_falls_back(tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="viralscan")
    bad = GROUPS.astype(object).copy()
    bad[3] = None
    kw = _inputs(tmp_path, sample_id=bad)
    prov = run_hostresponse(**kw)
    assert prov["cv_mode"] == "cell" and "missing" in prov["cv_auto_fallback"]
    assert [m for m in _msgs(caplog) if "falling back to cell" in m.lower()]
    with pytest.raises(ValueError, match="missing"):
        run_hostresponse(**kw, cv_mode="group", groups_column="sample_id")


def test_explicit_cell_ignores_columns_and_is_unchanged(tmp_path):
    kw = _inputs(tmp_path, sample_id=GROUPS)
    prov = run_hostresponse(**kw, cv_mode="cell")
    assert prov["cv_mode"] == "cell" and prov["cv_mode_requested"] == "cell"
    a = _metrics(kw)
    run_hostresponse(**kw, cv_mode="cell", groups_column="sample_id")
    pd.testing.assert_frame_equal(a, _metrics(kw))


def test_ebv_matched_style_call_without_cv_mode_is_unchanged(tmp_path):
    """One sample, no donor column, no cv_mode (as scripts/hostresponse_ebv_matched.py)."""
    kw = _inputs(tmp_path)
    prov = run_hostresponse(**kw)
    out = Path(kw["out_dir"])
    assert prov["cv_mode"] == "cell"
    auto = {f.name: f.read_bytes() for f in out.glob("*.csv")}
    run_hostresponse(**kw, cv_mode="cell")
    assert auto and auto == {f.name: f.read_bytes() for f in out.glob("*.csv")}
    assert not list(out.glob("*_cv_*.tsv"))


def test_validate_evaluation_accepts_auto():
    validate_evaluation(
        cv_mode="auto",
        groups_column=None,
        permutations=0,
        permutation_unit=None,
        permutation_block=None,
    )
    with pytest.raises(ValueError, match="cv_mode must be one of"):
        validate_evaluation(
            cv_mode="bogus",
            groups_column=None,
            permutations=0,
            permutation_unit=None,
            permutation_block=None,
        )
