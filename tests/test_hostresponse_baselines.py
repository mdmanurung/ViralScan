"""HR-04: same-fold depth / cell-type baselines with training-only learned objects."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tests.test_hostresponse_groups import N, _design, _write_inputs
from viralscan.scripts import hostresponse as hr
from viralscan.scripts.hostresponse import (
    _balanced_split,
    _baseline_design,
    _cell_mode_baselines,
    _depth_alone_auc,
    baseline_comparison,
    evaluate_group_cv,
    make_group_folds,
    run_hostresponse,
)

ALL_BASELINES = ("depth_only", "cell_type_only", "depth_plus_cell_type")


def _cell_types():
    """Three cell types cycling inside every group (so types are shared across groups)."""
    return np.array([["T", "B", "M"][i % 3] for i in range(N)])


def _run(X, y, groups, depth, ids, **kw):
    kw.setdefault("seeds", [0, 1])
    kw.setdefault("n_folds", 3)
    kw.setdefault("use_hvg", False)
    kw.setdefault("top_depth_frac", 1.0)
    return evaluate_group_cv(X, y, groups, depth, ids, **kw)


class TestBaselineDesign:
    def test_levels_and_scaler_come_from_training_rows_only(self):
        logd_tr = np.array([0.0, 2.0, 4.0])
        logd_te = np.array([2.0, 100.0])
        ct_tr = np.array(["a", "b", "a"])
        ct_te = np.array(["c", "a"])  # "c" is unseen in training
        tr, te = _baseline_design("depth_plus_cell_type", logd_tr, logd_te, ct_tr, ct_te)
        sd = np.sqrt(8 / 3)
        np.testing.assert_allclose(tr[:, 0], (logd_tr - 2.0) / sd)
        np.testing.assert_allclose(te[:, 0], (logd_te - 2.0) / sd)  # train mean/sd, not test
        np.testing.assert_array_equal(tr[:, 1:], [[1, 0], [0, 1], [1, 0]])  # levels: a, b
        np.testing.assert_array_equal(te[:, 1:], [[0, 0], [1, 0]])  # unseen "c" -> all zero

    def test_depth_only_has_one_column_and_cell_type_only_has_no_depth(self):
        tr, _ = _baseline_design("depth_only", np.arange(4.0), np.arange(2.0))
        assert tr.shape == (4, 1)
        tr, _ = _baseline_design(
            "cell_type_only", np.arange(3.0), np.arange(2.0), ["a", "b", "a"], ["a", "b"]
        )
        assert tr.shape == (3, 2)

    def test_unknown_model_rejected(self):
        with pytest.raises(ValueError, match="unknown baseline"):
            _baseline_design("nope", np.arange(3.0), np.arange(3.0))


class TestSameFoldBaselines:
    def test_every_model_uses_identical_folds_and_cohorts(self):
        X, y, groups, depth, ids = _design()
        res = _run(X, y, groups, depth, ids, cell_type=_cell_types(), baselines=ALL_BASELINES)
        fm = res.fold_metrics
        assert set(fm["model"]) == {"expression", *ALL_BASELINES}
        cols = ["seed", "fold", "n_train_cells", "n_test_cells", "n_train_groups", "n_test_groups"]
        ref = fm[fm["model"] == "expression"][cols].reset_index(drop=True)
        for m in ALL_BASELINES:
            pd.testing.assert_frame_equal(fm[fm["model"] == m][cols].reset_index(drop=True), ref)

    def test_same_fold_holds_with_depth_filtering(self):
        X, y, groups, depth, ids = _design()
        res = _run(X, y, groups, depth, ids, top_depth_frac=0.5, baselines=("depth_only",))
        fm = res.fold_metrics
        a = fm[fm["model"] == "expression"][["seed", "fold", "n_train_cells", "n_test_cells"]]
        b = fm[fm["model"] == "depth_only"][["seed", "fold", "n_train_cells", "n_test_cells"]]
        pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))

    def test_deterministic_and_row_order_invariant(self):
        X, y, groups, depth, ids = _design()
        ct = _cell_types()
        a = _run(X, y, groups, depth, ids, cell_type=ct, baselines=ALL_BASELINES)
        perm = np.random.default_rng(2).permutation(N)
        b = _run(
            X[perm],
            y[perm],
            groups[perm],
            depth[perm],
            ids[perm],
            cell_type=ct[perm],
            baselines=ALL_BASELINES,
        )
        pd.testing.assert_frame_equal(a.fold_metrics, b.fold_metrics)

    def test_cell_type_baseline_exposes_composition_prediction(self):
        # Label is fully determined by cell type; expression is pure noise.
        X, _, groups, depth, ids = _design()
        ct = _cell_types()
        y = (ct == "T").astype(int)
        X = np.random.default_rng(1).standard_normal(X.shape).astype(np.float32)
        res = _run(X, y, groups, depth, ids, cell_type=ct, baselines=ALL_BASELINES)
        auc = {
            m: res.summary.query("model == @m and metric == 'auc'")["mean"].iloc[0]
            for m in res.summary["model"].unique()
        }
        assert auc["cell_type_only"] == pytest.approx(1.0)
        assert auc["expression"] < 0.8
        comp = baseline_comparison(auc["expression"], {m: auc[m] for m in ALL_BASELINES})
        assert comp["expression_minus_best_baseline_auc"] < 0

    def test_depth_baseline_flags_depth_confounded_label(self):
        X, y, groups, _, ids = _design()
        depth = 1000.0 + 500.0 * y + np.random.default_rng(4).uniform(0, 100, size=N)
        res = _run(X, y, groups, depth, ids, baselines=("depth_only",))
        s = res.summary
        assert s.query("model == 'depth_only' and metric == 'auc'")["mean"].iloc[0] > 0.99

    def test_unseen_test_cell_type_is_handled(self):
        # Cell type "RARE" exists only in group g0, so it is unseen whenever g0 is held out.
        X, y, groups, depth, ids = _design()
        ct = _cell_types().astype(object)
        ct[groups == "g0"] = "RARE"
        res = _run(X, y, groups, depth, ids, cell_type=ct, baselines=("cell_type_only",))
        assert res.status == "ok"
        ok = res.fold_metrics.query("model == 'cell_type_only'")
        assert (ok["status"] == "ok").all() and ok["auc"].notna().all()

    def test_argument_validation(self):
        X, y, groups, depth, ids = _design()
        with pytest.raises(ValueError, match="cell_type"):
            _run(X, y, groups, depth, ids, baselines=("cell_type_only",))
        with pytest.raises(ValueError, match="unknown baseline"):
            _run(X, y, groups, depth, ids, baselines=("bogus",))
        with pytest.raises(ValueError, match="panel_iters"):
            _run(X, y, groups, depth, ids, baselines=("panel_depth",))


class TestInFoldPanel:
    def test_panel_is_selected_on_training_cells_only(self, monkeypatch):
        X, y, groups, depth, ids = _design(n_genes=15)
        plan = make_group_folds(y, groups, ids, n_folds=3, seed=0)
        held = plan.fold == 0
        X2 = X.copy()
        X2[held] = np.random.default_rng(9).standard_normal((int(held.sum()), 15)) * 40

        calls: list = []
        orig = hr._run_stability_selection

        def spy(x, vp, n_iter, seed=42, alpha_min=0.2):
            out = orig(x, vp, n_iter, seed=seed, alpha_min=alpha_min)
            calls.append((x.shape[0], out.copy()))
            return out

        monkeypatch.setattr(hr, "_run_stability_selection", spy)
        kw = dict(
            seeds=[0], n_folds=3, baselines=("panel_depth",), panel_iters=8, stab_min_prob=0.0
        )
        a = _run(X, y, groups, depth, ids, **kw)
        calls_a, calls[:] = list(calls), []
        _run(X2, y, groups, depth, ids, **kw)
        n_held = [int((plan.fold == k).sum()) for k in range(3)]
        assert [c[0] for c in calls_a] == [N - h for h in n_held]  # training rows only
        np.testing.assert_array_equal(
            calls_a[0][1], calls[0][1]
        )  # fold-0 panel ignores held-out cells
        assert (a.fold_metrics.query("model == 'panel_depth'")["status"] == "ok").all()

    def test_no_stable_panel_is_not_estimable_not_dropped(self):
        X, y, groups, depth, ids = _design(n_genes=6)
        X = np.random.default_rng(0).standard_normal(X.shape).astype(np.float32)
        res = _run(
            X, y, groups, depth, ids, baselines=("panel_depth",), panel_iters=3, stab_min_prob=1.01
        )
        pm = res.fold_metrics.query("model == 'panel_depth'")
        assert len(pm) == 2 * 3
        assert (pm["status"] == "not_estimable").all()
        assert set(pm["reason"]) == {"no_stable_panel_in_training_fold"}
        assert res.status == "ok"  # nuisance baseline failure does not invalidate the design


class TestCellModeBaselines:
    def _data(self):
        rng = np.random.default_rng(0)
        n = 120
        depth = rng.integers(500, 5000, size=n).astype(float)
        ct = np.array([["T", "B"][i % 2] for i in range(n)])
        return np.arange(30), np.arange(30, n), depth, ct

    def test_depth_only_matches_legacy_depth_alone_exactly(self):
        pos, neg, depth, ct = self._data()
        legacy = _depth_alone_auc(pos, neg, depth, [0, 1, 10])
        new = _cell_mode_baselines(pos, neg, depth, ct, [0, 1, 10], models=("depth_only",))
        assert new["depth_only"] == pytest.approx(legacy)

    def test_splits_are_the_expression_model_splits(self):
        pos, neg, depth, _ = self._data()
        X = np.random.default_rng(1).standard_normal((len(depth), 5)).astype(np.float32)
        x_tr, y_tr, x_pos, x_neg = _balanced_split(pos, neg, depth, X, 10)
        idx = np.arange(len(depth)).reshape(-1, 1)
        i_tr, _, i_pos, i_neg = _balanced_split(pos, neg, depth, idx, 10)
        np.testing.assert_array_equal(X[i_tr[:, 0]], x_tr)
        np.testing.assert_array_equal(X[i_pos[:, 0]], x_pos)
        np.testing.assert_array_equal(X[i_neg[:, 0]], x_neg)

    def test_cell_type_models_none_without_cell_type(self):
        pos, neg, depth, _ = self._data()
        out = _cell_mode_baselines(pos, neg, depth, None, [0], models=ALL_BASELINES)
        assert out["depth_only"] is not None
        assert out["cell_type_only"] is None and out["depth_plus_cell_type"] is None


class TestBaselineComparison:
    def test_pinned_values(self):
        c = baseline_comparison(
            0.9, {"depth_only": 0.6, "cell_type_only": 0.7, "depth_plus_cell_type": 0.75}
        )
        assert c == {
            "expression_auc": 0.9,
            "depth_auc": 0.6,
            "cell_type_auc": 0.7,
            "depth_cell_type_auc": 0.75,
            "expression_minus_best_baseline_auc": pytest.approx(0.15),
        }

    def test_missing_baselines_are_none_not_zero(self):
        c = baseline_comparison(0.8, {"depth_only": float("nan")})
        assert c["depth_auc"] is None and c["expression_minus_best_baseline_auc"] is None


class TestRunHostresponseBaselines:
    def _with_celltype(self, tmp_path):
        kw = _write_inputs(tmp_path)
        h = ad.read_h5ad(kw["host_h5ad"])
        h.obs["cell_type"] = _cell_types()
        h.write_h5ad(kw["host_h5ad"])
        return kw

    @pytest.mark.parametrize("mode", ["cell", "group"])
    def test_metrics_row_has_baseline_columns(self, tmp_path, mode):
        kw = self._with_celltype(tmp_path)
        extra = {"groups_column": "sample_id", "cv_folds": 3} if mode == "group" else {}
        prov = run_hostresponse(**kw, cv_mode=mode, cell_type_column="cell_type", **extra)
        row = pd.read_csv(Path(kw["out_dir"]) / "hostresponse_metrics.csv").iloc[0]
        for col in (
            "expression_auc",
            "depth_auc",
            "cell_type_auc",
            "depth_cell_type_auc",
            "expression_minus_best_baseline_auc",
        ):
            assert col in row and np.isfinite(row[col])
        assert prov["cell_type_column"] == "cell_type"
        assert "cell_type_only" in prov["baseline_models"]

    def test_group_mode_panel_in_fold_reports_depth_adjusted_auc(self, tmp_path):
        kw = _write_inputs(tmp_path)
        kw["n_stab_iter"] = 6
        run_hostresponse(
            **kw,
            cv_mode="group",
            groups_column="sample_id",
            cv_folds=3,
            panel_in_fold=True,
            stab_min_prob=0.0,
        )
        row = pd.read_csv(Path(kw["out_dir"]) / "hostresponse_metrics.csv").iloc[0]
        assert np.isfinite(row["model_auc_depth_adjusted_mean"])
        cv = pd.read_csv(Path(kw["out_dir"]) / "VIRUS_A_cv_metrics.tsv", sep="\t")
        assert "panel_depth" in set(cv["model"])

    def test_default_run_without_cell_type_has_no_cell_type_columns(self, tmp_path):
        kw = _write_inputs(tmp_path)
        run_hostresponse(**kw)
        row = pd.read_csv(Path(kw["out_dir"]) / "hostresponse_metrics.csv").iloc[0]
        assert "cell_type_auc" not in row and "depth_alone_auc_mean" in row

    def test_missing_cell_type_column_fails_before_cleanup(self, tmp_path):
        kw = _write_inputs(tmp_path)
        out = Path(kw["out_dir"])
        out.mkdir()
        sentinel = out / "OLD_gene_weights.csv"
        sentinel.write_text("keep")
        with pytest.raises(ValueError, match="cell_type"):
            run_hostresponse(**kw, cell_type_column="cell_type")
        assert sentinel.exists()
