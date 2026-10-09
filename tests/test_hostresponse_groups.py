"""HR-01: group-disjoint host-response evaluation (contract tests, synthetic data only)."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.scripts import hostresponse as hr
from viralscan.scripts.hostresponse import (
    FoldDesignError,
    evaluate_group_cv,
    make_group_folds,
    run_hostresponse,
    summarize_fold_metrics,
    validate_obs_column,
)

N_GROUPS = 6
PER_GROUP = 10
N = N_GROUPS * PER_GROUP


def _design(seed=0, n_genes=12):
    """6 groups x 10 cells; 4 positives per group; gene 0 marks positives."""
    rng = np.random.default_rng(seed)
    ids = np.array([f"c{i:03d}" for i in range(N)])
    groups = np.repeat([f"g{k}" for k in range(N_GROUPS)], PER_GROUP)
    y = np.tile([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], N_GROUPS)
    X = rng.standard_normal((N, n_genes)).astype(np.float32)
    X[:, 0] += 3.0 * y
    depth = rng.integers(1000, 5000, size=N).astype(float)
    return X, y, groups, depth, ids


class TestFoldPlan:
    def test_no_group_in_train_and_test(self):
        _, y, groups, _, ids = _design()
        for seed in (0, 1, 42):
            plan = make_group_folds(y, groups, ids, n_folds=3, seed=seed)
            assert plan.n_folds == 3
            for g in set(groups):
                assert len(set(plan.fold[groups == g])) == 1, f"group {g} spans folds"
            for k in range(plan.n_folds):
                test_groups = set(groups[plan.fold == k])
                train_groups = set(groups[plan.fold != k])
                assert test_groups and not (test_groups & train_groups)

    def test_every_fold_has_both_classes_on_both_sides(self):
        _, y, groups, _, ids = _design()
        plan = make_group_folds(y, groups, ids, n_folds=3, seed=0)
        for k in range(3):
            assert set(y[plan.fold == k]) == {0, 1}
            assert set(y[plan.fold != k]) == {0, 1}

    def test_row_order_invariant_membership(self):
        _, y, groups, _, ids = _design()
        base = make_group_folds(y, groups, ids, n_folds=3, seed=7)
        base_map = dict(zip(ids, base.fold))
        perm = np.random.default_rng(3).permutation(N)
        shuffled = make_group_folds(y[perm], groups[perm], ids[perm], n_folds=3, seed=7)
        assert dict(zip(ids[perm], shuffled.fold)) == base_map

    def test_deterministic_per_seed(self):
        _, y, groups, _, ids = _design()
        a = make_group_folds(y, groups, ids, n_folds=3, seed=5)
        b = make_group_folds(y, groups, ids, n_folds=3, seed=5)
        assert (a.fold == b.fold).all()

    def test_n_splits_capped_by_group_count(self):
        _, y, groups, _, ids = _design()
        plan = make_group_folds(y, groups, ids, n_folds=50, seed=0)
        assert plan.n_folds == N_GROUPS

    def test_single_group_fails_with_reason(self):
        _, y, _, _, ids = _design()
        with pytest.raises(FoldDesignError) as exc:
            make_group_folds(y, np.array(["only"] * N), ids, n_folds=5)
        assert exc.value.reason == "insufficient_groups"

    def test_single_class_fails_with_reason(self):
        _, _, groups, _, ids = _design()
        with pytest.raises(FoldDesignError) as exc:
            make_group_folds(np.zeros(N, dtype=int), groups, ids, n_folds=3)
        assert exc.value.reason == "single_class"

    def test_class_confined_to_one_group_is_invalid_not_dropped(self):
        # Positives live only in g0: whenever g0 is held out the training side has none.
        _, _, groups, _, ids = _design()
        y = (groups == "g0").astype(int) * np.tile([1, 1, 0, 0, 0, 0, 0, 0, 0, 0], N_GROUPS)
        with pytest.raises(FoldDesignError) as exc:
            make_group_folds(y, groups, ids, n_folds=3)
        assert exc.value.reason == "invalid_folds"

    def test_duplicate_cell_ids_fail(self):
        _, y, groups, _, ids = _design()
        ids = ids.copy()
        ids[1] = ids[0]
        with pytest.raises(FoldDesignError) as exc:
            make_group_folds(y, groups, ids)
        assert exc.value.reason == "duplicate_cell_ids"

    def test_groupkfold_fallback_is_recorded(self, monkeypatch):
        import sklearn.model_selection as ms

        class _Broken:
            def __init__(self, *a, **k):
                pass

            def split(self, *a, **k):
                raise ValueError("cannot stratify")

        monkeypatch.setattr(ms, "StratifiedGroupKFold", _Broken)
        _, y, groups, _, ids = _design()
        plan = make_group_folds(y, groups, ids, n_folds=3, seed=0)
        assert plan.method == "GroupKFold"
        assert plan.stratified is False
        assert "StratifiedGroupKFold unavailable" in plan.note


class TestEvaluateGroupCV:
    def _run(self, X, y, groups, depth, ids, **kw):
        kw.setdefault("seeds", [0, 1])
        kw.setdefault("n_folds", 3)
        kw.setdefault("use_hvg", False)
        kw.setdefault("top_depth_frac", 1.0)
        return evaluate_group_cv(X, y, groups, depth, ids, **kw)

    def test_fold_table_reports_group_counts_and_metrics(self):
        X, y, groups, depth, ids = _design()
        res = self._run(X, y, groups, depth, ids)
        assert res.status == "ok"
        ok = res.fold_metrics[res.fold_metrics["status"] == "ok"]
        assert len(ok) == 2 * 3  # seeds x folds
        # 6 groups, 3 folds -> 2 test groups and 4 train groups per fold, 60 cells in total.
        assert set(ok["n_test_groups"]) == {2}
        assert set(ok["n_train_groups"]) == {4}
        assert ((ok["n_train_cells"] + ok["n_test_cells"]) == N).all()
        assert ok["auc"].between(0, 1).all()
        # Planted signal in gene 0 must generalise to unseen groups.
        assert res.statistic is not None and res.statistic["auc"] > 0.9

    def test_no_leakage_in_saved_membership(self):
        X, y, groups, depth, ids = _design()
        res = self._run(X, y, groups, depth, ids)
        for _, sub in res.folds.groupby("seed"):
            assert (sub.groupby("group")["fold"].nunique() == 1).all()
            assert len(sub) == N

    def test_row_order_does_not_change_results(self):
        X, y, groups, depth, ids = _design()
        a = self._run(X, y, groups, depth, ids, use_hvg=True)
        perm = np.random.default_rng(11).permutation(N)
        b = self._run(X[perm], y[perm], groups[perm], depth[perm], ids[perm], use_hvg=True)
        pd.testing.assert_frame_equal(a.fold_metrics, b.fold_metrics)
        pd.testing.assert_frame_equal(a.folds, b.folds)
        assert a.fold_features == b.fold_features

    def test_deterministic_under_seed(self):
        X, y, groups, depth, ids = _design()
        a = self._run(X, y, groups, depth, ids)
        b = self._run(X, y, groups, depth, ids)
        pd.testing.assert_frame_equal(a.fold_metrics, b.fold_metrics)

    def test_hvg_selection_never_sees_held_out_cells(self, monkeypatch):
        """HVG is fitted on training cells only: rewriting held-out cells must not move it."""
        X, y, groups, depth, ids = _design(n_genes=40)
        plan = make_group_folds(y, groups, ids, n_folds=3, seed=0)
        held = plan.fold == 0
        X2 = X.copy()
        X2[held] = np.random.default_rng(99).standard_normal((int(held.sum()), 40)) * 50
        kw = dict(seeds=[0], n_folds=3, use_hvg=True, top_depth_frac=1.0)
        a = evaluate_group_cv(X, y, groups, depth, ids, **kw)
        b = evaluate_group_cv(X2, y, groups, depth, ids, **kw)
        assert a.fold_features[(0, 0)] == b.fold_features[(0, 0)]

        seen = []
        orig = hr._hvg_mask

        def spy(x_train):
            seen.append(x_train.shape[0])
            return orig(x_train)

        monkeypatch.setattr(hr, "_hvg_mask", spy)
        evaluate_group_cv(X, y, groups, depth, ids, **kw)
        n_test = [int((plan.fold == k).sum()) for k in range(3)]
        assert seen == [N - t for t in n_test]

    def test_invalid_fold_classes_after_filtering_fail_explicitly(self):
        X, y, groups, depth, ids = _design()
        depth = np.where(y == 1, 9000.0, 10.0) + np.arange(N)  # classes fully separated in depth
        res = self._run(X, y, groups, depth, ids, depth_match=True)
        assert res.status == "failed"
        assert res.statistic is None
        assert {f["reason"] for f in res.failures} == {"invalid_fold_classes"}
        assert (res.fold_metrics["status"] == "failed").all()

    def test_insufficient_groups_recorded_as_failure(self):
        X, y, _, depth, ids = _design()
        res = self._run(X, y, np.array(["one"] * N), depth, ids)
        assert res.status == "failed"
        assert res.failures[0]["reason"] == "insufficient_groups"

    def test_single_class_group_auc_is_not_estimable(self):
        X, y, groups, depth, ids = _design()
        y = y.copy()
        y[groups == "g0"] = 0  # g0 has no positives; others unchanged
        res = self._run(X, y, groups, depth, ids)
        g0 = res.group_auc[res.group_auc["group"] == "g0"]
        assert (g0["status"] == "not_estimable").all() and g0["auc"].isna().all()
        assert res.statistic["n_groups_auc_estimable"] == N_GROUPS - 1

    def test_summary_columns_and_pinned_values(self):
        fm = pd.DataFrame(
            {
                "model": ["expression"] * 4,
                "seed": [0] * 4,
                "fold": [0, 1, 2, 3],
                "status": ["ok", "ok", "ok", "failed"],
                "auc": [0.6, 0.8, 1.0, np.nan],
                "sensitivity": [1.0] * 4,
                "specificity": [1.0] * 4,
                "balanced_accuracy": [1.0] * 4,
                "mcc": [1.0] * 4,
            }
        )
        s = summarize_fold_metrics(fm)
        r = s[(s["model"] == "expression") & (s["metric"] == "auc")].iloc[0]
        assert list(s.columns) == [
            "model",
            "metric",
            "mean",
            "sd",
            "median",
            "min",
            "max",
            "n_valid_folds",
        ]
        assert r["mean"] == pytest.approx(0.8)
        assert r["sd"] == pytest.approx(np.std([0.6, 0.8, 1.0]))
        assert (r["median"], r["min"], r["max"], r["n_valid_folds"]) == (0.8, 0.6, 1.0, 3)

    def test_provenance_is_json_serialisable(self):
        import json

        X, y, groups, depth, ids = _design()
        prov = self._run(X, y, groups, depth, ids).provenance()
        json.dumps(prov)
        assert set(prov["fold_digests"]) == {0, 1}


class TestMetadataValidation:
    def test_missing_column(self):
        obs = pd.DataFrame({"a": ["x"]}, index=["c1"])
        with pytest.raises(ValueError, match="not found"):
            validate_obs_column(obs, "sample_id")

    def test_missing_values(self):
        obs = pd.DataFrame({"s": ["a", None, "b"]}, index=["c1", "c2", "c3"])
        with pytest.raises(ValueError, match="missing"):
            validate_obs_column(obs, "s")

    def test_blank_values(self):
        obs = pd.DataFrame({"s": ["a", "  ", "b"]}, index=["c1", "c2", "c3"])
        with pytest.raises(ValueError, match="blank"):
            validate_obs_column(obs, "s")

    def test_returns_strings(self):
        obs = pd.DataFrame({"s": [1, 2]}, index=["c1", "c2"])
        assert validate_obs_column(obs, "s").tolist() == ["1", "2"]


# ── end-to-end through run_hostresponse ───────────────────────────────────────


def _write_inputs(tmp_path, *, group_col="sample_id", dup=False, counts_layer=False):
    X, y, groups, _, ids = _design(n_genes=20)
    rng = np.random.default_rng(5)
    raw = rng.negative_binomial(5, 0.3, size=X.shape).astype(np.float32)
    raw[:, 0] += (y * 30).astype(np.float32)
    obs = pd.DataFrame(index=list(ids))
    if group_col:
        obs[group_col] = groups
    host = ad.AnnData(
        X=sp.csr_matrix(raw), obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(20)])
    )
    if counts_layer:
        host.layers["counts"] = sp.csr_matrix(raw)
    if dup:
        host = host[list(ids[:-1]) + [ids[0]]].copy()
    counts = np.zeros((N, 1), dtype=np.float32)
    counts[y == 1, 0] = 5
    virus = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=list(ids)),
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


class TestRunHostresponseGroupMode:
    def test_group_mode_writes_distinct_group_outputs(self, tmp_path):
        kw = _write_inputs(tmp_path)
        prov = run_hostresponse(**kw, cv_mode="group", groups_column="sample_id", cv_folds=3)
        out = Path(kw["out_dir"])
        assert prov["cv_mode"] == "group" and prov["group_column"] == "sample_id"
        assert prov["n_folds"] == 3 and prov["seeds"] == [0, 1] and prov["n_groups"] == N_GROUPS
        assert len(prov["inputs"]["host_h5ad"]["sha256"]) == 64
        assert prov["raw_depth_source"] == "X"
        row = pd.read_csv(out / "hostresponse_metrics.csv").iloc[0]
        assert row["cv_mode"] == "group" and row["cv_n_groups"] == N_GROUPS
        folds = pd.read_csv(out / "VIRUS_A_cv_folds.tsv", sep="\t")
        for _, sub in folds.groupby("seed"):
            assert (sub.groupby("group")["fold"].nunique() == 1).all()
        assert (out / "VIRUS_A_cv_metrics.tsv").is_file()
        assert (out / "VIRUS_A_cv_summary.tsv").is_file()
        assert (out / "hostresponse_status.tsv").is_file()

    def test_cell_mode_is_default_and_unchanged(self, tmp_path):
        kw = _write_inputs(tmp_path)
        run_hostresponse(**kw)
        out = Path(kw["out_dir"])
        a = pd.read_csv(out / "hostresponse_metrics.csv")
        w_a = pd.read_csv(out / "VIRUS_A_gene_weights.csv")
        assert a.loc[0, "cv_mode"] == "cell"
        assert not list(out.glob("*_cv_*.tsv"))
        run_hostresponse(**kw, cv_mode="cell", groups_column="sample_id")
        pd.testing.assert_frame_equal(a, pd.read_csv(out / "hostresponse_metrics.csv"))
        pd.testing.assert_frame_equal(w_a, pd.read_csv(out / "VIRUS_A_gene_weights.csv"))

    def test_missing_group_column_fails_without_fallback_or_cleanup(self, tmp_path):
        kw = _write_inputs(tmp_path, group_col=None)
        out = Path(kw["out_dir"])
        out.mkdir()
        sentinel = out / "OLD_gene_weights.csv"
        sentinel.write_text("keep")
        with pytest.raises(ValueError, match="sample_id"):
            run_hostresponse(**kw, cv_mode="group", groups_column="sample_id")
        assert sentinel.exists(), "validation must run before stale outputs are cleared"
        assert not (out / "hostresponse_metrics.csv").exists()

    def test_group_mode_requires_group_column_argument(self, tmp_path):
        kw = _write_inputs(tmp_path)
        with pytest.raises(ValueError, match="groups_column"):
            run_hostresponse(**kw, cv_mode="group")

    def test_insufficient_groups_reports_status_not_silent_cell_cv(self, tmp_path):
        kw = _write_inputs(tmp_path)
        h = ad.read_h5ad(kw["host_h5ad"])
        h.obs["sample_id"] = "one"
        h.write_h5ad(kw["host_h5ad"])
        prov = run_hostresponse(**kw, cv_mode="group", groups_column="sample_id")
        out = Path(kw["out_dir"])
        status = pd.read_csv(out / "hostresponse_status.tsv", sep="\t")
        assert status.loc[0, "status"] == "group_cv_failed"
        assert "insufficient_groups" in status.loc[0, "reason"]
        assert prov["status"] == "no_eligible_virus"
        assert not (out / "hostresponse_metrics.csv").exists()

    def test_duplicate_barcodes_fail_before_cleanup(self, tmp_path):
        kw = _write_inputs(tmp_path, dup=True)
        out = Path(kw["out_dir"])
        out.mkdir()
        sentinel = out / "OLD_stability.csv"
        sentinel.write_text("keep")
        with pytest.raises(ValueError, match="duplicate barcodes"):
            run_hostresponse(**kw)
        assert sentinel.exists()

    def test_counts_layer_is_the_raw_depth_source(self, tmp_path):
        kw = _write_inputs(tmp_path, counts_layer=True)
        prov = run_hostresponse(**kw)
        assert prov["raw_depth_source"] == "layers['counts']"

    def test_no_eligible_run_removes_stale_metrics(self, tmp_path):
        kw = _write_inputs(tmp_path)
        out = Path(kw["out_dir"])
        out.mkdir()
        (out / "hostresponse_metrics.csv").write_text("stale")
        run_hostresponse(**{**kw, "detection_threshold": 10_000})
        assert not (out / "hostresponse_metrics.csv").exists()
        assert (out / "hostresponse_status.tsv").is_file()


class TestRawDepthResolution:
    def test_negative_values_rejected(self):
        a = ad.AnnData(X=np.array([[1.0, -1.0], [2.0, 3.0]], dtype=np.float32))
        with pytest.raises(ValueError, match="negative"):
            hr._resolve_raw_depth(a)

    def test_log1p_without_counts_is_labelled_not_raw(self):
        a = ad.AnnData(X=np.ones((3, 2), dtype=np.float32))
        a.uns["log1p"] = {}
        assert hr._resolve_raw_depth(a)[1] == "X_log1p_not_raw"
