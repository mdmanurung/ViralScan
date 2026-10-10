"""HR-02: declared structured permutation null (contract tests, synthetic data only)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tests.test_hostresponse_groups import N_GROUPS, PER_GROUP, N, _design, _write_inputs
from viralscan.scripts import hostresponse as hr
from viralscan.scripts.hostresponse import (
    FoldDesignError,
    empirical_p_value,
    group_cv_permutation_null,
    permute_labels,
    run_hostresponse,
    run_structured_null,
    validate_null_spec,
)


def _varied_design():
    """Like _design but group k has k+2 positives, so label vectors differ per group."""
    X, _, groups, depth, ids = _design()
    y = np.zeros(N, dtype=int)
    for k in range(N_GROUPS):
        y[k * PER_GROUP : k * PER_GROUP + k + 2] = 1
    X = X.copy()
    X[:, 0] += 3.0 * y
    return X, y, groups, depth, ids


def _stub(aucs):
    """Procedure stub returning a fixed AUC sequence (observed first)."""
    it = iter(aucs)
    calls: list = []

    def procedure(y):
        calls.append(np.asarray(y).copy())
        a = next(it)
        if a is None:
            raise FoldDesignError("invalid_folds", "stub failure")
        return {"auc": a, "balanced_accuracy": 0.5, "mcc": 0.0, "n_groups": 6, "n_cells": 60}

    return procedure, calls


class TestEmpiricalP:
    def test_formula_is_pinned(self):
        # 2 of 3 null AUCs are >= 0.6 (ties count): (1 + 2) / (1 + 3)
        assert empirical_p_value([0.5, 0.6, 0.7], 0.6, 3) == pytest.approx(0.75)

    def test_never_zero(self):
        assert empirical_p_value([0.1] * 9, 0.9, 9) == pytest.approx(1 / 10)

    def test_all_null_at_least_observed(self):
        assert empirical_p_value([0.9, 0.9], 0.5, 2) == 1.0


class TestSpecValidation:
    def test_disabled_by_default_accepts_anything(self):
        validate_null_spec(None, False, 0)

    def test_unit_must_be_declared(self):
        with pytest.raises(ValueError, match="declared explicitly"):
            validate_null_spec(None, False, 10)
        with pytest.raises(ValueError, match="declared explicitly"):
            validate_null_spec("whatever", True, 10)

    def test_block_unit_needs_block_and_group_unit_forbids_it(self):
        with pytest.raises(ValueError, match="requires a permutation block"):
            validate_null_spec("cell_within_block", False, 10)
        with pytest.raises(ValueError, match="only applies"):
            validate_null_spec("group", True, 10)

    def test_negative_rejected(self):
        with pytest.raises(ValueError, match=">= 0"):
            validate_null_spec("group", False, -1)


class TestPermuteWithinBlock:
    def test_preserves_block_structure_and_prevalence(self):
        _, y, groups, _, ids = _varied_design()
        yp, why = permute_labels(y, ids, groups, "cell_within_block", block=groups, seed=1)
        assert why is None and (yp != y).any()
        for g in set(groups):
            assert yp[groups == g].sum() == y[groups == g].sum()  # per-block prevalence fixed
        assert yp.sum() == y.sum()

    def test_deterministic_per_seed_and_seed_sensitive(self):
        _, y, groups, _, ids = _varied_design()
        a = permute_labels(y, ids, groups, "cell_within_block", block=groups, seed=5)[0]
        b = permute_labels(y, ids, groups, "cell_within_block", block=groups, seed=5)[0]
        c = permute_labels(y, ids, groups, "cell_within_block", block=groups, seed=6)[0]
        assert (a == b).all() and (a != c).any()

    def test_row_order_invariant(self):
        _, y, groups, _, ids = _varied_design()
        base = permute_labels(y, ids, groups, "cell_within_block", block=groups, seed=3)[0]
        perm = np.random.default_rng(0).permutation(N)
        shuf = permute_labels(
            y[perm], ids[perm], groups[perm], "cell_within_block", block=groups[perm], seed=3
        )[0]
        assert dict(zip(ids, base)) == dict(zip(ids[perm], shuf))

    def test_constant_blocks_stay_and_all_constant_is_not_estimable(self):
        _, y, groups, _, ids = _varied_design()
        y2 = y.copy()
        y2[groups == "g0"] = 0  # constant block cannot change
        yp, _ = permute_labels(y2, ids, groups, "cell_within_block", block=groups, seed=2)
        assert (yp[groups == "g0"] == 0).all()
        yp, why = permute_labels(
            np.zeros(N, dtype=int), ids, groups, "cell_within_block", block=groups
        )
        assert yp is None and why == "all_blocks_constant"

    def test_missing_block_rejected(self):
        _, y, groups, _, ids = _varied_design()
        blk = groups.astype(object)
        blk[3] = None
        with pytest.raises(ValueError, match="no missing"):
            permute_labels(y, ids, groups, "cell_within_block", block=blk)
        with pytest.raises(ValueError, match="requires block"):
            permute_labels(y, ids, groups, "cell_within_block", block=None)


class TestPermuteGroups:
    def test_exchanges_whole_label_vectors_between_equal_size_groups(self):
        _, y, groups, _, ids = _varied_design()
        yp, why = permute_labels(y, ids, groups, "group", seed=4)
        assert why is None and (yp != y).any()
        originals = {tuple(y[groups == g]) for g in set(groups)}
        for g in set(groups):
            assert tuple(yp[groups == g]) in originals  # each group got a donor's vector
        assert sorted(tuple(yp[groups == g]) for g in set(groups)) == sorted(originals)

    def test_deterministic_and_row_order_invariant(self):
        _, y, groups, _, ids = _varied_design()
        a = permute_labels(y, ids, groups, "group", seed=4)[0]
        perm = np.random.default_rng(1).permutation(N)
        b = permute_labels(y[perm], ids[perm], groups[perm], "group", seed=4)[0]
        assert dict(zip(ids, a)) == dict(zip(ids[perm], b))

    def test_identical_label_vectors_are_not_exchangeable(self):
        _, y, groups, _, ids = _design()  # every group has the same label vector
        yp, why = permute_labels(y, ids, groups, "group", seed=0)
        assert yp is None and why == "no_exchangeable_groups"

    def test_groups_of_unique_size_are_not_exchangeable(self):
        _, y, groups, _, ids = _varied_design()
        keep = np.ones(N, dtype=bool)
        for k in range(N_GROUPS):  # shrink group k to 4 + k cells -> all sizes distinct
            keep[np.where(groups == f"g{k}")[0][4 + k :]] = False
        yp, why = permute_labels(y[keep], ids[keep], groups[keep], "group", seed=0)
        assert yp is None and why == "no_exchangeable_groups"


class TestRunStructuredNull:
    kw = dict(unit="cell_within_block", n_permutations=4, seed=11)

    def _run(self, procedure, **over):
        _, y, groups, _, ids = _varied_design()
        return run_structured_null(
            procedure, y, ids, groups, block=groups, virus="V", **{**self.kw, **over}
        )

    def test_disabled_by_default_never_calls_procedure(self):
        proc, calls = _stub([])
        res = self._run(proc, n_permutations=0)
        assert res.summary["status"] == "disabled" and not calls and res.replicates.empty

    def test_complete_null_reports_pinned_p_and_columns(self):
        proc, calls = _stub([0.6, 0.5, 0.6, 0.7, 0.4])  # observed, then 4 null AUCs
        res = self._run(proc)
        s = res.summary
        assert s["status"] == "ok" and s["n_valid"] == 4 and s["n_failed"] == 0
        assert s["observed_group_auc"] == 0.6
        assert s["empirical_p_auc"] == pytest.approx((1 + 2) / (1 + 4))
        assert s["null_auc_mean"] == pytest.approx(0.55)
        assert s["null_auc_sd"] == pytest.approx(np.std([0.5, 0.6, 0.7, 0.4]))
        assert s["observed_minus_null_auc"] == pytest.approx(0.05)
        assert s["permutation_unit"] == "cell_within_block" and s["assumption"]
        assert list(res.replicates.columns)[:9] == [
            "virus", "permutation_id", "seed", "permutation_unit", "auc",
            "balanced_accuracy", "mcc", "n_groups", "n_cells",
        ]  # fmt: skip
        assert list(res.replicates["permutation_id"]) == [1, 2, 3, 4]
        assert len(calls) == 5

    def test_every_replicate_receives_its_own_deterministic_permuted_labels(self):
        proc, calls = _stub([0.6, 0.5, 0.5, 0.5, 0.5])
        res = self._run(proc)
        _, y, groups, _, ids = _varied_design()
        assert (calls[0] == y).all()  # observed
        for r in range(1, 5):
            seed = int(res.replicates.loc[r - 1, "seed"])
            expect = permute_labels(y, ids, groups, "cell_within_block", block=groups, seed=seed)[0]
            assert (calls[r] == expect).all() and (calls[r] != y).any()
        proc2, calls2 = _stub([0.6, 0.5, 0.5, 0.5, 0.5])
        res2 = self._run(proc2)
        pd.testing.assert_frame_equal(res.replicates, res2.replicates)
        assert all((a == b).all() for a, b in zip(calls, calls2))

    def test_failed_replicate_is_retained_and_null_not_estimable(self):
        proc, _ = _stub([0.6, 0.5, None, 0.7, 0.4])
        res = self._run(proc)
        s = res.summary
        assert s["status"] == "not_estimable" and s["empirical_p_auc"] is None
        assert s["null_auc_mean"] is None and s["observed_minus_null_auc"] is None
        assert (s["n_requested"], s["n_valid"], s["n_failed"]) == (4, 3, 1)
        assert len(res.replicates) == 4  # nothing dropped
        failed = res.replicates[res.replicates["status"] == "failed"]
        assert list(failed["permutation_id"]) == [2] and "stub failure" in failed["reason"].iloc[0]

    def test_invalid_observed_is_not_estimable(self):
        proc, calls = _stub([None])
        res = self._run(proc)
        assert res.summary["status"] == "not_estimable"
        assert "observed_invalid" in res.summary["reason"] and len(calls) == 1

    def test_non_exchangeable_structure_is_a_failed_replicate(self):
        proc, _ = _stub([0.6] + [0.5] * 4)
        _, y, groups, _, ids = _varied_design()
        res = run_structured_null(
            proc,
            np.zeros(N, dtype=int),
            ids,
            groups,
            block=groups,
            observed={"auc": 0.6, **{}},
            **self.kw,
        )
        assert res.summary["status"] == "not_estimable"
        assert set(res.replicates["reason"]) == {"all_blocks_constant"}

    def test_supplied_observed_skips_recompute(self):
        proc, calls = _stub([0.5] * 4)
        res = self._run(proc, observed={"auc": 0.6})
        assert len(calls) == 4 and res.summary["status"] == "ok"

    def test_undeclared_unit_rejected(self):
        proc, _ = _stub([])
        with pytest.raises(ValueError, match="declared explicitly"):
            self._run(proc, unit=None)


class TestFullRefitPerReplicate:
    def _kw(self):
        return dict(seeds=[0], n_folds=3, use_hvg=True, top_depth_frac=0.5)

    def test_hvg_and_cohorts_are_recomputed_from_permuted_labels(self, monkeypatch):
        X, y, groups, depth, ids = _varied_design()
        hvg_rows: list = []
        side_labels: list = []
        orig_hvg, orig_side = hr._hvg_mask, hr._side_cohort

        def hvg_spy(x_train):
            hvg_rows.append(x_train.shape[0])
            return orig_hvg(x_train)

        def side_spy(idx, y_in, *a, **k):
            side_labels.append(np.asarray(y_in).copy())
            return orig_side(idx, y_in, *a, **k)

        monkeypatch.setattr(hr, "_hvg_mask", hvg_spy)
        monkeypatch.setattr(hr, "_side_cohort", side_spy)
        res = group_cv_permutation_null(
            X, y, groups, depth, ids, unit="cell_within_block", block=groups, n_permutations=2, seed=3,
            **self._kw(),
        )  # fmt: skip
        # observed (3 folds) + 2 replicates x 3 folds: HVG is fitted again for every replicate fold
        assert len(hvg_rows) == 3 * 3
        assert all(r < N for r in hvg_rows)  # always on training rows only
        # cohort selection saw the permuted label vector, not the observed one
        distinct = {tuple(v) for v in side_labels}
        assert tuple(y) in distinct and len(distinct) == 3
        assert res.summary["status"] == "ok"

    def test_observed_reuse_is_equivalent_to_recompute(self):
        X, y, groups, depth, ids = _varied_design()
        a = group_cv_permutation_null(
            X, y, groups, depth, ids, unit="cell_within_block", block=groups, n_permutations=2,
            seed=3, **self._kw(),
        )  # fmt: skip
        obs = {
            "auc": a.summary["observed_group_auc"], "balanced_accuracy": 0, "mcc": 0,
            "n_groups": N_GROUPS, "n_cells": N,
        }  # fmt: skip
        b = group_cv_permutation_null(
            X, y, groups, depth, ids, unit="cell_within_block", block=groups, n_permutations=2,
            seed=3, observed=obs, **self._kw(),
        )  # fmt: skip
        pd.testing.assert_frame_equal(a.replicates, b.replicates)

    def test_invalid_replicate_design_is_a_failed_replicate(self):
        # Depth separates the classes after permutation only by chance; force invalid folds
        # by requiring more class cells per side than any fold can offer.
        X, y, groups, depth, ids = _varied_design()
        res = group_cv_permutation_null(
            X, y, groups, depth, ids, unit="cell_within_block", block=groups, n_permutations=2,
            seed=1, observed={"auc": 0.9}, seeds=[0], n_folds=3, use_hvg=False, min_class_cells=999,
        )  # fmt: skip
        assert res.summary["status"] == "not_estimable"
        assert (res.replicates["status"] == "failed").all()
        assert res.summary["empirical_p_auc"] is None

    def test_signal_gives_small_p_and_noise_does_not(self):
        X, y, groups, depth, ids = _varied_design()
        kw = dict(unit="cell_within_block", block=groups, n_permutations=19, seed=2,
                  seeds=[0], n_folds=3, use_hvg=False, top_depth_frac=1.0)  # fmt: skip
        sig = group_cv_permutation_null(X, y, groups, depth, ids, **kw)
        assert sig.summary["status"] == "ok" and sig.summary["observed_group_auc"] > 0.9
        assert sig.summary["empirical_p_auc"] <= 0.1
        Xn = np.random.default_rng(0).standard_normal(X.shape).astype(np.float32)
        noise = group_cv_permutation_null(Xn, y, groups, depth, ids, **kw)
        assert noise.summary["empirical_p_auc"] > sig.summary["empirical_p_auc"]


class TestRunHostresponseNull:
    def test_default_has_no_permutation_outputs(self, tmp_path):
        kw = _write_inputs(tmp_path)
        prov = run_hostresponse(**kw, cv_mode="group", groups_column="sample_id", cv_folds=3)
        out = Path(kw["out_dir"])
        assert not list(out.glob("*_permutation_metrics.tsv"))
        assert prov["permutation_count"] == 0 and prov["permutation_unit"] is None

    def test_within_block_null_end_to_end(self, tmp_path):
        kw = _write_inputs(tmp_path)
        prov = run_hostresponse(
            **kw, cv_mode="group", groups_column="sample_id", cv_folds=3,
            permutations=5, permutation_unit="cell_within_block", permutation_block="sample_id",
        )  # fmt: skip
        out = Path(kw["out_dir"])
        reps = pd.read_csv(out / "VIRUS_A_permutation_metrics.tsv", sep="\t")
        assert list(reps["permutation_id"]) == [1, 2, 3, 4, 5]
        assert set(reps["permutation_unit"]) == {"cell_within_block"}
        row = pd.read_csv(out / "hostresponse_metrics.csv").iloc[0]
        assert row["permutation_status"] == "ok" and row["permutation_n_valid"] == 5
        assert (row["empirical_p_auc"] * 6) == pytest.approx(round(row["empirical_p_auc"] * 6))
        assert 1 <= round(row["empirical_p_auc"] * 6) <= 6  # (1 + k) / (1 + 5)
        assert row["observed_group_auc"] == pytest.approx(row["auc_mean"])  # same frozen statistic
        null = prov["viruses"]["VIRUS_A"]["null"]
        assert null["n_requested"] == 5 and null["assumption"]
        assert prov["permutation_block"] == "sample_id"

    def test_unexchangeable_group_null_reports_not_estimable(self, tmp_path):
        kw = _write_inputs(tmp_path)  # identical label vector in every group
        run_hostresponse(
            **kw, cv_mode="group", groups_column="sample_id", cv_folds=3,
            permutations=3, permutation_unit="group",
        )  # fmt: skip
        out = Path(kw["out_dir"])
        row = pd.read_csv(out / "hostresponse_metrics.csv").iloc[0]
        assert row["permutation_status"] == "not_estimable"
        assert "empirical_p_auc" not in row or pd.isna(row["empirical_p_auc"])
        reps = pd.read_csv(out / "VIRUS_A_permutation_metrics.tsv", sep="\t")
        assert (reps["status"] == "failed").all() and len(reps) == 3

    def test_invalid_specs_fail_before_cleanup(self, tmp_path):
        kw = _write_inputs(tmp_path)
        out = Path(kw["out_dir"])
        out.mkdir()
        sentinel = out / "OLD_gene_weights.csv"
        sentinel.write_text("keep")
        for bad in (
            dict(cv_mode="group", groups_column="sample_id", permutations=3),  # unit undeclared
            dict(permutations=3, permutation_unit="group"),  # cell mode
            dict(
                cv_mode="group",
                groups_column="sample_id",
                permutations=3,
                permutation_unit="cell_within_block",
                permutation_block="nope",
            ),  # missing block column  # fmt: skip
        ):
            with pytest.raises(ValueError):
                run_hostresponse(**kw, **bad)
            assert sentinel.exists()

    def test_stale_permutation_file_is_cleared(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        (out / "EBV_permutation_metrics.tsv").write_text("old")
        assert hr.clear_stale_virus_outputs(out)
        assert not (out / "EBV_permutation_metrics.tsv").exists()
