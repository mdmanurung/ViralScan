"""Independent scientific-boundary checks for the isolated Track A analysis."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")
pytest.importorskip("sklearn")
pytestmark = pytest.mark.research

ANALYSIS = Path(__file__).resolve().parents[1] / "analysis/trackA_host_programs"
sys.path.insert(0, str(ANALYSIS))
import pilot  # noqa: E402
import tracka_analysis as analysis  # noqa: E402
import tracka_data as data  # noqa: E402
from setup_runtime import apply_source_fixes  # noqa: E402

sys.path.pop(0)


def test_convergence_retry_preserves_original_pilot():
    original, retry = pilot.config(), pilot.config(convergence_check=True)
    assert original["maxiter"] == 100 and retry["maxiter"] == 500
    assert original["output_root"] != retry["output_root"]
    assert all(
        original[key] == retry[key] for key in original if key not in {"maxiter", "output_root"}
    )


@pytest.mark.parametrize(
    "method,stage,budget,details_present,native_status,expected",
    [
        ("gbcd", "initialization", 25, True, "not_exposed", None),
        ("gbcd", "terminal_or_unknown", 500, True, "not_exposed", "iteration_limit"),
        ("gbcd", "initialization", 500, True, "not_exposed", "iteration_limit"),
        ("ebcd", "initialization", 25, True, "not_exposed", "iteration_limit"),
        ("gbcd", "initialization", 25, False, "not_exposed", "iteration_limit"),
        ("gbcd", "initialization", 25, True, "iteration_limit", "iteration_limit"),
    ],
)
def test_only_recorded_gbcd_initialization_limits_are_nonterminal(
    method, stage, budget, details_present, native_status, expected
):
    warning = "Maximum number of iterations reached."
    diagnostics = {
        "method": method,
        "diagnostics": {"convergence": native_status},
        "warnings": warning,
    }
    if details_present:
        diagnostics["warning_details"] = [{"message": warning, "stage": stage, "maxiter": budget}]
    assert analysis.fit_rejection(diagnostics) == expected
    if expected is None:
        diagnostics["warnings"] = [warning, warning]
        assert analysis.fit_rejection(diagnostics) == "iteration_limit"
        diagnostics["warnings"] = warning
        diagnostics["warning_details"][0]["message"] = "different warning"
        assert analysis.fit_rejection(diagnostics) == "iteration_limit"


@pytest.mark.parametrize(
    "warnings,expected,initialization_recorded",
    [
        (["Maximum number of iterations reached."], "iteration_limit", False),
        ("Maximum number of iterations reached.", "iteration_limit", False),
        (["unrelated diagnostic warning"], None, False),
        ([], None, False),
        (["Maximum number of iterations reached."], None, True),
    ],
)
def test_gbcd_warning_rejects_fresh_and_cached_fits(
    tmp_path, monkeypatch, warnings, expected, initialization_recorded
):
    diagnostics = {
        "diagnostics": {"convergence": "not_exposed_by_upstream_gbcd"},
        "warnings": warnings,
    }
    if initialization_recorded:
        diagnostics.update(
            method="gbcd",
            warning_details=[{"message": warnings[0], "stage": "initialization", "maxiter": 25}],
        )
    assert analysis.fit_rejection(diagnostics) == expected
    cfg = {
        "maxiter": 100,
        "divas_nsim": 40,
        "r_libraries": [],
        "rscript": "unused-mocked-subprocess",
    }
    genes = np.array(["g1", "g2", "g3"])
    cells = pd.DataFrame({"cell_id": ["c1", "c2"], "dataset": ["a", "b"]})

    def native_fit(*args, **kwargs):
        pd.DataFrame({"gene_id": genes, "factor1": [1.0, 2.0, 3.0]}).to_csv(
            tmp_path / "signatures.tsv", sep="\t", index=False
        )
        data.write_json(tmp_path / "diagnostics.json", diagnostics)
        return type("Result", (), {"returncode": 0})()

    monkeypatch.setattr(analysis.subprocess, "run", native_fit)
    monkeypatch.setattr(analysis, "sha256", lambda path: "test-only-source-hash")
    args = (cfg, "gbcd", np.ones((2, 3)), cells, genes, 3, 1, tmp_path)
    F = analysis.fit_method(*args)
    assert (F is None) == (expected is not None)
    status = json.loads((tmp_path / "status.json").read_text())
    assert status["status"] == ("method_unqualified" if expected else "ok")
    if expected:
        # Old receipts can say ok; recheck their native diagnostics on cache reads.
        status.update(status="ok", signature_sha256="test-only-source-hash")
        data.write_json(tmp_path / "status.json", status)
    cached = analysis.fit_method(*args)
    assert (cached is None) == (expected is not None)


def test_upstream_rank_one_correction_is_exact_and_idempotent(tmp_path):
    path = tmp_path / "signal.R"
    original = "before\nVBar <- VHat[, validPC]\nUBar <- UHat[, validPC]\nafter\n"
    path.write_text(original)
    pin = {
        "source_fixes": [
            {
                "path": "signal.R",
                "before": f"{x}Bar <- {x}Hat[, validPC]",
                "after": f"{x}Bar <- {x}Hat[, validPC, drop = FALSE]",
            }
            for x in ["V", "U"]
        ]
    }
    apply_source_fixes(tmp_path, pin)
    corrected = path.read_text()
    assert corrected == original.replace("[, validPC]", "[, validPC, drop = FALSE]")
    apply_source_fixes(tmp_path, pin)
    assert path.read_text() == corrected
    path.write_text("different upstream source\n")
    with pytest.raises(ValueError, match="differs"):
        apply_source_fixes(tmp_path, pin)


def test_frozen_barcode_join_and_integer_counts(tmp_path):
    raw = sparse.csr_matrix([[3, 0, 1], [0, 7, 2], [5, 6, 0]])
    a = ad.AnnData(
        X=raw,
        obs=pd.DataFrame(index=["c", "a", "b"]),
        var=pd.DataFrame(
            {"gene_id": ["g1", "g2", "v"], "is_viral": [False, False, True]},
            index=["g1", "g2", "v"],
        ),
    )
    for layer in ["counts_original", "counts_multimap_host_conservative", "counts_unique_viral"]:
        a.layers[layer] = raw.copy()
    path = tmp_path / "counts.h5ad"
    a.write_h5ad(path)
    before = data.sha256(path)
    got, _, _, _, barcodes = data.read_called_matrix(path, ["b", "c"])
    np.testing.assert_array_equal(got.toarray(), raw[[0, 2]].toarray())
    assert barcodes.tolist() == ["c", "b"]
    assert data.sha256(path) == before
    with pytest.raises(ValueError, match="exactly once"):
        data.read_called_matrix(path, ["absent"])
    with pytest.raises(ValueError, match="exactly once"):
        data.read_called_matrix(path, ["b", "b"])
    with pytest.raises(ValueError, match="integer"):
        data.validate_counts(sparse.csr_matrix([[0.5]]))


def test_activity_label_is_relative_burden_not_infection():
    counts = np.array([20, 2, 0, 10], dtype=float)
    depth = np.array([10000, 10, 100, 100], dtype=float)
    # Raw positives are 0 and 3; relative activity positives are 1 and 3.
    assert data.cpm_labels(counts, depth).tolist() == [0, 1, 0, 1]
    assert not data.cpm_labels(np.zeros(4), depth).any()


def test_gene_selection_uses_only_supplied_training_rows():
    X = sparse.csr_matrix([[1.0, 0.0, 1.0], [2.0, 0.0, 1.0], [3.0, 0.0, 1.0], [4.0, 0.0, 1.0]])
    cells = pd.DataFrame({"dataset": ["a"] * 4})
    # Add two independent variable genes so the minimum feasible input passes.
    X = sparse.hstack([X[:, 0], X[:, 0].multiply(2), X[:, 0].multiply(3)], format="csr")
    selected = analysis.training_genes(X, cells, {"min_detected_cells": 1, "n_genes": 1})
    assert selected.tolist() == [2]


def test_nested_splits_never_separate_repeated_units():
    cells = pd.DataFrame(
        {
            "dataset": ["a"] * 6 + ["b"] * 6,
            "biological_unit": [
                "a1",
                "a1",
                "a2",
                "a2",
                "a3",
                "a3",
                "b1",
                "b1",
                "b2",
                "b2",
                "b3",
                "b3",
            ],
        }
    )
    train, validation = analysis.inner_split(cells)
    assert not set(cells.iloc[train].biological_unit) & set(cells.iloc[validation].biological_unit)
    assert set(cells.iloc[train].dataset) == {"a", "b"}
    meta = pd.DataFrame(
        {
            "dataset": ["a"] * 4 + ["b"] * 3,
            "library": ["a1", "a1_repeat", "a2", "a3", "b1", "b2", "b3"],
            "biological_unit": ["A1", "A1", "A2", "A3", "B1", "B2", "B3"],
            "unit_known": [True] * 7,
            "included": [True] * 7,
        }
    )
    for fold in analysis.outer_folds(meta):
        train_units = set(meta[meta.library.isin(fold["train"])].biological_unit)
        test_units = set(meta[meta.library.isin(fold["test"])].biological_unit)
        assert not train_units & test_units


def test_reconstruction_scores_unobserved_genes():
    F = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])
    Y = np.array([[2.0], [7.0], [1.0]]) @ F.T
    assert analysis.reconstruction_error(Y, F) < 1e-9
    bad = F.copy()
    bad[[0, 4]] = bad[[4, 0]]
    assert analysis.reconstruction_error(Y, bad) > 0.1


def test_batched_projection_keeps_all_rows_and_handles_empty_libraries():
    rng = np.random.default_rng(7)
    X = sparse.csr_matrix(rng.normal(size=(4097, 3)))
    F = np.eye(3)[:, :2]
    result = analysis.project_batches(X, np.arange(3), F)
    np.testing.assert_allclose(result, X[:, :2].toarray() / (1 + 1e-6))
    assert analysis.project_batches(X[:0], np.arange(3), F).shape == (0, 2)


def test_stability_handles_sign_permutation_and_missing_seeds():
    F = np.eye(3)[:, :2]
    candidate = -F[:, ::-1]
    _, support = analysis.signature_stability(
        F, {1: candidate}, {"seeds": [1, 2], "signature_cosine": 0.7}
    )
    np.testing.assert_array_equal(support, [0.5, 0.5])


def test_cells_do_not_inflate_biological_replication():
    cfg = {"min_class_cells": 2, "min_biological_units": 3, "bootstrap_replicates": 100}
    cells = pd.DataFrame(
        {
            "dataset": ["a"] * 80,
            "library": np.repeat(["l1", "l2"], 40),
            "biological_unit": ["same_donor"] * 80,
            "unit_known": [True] * 80,
            "cell_type": ["T"] * 80,
            "target_resolved": [True] * 80,
            "host_depth": np.tile(np.arange(40) + 100, 2),
            "n_genes": [20] * 80,
            "mito_fraction": [0.01] * 80,
            "viral_cpm": np.tile(np.arange(40) ** 2, 2),
            "viral_unique_molecules": np.tile(np.arange(40) ** 2, 2),
            "activity_label": np.tile([0, 1], 40),
        }
    )
    scores = np.log1p(cells.viral_cpm.to_numpy())[:, None]
    result = analysis.associations(scores, cells, np.array([1.0]), cfg)
    assert result.n_units.max() == 1
    assert not result.inference_eligible.any()
    assert result.pvalue.isna().all()


def test_depth_only_signal_disappears_after_adjustment():
    rng = np.random.default_rng(1)
    depth = np.exp(rng.uniform(3, 10, 200))
    cells = pd.DataFrame(
        {"host_depth": depth, "mito_fraction": rng.uniform(0, 0.1, 200), "n_genes": np.sqrt(depth)}
    )
    # Both host score and viral burden are deterministic functions of depth.
    result = analysis.adjusted_slope(np.log1p(depth), depth, cells)
    assert result is None


def test_study_confounded_signal_has_no_within_library_evidence():
    study = np.repeat(["study_a", "study_b"], 120)
    label = (study == "study_b").astype(int)
    cells = pd.DataFrame(
        {
            "dataset": study,
            "library": np.repeat([f"l{i}" for i in range(6)], 40),
            "biological_unit": np.repeat([f"donor{i}" for i in range(6)], 40),
            "unit_known": True,
            "target_resolved": True,
            "cell_type": "T",
            "host_depth": np.tile(np.arange(40) + 100, 6),
            "n_genes": 20,
            "mito_fraction": 0.01,
            "viral_cpm": label * 1000,
            "viral_unique_molecules": label * 20,
            "activity_label": label,
        }
    )
    scores = label[:, None] + np.random.default_rng(1).normal(0, 0.01, (240, 1))
    result = analysis.associations(
        scores,
        cells,
        np.array([1.0]),
        dict(min_class_cells=5, min_biological_units=3, bootstrap_replicates=100),
    )
    assert result.n_units.max() == 0
    assert not result.inference_eligible.any()
    assert result.pvalue.isna().all()


def test_fdr_and_single_class_denominators():
    np.testing.assert_allclose(analysis.bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
    with pytest.raises(analysis.Unevaluable, match="one class"):
        analysis.metrics(np.zeros(5), np.linspace(0.1, 0.9, 5))


def test_descriptive_sharing_is_not_formal_support_and_requires_concordance():
    rows = [
        dict(
            dataset=study,
            method="gbcd",
            scope="pooled",
            program=0,
            cell_type="unannotated",
            sensitivity=sensitivity,
            effect=0.4,
            n_known_units=2,
            qvalue=np.nan,
            inference_eligible=False,
            stability_support=1.0,
        )
        for study in ["a", "b"]
        for sensitivity in ["primary", "unique_viral", "depth_matched"]
    ]
    table = pd.DataFrame(rows)
    result = analysis.sharing_classification(table, {"fdr": 0.05})[0]
    assert result["classification"] == "descriptive_cross_context"
    assert result["descriptive_studies"] == "a;b"
    assert result["descriptive_annotation_limited"]
    assert not result["passes_unit_inference_gates"]
    assert not result["virus_specificity_established"]
    formal_table = table.copy()
    formal_table["cell_type"] = "T"
    formal_table["n_known_units"] = 6
    formal_table["inference_eligible"] = True
    formal_table["qvalue"] = 0.01
    formal = analysis.sharing_classification(formal_table, {"fdr": 0.05})[0]
    assert formal["classification"] == "candidate_shared"
    assert formal["passes_unit_inference_gates"]
    table.loc[(table.dataset == "b") & table.sensitivity.eq("unique_viral"), "effect"] = -0.4
    assert (
        analysis.sharing_classification(table, {"fdr": 0.05})[0]["classification"]
        == "context_associated_or_unresolved"
    )
    table.loc[table.dataset == "b", "effect"] = 0.4
    table.loc[table.dataset == "b", "n_known_units"] = 1
    assert (
        analysis.sharing_classification(table, {"fdr": 0.05})[0]["classification"]
        == "context_associated_or_unresolved"
    )


def test_signflip_power_bound_prevents_unattainable_shared_claims():
    meta = pd.DataFrame(
        [
            dict(dataset=study, included=True, unit_known=True, biological_unit=f"{study}{u}")
            for study, n in [("a", 5), ("b", 15), ("c", 2)]
            for u in range(n)
        ]
    )
    result = analysis.inference_power(meta, {"targets": {"a": [], "b": [], "c": []}, "fdr": 0.05})
    assert [r["optimistic_min_two_sided_signflip_p"] for r in result] == [1 / 16, 1 / 16384, 1 / 2]
    assert sum(r["can_reach_raw_threshold"] for r in result) == 1


def test_divas_single_study_is_explicitly_unevaluable(tmp_path):
    cells = pd.DataFrame({"dataset": ["one"] * 3})
    with pytest.raises(analysis.Unevaluable, match="two study blocks"):
        analysis.learn_programs({}, "divas", None, cells, None, tmp_path)


def test_partial_report_retains_missing_and_failed_denominators(tmp_path, monkeypatch):
    meta = pd.DataFrame(
        {
            "dataset": ["a"],
            "library": ["l1"],
            "biological_unit": ["u1"],
            "unit_known": [True],
            "included": [True],
        }
    )
    meta.to_csv(tmp_path / "libraries.tsv", sep="\t", index=False)
    failed = tmp_path / "discovery/pooled/gbcd"
    failed.mkdir(parents=True)
    (failed / "status.json").write_text(json.dumps({"status": "method_failed"}))
    # Partial outputs from the failed task must never be consumed as valid results.
    (failed / "associations.tsv").write_text("corrupt_partial_output\n")
    monkeypatch.setattr(analysis, "run_root", lambda cfg: tmp_path)
    cfg = {"targets": {"a": []}, "methods": ["gbcd"], "fdr": 0.05}
    analysis.report(cfg)
    receipt = json.loads((tmp_path / "report.json").read_text())
    assert receipt["status"] == "incomplete_missing_tasks"
    assert receipt["tasks"] == 4
    assert receipt["failed"] == 1
    assert receipt["not_run"] == 3
    assert receipt["candidate_shared"] == 0


def test_zero_signature_columns_retain_native_program_mapping():
    F = np.array([[1.0, 0.0, 3.0], [0.0, 0.0, -1.0]])
    normalized = analysis.normalize_signatures(F)
    assert normalized.shape == F.shape
    assert not normalized[:, 1].any()
    with pytest.raises(analysis.Unevaluable):
        analysis.normalize_signatures(np.zeros((3, 2)))


def test_prediction_uncertainty_requires_independent_units():
    cells = pd.DataFrame(
        {
            "unit_known": [True] * 100,
            "biological_unit": ["one"] * 100,
            "activity_label": [0, 1] * 50,
        }
    )
    got = analysis.prediction_uncertainty(
        cells, {"QC_only": np.full(100, 0.5)}, {"min_biological_units": 3}
    )
    assert got["status"] == "insufficient_known_independent_test_units"


def test_independent_mock_contrasts_do_not_require_viral_calls():
    rng = np.random.default_rng(4)
    rows, scores = [], []
    for replicate in range(3):
        for hour in [0, 1, 3, 5]:
            for _ in range(20):
                rows.append(
                    dict(
                        dataset="hsv1",
                        cell_type="fibroblast",
                        biological_unit=f"rep{replicate}",
                        unit_known=True,
                        independent_status="uninfected" if hour == 0 else "infected",
                        time=str(hour),
                        host_depth=rng.uniform(1000, 2000),
                        mito_fraction=rng.uniform(0, 0.1),
                        n_genes=rng.uniform(100, 200),
                        target_resolved=False,
                    )
                )
                scores.append([hour + rng.normal(0, 0.1)])
    result = analysis.independent_exposure_associations(
        np.array(scores),
        pd.DataFrame(rows),
        np.array([1.0]),
        dict(min_class_cells=5, min_biological_units=3, bootstrap_replicates=100),
    )
    assert set(result.sensitivity) == {
        "independent_exposure_1",
        "independent_exposure_3",
        "independent_exposure_5",
        "independent_exposure_trend",
    }
    assert (result.n_units == 3).all()
    assert result.inference_eligible.all()
    assert (result.effect > 0.9).all()


def test_nested_prediction_never_fits_heldout_study(tmp_path, monkeypatch):
    """Test Python orchestration with a fixed basis; real engines have separate oracles."""
    cfg = dict(
        output_root=str(tmp_path),
        targets={"a": [], "b": [], "c": []},
        methods=["gbcd"],
        ranks=[2],
        seeds=[1],
        n_genes=3,
        max_cells_per_library=10,
        max_cells_per_study=20,
        min_detected_cells=1,
        min_class_cells=3,
        min_biological_units=3,
        signature_cosine=0.7,
        classifier_C=[0.1, 1],
        bootstrap_replicates=10,
    )
    rng = np.random.default_rng(10)
    meta = []
    for study in ["a", "b", "c"]:
        for unit in [1, 2]:
            library = f"{study}{unit}"
            directory = tmp_path / "prepared" / library
            directory.mkdir(parents=True)
            y = np.tile([0, 1], 10)
            X = rng.integers(1, 5, (len(y), 20))
            X[:, 0] = np.where(y, 25, 2)
            X[:, 1] = np.where(y, 2, 25)
            X[:, 19] = 1 if study != "c" else 1000000
            sparse.save_npz(directory / "host_counts.npz", sparse.csr_matrix(X))
            pd.DataFrame({"gene_id": [f"g{i}" for i in range(20)]}).to_csv(
                directory / "genes.tsv", sep="\t", index=False
            )
            cells = pd.DataFrame(
                {
                    "cell_id": [f"{library}:cell{i}" for i in range(len(y))],
                    "dataset": study,
                    "library": library,
                    "biological_unit": library,
                    "unit_known": True,
                    "target_resolved": True,
                    "cell_type": "synthetic",
                    "host_depth": X.sum(axis=1),
                    "mito_fraction": 0.01,
                    "n_genes": 20,
                    "activity_label": y,
                    "viral_molecules": y * 20,
                }
            )
            cells.to_csv(directory / "cells.tsv.gz", sep="\t", index=False)
            meta.append(
                dict(
                    dataset=study,
                    library=library,
                    biological_unit=library,
                    unit_known=True,
                    included=True,
                )
            )
    pd.DataFrame(meta).to_csv(tmp_path / "libraries.tsv", sep="\t", index=False)
    seen = []

    def fixed_basis(cfg, method, Y, cells, genes, rank, seed, directory):
        seen.extend(cells.cell_id.tolist())
        assert "g19" not in genes  # Held-out-only variance must never select a gene.
        directory.mkdir(parents=True, exist_ok=True)
        F = np.eye(len(genes), 2)
        pd.DataFrame(F, index=pd.Index(genes, name="gene_id")).to_csv(
            directory / "signatures.tsv", sep="\t"
        )
        return F

    monkeypatch.setattr(analysis, "run_root", lambda cfg: tmp_path)
    monkeypatch.setattr(analysis, "verify_runtime", lambda cfg: None)
    monkeypatch.setattr(analysis, "fit_method", fixed_basis)
    analysis.predict_fold(cfg, "gbcd", "study_c")
    out = tmp_path / "prediction/study_c/gbcd"
    assert json.loads((out / "status.json").read_text())["status"] == "ok"
    assert seen and not any(cell.startswith("c") for cell in seen)
    predictions = pd.read_csv(out / "predictions.tsv.gz", sep="\t")
    assert set(predictions.dataset) == {"c"}
    assert len(predictions) == 40
    assert all(
        np.isfinite(predictions[name]).all()
        for name in [
            "QC_only",
            "host_programs",
            "host_programs_plus_QC",
            "scrambled_host_programs_plus_QC",
        ]
    )
    assert (
        json.loads((out / "metrics.json").read_text())["uncertainty"]["status"]
        == "insufficient_known_independent_test_units"
    )
