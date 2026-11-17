"""CONTRACT-DIAG-01: EM non-convergence reaches the user (log and summary.txt), not only adata.uns."""

import logging

import numpy as np
import pytest

from viralscan.multimapping import em_gene_abundances, em_nonconvergence_note


def test_no_note_for_converged_or_non_em_runs() -> None:
    assert em_nonconvergence_note(None) is None
    assert em_nonconvergence_note({}) is None
    assert em_nonconvergence_note({"method": "equal", "model_scope": "deterministic"}) is None
    converged = {
        "method": "em-cell",
        "global_model": {"converged": True, "iterations": 7},
        "cell_models": {"n_cells": 5, "n_converged": 5},
    }
    assert em_nonconvergence_note(converged) is None


def test_note_names_the_global_model_and_the_remedy() -> None:
    note = em_nonconvergence_note(
        {
            "method": "em-global",
            "tolerance": 1e-6,
            "max_iter": 100,
            "global_model": {"converged": False, "iterations": 100},
        }
    )
    assert note is not None
    assert "global model stopped after 100 iterations" in note
    assert "--multimap-em-max-iter" in note and "--multimap-em-tol" in note


def test_note_counts_cell_models_that_hit_the_limit() -> None:
    note = em_nonconvergence_note(
        {
            "method": "em-cell",
            "tolerance": 1e-6,
            "max_iter": 100,
            "global_model": {"converged": True, "iterations": 3},
            "cell_models": {"n_cells": 40, "n_converged": 37},
        }
    )
    assert note is not None and "3 of 40 cell models" in note
    assert "global model" not in note


def test_a_real_one_iteration_fit_produces_the_note() -> None:
    diagnostics: dict = {}
    em_gene_abundances(
        {(0, 1): 10.0, (1, 2): 10.0},
        np.array([1.0, 0.0, 1.0]),
        1.0,
        1,
        1e-12,
        diagnostics=diagnostics,
    )
    assert diagnostics["converged"] is False
    note = em_nonconvergence_note(
        {"method": "em-global", "tolerance": 1e-12, "max_iter": 1, "global_model": diagnostics}
    )
    assert note is not None and "1 iterations" in note


def test_the_h5ad_contract_requires_the_diagnostics() -> None:
    import json
    from importlib import resources

    contract = json.loads(
        resources.files("viralscan.schemas.v3").joinpath("h5ad_contract.json").read_text("utf-8")
    )
    assert "multimap_diagnostics" in contract["required_uns"]


@pytest.fixture(autouse=True)
def _quiet(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.WARNING)
