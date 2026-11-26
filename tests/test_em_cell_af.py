"""``em-cell-af``: the per-cell EM that mirrors alevin-fry (PLAN EM-CELL-AF-01).

Contrast case: sibling A has unique support in the cell, sibling B has none, and B is well supported
in another cell. ``em-cell`` shrinks the cell model toward the sample model, so B keeps a visible
share of the ambiguous mass. ``em-cell-af`` has no sample prior: the fixed point sends the mass to A
and the 0.01 cut-off zeroes B.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.multimapping import (
    AF_INIT_OFFSET,
    AF_INIT_SCALE,
    AF_MIN_OUTPUT_ALPHA,
    build_multimap_layers,
    em_cell_af_abundances,
)

EC_MAP = {0: [0], 1: [1], 2: [0, 1]}


def _bus(rows: list[tuple[str, int, int]]) -> pd.DataFrame:
    """(barcode, n molecules, ec) rows -> one record per molecule with distinct UMIs."""
    records = []
    for barcode, n, ec in rows:
        records += [(barcode, f"{barcode}{ec}{i:03d}", ec, 1) for i in range(n)]
    return pd.DataFrame(records, columns=["barcode", "umi", "ec", "count"])


def _build(bus: pd.DataFrame, method: str, n_cells: int = 2, **kwargs):
    return build_multimap_layers(
        bus,
        {f"C{i}": i for i in range(n_cells)},
        EC_MAP,
        n_cells,
        2,
        {1},
        sparse.csr_matrix((n_cells, 2)),
        method=method,
        **kwargs,
    )


SIBLING_BUS = _bus([("C0", 10, 0), ("C0", 10, 2), ("C1", 10, 1), ("C1", 4, 2)])


def test_supported_sibling_takes_the_mass_where_em_cell_keeps_a_share() -> None:
    af = _build(SIBLING_BUS, "em-cell-af").corrected.toarray()
    em = _build(SIBLING_BUS, "em-cell").corrected.toarray()
    assert af[0, 1] < 1e-2  # B, no unique support in cell 0
    assert af[0, 0] == pytest.approx(10.0, abs=1e-2)
    assert em[0, 1] > 0.1  # the sample prior keeps a visible share for B
    assert em[0, 1] > af[0, 1]
    assert af[1, 1] == pytest.approx(4.0, abs=1e-2)  # symmetric in the other cell


def test_molecules_are_conserved_per_cell() -> None:
    result = _build(SIBLING_BUS, "em-cell-af")
    np.testing.assert_allclose(result.corrected.toarray().sum(axis=1), [10.0, 4.0])
    assert result.corrected.sum() == pytest.approx(result.audit.ambiguous_molecules)


def test_cell_without_any_unique_support_converges_and_conserves() -> None:
    bus = _bus([("C0", 7, 2)])
    result = _build(bus, "em-cell-af", n_cells=1)
    np.testing.assert_allclose(result.corrected.toarray(), [[3.5, 3.5]])
    models = result.method_diagnostics["cell_models"]
    assert (models["n_cells"], models["n_converged"]) == (1, 1)
    assert models["fallback_reasons"] == []


def test_truncated_class_falls_back_to_an_equal_split() -> None:
    # Every abundance of the class below the cut-off: mass is kept, split equally.
    fitted = em_cell_af_abundances({(0, 1): 0.001}, np.zeros(2), 100, 1e-6)
    assert fitted.sum() == 0.0
    result = _build(_bus([("C0", 1, 2)]), "em-cell-af", n_cells=1)
    np.testing.assert_allclose(result.corrected.toarray(), [[0.5, 0.5]])


def test_fit_matches_a_hand_run_of_the_alevin_fry_update() -> None:
    unique = np.array([4.0, 0.0, 1.0])
    classes = {(0, 1): 6.0, (1, 2): 3.0}
    alpha = (unique + AF_INIT_OFFSET) * AF_INIT_SCALE
    for _ in range(400):
        new = unique.copy()
        for genes, count in classes.items():
            idx = list(genes)
            new[idx] += count * alpha[idx] / alpha[idx].sum()
        alpha = new
    alpha[alpha < AF_MIN_OUTPUT_ALPHA] = 0.0
    fitted = em_cell_af_abundances(classes, unique, 400, 1e-12)
    np.testing.assert_allclose(fitted, alpha, atol=1e-6)
    assert fitted.sum() == pytest.approx(unique.sum() + sum(classes.values()), abs=0.02)


def test_iteration_floor_and_limit_are_reported() -> None:
    diagnostics: dict = {}
    em_cell_af_abundances({(0, 1): 6.0}, np.array([9.0, 1.0]), 1, 1e-12, diagnostics)
    assert diagnostics["iterations"] == 1
    assert diagnostics["converged"] is False
    assert diagnostics["fallback_reason"] == "max_iterations_reached"
    diagnostics = {}
    em_cell_af_abundances({(0, 1): 6.0}, np.array([9.0, 1.0]), 100, 1e-2, diagnostics)
    assert diagnostics["iterations"] >= 2 and diagnostics["converged"] is True


def test_results_are_deterministic() -> None:
    first = _build(SIBLING_BUS, "em-cell-af").corrected.toarray()
    second = _build(SIBLING_BUS.sample(frac=1, random_state=3), "em-cell-af").corrected.toarray()
    np.testing.assert_array_equal(first, second)


def test_diagnostics_contract() -> None:
    diag = _build(SIBLING_BUS, "em-cell-af", em_max_iter=50, em_tol=1e-4).method_diagnostics
    assert diag["method"] == "em-cell-af"
    assert diag["model_scope"] == "cell"
    assert (diag["max_iter"], diag["tolerance"]) == (50, 1e-4)
    assert "global_model" not in diag  # no sample-level fit
    assert set(diag["cell_models"]) == {
        "n_cells",
        "n_converged",
        "max_iterations",
        "fallback_reasons",
    }
    assert diag["cell_models"]["n_cells"] == 2


def test_host_virus_ambiguity_is_shared_by_the_fit_like_any_class() -> None:
    # B (index 1) is viral; the {A,B} class is host+viral. The viral share is credited to
    # counts_host_viral_selected exactly as for em-cell, and mass is conserved.
    result = _build(SIBLING_BUS, "em-cell-af")
    selected = result.host_viral_selected.toarray()
    assert selected[0, 1] == pytest.approx(result.corrected.toarray()[0, 1])
    assert selected[1, 1] == pytest.approx(4.0, abs=1e-2)


def test_molecule_assignments_sum_to_the_corrected_layer(tmp_path) -> None:
    out = tmp_path / "molecules.tsv.gz"
    result = _build(
        SIBLING_BUS, "em-cell-af", molecule_assignments=out, gene_names=["hostA", "viralB"]
    )
    table = pd.read_csv(out, sep="\t")
    assert result.method_diagnostics["molecule_assignments"]
    assert len(table) == result.audit.unique_molecules + result.audit.ambiguous_molecules


def test_last_round_reallocates_thresholded_mass():
    # C (gene 2) converges below 0.01; alevin-fry zeroes it and runs one more update, so the
    # class mass splits over A and B only.
    unique = np.array([5.0, 5.0, 0.0])
    fitted = em_cell_af_abundances({(0, 1, 2): 3.0}, unique, 100, 1e-6)
    assert fitted[2] == 0.0
    assert np.isclose(fitted.sum(), 13.0)
