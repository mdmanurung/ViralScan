"""REVIEW-FIX-A: em-cell per-EC fallback, ambiguity fraction, swap diagnostics, off-list audit."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan import menu
from viralscan.multimapping import build_multimap_layers
from viralscan.scripts import detection as D
from viralscan.scripts.multimap import load_barcodes

# genes: 0 host H, 1 viral V, 2 host X, 3 host Y
ECS = {0: [0, 1], 1: [2, 3], 2: [1], 3: [2]}


def _frame(rows):
    return pd.DataFrame(rows, columns=["barcode", "umi", "ec", "count"])


def test_em_cell_zero_support_molecule_stays_equal_split_beside_other_support():
    rows = [(f"C{c}", f"V{u:03d}", 2, 1) for c in range(5) for u in range(30)]
    rows += [("CB", "M", 0, 1), ("CB", "X", 3, 1), ("CB", "Z", 1, 1)]
    bc = {f"C{c}": c for c in range(5)} | {"CB": 5}
    layers = build_multimap_layers(
        _frame(rows), bc, ECS, 6, 4, {1}, sparse.csr_matrix((6, 4)), method="em-cell"
    )
    row = layers.corrected[5].toarray()[0]
    # H/V molecule has no local support: 0.5/0.5. X/Y molecule has support at X.
    assert row[0] == pytest.approx(0.5)
    assert row[1] == pytest.approx(0.5)
    assert row.sum() == pytest.approx(2.0)


def _ad(x, hv, unique):
    a = ad.AnnData(
        X=sparse.csr_matrix(np.asarray(x, float)),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(len(x))]),
        var=pd.DataFrame(index=["h1", "v1"]),
    )
    a.layers["counts_host_viral_ambiguous"] = sparse.csr_matrix(np.asarray(hv, float))
    a.layers["counts_unique"] = sparse.csr_matrix(np.asarray(unique, float))
    return a


@pytest.mark.parametrize(
    ("x", "mixed", "unique", "expected"),
    [
        ([[9, 1], [5, 0]], [[0, 9], [0, 0]], [[0, 1], [0, 0]], 0.9),
        ([[10, 0], [5, 0]], [[0, 10], [0, 0]], [[0, 0], [0, 0]], 1.0),
    ],
)
def test_ambig_fraction_is_mixed_over_unique_plus_mixed(x, mixed, unique, expected):
    stats, _ = D.compute_stats(_ad(x, mixed, unique), ["v1"], {"V": ["v1"]}, ["v1"])
    assert stats["V"]["host_viral_ambig_fraction"] == expected


def test_headline_unique_reads_alias_only_layer():
    a = _ad([[0, 5], [0, 0]], [[0, 0], [0, 0]], [[0, 5], [0, 0]])
    a.layers["counts_original"] = a.layers.pop("counts_unique")
    assert D._headline_totals(a, {"V": ["v1"]})["unique"] == 5


def test_load_barcodes_collision_raises(tmp_path):
    path = tmp_path / "b.txt"
    path.write_text("AAAC\nAAAC-1\nGGGT\n")
    with pytest.raises(ValueError, match="collide"):
        load_barcodes(str(path))


def test_offlist_barcode_molecules_are_counted():
    # two molecules for ZZ (one with two records), one on-list
    df = _frame([("A", "U1", 0, 1), ("ZZ", "U2", 0, 1), ("ZZ", "U2", 0, 1), ("ZZ", "U3", 0, 1)])
    layers = build_multimap_layers(df, {"A": 0}, {0: [0]}, 1, 2, set(), sparse.csr_matrix((1, 2)))
    assert layers.audit.input_molecules == 1
    assert layers.audit.offlist_barcode_molecules == 2


def test_swap_updates_diagnostics_method_and_scope(tmp_path):
    n_cells, n_genes = 2, 2
    z = sparse.csr_matrix(np.zeros((n_cells, n_genes)))
    a = ad.AnnData(
        X=z,
        obs=pd.DataFrame(index=["a", "b"]),
        var=pd.DataFrame(index=["h", "v"]),
    )
    for name in (
        "counts_unique",
        "counts_multimap_equal",
        "counts_multimap_unique_weighted",
        "counts_host_viral_selected_equal",
    ):
        a.layers[name] = z
    a.uns["count_schema_version"] = "3.0.0"
    a.uns["multimap_method"] = "em-cell"
    a.uns["multimap_diagnostics"] = {
        "method": "em-cell",
        "model_scope": "cell",
        "max_iter": 100,
        "tolerance": 1e-6,
        "global_model": {"converged": True},
        "cell_models": {"n_cells": 2},
    }
    path = tmp_path / "adata_multimap.h5ad"
    a.write_h5ad(path)
    assert menu._swap_multimap_layer(path, "equal")
    diag = ad.read_h5ad(path).uns["multimap_diagnostics"]
    assert diag["method"] == "equal"
    assert diag["model_scope"] == "deterministic"
    assert "cell_models" not in diag
    assert "global_model" not in diag
