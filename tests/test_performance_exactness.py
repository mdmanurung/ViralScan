"""Exact contracts against oracles frozen before the performance edits."""

import copy
import runpy
import struct
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan import multimapping as mm
from viralscan.defaults import MULTIMAP_METHODS
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection, umap

ORIGINAL = runpy.run_path(str(Path(__file__).parent / "fixtures/performance_original.py"))


def outcome(fn, *args):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            value = fn(*args)
            if isinstance(value, tuple):
                result = ([type(v) for v in value], struct.pack("ddd", *value))
            elif isinstance(value, float):
                result = (type(value), struct.pack("d", value))
            else:
                result = (type(value), value)
        except Exception as exc:
            result = (type(exc), str(exc))
        return result, [(type(w.message), str(w.message)) for w in caught]


@pytest.mark.parametrize("digits", [-2, 0, 6, 12])
def test_scalar_bits_types_errors_and_numpy_error_policy(digits):
    values = [
        0.0,
        -0.0,
        1 / 3,
        np.nextafter(0.0, 1.0),
        2.0**53,
        2.0**63,
        1e100,
        1e308,
        np.nan,
        np.inf,
        -np.inf,
        None,
        "bad",
        "2.5",
    ]
    for n in (0, 1, 10, 10000, 2**53):
        for sign in (-1, 1):
            edge = n + sign * (1e-8 + 1e-5 * abs(n))
            values.extend([np.nextafter(edge, -np.inf), edge, np.nextafter(edge, np.inf)])
    for policy in ("warn", "raise"):
        with np.errstate(all=policy):
            for value in values:
                assert outcome(detection._count_value, value, digits) == outcome(
                    ORIGINAL["_count_value"], value, digits
                )


@pytest.mark.parametrize("dtype", [np.bool_, np.int8, np.uint8, np.int64, np.float32])
@pytest.mark.parametrize("permutations", [-2, -1, 0, 1, 31])
def test_enrichment_bits_rng_and_ties(dtype, permutations):
    rng = np.random.default_rng(42)
    coords = rng.normal(size=(30, 2))
    coords[1:4] = coords[0]
    labels = np.zeros(30, dtype=dtype)
    labels[::3] = 1
    before = (coords.tobytes(), labels.tobytes())
    old_rng, new_rng = np.random.default_rng(18), np.random.default_rng(18)
    expected = outcome(
        ORIGINAL["viral_neighbor_enrichment"], coords, labels, 9, permutations, old_rng
    )
    actual = outcome(umap.viral_neighbor_enrichment, coords, labels, 9, permutations, new_rng)
    assert actual == expected
    assert old_rng.bit_generator.state == new_rng.bit_generator.state
    assert before == (coords.tobytes(), labels.tobytes())


@pytest.mark.parametrize(
    "labels",
    [
        np.array([1, 0, 2, 1]),
        np.array([1.0, 0.0, np.nan, 0.0]),
        np.array([[1], [0], [0], [1]]),
        [1, 0, 0, 1],
    ],
)
def test_enrichment_unsupported_input_keeps_original_behavior(labels):
    args = (np.zeros((4, 2)), labels, 3, 17, 0)
    assert outcome(umap.viral_neighbor_enrichment, *args) == outcome(
        ORIGINAL["viral_neighbor_enrichment"], *args
    )


def matrix_buffers(matrix):
    if sparse.issparse(matrix):
        return (
            matrix.dtype.str,
            matrix.data.tobytes(),
            matrix.indices.tobytes(),
            matrix.indptr.tobytes(),
        )
    return (matrix.dtype.str, matrix.tobytes())


def test_evidence_cache_budget_falls_back_without_submatrix_allocation(monkeypatch):
    a = ad.AnnData(
        sparse.csr_matrix(np.ones((3, 32))), var=pd.DataFrame(index=[f"g{i}" for i in range(32)])
    )
    a.layers["counts_corrected"] = a.X.copy()
    groups = {"v": ["g0", "g1"]}
    expected = ORIGINAL["summarize_multimap_evidence"](a, groups, RunConfig())
    calls = []
    getitem = sparse.csr_matrix.__getitem__

    def track(matrix, key):
        calls.append(key)
        return getitem(matrix, key)

    monkeypatch.setattr(np, "count_nonzero", lambda *args, **kwargs: 10**12)
    monkeypatch.setattr(sparse.csr_matrix, "__getitem__", track)
    actual = mm.summarize_multimap_evidence(a, groups, RunConfig())
    pd.testing.assert_frame_equal(expected, actual, check_exact=True)
    assert calls and all(isinstance(key[1], int) for key in calls)


def test_evidence_generator_groups_remain_single_pass():
    a = ad.AnnData(sparse.csr_matrix([[1.0, 2.0]]), var=pd.DataFrame(index=["a", "b"]))
    a.layers["counts_corrected"] = a.X.copy()
    expected = ORIGINAL["summarize_multimap_evidence"](
        a, {"v": (g for g in ["a", "b"])}, RunConfig()
    )
    actual = mm.summarize_multimap_evidence(a, {"v": (g for g in ["a", "b"])}, RunConfig())
    pd.testing.assert_frame_equal(expected, actual, check_exact=True)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int64])
@pytest.mark.parametrize("layout", ["dense", "csr", "csc", "unsorted", "duplicates"])
@pytest.mark.parametrize("method", MULTIMAP_METHODS)
def test_evidence_exact_columns_bytes_and_input_immutability(tmp_path, dtype, layout, method):
    values = np.random.default_rng(19).integers(0, 4, size=(65, 32)).astype(dtype)
    if np.issubdtype(dtype, np.floating):
        values /= 3
        values[:, 0] = np.array([1.0] + [1e-16] * 64, dtype=dtype)
    matrix = sparse.csr_matrix(values)
    matrix.data[::17] = 0
    if layout == "dense":
        matrix = matrix.toarray()
    elif layout == "csc":
        matrix = matrix.tocsc()
    elif layout == "unsorted":
        for row in range(matrix.shape[0]):
            lo, hi = matrix.indptr[row : row + 2]
            matrix.indices[lo:hi] = matrix.indices[lo:hi][::-1]
            matrix.data[lo:hi] = matrix.data[lo:hi][::-1]
        matrix.has_sorted_indices = False
    elif layout == "duplicates":
        matrix = sparse.csr_matrix(
            (
                np.tile(np.array([1.0, 1e-16, 0.0], dtype=dtype), 65),
                np.tile([0, 0, 1], 65),
                np.arange(0, 196, 3),
            ),
            shape=(65, 32),
        )
    a = ad.AnnData(matrix, var=pd.DataFrame(index=[f"g{i}" for i in range(32)]))
    for name in (
        "counts_unique_viral",
        "counts_corrected",
        "counts_host_viral_ambiguous",
        "counts_host_viral_selected",
        "counts_viral_ambiguous_upper",
    ):
        a.layers[name] = matrix.copy()
    groups = {"z": ["g0", "unknown", "g1", "g0"], "a": ["g3", "g2"], "empty": []}
    config = RunConfig(multimap_method=method)
    before = {name: matrix_buffers(x) for name, x in a.layers.items()}
    old = ORIGINAL["summarize_multimap_evidence"](copy.deepcopy(a), groups, config)
    new = mm.summarize_multimap_evidence(a, groups, config)
    pd.testing.assert_frame_equal(old, new, check_exact=True)
    for column in old:
        if old[column].dtype.kind != "O":
            assert old[column].to_numpy().tobytes() == new[column].to_numpy().tobytes()
    mm.write_multimap_evidence(old, str(tmp_path / "old"))
    mm.write_multimap_evidence(new, str(tmp_path / "new"))
    assert (tmp_path / "old/results/multimap_evidence.tsv").read_bytes() == (
        tmp_path / "new/results/multimap_evidence.tsv"
    ).read_bytes()
    assert before == {name: matrix_buffers(x) for name, x in a.layers.items()}
