"""Strict B21 oracle for dense-then-cast host features; no fitted-model changes."""

import importlib.util
import logging
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan import utils
from viralscan.scripts import hostresponse as candidate

SPEC = importlib.util.spec_from_file_location(
    "frozen_b21", Path(__file__).parent / "_frozen_oracles/frozen_host_features_d056986.py"
)
original = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(original)


def _adata(matrix):
    return ad.AnnData(matrix, var=pd.DataFrame(index=[f"g{i}" for i in range(matrix.shape[1])]))


def _capture(fn, adata, use_hvg, caplog):
    caplog.clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            array, names = fn(adata, use_hvg)
            result = (array.dtype.str, array.shape, array.strides, array.tobytes(), names)
            error = None
        except Exception as exc:  # compare the original exception type and message
            result = None
            error = (type(exc), str(exc))
    return (
        result,
        error,
        [(w.category, str(w.message)) for w in caught],
        [(r.levelname, r.getMessage()) for r in caplog.records],
    )


@pytest.mark.parametrize("layout", ["csr", "csc", "dense", "fortran", "duplicate"])
@pytest.mark.parametrize("dtype", [np.float64, np.float32, np.int64, np.uint64, np.complex128])
@pytest.mark.parametrize("shape", [(0, 3), (7, 0), (9, 7)])
def test_frozen_feature_bits_layout_warnings_and_input_immutability(
    monkeypatch, caplog, layout, dtype, shape
):
    caplog.set_level(logging.INFO, logger="viralscan.scripts.hostresponse")
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 64)
    values = np.arange(np.prod(shape)).reshape(shape).astype(dtype)
    if values.size and np.dtype(dtype).kind in "fc":
        values = values / 7
        values[0, 0] = -0.0
    if layout in ("csr", "duplicate"):
        matrix = sparse.csr_matrix(values)
        if layout == "duplicate" and shape[0] and shape[1]:
            matrix = sparse.csr_matrix(
                (
                    np.array([1e16, 1.0, -1e16]),
                    np.array([0, 0, 0]),
                    np.array([0, 3] + [3] * (shape[0] - 1)),
                ),
                shape=shape,
            )
    elif layout == "csc":
        matrix = sparse.csc_matrix(values)
    else:
        matrix = np.asfortranarray(values) if layout == "fortran" else values
    a, b = _adata(matrix.copy()), _adata(matrix.copy())
    before = (
        (b.X.data.copy(), b.X.indices.copy(), b.X.indptr.copy())
        if sparse.issparse(b.X)
        else b.X.copy()
    )
    assert _capture(original._select_features, a, False, caplog) == _capture(
        candidate._select_features, b, False, caplog
    )
    if sparse.issparse(b.X):
        assert all(
            x.tobytes() == y.tobytes()
            for x, y in zip(before, (b.X.data, b.X.indices, b.X.indptr), strict=True)
        )
    else:
        assert before.tobytes() == b.X.tobytes()


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, 1e300, -1e300, 1e-45, 1e-300])
@pytest.mark.parametrize("policy", ["ignore", "warn", "raise"])
def test_exceptional_cast_policies_match(monkeypatch, caplog, value, policy):
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 32)
    values = np.full((7, 3), value, dtype=np.float64)
    with np.errstate(all=policy):
        assert _capture(
            original._select_features, _adata(sparse.csr_matrix(values)), False, caplog
        ) == _capture(candidate._select_features, _adata(sparse.csr_matrix(values)), False, caplog)


def test_blocked_path_never_densifies_all_rows(monkeypatch):
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 64)
    source = sparse.csr_matrix(np.arange(63, dtype=float).reshape(9, 7))
    seen = []
    toarray = sparse.csr_matrix.toarray

    def observed(self, *args, **kwargs):
        seen.append(self.shape)
        return toarray(self, *args, **kwargs)

    monkeypatch.setattr(sparse.csr_matrix, "toarray", observed)
    result, names = candidate._select_features(_adata(source), False)
    assert result.dtype == np.float32 and result.flags.c_contiguous
    assert names == [f"g{i}" for i in range(7)]
    assert seen == [(1, 7)] * 9
    assert result.tobytes() == np.arange(63, dtype=float).reshape(9, 7).astype(np.float32).tobytes()


def test_hvg_feature_selection_and_mutations_match(monkeypatch, caplog):
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 64)
    values = np.random.default_rng(410).gamma(2.0, 1.0, size=(100, 30))
    a, b = _adata(sparse.csr_matrix(values)), _adata(sparse.csr_matrix(values))
    assert _capture(original._select_features, a, True, caplog) == _capture(
        candidate._select_features, b, True, caplog
    )
    pd.testing.assert_frame_equal(a.var, b.var, check_exact=True)


@pytest.mark.parametrize("mutation", ["dtype", "ulp", "layout", "names", "error"])
def test_comparator_rejects_output_mutations(monkeypatch, caplog, mutation):
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 64)
    source = _adata(sparse.csr_matrix(np.arange(12, dtype=float).reshape(4, 3) / 7))
    expected = _capture(original._select_features, source.copy(), False, caplog)
    assert expected == _capture(original._select_features, source.copy(), False, caplog)

    def mutated(adata, use_hvg):
        array, names = original._select_features(adata, use_hvg)
        if mutation == "dtype":
            array = array.astype(np.float64)
        elif mutation == "ulp":
            array[0, 0] = np.nextafter(array[0, 0], np.float32(1))
        elif mutation == "layout":
            array = np.asfortranarray(array)
        elif mutation == "names":
            names.reverse()
        else:
            raise ValueError("mutated conversion error")
        return array, names

    assert expected != _capture(mutated, source.copy(), False, caplog)
