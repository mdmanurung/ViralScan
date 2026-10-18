# ruff: noqa: B023  (closures below are called within the loop iteration that binds them)
"""Strict (bitwise) oracle for the B06 row-blocked helper.

The oracle is the baseline expression, reproduced inline, run on the same input.
Comparison is on raw bytes (``tobytes``), never ``allclose``.
"""

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan import utils


def _oracle(adata, matrix, gene_ids):
    sub = utils.matrix_for_genes(adata, matrix, [g for g in gene_ids if g in adata.var_names])
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    return np.asarray(sub, dtype=float).sum(axis=1).reshape(-1)


def _adata(n_rows, n_cols):
    names = [f"g{i}" for i in range(n_cols)]
    return ad.AnnData(
        np.zeros((n_rows, n_cols), dtype=np.float32), var=pd.DataFrame(index=names)
    ), names


def _values(n_rows, n_cols, seed):
    rng = np.random.default_rng(seed)
    v = rng.integers(0, 4, size=(n_rows, n_cols)).astype(np.float64)
    v *= rng.random((n_rows, n_cols)) < 0.6
    v = v / 3.0 + (rng.random((n_rows, n_cols)) < 0.1) * 1e-16
    if n_rows and n_cols:
        v[0, 0] = 1e16
        if n_cols > 1:
            v[0, 1] = 1.0
        v[-1, -1] = -0.0
    return v


def _run_capture(fn):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            out = fn()
            err = None
        except Exception as exc:  # noqa: BLE001 - compare type and message
            out, err = None, (type(exc), str(exc))
    return out, err, [(w.category, str(w.message)) for w in caught]


def _same(a, b):
    assert a.dtype == b.dtype and a.shape == b.shape
    assert a.tobytes() == b.tobytes()


@pytest.fixture
def spy(monkeypatch):
    calls = {"blocks": 0}
    real = utils._sparse_row_block

    def counting(n_rows, n_cols):
        r = real(n_rows, n_cols)
        if r is not None:
            calls["blocks"] += 1
        return r

    monkeypatch.setattr(utils, "_sparse_row_block", counting)
    return calls


WIDTHS = [1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 257]


@pytest.mark.parametrize("dtype", [np.float64, np.float32, np.int64, np.bool_])
@pytest.mark.parametrize("n_cols", WIDTHS)
def test_csr_blocks_bitwise(monkeypatch, spy, dtype, n_cols):
    n_rows = 23
    vals = _values(n_rows, n_cols, n_cols)
    if dtype is np.bool_:
        vals = vals > 0
    elif np.issubdtype(dtype, np.integer):
        vals = np.rint(vals * 3)  # keep 1e16 exactly representable
    mat = sparse.csr_matrix(np.asarray(vals, dtype=dtype))
    adata, names = _adata(n_rows, n_cols)
    expect = _oracle(adata, mat, names)
    for block in (1, 2, 3, 5, 7, 22, 23, 24, 1000):
        monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", block * n_cols * 8)
        got = utils.virus_molecules_per_cell(adata, mat, names)
        _same(got, expect)
    assert spy["blocks"] >= 1  # blocked path really ran (block < n_rows cases)


@pytest.mark.parametrize("n_rows", [1, 2, 9, 4097])
@pytest.mark.parametrize("n_cols", [1, 5, 8, 9, 20, 64, 70])
def test_row_counts(monkeypatch, n_rows, n_cols):
    vals = _values(n_rows, n_cols, 3)
    mat = sparse.csr_matrix(vals)
    adata, names = _adata(n_rows, n_cols)
    expect = _oracle(adata, mat, names)
    for block in (1, 2, 3, 64, n_rows - 1 if n_rows > 1 else 1, n_rows, n_rows + 1):
        monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", max(block, 1) * n_cols * 8)
        _same(utils.virus_molecules_per_cell(adata, mat, names), expect)


@pytest.mark.parametrize(
    "genes",
    [["g1", "g1", "g3"], ["g9", "g0", "g3", "g2"], ["nope", "g0", "g7"], [], ["nope"]],
)
def test_selection_variants(monkeypatch, genes):
    vals = _values(17, 10, 5)
    mat = sparse.csr_matrix(vals)
    adata, _ = _adata(17, 10)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 8)  # tiny: one row per block
    r1 = _run_capture(lambda: utils.virus_molecules_per_cell(adata, mat, genes))
    monkeypatch.undo()
    r2 = _run_capture(lambda: _oracle(adata, mat, genes))
    assert r1[1] == r2[1] and r1[2] == r2[2]
    if r2[0] is not None:
        _same(r1[0], r2[0])


def test_other_layouts_unchanged(monkeypatch, spy):
    vals = _values(20, 12, 8)
    adata, names = _adata(20, 12)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 8)
    for m in (
        sparse.csc_matrix(vals),
        sparse.coo_matrix(vals).tocsr().tocsc(),
        vals,
        vals.T.copy().T,
    ):
        _same(utils.virus_molecules_per_cell(adata, m, names), _oracle(adata, m, names))
    assert spy["blocks"] == 0  # csc/dense never block


def test_unsorted_and_duplicate_csr(monkeypatch):
    # duplicates and unsorted indices inside a row
    data = np.array([1e16, 1.0, 1e-16, 0.5, 0.25, 3.0, 2.0], dtype=np.float64)
    indices = np.array([4, 0, 4, 2, 2, 1, 1], dtype=np.int32)
    indptr = np.array([0, 3, 5, 5, 7], dtype=np.int32)
    mat = sparse.csr_matrix((data, indices, indptr), shape=(4, 5))
    adata, names = _adata(4, 5)
    expect = _oracle(adata, mat, names)
    for block in (1, 2, 3):
        monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", block * 5 * 8)
        _same(utils.virus_molecules_per_cell(adata, mat, names), expect)


def test_csr_array_type(monkeypatch):
    vals = _values(11, 9, 1)
    mat = sparse.csr_array(vals)
    adata, names = _adata(11, 9)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 2 * 9 * 8)
    _same(utils.virus_molecules_per_cell(adata, mat, names), _oracle(adata, mat, names))


def _nonfinite_cases():
    base = _values(13, 6, 2)
    for label, mut in (
        ("nan_row0", lambda v: v.__setitem__((0, 1), np.nan)),
        ("nan_last", lambda v: v.__setitem__((12, 5), np.nan)),
        ("inf_mid", lambda v: v.__setitem__((6, 2), np.inf)),
        (
            "inf_minus_inf_row",
            lambda v: (v.__setitem__((4, 0), np.inf), v.__setitem__((4, 1), -np.inf)),
        ),
        (
            "overflow_finite",
            lambda v: (v.__setitem__((3, 0), 1.7e308), v.__setitem__((3, 1), 1.7e308)),
        ),
        (
            "overflow_two_rows",
            lambda v: (
                v.__setitem__((1, 0), 1.7e308),
                v.__setitem__((1, 1), 1.7e308),
                v.__setitem__((11, 0), -1.7e308),
                v.__setitem__((11, 1), -1.7e308),
            ),
        ),
        ("many_nan", lambda v: v.__setitem__((slice(None), 3), np.nan)),
    ):
        v = base.copy()
        mut(v)
        yield label, v


@pytest.mark.parametrize("label,vals", list(_nonfinite_cases()))
@pytest.mark.parametrize("errstate", [None, "raise"])
def test_nonfinite_warnings_and_errors(monkeypatch, label, vals, errstate):
    mat = sparse.csr_matrix(vals)
    adata, names = _adata(*vals.shape)

    def run(fn):
        if errstate is None:
            return _run_capture(fn)
        with np.errstate(all=errstate):
            return _run_capture(fn)

    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 2 * vals.shape[1] * 8)
    got = run(lambda: utils.virus_molecules_per_cell(adata, mat, names))
    expect = run(lambda: _oracle(adata, mat, names))
    assert got[1] == expect[1]
    assert got[2] == expect[2]
    if expect[0] is not None:
        assert got[0].tobytes() == expect[0].tobytes()  # NaN payload/sign included


def test_underflow_errstate_untouched(monkeypatch):
    vals = np.full((6, 4), 5e-324)
    mat = sparse.csr_matrix(vals)
    adata, names = _adata(6, 4)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 2 * 4 * 8)
    with np.errstate(all="raise"):
        got = _run_capture(lambda: utils.virus_molecules_per_cell(adata, mat, names))
        expect = _run_capture(lambda: _oracle(adata, mat, names))
    assert got[1] == expect[1] and got[2] == expect[2]
    if expect[0] is not None:
        _same(got[0], expect[0])


def test_zero_rows_columns_and_shape(monkeypatch):
    for shape in ((0, 5), (5, 0), (0, 0)):
        mat = sparse.csr_matrix(shape, dtype=np.float32)
        adata, names = _adata(*shape)
        monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 8)
        got = utils.virus_molecules_per_cell(adata, mat, names)
        _same(got, _oracle(adata, mat, names))


def test_helper_budget_contract():
    f = utils._sparse_row_block
    assert f(10, 0) is None
    assert f(0, 5) is None
    old = utils._ROW_BLOCK_BYTES
    try:
        utils._ROW_BLOCK_BYTES = 8  # smaller than one row of 100 cols -> single-row blocks
        assert f(10, 100) == 1
        utils._ROW_BLOCK_BYTES = 10 * 100 * 8
        assert f(10, 100) is None
        assert f(11, 100) == 10
    finally:
        utils._ROW_BLOCK_BYTES = old


# ---- dtype guard: complex/longdouble/etc. must keep baseline warning multiplicity ----


def _dtype_vals(dtype, n_rows=7, n_cols=5):
    v = _values(n_rows, n_cols, 11)
    kind = np.dtype(dtype).kind
    if kind == "c":
        return (v + 1j * (v[::-1] + 0.5)).astype(dtype)
    if kind in "iu":
        return (np.rint(np.abs(v) * 3) % 100).astype(dtype)
    return v.astype(dtype)


@pytest.mark.parametrize("filt", ["always", "default", "error"])
@pytest.mark.parametrize("block", [1, 2, 3, 6, 7, 8, 100])
@pytest.mark.parametrize(
    "dtype",
    [
        np.complex64,
        np.complex128,
        np.longdouble,
        np.float32,
        np.float64,
        np.int8,
        np.uint8,
        np.int32,
        np.uint64,
        np.int64,
        np.bool_,
        np.clongdouble,
    ],
)
def test_dtype_guard_warnings_and_bytes(monkeypatch, spy, dtype, block, filt):
    vals = _dtype_vals(dtype)
    if dtype is np.bool_:
        vals = vals != 0
    mat = sparse.csr_matrix(vals)
    adata, names = _adata(*vals.shape)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", block * vals.shape[1] * 8)

    def capture(fn):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter(filt)
            try:
                out, err = fn(), None
            except Exception as exc:  # noqa: BLE001
                out, err = None, (type(exc), str(exc))
        return out, err, [(w.category, str(w.message)) for w in caught]

    got = capture(lambda: utils.virus_molecules_per_cell(adata, mat, names))
    expect = capture(lambda: _oracle(adata, mat, names))
    assert got[1] == expect[1] and got[2] == expect[2]
    if expect[0] is not None:
        _same(got[0], expect[0])
    if np.dtype(dtype).kind in "cO" or np.dtype(dtype).itemsize > 8:
        assert spy["blocks"] == 0  # non-proven dtype never takes the blocked path


def test_complex_emits_single_warning_like_baseline(monkeypatch):
    # regression: blocked path used to emit one ComplexWarning per block
    vals = (np.arange(12).reshape(4, 3) + 1j).astype(np.complex128)
    mat = sparse.csr_matrix(vals)
    adata, names = _adata(4, 3)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 3 * 8)  # one row per block
    got = _run_capture(lambda: utils.virus_molecules_per_cell(adata, mat, names))
    expect = _run_capture(lambda: _oracle(adata, mat, names))
    assert len(expect[2]) == 1 and got[2] == expect[2]
    _same(got[0], expect[0])


def test_longdouble_overflow_warning_multiplicity(monkeypatch):
    # finite longdouble beyond float64 range: cast overflows; baseline warns once
    vals = np.full((6, 3), 1e4000, dtype=np.longdouble)
    mat = sparse.csr_matrix(vals)
    adata, names = _adata(6, 3)
    monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", 3 * 8)
    for es in (None, "raise"):

        def run(fn):
            if es is None:
                return _run_capture(fn)
            with np.errstate(all=es):
                return _run_capture(fn)

        got = run(lambda: utils.virus_molecules_per_cell(adata, mat, names))
        expect = run(lambda: _oracle(adata, mat, names))
        assert got[1] == expect[1] and got[2] == expect[2]
        if expect[0] is not None:
            _same(got[0], expect[0])


def test_large_integer_values_no_warning(monkeypatch):
    # int64/uint64 above 2**53 round to float64 silently, per element, block-independent
    vals = np.array(
        [[2**62 + 1, 3, 2**63 - 1], [5, 2**61 + 7, 1], [0, 0, 0], [2**53 + 1, 2**53 + 3, 9]],
        dtype=np.int64,
    )
    adata, names = _adata(4, 3)
    for dt in (np.int64, np.uint64):
        mat = sparse.csr_matrix(vals.astype(dt))
        for block in (1, 2, 3, 4, 5):
            monkeypatch.setattr(utils, "_ROW_BLOCK_BYTES", block * 3 * 8)
            got = _run_capture(lambda: utils.virus_molecules_per_cell(adata, mat, names))
            expect = _run_capture(lambda: _oracle(adata, mat, names))
            assert got[1] == expect[1] and got[2] == expect[2] == []
            _same(got[0], expect[0])
