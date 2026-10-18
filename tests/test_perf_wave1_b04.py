# ruff: noqa: B023  (closures below are called within the loop iteration that binds them)
"""B04 strict oracle: detect_genes with the narrow CSR view vs the frozen original loop.

The oracle re-implements the baseline per-gene logic verbatim (full-matrix column access
through the unchanged ``_gene_counts_from_matrix``). Acceptance is exact: same dict order,
same keys, same scalar type, same dtype and same IEEE bits; same error type/message.
"""

import struct

import numpy as np
import pytest
from scipy import sparse

from viralscan.scripts import detection


def oracle(var_names, matrix, viral_accessions, threshold=1, groups=None):
    names = list(var_names)
    if len(set(names)) != len(names):
        dupes = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"detect_genes: duplicate var_names {dupes[:5]!r}; make them unique first")
    pos = {g: i for i, g in enumerate(names)}
    totals = {
        g: detection._gene_counts_from_matrix(matrix, pos[g]).sum()
        for g in viral_accessions
        if g in pos
    }
    bad = sorted(g for g, t in totals.items() if not np.isfinite(t))
    if bad:
        raise ValueError(f"detect_genes: non-finite counts (NaN/inf) in viral gene(s) {bad[:5]!r}")
    if groups is None:
        groups = {g: [g] for g in totals}
    found = {}
    for genes in groups.values():
        member = {g: totals[g] for g in genes if g in totals}
        if sum(member.values()) >= threshold:
            found.update({g: c for g, c in member.items() if c > 0})
    return found


def sig(value):
    if isinstance(value, dict):
        return [(k, sig(v)) for k, v in value.items()]  # keep insertion order
    if isinstance(value, np.generic):
        return (type(value).__module__, type(value).__name__, value.dtype.str, value.tobytes())
    if isinstance(value, float):
        return ("float", struct.pack("<d", value))
    return (type(value).__name__, value)


def run(fn):
    try:
        return ("ok", sig(fn()))
    except Exception as exc:  # noqa: BLE001
        return ("err", type(exc).__name__, str(exc))


def layouts(rows):
    out = {
        "dense": np.array(rows),
        "csr": sparse.csr_matrix(rows),
        "csc": sparse.csc_matrix(rows),
        "coo": sparse.coo_matrix(rows),
        "csr_array": sparse.csr_array(rows),
        "csc_array": sparse.csc_array(rows),
        "csr_f32": sparse.csr_matrix(np.asarray(rows, dtype=np.float32)),
        "csr_i64": sparse.csr_matrix(np.nan_to_num(np.asarray(rows)).astype(np.int64)),
        "csr_i32": sparse.csr_matrix(np.nan_to_num(np.asarray(rows)).astype(np.int32)),
    }
    unsorted = sparse.csr_matrix(rows)
    for a, b in zip(unsorted.indptr[:-1], unsorted.indptr[1:], strict=True):
        unsorted.indices[a:b] = unsorted.indices[a:b][::-1].copy()
        unsorted.data[a:b] = unsorted.data[a:b][::-1].copy()
    unsorted.has_sorted_indices = False
    out["csr_unsorted"] = unsorted
    return out


ROWS = np.array(
    [
        [1e16, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0.0] * 12,
        [1.0, 1e-16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
        [-2.5, 3.0, -0.0, 0, 5.0, 0, 0, 0, 0, 0, 0, 0],
        [1 / 3, 2 / 3, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
)
NAMES = [f"g{i}" for i in range(12)]
ACC_CASES = {
    "subset": ["g0", "g1", "g4", "g9", "absent", "between"],
    "repeated": ["g4", "g0", "g4", "g9", "g0"],
    "reverse": ["g9", "g4", "g1", "g0"],
    "none": [],
    "all_columns": list(NAMES),
    "all_plus_missing": list(NAMES) + ["nope"],
    "tuple": ("g1", "g0"),
}
GROUPS = {"A": ["g0", "g1"], "B": ["g2", "g3", "g4"], "C": ["g9", "absent"]}


@pytest.mark.parametrize("layout", list(layouts(ROWS)))
@pytest.mark.parametrize("acc", list(ACC_CASES))
@pytest.mark.parametrize("thr", [1, 2, 6, 1e300])
@pytest.mark.parametrize("grouped", [False, True])
def test_matches_oracle(layout, acc, thr, grouped):
    m = layouts(ROWS)[layout]
    g = GROUPS if grouped else None
    a = ACC_CASES[acc]
    assert run(lambda: detection.detect_genes(NAMES, m, a, thr, g)) == run(
        lambda: oracle(NAMES, m, a, thr, g)
    )


@pytest.mark.parametrize("layout", list(layouts(ROWS)))
def test_nonfinite_errors(layout):
    bad = ROWS.copy()
    bad[0, 0] = np.nan
    bad[1, 4] = np.inf
    bad[2, 9] = -np.inf
    m = layouts(bad)[layout] if layout not in ("csr_i64", "csr_i32") else layouts(ROWS)[layout]
    a = ["g9", "g4", "g0", "g1"]
    assert run(lambda: detection.detect_genes(NAMES, m, a, 1, None)) == run(
        lambda: oracle(NAMES, m, a, 1, None)
    )


def test_duplicate_names_checked_before_access():
    names = NAMES[:3] + ["g0"] + NAMES[4:]
    m = sparse.csr_matrix(ROWS)
    assert run(lambda: detection.detect_genes(names, m, ["g0"], 1, None)) == run(
        lambda: oracle(names, m, ["g0"], 1, None)
    )


def test_width_mismatch_and_generators_and_unhashable():
    m = sparse.csr_matrix(ROWS)
    short = NAMES[:5]  # fewer names than columns
    assert run(lambda: detection.detect_genes(short, m, ["g1", "g4"])) == run(
        lambda: oracle(short, m, ["g1", "g4"])
    )
    long = NAMES + ["x", "y"]  # more names than columns -> IndexError at access
    assert run(lambda: detection.detect_genes(long, m, ["g1", "y"])) == run(
        lambda: oracle(long, m, ["g1", "y"])
    )
    assert run(lambda: detection.detect_genes(NAMES, m, (g for g in ["g0", "g4", "g0"]))) == run(
        lambda: oracle(NAMES, m, (g for g in ["g0", "g4", "g0"]))
    )
    assert run(lambda: detection.detect_genes(NAMES, m, ["g0", ["g1"]])) == run(
        lambda: oracle(NAMES, m, ["g0", ["g1"]])
    )


def test_duplicate_coordinates_csr():
    data = np.array([1e16, 1.0, -1e16, 0.1, 0.2, 0.3, 5.0])
    indices = np.array([3, 3, 3, 7, 7, 7, 3], dtype=np.int32)
    indptr = np.array([0, 3, 3, 6, 7], dtype=np.int32)
    m = sparse.csr_matrix((data, indices, indptr), shape=(4, 12))
    for acc in (["g3", "g7"], ["g7"], ["g3", "g0"]):
        assert run(lambda: detection.detect_genes(NAMES, m, acc, 1, None)) == run(
            lambda: oracle(NAMES, m, acc, 1, None)
        )


def test_input_not_mutated():
    m = layouts(ROWS)["csr_unsorted"]
    before = (m.data.copy(), m.indices.copy(), m.indptr.copy(), m.has_sorted_indices, m.dtype)
    detection.detect_genes(NAMES, m, ["g4", "g0", "g9"], 1, None)
    after = (m.data, m.indices, m.indptr, m.has_sorted_indices, m.dtype)
    for x, y in zip(before[:3], after[:3], strict=True):
        assert x.tobytes() == y.tobytes()
    assert before[3:] == after[3:]


def test_narrow_view_used_and_fallbacks(monkeypatch):
    seen = []
    orig = detection._gene_counts_from_matrix

    def spy(matrix, idx):
        seen.append(getattr(matrix, "shape", None))
        return orig(matrix, idx)

    monkeypatch.setattr(detection, "_gene_counts_from_matrix", spy)
    m = sparse.csr_matrix(ROWS)
    detection.detect_genes(NAMES, m, ["g4", "g0", "g4"])
    assert seen == [(6, 2)] * 3
    for other in (sparse.csc_matrix(ROWS), ROWS):
        seen.clear()
        detection.detect_genes(NAMES, other, ["g4", "g0"])
        assert seen == [(6, 12)] * 2
    seen.clear()
    detection.detect_genes(NAMES, m, list(NAMES))  # all columns -> original matrix
    assert seen == [(6, 12)] * 12
    seen.clear()
    detection.detect_genes(NAMES, m, (g for g in ["g4"]))  # generator -> original matrix
    assert seen == [(6, 12)]


@pytest.mark.parametrize("seed", range(25))
def test_random_fuzz(seed):
    rng = np.random.default_rng(seed)
    n, p = int(rng.integers(1, 80)), int(rng.integers(2, 60))
    dense = rng.standard_normal((n, p)) * 10.0 ** rng.integers(-8, 17, size=(n, p))
    dense[rng.random((n, p)) < 0.8] = 0.0
    names = [f"v{i}" for i in range(p)]
    acc = [names[i] for i in rng.integers(0, p, size=int(rng.integers(0, 2 * p)))] + ["zz"]
    groups = {f"k{j}": names[j::5] for j in range(5)}
    for layout, m in layouts(dense).items():
        if layout in ("csr_i64", "csr_i32"):
            continue
        for g in (None, groups):
            assert run(lambda: detection.detect_genes(names, m, acc, 1e-3, g)) == run(
                lambda: oracle(names, m, acc, 1e-3, g)
            ), (seed, layout)
