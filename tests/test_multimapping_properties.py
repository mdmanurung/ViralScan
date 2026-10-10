"""Property tests for the molecule counting contract over random fixtures (PLAN DEF-06).

Every method must satisfy the same conservation laws, checked against an oracle that
re-derives the contract from the raw records without using ``viralscan.multimapping``:

* each (cell, UMI) molecule is exactly one of unique, ambiguous or unresolved;
* allocated mass equals the number of ambiguous molecules, per cell and in total;
* ``X = unique + allocated`` so every resolved molecule counts exactly once;
* values are finite and non-negative, and mass only reaches compatible genes;
* ``sibling-weighted`` moves mass only inside one sibling group;
* results do not depend on record order, read multiplicity, off-list barcodes or the
  streaming buffer size.

Seeded ``numpy`` loops rather than ``hypothesis`` (not a dependency).
"""

from collections import defaultdict
from itertools import product

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.defaults import MULTIMAP_METHODS
from viralscan.multimapping import build_multimap_layers

SEEDS = range(40)
PSEUDOCOUNT = 1.0


def random_fixture(seed: int):
    """Random BUS records, EC map, indexed viral genes and sibling groups."""
    rng = np.random.default_rng(seed)
    n_cells = int(rng.integers(2, 6))
    n_genes = int(rng.integers(5, 10))
    n_ec = int(rng.integers(4, 10))
    ecs = {
        ec: [int(g) for g in rng.integers(0, n_genes, size=int(rng.integers(1, 4)))]
        for ec in range(n_ec)
    }
    viral = {int(g) for g in rng.choice(n_genes, size=int(rng.integers(1, n_genes)), replace=False)}
    viral_sorted = sorted(viral)
    n_groups = int(rng.integers(1, 3))
    members = [g for g in viral_sorted if rng.random() < 0.7]
    groups = {g: f"S{int(rng.integers(0, n_groups))}" for g in members}
    umis = ["".join(p) for p in product("ACG", repeat=2)]  # few UMIs: collisions are common
    rows = []
    for cell in range(n_cells + 1):  # the last barcode is off-list
        for _ in range(int(rng.integers(3, 14))):
            rows.append(
                (
                    f"B{cell}",
                    umis[int(rng.integers(len(umis)))],
                    int(rng.integers(n_ec)),
                    int(rng.integers(1, 4)),
                )
            )
    frame = pd.DataFrame(rows, columns=["barcode", "umi", "ec", "count"])
    barcode_to_idx = {f"B{c}": c for c in range(n_cells)}
    return frame, barcode_to_idx, ecs, n_cells, n_genes, viral, groups


def build(fixture, *, method, frame=None, groups=True, **kwargs):
    base, barcode_to_idx, ecs, n_cells, n_genes, viral, sibling = fixture
    return build_multimap_layers(
        base if frame is None else frame,
        barcode_to_idx,
        ecs,
        n_cells,
        n_genes,
        viral,
        sparse.csr_matrix((n_cells, n_genes)),
        method=method,
        pseudocount=PSEUDOCOUNT,
        sibling_groups=sibling if groups else None,
        **kwargs,
    )


class Oracle:
    """The contract, re-derived from the records: intersect the ECs of each CB-UMI."""

    def __init__(self, fixture):
        frame, barcode_to_idx, ecs, n_cells, n_genes, _, _ = fixture
        by_molecule = defaultdict(list)
        for row in frame.itertuples(index=False):
            if row.barcode in barcode_to_idx:
                by_molecule[(barcode_to_idx[row.barcode], row.umi)].append(row.ec)
        self.n_cells, self.n_genes = n_cells, n_genes
        self.unique = np.zeros((n_cells, n_genes))
        self.ambiguous: list[tuple[int, frozenset[int]]] = []
        self.unresolved = 0
        for (cell, _umi), ec_ids in by_molecule.items():
            genes = set.intersection(*(set(ecs[e]) for e in ec_ids))
            if not genes:
                self.unresolved += 1
            elif len(genes) == 1:
                self.unique[cell, next(iter(genes))] += 1
            else:
                self.ambiguous.append((cell, frozenset(genes)))
        self.n_input = len(by_molecule)
        self.n_unique = int(self.unique.sum())

    def compatible(self) -> np.ndarray:
        """Cell x gene mask of genes some ambiguous molecule of that cell is compatible with."""
        mask = np.zeros((self.n_cells, self.n_genes), dtype=bool)
        for cell, genes in self.ambiguous:
            mask[cell, list(genes)] = True
        return mask

    def sibling_weighted(self, sibling_groups) -> np.ndarray:
        """Equal split, except weights unique+pseudocount when all genes share one group."""
        out = np.zeros((self.n_cells, self.n_genes))
        for cell, genes in self.ambiguous:
            labels = {sibling_groups.get(g) for g in genes}
            ordered = sorted(genes)
            if len(labels) == 1 and None not in labels:
                weights = np.array([self.unique[cell, g] + PSEUDOCOUNT for g in ordered])
            else:
                weights = np.ones(len(ordered))
            out[cell, ordered] += weights / weights.sum()
        return out


@pytest.mark.parametrize("method", MULTIMAP_METHODS)
@pytest.mark.parametrize("seed", SEEDS)
def test_conservation_laws_hold_for_every_method(method, seed):
    fixture = random_fixture(seed)
    oracle = Oracle(fixture)
    result = build(fixture, method=method)

    # each molecule is exactly one of unique / ambiguous / unresolved
    audit = result.audit
    assert audit.input_molecules == oracle.n_input
    assert audit.unique_molecules == oracle.n_unique
    assert audit.ambiguous_molecules == len(oracle.ambiguous)
    assert audit.unresolved_molecules == oracle.unresolved
    assert (
        audit.unique_molecules + audit.ambiguous_molecules + audit.unresolved_molecules
        == audit.input_molecules
    )

    unique = result.unique.toarray()
    allocated = result.corrected.toarray()
    np.testing.assert_allclose(unique, oracle.unique)

    # allocated mass equals the ambiguous count, per cell and in total
    per_cell = np.zeros(oracle.n_cells)
    for cell, _ in oracle.ambiguous:
        per_cell[cell] += 1
    np.testing.assert_allclose(allocated.sum(axis=1), per_cell, atol=1e-9)
    assert allocated.sum() == pytest.approx(len(oracle.ambiguous), abs=1e-9)

    # X = unique + allocated: every resolved molecule counts exactly once
    x = unique + allocated
    assert x.sum() == pytest.approx(audit.resolved_molecules, abs=1e-9)

    # finite, non-negative, and mass only reaches genes the molecule is compatible with
    for matrix in (result.unique, result.corrected):
        assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all()
    assert not (allocated[~oracle.compatible()] > 1e-12).any()


@pytest.mark.parametrize("seed", SEEDS)
def test_sibling_weighted_matches_the_oracle_and_stays_inside_groups(seed):
    fixture = random_fixture(seed)
    oracle = Oracle(fixture)
    result = build(fixture, method="sibling-weighted")
    expected = oracle.sibling_weighted(fixture[-1])
    np.testing.assert_allclose(result.corrected.toarray(), expected, atol=1e-9)
    np.testing.assert_allclose(result.sibling_weighted.toarray(), expected, atol=1e-9)

    # no sibling groups: every molecule keeps its equal allocation
    empty = build_multimap_layers(
        fixture[0],
        fixture[1],
        fixture[2],
        fixture[3],
        fixture[4],
        fixture[5],
        sparse.csr_matrix((fixture[3], fixture[4])),
        method="sibling-weighted",
        pseudocount=PSEUDOCOUNT,
        sibling_groups={},
    )
    np.testing.assert_allclose(empty.corrected.toarray(), empty.equal.toarray(), atol=1e-12)


@pytest.mark.parametrize("method", MULTIMAP_METHODS)
@pytest.mark.parametrize("seed", range(12))
def test_results_do_not_depend_on_record_order_or_read_multiplicity(method, seed):
    fixture = random_fixture(seed)
    reference = build(fixture, method=method)
    frame = fixture[0].sample(frac=1, random_state=seed)
    frame = frame.assign(count=np.random.default_rng(seed).integers(1, 9, size=len(frame)))
    shuffled = build(fixture, method=method, frame=frame)
    np.testing.assert_allclose(
        shuffled.corrected.toarray(), reference.corrected.toarray(), atol=1e-9
    )
    np.testing.assert_allclose(shuffled.unique.toarray(), reference.unique.toarray())
    assert shuffled.audit.unique_molecules == reference.audit.unique_molecules
    assert shuffled.audit.ambiguous_molecules == reference.audit.ambiguous_molecules
    assert shuffled.audit.unresolved_molecules == reference.audit.unresolved_molecules


@pytest.mark.parametrize("method", MULTIMAP_METHODS)
@pytest.mark.parametrize("seed", range(6))
def test_streaming_buffer_size_does_not_change_the_result(method, seed, tmp_path):
    fixture = random_fixture(seed)
    path = tmp_path / "bus.txt"
    fixture[0].sort_values(["barcode", "umi", "ec"]).to_csv(
        path, sep="\t", header=False, index=False
    )
    reference = build(fixture, method=method)
    for buffer_size in (2, 8192):
        streamed = build(fixture, method=method, frame=path, bus_buffer_size=buffer_size)
        np.testing.assert_allclose(
            streamed.corrected.toarray(), reference.corrected.toarray(), atol=1e-9
        )
        assert streamed.audit == reference.audit


# ── EMC-01: em-cell keeps no dense per-cell abundance vector ─────────────────


def dense_em_cell_reference(fixture, em_max_iter=100, em_tol=1e-6):
    """The pre-EMC-01 em-cell allocation, written densely from the public EM functions."""
    from collections import Counter

    from viralscan.multimapping import em_cell_abundances, em_gene_abundances

    oracle = Oracle(fixture)
    pooled = Counter(tuple(sorted(genes)) for _, genes in oracle.ambiguous)
    theta = em_gene_abundances(
        {k: float(v) for k, v in pooled.items()},
        oracle.unique.sum(axis=0),
        PSEUDOCOUNT,
        em_max_iter,
        em_tol,
    )
    out = np.zeros((oracle.n_cells, oracle.n_genes))
    by_cell = defaultdict(list)
    for cell, genes in oracle.ambiguous:
        by_cell[cell].append(tuple(sorted(genes)))
    for cell, molecules in by_cell.items():
        counts = {k: float(v) for k, v in Counter(molecules).items()}
        fitted = em_cell_abundances(
            counts, oracle.unique[cell], theta, PSEUDOCOUNT, em_max_iter, em_tol
        )
        for genes in molecules:
            idx = np.asarray(genes)
            w = fitted[idx]
            out[cell, idx] += w / w.sum() if w.sum() > 0 else 1 / len(idx)
    return out


@pytest.mark.parametrize("seed", SEEDS)
def test_em_cell_matches_the_dense_reference(seed):
    fixture = random_fixture(seed)
    result = build(fixture, method="em-cell")
    np.testing.assert_allclose(
        result.corrected.toarray(), dense_em_cell_reference(fixture), rtol=0, atol=1e-12
    )


def test_em_cell_memory_does_not_scale_with_cells_times_genes():
    """300 cells x 100k genes: a dense theta per cell alone would be 240 MB."""
    import tracemalloc

    n_cells, n_genes = 300, 100_000
    rows = []
    for cell in range(n_cells):
        rows += [(f"B{cell}", "AA", 0, 1), (f"B{cell}", "CC", 1, 1), (f"B{cell}", "GG", 2, 1)]
    frame = pd.DataFrame(rows, columns=["barcode", "umi", "ec", "count"])
    ecs = {0: [1, 2], 1: [1], 2: [2, 3, 4]}
    barcode_to_idx = {f"B{c}": c for c in range(n_cells)}
    tracemalloc.start()
    result = build_multimap_layers(
        frame,
        barcode_to_idx,
        ecs,
        n_cells,
        n_genes,
        {1, 2, 3, 4},
        sparse.csr_matrix((n_cells, n_genes)),
        method="em-cell",
    )
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert result.corrected.sum() == pytest.approx(2 * n_cells)
    assert peak < 60e6, f"peak {peak / 1e6:.0f} MB suggests a dense per-cell vector is kept"
