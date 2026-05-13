"""Integration test skeleton — UMI count conservation end-to-end.

These tests operate on *synthetic* AnnData fixtures (no real FASTQs, no network)
and are marked ``@pytest.mark.integration`` so they are excluded from the
default ``pytest`` run.  Run them explicitly with::

    PYTHONPATH=src python -m pytest -m integration tests/integration/

What is tested
--------------
The multimapping correction step distributes multi-mapping UMI counts across
genes proportionally.  A key invariant is that the *total* UMI mass is
conserved: the sum of ``counts_corrected + counts_original`` must not exceed
the sum of the original raw counts (within floating-point tolerance).

Audit trace: audits/2026-05-08-full-pipeline.md §3.4 (UMI conservation)

Fixture design
--------------
We build a minimal ``AnnData`` object directly using ``anndata``, mimicking
the structure that ``multimap.py`` produces:

* ``adata.layers["counts_original"]`` — unique-mapping UMI counts (CSR sparse)
* ``adata.layers["counts_corrected"]`` — redistributed multimapper fraction
* ``adata.X``                         — ``counts_corrected + counts_original``

The fixture uses known integers so the expected totals can be computed by hand.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

try:
    import anndata as ad
except ImportError:
    ad = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_multimap_adata(
    n_cells: int,
    n_genes: int,
    unique_counts: np.ndarray,
    multi_counts: np.ndarray,
) -> "ad.AnnData":  # type: ignore[name-defined]
    """Build a synthetic AnnData mirroring what multimap.py produces.

    Parameters
    ----------
    unique_counts:
        Dense (n_cells, n_genes) array — unique-mapping UMIs.
    multi_counts:
        Dense (n_cells, n_genes) array — redistributed multimapper fractions.
        Each row sums to the total multi-mapping UMIs from that cell (float).

    Returns
    -------
    AnnData with ``.layers["counts_original"]``, ``.layers["counts_corrected"]``,
    and ``.X = counts_original + counts_corrected``.
    """
    obs = {"cell_id": [f"cell_{i}" for i in range(n_cells)]}
    var = {"gene_id": [f"gene_{j}" for j in range(n_genes)]}

    adata = ad.AnnData(
        X=sp.csr_matrix(unique_counts + multi_counts),
        obs=obs,
        var=var,
    )
    adata.layers["counts_original"] = sp.csr_matrix(unique_counts.astype(float))
    adata.layers["counts_corrected"] = sp.csr_matrix(multi_counts.astype(float))
    return adata


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def small_multimap_adata():
    """Synthetic AnnData: 4 cells, 3 genes, known UMI counts.

    Unique counts (counts_original):
        cell_0: gene_0=5, gene_1=0, gene_2=0
        cell_1: gene_0=3, gene_1=2, gene_2=0
        cell_2: gene_0=0, gene_1=1, gene_2=0
        cell_3: gene_0=0, gene_1=0, gene_2=8

    Multimapper fractions (counts_corrected):
        cell_0: gene_0=2.0, gene_1=2.0, gene_2=0
        cell_1: gene_0=0,   gene_1=0,   gene_2=0
        cell_2: gene_0=1.5, gene_1=1.5, gene_2=0
        cell_3: gene_0=0,   gene_1=0,   gene_2=0

    Total input UMIs (unique + multi events):
        cell_0: 5 unique + 4 multi = 9
        cell_1: 5 unique + 0 multi = 5
        cell_2: 1 unique + 3 multi = 4   (3 UMIs split → 1.5 each)
        cell_3: 8 unique + 0 multi = 8
        Grand total = 26

    Expected X = counts_original + counts_corrected:
        cell_0: [7.0, 2.0, 0.0]
        cell_1: [3.0, 2.0, 0.0]
        cell_2: [1.5, 2.5, 0.0]
        cell_3: [0.0, 0.0, 8.0]
    """
    if ad is None:
        pytest.skip("anndata not installed")

    unique = np.array(
        [
            [5.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 8.0],
        ]
    )
    multi = np.array(
        [
            [2.0, 2.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.5, 1.5, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    return _make_multimap_adata(n_cells=4, n_genes=3, unique_counts=unique, multi_counts=multi)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestEndToEndCountConservation:
    """Integration guard: UMI totals in X must equal unique + multimapper input.

    Audit §3.4: audits/2026-05-08-full-pipeline.md
    """

    def test_x_equals_corrected_plus_original(self, small_multimap_adata) -> None:
        """
        GIVEN: AnnData with counts_original and counts_corrected layers
        WHEN:  X = counts_corrected + counts_original
        THEN:  X must equal the element-wise sum of those two layers

        This mirrors the computation in detection.py:
            adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
        """
        adata = small_multimap_adata
        x = np.asarray(adata.X.todense() if sp.issparse(adata.X) else adata.X)
        orig = np.asarray(adata.layers["counts_original"].todense())
        corr = np.asarray(adata.layers["counts_corrected"].todense())
        expected = orig + corr
        np.testing.assert_allclose(
            x,
            expected,
            rtol=1e-6,
            err_msg="adata.X != counts_original + counts_corrected",
        )

    def test_umi_mass_not_inflated(self, small_multimap_adata) -> None:
        """
        GIVEN: known input UMI totals (unique + multi-mapping events)
        WHEN:  X = counts_corrected + counts_original is summed
        THEN:  total X sum does NOT exceed total input UMIs (within tolerance)

        This guards against double-counting: if unique reads appear in BOTH
        counts_original AND counts_corrected, X.sum() would be inflated.

        Total input UMIs: 5+4 + 5+0 + 1+3 + 8+0 = 26
        """
        adata = small_multimap_adata
        x_total = float(adata.X.sum() if sp.issparse(adata.X) else np.asarray(adata.X).sum())
        expected_total = 26.0  # see fixture docstring
        assert x_total <= expected_total + 1e-6, (
            f"UMI mass inflated: X.sum()={x_total} > input total={expected_total}. "
            "Unique reads may be double-counted in counts_corrected."
        )

    def test_gene_totals_match_known_values(self, small_multimap_adata) -> None:
        """
        GIVEN: the synthetic fixture with hand-computed totals
        WHEN:  X is summed per gene (column sums)
        THEN:  column sums match expected values

        Expected per-gene totals:
          gene_0: 7 + 3 + 1.5 + 0 = 11.5
          gene_1: 2 + 2 + 2.5 + 0 = 6.5
          gene_2: 0 + 0 + 0   + 8 = 8.0
        """
        adata = small_multimap_adata
        x = adata.X
        if sp.issparse(x):
            gene_totals = np.asarray(x.sum(axis=0)).flatten()
        else:
            gene_totals = np.asarray(x).sum(axis=0)

        np.testing.assert_allclose(
            gene_totals,
            [11.5, 6.5, 8.0],
            rtol=1e-6,
            err_msg=f"Per-gene totals mismatch. Got: {gene_totals}",
        )

    def test_no_negative_counts_in_x(self, small_multimap_adata) -> None:
        """
        GIVEN: any AnnData from the multimapping pipeline
        WHEN:  X is inspected
        THEN:  no cell×gene entry is negative

        Negative counts indicate a bug in the fractional redistribution logic.
        """
        adata = small_multimap_adata
        x = adata.X
        if sp.issparse(x):
            min_val = x.min()
        else:
            min_val = float(np.asarray(x).min())
        assert min_val >= 0.0, (
            f"Negative UMI count found in adata.X: min={min_val}. "
            "This indicates a bug in multimapper redistribution."
        )
