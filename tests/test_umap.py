"""Tests for UMAP-related statistical functions in scripts/umap.py.

``umap.py`` is a Snakemake script (references ``snakemake.*`` at module level)
and depends on scanpy, so we cannot import it directly in the test environment.
Instead we test the core statistical logic — ``viral_neighbor_enrichment`` —
via a standalone re-implementation that mirrors the *fixed* version of the
function.  The re-implementation pattern is the same used by test_analysis.py
and test_multimap.py.

Audit findings covered:
  §3.1 — ``viral_neighbor_enrichment`` and sc PCA/UMAP calls lack random seeds,
          making UMAP coordinates and permutation p-values non-reproducible.

Regression for: audits/2026-05-08-full-pipeline.md §3.1
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.neighbors import NearestNeighbors


# ---------------------------------------------------------------------------
# Standalone re-implementation of viral_neighbor_enrichment
# (mirrors the *fixed* code with random_state parameter)
# ---------------------------------------------------------------------------


def _viral_neighbor_enrichment(
    coords: np.ndarray,
    labels: np.ndarray,
    k: int,
    n_permutations: int = 1000,
    random_state: int = 0,
) -> tuple[float, float, float]:
    """Permutation test for viral-cell spatial clustering in UMAP space.

    Fixed version: uses ``np.random.default_rng(random_state)`` instead of
    the global ``np.random.permutation`` to guarantee reproducibility.
    """
    rng = np.random.default_rng(random_state)

    viral_cells = np.where(labels == 1)[0]
    if len(viral_cells) == 0:
        return 0.0, 0.0, 1.0

    nbrs = NearestNeighbors(n_neighbors=k).fit(coords)
    _distances, indices = nbrs.kneighbors(coords)

    counts = []
    for i in viral_cells:
        neighbor_idx = indices[i][1:]
        counts.append(np.mean(labels[neighbor_idx]))
    observed = float(np.mean(counts))

    permuted = []
    for _ in range(n_permutations):
        shuffled = rng.permutation(labels)
        counts_perm = []
        for i in viral_cells:
            neighbor_idx = indices[i][1:]
            counts_perm.append(np.mean(shuffled[neighbor_idx]))
        permuted.append(np.mean(counts_perm))

    expected = float(np.mean(permuted))
    p_value = (np.sum(np.array(permuted) >= observed) + 1) / (n_permutations + 1)
    return observed, expected, float(p_value)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def clustered_viral_data() -> tuple[np.ndarray, np.ndarray]:
    """100 cells: first 10 viral-positive and spatially clustered near (0,0)."""
    rng = np.random.default_rng(42)
    viral_coords = rng.normal(loc=[0.0, 0.0], scale=0.5, size=(10, 2))
    non_viral_coords = rng.normal(loc=[5.0, 5.0], scale=2.0, size=(90, 2))
    coords = np.vstack([viral_coords, non_viral_coords])
    labels = np.array([1] * 10 + [0] * 90)
    return coords, labels


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestViralNeighborEnrichmentReproducibility:
    """Audit §3.1: enrichment p-value must be reproducible with same random_state.

    Without a seeded RNG the permutation p-value changes between identical runs,
    making it possible for the same data to flip across the 0.05 significance
    threshold between runs.
    """

    def test_same_seed_produces_identical_p_values(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """
        GIVEN: identical coords, labels, k, and random_state=0
        WHEN:  _viral_neighbor_enrichment is called twice
        THEN:  both calls return exactly the same p-value

        Regression for: audits/2026-05-08-full-pipeline.md §3.1
        """
        coords, labels = clustered_viral_data
        _, _, p1 = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=200, random_state=0
        )
        _, _, p2 = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=200, random_state=0
        )
        assert p1 == p2, (
            f"Same random_state=0 yielded different p-values: {p1} vs {p2}. "
            "The permutation loop must use a seeded RNG, not np.random.permutation."
        )

    def test_same_seed_produces_identical_observed_enrichment(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """observed enrichment score must also be identical across identical calls."""
        coords, labels = clustered_viral_data
        obs1, _, _ = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=50, random_state=7
        )
        obs2, _, _ = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=50, random_state=7
        )
        assert obs1 == obs2, "Observed enrichment differs — kNN structure is non-deterministic"

    def test_no_viral_cells_returns_p_value_one(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """With no viral-positive cells, p_value must be 1.0 without crashing."""
        coords, _ = clustered_viral_data
        labels = np.zeros(len(coords), dtype=int)
        observed, expected, p_value = _viral_neighbor_enrichment(
            coords, labels, k=10, random_state=0
        )
        assert p_value == 1.0, f"Expected p=1.0 with no viral cells, got {p_value}"
        assert observed == 0.0
        assert expected == 0.0

    def test_p_value_bounded_in_unit_interval(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """p-value must always lie in [0, 1]."""
        coords, labels = clustered_viral_data
        _, _, p = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=100, random_state=0
        )
        assert 0.0 <= p <= 1.0, f"p-value {p} outside [0, 1]"

    def test_clustered_viral_cells_yield_significant_p_value(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """Spatially clustered viral cells must produce a low enrichment p-value."""
        coords, labels = clustered_viral_data
        _, _, p = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=500, random_state=0
        )
        assert p < 0.05, (
            f"Strongly clustered viral cells produced p={p:.4f} >= 0.05. "
            "The enrichment test may not be functioning correctly."
        )
