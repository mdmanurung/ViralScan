"""Tests for viral gene detection threshold logic in scripts/detection.py.

``detection.py`` is a Snakemake script that references ``snakemake.*`` at module
level.  We test the core logic via standalone re-implementations, following the
established pattern from test_analysis.py.

Audit findings covered:
  §2.1 — ``preprocessing()`` uses ``total_count >= threshold``.  We guard against
          regressions where the comparison is changed (e.g. to ``>`` which would
          silently break at threshold=1) or where threshold is read incorrectly.

Regression for: audits/2026-05-08-full-pipeline.md §2.1
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.enrichment import _bh_adjust, cell_type_enrichment


# ---------------------------------------------------------------------------
# Standalone re-implementation of the detection filtering logic
# (mirrors the relevant section of preprocessing() in detection.py)
# ---------------------------------------------------------------------------


def _detect_genes(
    var_names: list[str],
    counts_matrix,
    viral_accessions: set[str],
    threshold: int = 1,
) -> dict[str, float]:
    """Mirror of the detection loop in detection.py:preprocessing().

    Parameters
    ----------
    var_names:
        List of gene IDs (``adata.var_names``).
    counts_matrix:
        Dense or sparse array, shape (n_cells, n_genes).  Rows are cells,
        columns are genes indexed by ``var_names``.
    viral_accessions:
        Set of viral gene IDs to look for in ``var_names``.
    threshold:
        Minimum total UMI count across all cells for a viral gene to be
        reported as detected.  The audit confirmed the check is ``>=``.

    Returns
    -------
    dict mapping detected gene_id → total_count (float).
    """
    if hasattr(counts_matrix, "toarray"):
        dense = counts_matrix.toarray()
    else:
        dense = np.asarray(counts_matrix)

    found: dict[str, float] = {}
    for gene_id in viral_accessions:
        if gene_id not in var_names:
            continue
        idx = var_names.index(gene_id)
        total_count = float(dense[:, idx].sum())
        if total_count >= threshold:
            found[gene_id] = total_count
    return found


def _count_value(value: float, ndigits: int = 6):
    """Mirror detection.py helper: keep whole UMI counts tidy, preserve fractions."""
    value = float(value)
    rounded = round(value)
    if np.isclose(value, rounded):
        return int(rounded)
    return round(value, ndigits)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def simple_count_setup():
    """Three cells, three genes: host_A, viral_X, viral_Y.

    UMI counts:
        host_A  : [10, 20, 5]  — never viral
        viral_X : [0,  0,  5]  — total = 5
        viral_Y : [0,  0,  0]  — total = 0
    """
    var_names = ["host_A", "viral_X", "viral_Y"]
    counts = np.array(
        [
            [10, 0, 0],
            [20, 0, 0],
            [5, 5, 0],
        ],
        dtype=float,
    )
    viral_accessions = {"viral_X", "viral_Y"}
    return var_names, counts, viral_accessions


# ---------------------------------------------------------------------------
# §2.1 regression tests
# ---------------------------------------------------------------------------


class TestDetectionThreshold:
    """Audit §2.1: detection_threshold filters genes correctly using >=.

    Regression for: audits/2026-05-08-full-pipeline.md §2.1
    """

    def test_gene_with_zero_counts_never_detected(self, simple_count_setup) -> None:
        """
        GIVEN: a viral gene with total_count == 0
        WHEN:  threshold is 1 (default)
        THEN:  gene is NOT in found_genes

        Regression for: audits/2026-05-08-full-pipeline.md §2.1
        """
        var_names, counts, viral_accessions = simple_count_setup
        found = _detect_genes(var_names, counts, viral_accessions, threshold=1)
        assert "viral_Y" not in found, (
            "viral_Y has total_count=0 and must never be detected (would be FP)."
        )

    def test_gene_at_threshold_is_detected(self, simple_count_setup) -> None:
        """
        GIVEN: a viral gene with total_count == threshold
        WHEN:  the detection function is called
        THEN:  gene IS in found_genes ('>=' comparison, inclusive)

        This guards against an accidental change to strict '>' which would
        silently drop genes at exactly the threshold (common off-by-one).
        """
        var_names, counts, viral_accessions = simple_count_setup
        # viral_X has total_count=5; set threshold to 5 → still detected
        found = _detect_genes(var_names, counts, viral_accessions, threshold=5)
        assert "viral_X" in found, (
            "viral_X has total_count==threshold (5); must be detected (inclusive >=). "
            "Bug: comparison may have been changed from '>=' to '>'."
        )
        assert found["viral_X"] == 5.0

    def test_gene_below_threshold_excluded(self, simple_count_setup) -> None:
        """
        GIVEN: a viral gene with total_count == 5
        WHEN:  threshold is 6
        THEN:  gene is NOT in found_genes
        """
        var_names, counts, viral_accessions = simple_count_setup
        found = _detect_genes(var_names, counts, viral_accessions, threshold=6)
        assert "viral_X" not in found, (
            f"viral_X has total_count=5 < threshold=6; must not be detected. Got: {found}"
        )

    def test_host_genes_never_in_found_even_when_high_count(self, simple_count_setup) -> None:
        """
        GIVEN: a gene with high UMI counts that is NOT in viral_accessions
        WHEN:  detection is run
        THEN:  it is never in found_genes regardless of count

        Regression for: audits/2026-05-08-full-pipeline.md §2.1
        """
        var_names, counts, viral_accessions = simple_count_setup
        # host_A has total_count=35 but is not in viral_accessions
        found = _detect_genes(var_names, counts, viral_accessions, threshold=1)
        assert "host_A" not in found, (
            "host_A is not a viral accession and must never appear in found_genes."
        )

    def test_unknown_viral_accession_silently_skipped(self, simple_count_setup) -> None:
        """
        GIVEN: viral_accessions contains an ID not in var_names
        WHEN:  detection is run
        THEN:  the unknown ID is absent from found_genes (no KeyError)
        """
        var_names, counts, viral_accessions = simple_count_setup
        accessions_with_unknown = {"viral_X", "viral_UNKNOWN_123"}
        found = _detect_genes(var_names, counts, accessions_with_unknown, threshold=1)
        assert "viral_UNKNOWN_123" not in found
        assert "viral_X" in found

    def test_sparse_input_handled_identically(self, simple_count_setup) -> None:
        """
        GIVEN: counts_matrix is a scipy sparse matrix (csr)
        WHEN:  detection is run
        THEN:  results are identical to the dense matrix case
        """
        var_names, counts, viral_accessions = simple_count_setup
        sparse_counts = sp.csr_matrix(counts)
        found_dense = _detect_genes(var_names, counts, viral_accessions, threshold=1)
        found_sparse = _detect_genes(var_names, sparse_counts, viral_accessions, threshold=1)
        assert found_dense == found_sparse, (
            f"Dense and sparse inputs give different results: {found_dense} vs {found_sparse}"
        )

    def test_total_count_value_is_sum_across_all_cells(self) -> None:
        """
        GIVEN: viral_X has UMIs [1, 2, 3] in three cells
        WHEN:  detection is run
        THEN:  found_genes['viral_X'] == 6.0 (sum across all cells)
        """
        var_names = ["viral_X"]
        counts = np.array([[1], [2], [3]], dtype=float)
        found = _detect_genes(var_names, counts, {"viral_X"}, threshold=1)
        assert "viral_X" in found
        assert found["viral_X"] == 6.0, (
            f"Expected total_count=6.0 (sum across cells), got {found['viral_X']}"
        )


class TestCellTypeEnrichment:
    """Task 2: verify enrichment math and adjusted p-values schema.

    Uses the real ``cell_type_enrichment`` from ``viralscan.enrichment``
    (extracted from detection.py so it can be imported without Snakemake globals).
    """

    def test_enrichment_outputs_expected_columns(self, tmp_path) -> None:
        import anndata as ad

        var_names = ["virus_a", "host"]
        counts = np.array([[3.0, 0.0], [1.0, 5.0], [0.0, 2.0], [0.0, 4.0]])
        barcodes = ["BC1", "BC2", "BC3", "BC4"]
        adata = ad.AnnData(
            X=counts,
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=var_names),
        )

        labels_df = pd.DataFrame({"barcode": barcodes, "cell_type": ["T", "T", "B", "B"]})
        csv_path = tmp_path / "cell_types.csv"
        labels_df.to_csv(csv_path, index=False)

        cfg = {"cell_types": str(csv_path)}
        group_by_virus = {"VirusA": ["virus_a"]}

        result = cell_type_enrichment(adata, group_by_virus, cfg)
        assert len(result) == 2
        assert set(result.columns) == {
            "virus",
            "cell_type",
            "n_infected",
            "n_total",
            "pct",
            "OR",
            "pvalue",
            "padj",
        }

    def test_bh_adjust_is_monotonic_and_bounded(self) -> None:
        pvals = [0.001, 0.01, 0.2, 0.8]
        adj = _bh_adjust(pvals)
        assert len(adj) == len(pvals)
        assert np.all(adj >= 0.0)
        assert np.all(adj <= 1.0)
        assert np.all(np.diff(np.sort(adj)) >= 0.0)

    def test_cell_type_without_infected_cells_has_zero_pct(self, tmp_path) -> None:
        import anndata as ad

        var_names = ["virus_a"]
        counts = np.array([[2.0], [0.0], [0.0]])
        barcodes = ["BC1", "BC2", "BC3"]
        adata = ad.AnnData(
            X=counts,
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=var_names),
        )

        labels_df = pd.DataFrame({"barcode": barcodes, "cell_type": ["Mono", "B", "B"]})
        csv_path = tmp_path / "cell_types.csv"
        labels_df.to_csv(csv_path, index=False)

        cfg = {"cell_types": str(csv_path)}
        group_by_virus = {"VirusA": ["virus_a"]}

        result = cell_type_enrichment(adata, group_by_virus, cfg)
        b_row = result[result["cell_type"] == "B"].iloc[0]
        assert b_row["n_infected"] == 0
        assert b_row["pct"] == 0.0


class TestFractionalMultimapCounts:
    def test_whole_counts_stay_integer_like(self) -> None:
        assert _count_value(5.0) == 5

    def test_fractional_counts_are_not_truncated(self) -> None:
        assert _count_value(2.5) == 2.5

    def test_fractional_counts_are_rounded_not_floored(self) -> None:
        assert _count_value(1.0 / 3.0) == 0.333333
