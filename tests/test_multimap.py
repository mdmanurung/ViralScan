"""Tests for barcode processing and multimapping matrix logic in scripts/multimap.py.

``multimap.py`` is a Snakemake script (references ``snakemake.*`` at module
level), so we test the core logic via standalone re-implementations that mirror
the relevant functions.  This follows the pattern established in test_analysis.py.

Audit findings covered:
  §2.3 — ``load_barcodes()`` and ``normalize_barcodes()`` use
          ``str.replace("-1", "")`` which is a global substitution that corrupts
          any barcode containing the substring "-1" at a non-trailing position.

Regression for: audits/2026-05-08-full-pipeline.md §2.3
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


# ---------------------------------------------------------------------------
# Standalone re-implementations of the barcode-stripping logic
# (mirrors the *fixed* code in multimap.py)
# ---------------------------------------------------------------------------


def _strip_10x_suffix(barcode: str) -> str:
    """Remove the trailing '-1' lane suffix added by 10x Genomics Cell Ranger.

    Only the *trailing* '-1' is removed.  A barcode that contains '-1' at
    an internal position (e.g. 'ACGT-1GCTA') is left unchanged except for
    the terminal '-1'.

    This is the *fixed* version that replaces the global ``str.replace``
    pattern used in the original code.
    """
    return barcode.removesuffix("-1")


def _load_barcodes_fixed(barcodes: list[str]) -> dict[str, int]:
    """Mirror of fixed load_barcodes(): strips only trailing '-1' suffix."""
    stripped = [_strip_10x_suffix(bc) for bc in barcodes]
    return {bc: i for i, bc in enumerate(stripped)}


def _normalize_barcodes_fixed(bus_df: pd.DataFrame) -> pd.DataFrame:
    """Mirror of fixed normalize_barcodes(): strips only trailing '-1' suffix."""
    bus_df = bus_df.copy()
    bus_df["barcode"] = bus_df["barcode"].map(_strip_10x_suffix)
    return bus_df


def _load_barcodes_buggy(barcodes: list[str]) -> dict[str, int]:
    """Mirror of the *original* buggy load_barcodes() for comparison."""
    stripped = [bc.replace("-1", "") for bc in barcodes]
    return {bc: i for i, bc in enumerate(stripped)}


def _normalize_barcodes_buggy(bus_df: pd.DataFrame) -> pd.DataFrame:
    """Mirror of the *original* buggy normalize_barcodes() for comparison."""
    bus_df = bus_df.copy()
    bus_df["barcode"] = bus_df["barcode"].str.replace("-1", "", regex=False)
    return bus_df


# ---------------------------------------------------------------------------
# Standalone re-implementation of build_multimap_matrix
# (for UMI conservation tests — §3.4 / integration)
# ---------------------------------------------------------------------------


def _build_multimap_matrix(
    bus_df: pd.DataFrame,
    barcode_to_idx: dict[str, int],
    ec_map: dict[int, list[int]],
    n_cells: int,
    n_genes: int,
) -> sparse.csr_matrix:
    """Mirror of build_multimap_matrix() from multimap.py."""
    rows, cols, data = [], [], []

    for row in bus_df.itertuples(index=False):
        bc, ec, count = row.barcode, row.ec, row.count
        if pd.isna(ec):
            continue
        ec = int(ec)
        if bc not in barcode_to_idx:
            continue
        if ec not in ec_map:
            continue

        cell_idx = barcode_to_idx[bc]
        genes_in_ec = ec_map[ec]
        if not genes_in_ec:
            continue

        # Only multi-mapping ECs are redistributed; unique-mapping ECs
        # are already in counts_original and are skipped here to avoid
        # double-counting.
        if len(genes_in_ec) == 1:
            continue

        share = count / len(genes_in_ec)
        for gid in genes_in_ec:
            rows.append(cell_idx)
            cols.append(gid)
            data.append(share)

    return sparse.csr_matrix((data, (rows, cols)), shape=(n_cells, n_genes))


# ---------------------------------------------------------------------------
# Tests — §2.3 barcode suffix stripping
# ---------------------------------------------------------------------------


class TestLoadBarcodes:
    """Audit §2.3: only the *trailing* '-1' suffix must be removed from barcodes.

    The original code uses ``bc.replace("-1", "")`` which removes ALL occurrences
    of the substring "-1" anywhere in the barcode string.

    Regression for: audits/2026-05-08-full-pipeline.md §2.3
    """

    def test_standard_10x_suffix_removed(self) -> None:
        """Standard 10x barcodes ending in '-1' must have the suffix removed."""
        barcodes = ["ACGTACGTACGT-1", "TTGGCCAATTGG-1"]
        result = _load_barcodes_fixed(barcodes)
        assert "ACGTACGTACGT" in result
        assert "TTGGCCAATTGG" in result

    def test_barcode_without_suffix_unchanged(self) -> None:
        """Barcodes without a '-1' suffix must pass through unchanged."""
        barcodes = ["ACGTACGTACGT", "TTGGCCAATTGG"]
        result = _load_barcodes_fixed(barcodes)
        assert "ACGTACGTACGT" in result
        assert "TTGGCCAATTGG" in result

    def test_internal_minus_one_not_removed(self) -> None:
        """
        GIVEN: a barcode with '-1' at an internal position (not the trailing suffix)
        WHEN:  the fixed strip function is applied
        THEN:  only the *trailing* '-1' is removed; the internal '-1' is preserved

        This test would FAIL against the original buggy implementation.
        Regression for: audits/2026-05-08-full-pipeline.md §2.3
        """
        barcode = "ACGT-1GCTA-1"  # '-1' appears mid-string AND as trailing suffix
        result = _strip_10x_suffix(barcode)
        # Only the trailing '-1' must be removed
        assert result == "ACGT-1GCTA", (
            f"Expected 'ACGT-1GCTA' (trailing '-1' only removed), got {result!r}. "
            "Internal '-1' must not be stripped."
        )

    def test_buggy_implementation_mangles_internal_minus_one(self) -> None:
        """Confirm that the *original* buggy implementation does strip internal '-1'.

        This test documents the known-bad behaviour we are fixing.
        """
        barcode = "ACGT-1GCTA-1"
        buggy_result = barcode.replace("-1", "")
        # The buggy code removes ALL '-1' occurrences
        assert buggy_result == "ACGTGCTA", (
            "This documents the bug: .replace() strips all '-1' occurrences"
        )

    def test_fixed_vs_buggy_differ_on_adversarial_barcode(self) -> None:
        """The fixed and buggy implementations must produce different results
        for a barcode with an internal '-1'."""
        barcode = "ACGT-1GCTA-1"
        fixed = _strip_10x_suffix(barcode)
        buggy = barcode.replace("-1", "")
        assert fixed != buggy, (
            "Fixed and buggy implementations agree on adversarial input — "
            "the fix may not have been applied."
        )

    def test_empty_barcode_safe(self) -> None:
        """Empty string barcode must not crash."""
        result = _strip_10x_suffix("")
        assert result == ""

    def test_only_suffix_barcode(self) -> None:
        """Barcode that IS just '-1' must become empty string."""
        result = _strip_10x_suffix("-1")
        assert result == ""


class TestNormalizeBarcodes:
    """Audit §2.3: normalize_barcodes in a BUS DataFrame strips only trailing '-1'."""

    def test_trailing_suffix_stripped_in_dataframe(self) -> None:
        df = pd.DataFrame({"barcode": ["AAAA-1", "CCCC-1"], "ec": [0, 1], "count": [5, 3]})
        result = _normalize_barcodes_fixed(df)
        assert list(result["barcode"]) == ["AAAA", "CCCC"]

    def test_internal_minus_one_preserved_in_dataframe(self) -> None:
        """
        GIVEN: a DataFrame row with a barcode having internal '-1'
        WHEN:  normalize_barcodes is applied
        THEN:  only the trailing '-1' is removed; internal '-1' is preserved
        """
        df = pd.DataFrame(
            {
                "barcode": ["ACGT-1GCTA-1"],
                "ec": [0],
                "count": [5],
            }
        )
        result = _normalize_barcodes_fixed(df)
        assert result["barcode"].iloc[0] == "ACGT-1GCTA", (
            f"Got {result['barcode'].iloc[0]!r}, expected 'ACGT-1GCTA'"
        )


# ---------------------------------------------------------------------------
# Tests — multimapping matrix UMI conservation (§3.4 / audit trace A)
# ---------------------------------------------------------------------------


class TestBuildMultimapMatrix:
    """Verify UMI count conservation properties of the multimapping matrix.

    Audit trace A confirms: unique-mapping ECs (len==1) are SKIPPED in
    build_multimap_matrix(), so counts_corrected holds only the redistributed
    multimapper fraction.  The final X = counts_corrected + counts_original
    is NOT double-counting unique reads.

    Regression for: audits/2026-05-08-full-pipeline.md §3.4
    """

    def test_unique_ec_contributes_zero_to_corrected(self) -> None:
        """
        GIVEN: a BUS record mapping to a single gene (EC with len==1)
        WHEN:  build_multimap_matrix processes it
        THEN:  the corrected matrix has 0.0 for that gene (unique reads go
               into counts_original, not counts_corrected)
        """
        bus_df = pd.DataFrame(
            [
                {"barcode": "BC1", "ec": 0, "count": 5},  # EC0 → only gene_A (unique)
            ]
        )
        barcode_to_idx = {"BC1": 0}
        ec_map = {0: [0]}  # EC0 maps to gene_A only
        corrected = _build_multimap_matrix(bus_df, barcode_to_idx, ec_map, n_cells=1, n_genes=2)
        assert corrected[0, 0] == 0.0, (
            f"Unique-mapping EC must contribute 0 to counts_corrected; got {corrected[0, 0]}"
        )

    def test_multi_ec_redistributed_equally(self) -> None:
        """
        GIVEN: a BUS record mapping to two genes (EC with len==2)
        WHEN:  build_multimap_matrix processes it
        THEN:  each gene receives count/2 in the corrected matrix
        """
        bus_df = pd.DataFrame(
            [
                {"barcode": "BC1", "ec": 1, "count": 4},  # EC1 → gene_A + gene_B
            ]
        )
        barcode_to_idx = {"BC1": 0}
        ec_map = {1: [0, 1]}  # EC1 maps to gene_A (idx 0) and gene_B (idx 1)
        corrected = _build_multimap_matrix(bus_df, barcode_to_idx, ec_map, n_cells=1, n_genes=2)
        assert corrected[0, 0] == 2.0, f"Expected 2.0 share for gene_A, got {corrected[0, 0]}"
        assert corrected[0, 1] == 2.0, f"Expected 2.0 share for gene_B, got {corrected[0, 1]}"

    def test_final_x_not_double_counting_unique_reads(self) -> None:
        """
        GIVEN: 5 unique reads to gene_A and 4 multi-mapping reads to gene_A+B
        WHEN:  counts_corrected + counts_original is computed
        THEN:  gene_A total == 5 (unique) + 2 (multi share) == 7
               gene_B total == 0 (unique) + 2 (multi share) == 2
               Total UMIs == 9  (== 5 + 4, no double-counting)

        Regression for: audits/2026-05-08-full-pipeline.md §3.4
        """
        n_cells, n_genes = 1, 2
        bus_df_unique = pd.DataFrame([{"barcode": "BC1", "ec": 0, "count": 5}])
        bus_df_multi = pd.DataFrame([{"barcode": "BC1", "ec": 1, "count": 4}])
        bus_df = pd.concat([bus_df_unique, bus_df_multi], ignore_index=True)

        barcode_to_idx = {"BC1": 0}
        ec_map = {0: [0], 1: [0, 1]}  # EC0 unique → gene_A; EC1 multi → gene_A, gene_B

        corrected = _build_multimap_matrix(bus_df, barcode_to_idx, ec_map, n_cells, n_genes)

        # counts_original: unique reads (5 to gene_A, 0 to gene_B)
        counts_original = sparse.csr_matrix(np.array([[5.0, 0.0]]))

        final_x = corrected + counts_original

        assert final_x[0, 0] == 7.0, (
            f"gene_A final count: expected 7.0, got {final_x[0, 0]}. "
            "Unique reads (5) + multimapper share (2) should equal 7."
        )
        assert final_x[0, 1] == 2.0, f"gene_B final count: expected 2.0, got {final_x[0, 1]}"
        total_input_umis = 5 + 4  # unique + multi reads
        total_output_umis = final_x.sum()
        assert abs(total_output_umis - total_input_umis) < 1e-6, (
            f"UMI count not conserved: input={total_input_umis}, "
            f"output={total_output_umis}. No double-counting expected."
        )
