"""Tests for ambiguity-aware multimapper evidence logic."""

from __future__ import annotations

import pytest
import numpy as np
import pandas as pd
from scipy import sparse

from viralscan.defaults import DEFAULTS
from viralscan.multimapping import (
    MULTIMAP_EVIDENCE_COLUMNS,
    build_multimap_layers,
    em_gene_abundances,
    select_detection_matrix,
    should_write_multimap_evidence,
    summarize_multimap_evidence,
)
from viralscan.runconfig import RunConfig


def _toy_bus() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"barcode": "BC1", "ec": 0, "count": 10},  # unique host
            {"barcode": "BC1", "ec": 1, "count": 2},  # unique virus
            {"barcode": "BC1", "ec": 2, "count": 4},  # host + virus ambiguous
            {"barcode": "BC2", "ec": 2, "count": 2},  # host + virus ambiguous
            {"barcode": "BC2", "ec": 3, "count": 3},  # viral-only ambiguous
        ]
    )


def _toy_inputs():
    barcode_to_idx = {"BC1": 0, "BC2": 1}
    ec_map = {
        0: [0],
        1: [1],
        2: [0, 1],
        3: [1, 2],
    }
    viral_gene_indices = {1, 2}
    unique_counts = sparse.csr_matrix(
        np.array(
            [
                [10.0, 2.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        )
    )
    return _toy_bus(), barcode_to_idx, ec_map, viral_gene_indices, unique_counts


class TestBuildMultimapLayers:
    def test_equal_method_matches_current_equal_split(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df,
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            method="equal",
            pseudocount=1.0,
        )
        corrected = result.corrected.toarray()
        expected = np.array(
            [
                [2.0, 2.0, 0.0],
                [1.0, 2.5, 1.5],
            ]
        )
        np.testing.assert_allclose(corrected, expected)

    def test_unique_ec_contributes_zero_to_corrected(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df.iloc[:2],
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            method="equal",
            pseudocount=1.0,
        )
        assert result.corrected.sum() == 0

    def test_host_conservative_excludes_viral_share_from_host_virus_ec(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df,
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            method="host-conservative",
            pseudocount=1.0,
        )
        corrected = result.corrected.toarray()
        assert corrected[0, 1] == 0.0
        assert corrected[1, 1] == 1.5
        assert corrected[0, 0] == 2.0
        assert corrected[1, 0] == 1.0

    def test_default_method_is_host_conservative(self) -> None:
        assert DEFAULTS["multimap_method"] == "host-conservative"
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df,
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            pseudocount=1.0,
        )
        corrected = result.corrected.toarray()
        assert corrected[0, 1] == 0.0
        assert corrected[1, 1] == 1.5

    def test_unique_weighted_favors_high_unique_host_evidence(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df[bus_df["ec"] == 2],
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            method="unique-weighted",
            pseudocount=1.0,
        )
        corrected = result.corrected.toarray()
        assert corrected[0, 0] > corrected[0, 1]
        assert corrected[0, 1] > 0

    def test_mass_conserved_for_equal_and_unique_weighted(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        ambiguous_count_sum = 9.0
        for method in ("equal", "unique-weighted"):
            result = build_multimap_layers(
                bus_df,
                barcode_to_idx,
                ec_map,
                n_cells=2,
                n_genes=3,
                viral_gene_indices=viral_gene_indices,
                original_counts=unique_counts,
                method=method,
                pseudocount=1.0,
            )
            assert result.corrected.sum() == ambiguous_count_sum

    def test_duplicate_same_gene_ec_preserves_legacy_equal_split(self) -> None:
        bus_df = pd.DataFrame([{"barcode": "BC1", "ec": 0, "count": 4}])
        result = build_multimap_layers(
            bus_df,
            {"BC1": 0},
            {0: [1, 1]},
            n_cells=1,
            n_genes=2,
            viral_gene_indices={1},
            original_counts=sparse.csr_matrix(np.zeros((1, 2))),
            method="equal",
            pseudocount=1.0,
        )
        corrected = result.corrected.toarray()
        assert corrected[0, 1] == 4.0
        assert result.unique_viral.sum() == 0.0

    def test_mixed_duplicate_host_virus_ec_preserves_selected_mass_and_upper_bound(self) -> None:
        bus_df = pd.DataFrame([{"barcode": "BC1", "ec": 0, "count": 6}])
        result = build_multimap_layers(
            bus_df,
            {"BC1": 0},
            {0: [0, 1, 1]},
            n_cells=1,
            n_genes=2,
            viral_gene_indices={1},
            original_counts=sparse.csr_matrix(np.zeros((1, 2))),
            method="equal",
            pseudocount=1.0,
        )
        assert result.corrected[0, 0] == 2.0
        assert result.corrected[0, 1] == 4.0
        assert result.host_viral_ambiguous[0, 1] == 6.0
        assert result.host_viral_selected[0, 1] == 4.0
        assert result.viral_ambiguous_upper[0, 1] == 6.0


class TestMultimapEvidenceSummary:
    def test_evidence_columns_are_stable(self) -> None:
        empty = summarize_multimap_evidence(
            adata=None,
            group_by_virus={},
            config=RunConfig(multimap_method="equal", detection_threshold=1),
        )
        assert list(empty.columns) == MULTIMAP_EVIDENCE_COLUMNS

    def test_confidence_tiers(self) -> None:
        import anndata as ad

        adata = ad.AnnData(
            X=sparse.csr_matrix(np.zeros((2, 4))),
            obs=pd.DataFrame(index=["BC1", "BC2"]),
            var=pd.DataFrame(index=["virus_strong", "virus_ambiguous", "virus_low", "virus_none"]),
        )
        adata.layers["counts_unique_viral"] = sparse.csr_matrix(
            np.array([[2.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        )
        adata.layers["counts_corrected"] = sparse.csr_matrix(
            np.array([[0.0, 2.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        )
        adata.layers["counts_host_viral_ambiguous"] = sparse.csr_matrix(
            np.array([[0.0, 0.0, 2.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        )
        adata.layers["counts_host_viral_selected"] = sparse.csr_matrix(
            np.array([[0.0, 0.0, 2.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        )
        group_by_virus = {
            "StrongVirus": ["virus_strong"],
            "AmbiguousVirus": ["virus_ambiguous"],
            "LowVirus": ["virus_low"],
            "NoVirus": ["virus_none"],
        }
        result = summarize_multimap_evidence(
            adata,
            group_by_virus,
            RunConfig(multimap_method="equal", detection_threshold=2),
        )
        tiers = dict(zip(result["virus_name"], result["call_confidence"]))
        assert tiers["StrongVirus"] == "strong"
        assert tiers["AmbiguousVirus"] == "ambiguous"
        assert tiers["LowVirus"] == "low_confidence"
        assert tiers["NoVirus"] == "not_detected"

    def test_host_virus_only_equal_split_signal_is_low_confidence(self) -> None:
        import anndata as ad

        adata = ad.AnnData(
            X=sparse.csr_matrix(np.zeros((1, 1))),
            obs=pd.DataFrame(index=["BC1"]),
            var=pd.DataFrame(index=["virus_low"]),
        )
        adata.layers["counts_unique_viral"] = sparse.csr_matrix([[0.0]])
        adata.layers["counts_corrected"] = sparse.csr_matrix([[2.0]])
        adata.layers["counts_host_viral_ambiguous"] = sparse.csr_matrix([[4.0]])
        adata.layers["counts_host_viral_selected"] = sparse.csr_matrix([[2.0]])
        result = summarize_multimap_evidence(
            adata,
            {"LowVirus": ["virus_low"]},
            RunConfig(multimap_method="equal", detection_threshold=2),
        )
        assert result.loc[0, "call_confidence"] == "low_confidence"


class TestDetectionMatrixSelection:
    def test_legacy_primary_call_uses_combined_x(self) -> None:
        import anndata as ad

        adata = ad.AnnData(X=sparse.csr_matrix([[0.0, 2.0]]))
        adata.layers["counts_unique_viral"] = sparse.csr_matrix([[0.0, 0.0]])
        selected = select_detection_matrix(
            adata, RunConfig(multimapping=True, multimap_primary_call="legacy")
        )
        assert selected is adata.X

    def test_unique_only_primary_call_uses_unique_viral_layer(self) -> None:
        import anndata as ad

        adata = ad.AnnData(X=sparse.csr_matrix([[0.0, 2.0]]))
        unique = sparse.csr_matrix([[0.0, 0.0]])
        adata.layers["counts_unique_viral"] = unique
        selected = select_detection_matrix(
            adata, RunConfig(multimapping=True, multimap_primary_call="unique-only")
        )
        assert selected is unique

    def test_no_multimapping_does_not_write_multimap_evidence(self) -> None:
        assert should_write_multimap_evidence(RunConfig(multimapping=False)) is False

    def test_multimapping_writes_multimap_evidence(self) -> None:
        assert should_write_multimap_evidence(RunConfig(multimapping=True)) is True


class TestEMMultimapper:
    def test_em_abundances_squeeze_zero_unique_gene_toward_abundant_host(self) -> None:
        # Gene 0 (host) has strong unique support; gene 1 (viral) has none.
        # An ambiguous EC {0,1} should converge to put almost all mass on gene 0.
        theta = em_gene_abundances(
            ec_counts={(0, 1): 10.0},
            unique_per_gene=np.array([100.0, 0.0]),
            pseudocount=1.0,
            max_iter=100,
            tol=1e-9,
        )
        assert theta[0] > theta[1]
        # the viral gene's converged abundance is far below the equal-split value (5)
        assert theta[1] < 1.0

    def test_em_abundances_split_proportional_for_equal_unique(self) -> None:
        # Two genes with identical unique support split an ambiguous EC ~evenly.
        theta = em_gene_abundances(
            ec_counts={(0, 1): 8.0},
            unique_per_gene=np.array([5.0, 5.0]),
            pseudocount=1.0,
            max_iter=100,
            tol=1e-12,
        )
        np.testing.assert_allclose(theta[0], theta[1])

    def test_em_method_conserves_multimapper_mass(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df,
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            method="em",
        )
        # multi-gene ECs in the toy set: ec2 (count 4+2) + ec3 (count 3) = 9 total
        assert result.corrected.sum() == pytest.approx(9.0)

    def test_em_downweights_gene_with_no_unique_support_vs_equal(self) -> None:
        # Gene 2 (viral, zero unique support) shares ec3 with the better-supported
        # gene 1. EM should credit gene 2 LESS than the naive equal split, because
        # the iterated abundances favor the gene with evidence.
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        kwargs = dict(
            barcode_to_idx=barcode_to_idx,
            ec_map=ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
        )
        em = build_multimap_layers(bus_df, method="em", **kwargs).corrected.toarray()
        equal = build_multimap_layers(bus_df, method="equal", **kwargs).corrected.toarray()
        # column 2 = gene with zero unique support
        assert em[:, 2].sum() < equal[:, 2].sum()
