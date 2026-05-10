"""Tests for ambiguity-aware multimapper evidence logic."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from viralscan.multimapping import (
    MULTIMAP_EVIDENCE_COLUMNS,
    build_multimap_layers,
    select_detection_matrix,
    should_write_multimap_evidence,
    summarize_multimap_evidence,
)


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
            config={"multimap_method": "equal", "detection_threshold": 1},
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
            {"multimap_method": "equal", "detection_threshold": 2},
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
            {"multimap_method": "equal", "detection_threshold": 2},
        )
        assert result.loc[0, "call_confidence"] == "low_confidence"


class TestDetectionMatrixSelection:
    def test_legacy_primary_call_uses_combined_x(self) -> None:
        import anndata as ad

        adata = ad.AnnData(X=sparse.csr_matrix([[0.0, 2.0]]))
        adata.layers["counts_unique_viral"] = sparse.csr_matrix([[0.0, 0.0]])
        selected = select_detection_matrix(
            adata, {"multimapping": True, "multimap_primary_call": "legacy"}
        )
        assert selected is adata.X

    def test_unique_only_primary_call_uses_unique_viral_layer(self) -> None:
        import anndata as ad

        adata = ad.AnnData(X=sparse.csr_matrix([[0.0, 2.0]]))
        unique = sparse.csr_matrix([[0.0, 0.0]])
        adata.layers["counts_unique_viral"] = unique
        selected = select_detection_matrix(
            adata, {"multimapping": True, "multimap_primary_call": "unique-only"}
        )
        assert selected is unique

    def test_no_multimapping_does_not_write_multimap_evidence(self) -> None:
        assert should_write_multimap_evidence({"multimapping": False}) is False

    def test_multimapping_writes_multimap_evidence(self) -> None:
        assert should_write_multimap_evidence({"multimapping": True}) is True
