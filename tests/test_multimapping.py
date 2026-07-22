"""Tests for ambiguity-aware multimapper evidence logic."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.defaults import DEFAULTS
from viralscan.multimapping import (
    MULTIMAP_EVIDENCE_COLUMNS,
    build_multimap_layers,
    em_cell_abundances,
    em_gene_abundances,
    resolve_cb_umi_molecules,
    select_detection_matrix,
    should_write_multimap_evidence,
    summarize_multimap_evidence,
)
from viralscan.runconfig import RunConfig


def _toy_bus() -> pd.DataFrame:
    rows = []
    for barcode, ec, molecules in [
        ("BC1", 0, 10),
        ("BC1", 1, 2),
        ("BC1", 2, 4),
        ("BC2", 2, 2),
        ("BC2", 3, 3),
    ]:
        for i in range(molecules):
            rows.append({"barcode": barcode, "umi": f"U{ec}_{i}", "ec": ec, "count": 1})
    return pd.DataFrame(rows)


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
        # Mixed host-virus molecules retain their full mass on compatible host genes.
        assert corrected[0, 0] == 4.0
        assert corrected[1, 0] == 2.0

    def test_default_method_is_host_conservative(self) -> None:
        assert DEFAULTS["multimap_method"] == "host-conservative"
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        common = dict(
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            pseudocount=1.0,
        )
        default_result = build_multimap_layers(bus_df, barcode_to_idx, ec_map, **common)
        explicit = build_multimap_layers(
            bus_df, barcode_to_idx, ec_map, method="host-conservative", **common
        )
        # Calling with no method selects the default (host-conservative), which keeps
        # host-virus ambiguous mass off viral genes.
        np.testing.assert_allclose(default_result.corrected.toarray(), explicit.corrected.toarray())

    def test_unique_weighted_favors_high_unique_host_evidence(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        result = build_multimap_layers(
            bus_df,
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

    def test_unique_layer_is_built_from_same_resolved_molecule_stream(self) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        # An incompatible pre-v3 matrix must not influence v3 unique counts.
        incompatible = sparse.csr_matrix(np.full((2, 3), 999.0))
        result = build_multimap_layers(
            bus_df,
            barcode_to_idx,
            ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=incompatible,
            method="host-conservative",
        )
        assert result.unique.sum() == result.audit.unique_molecules == 12
        assert result.unique[0, 0] == 10
        assert result.unique[0, 1] == 2

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

    def test_same_gene_multitranscript_ec_is_unique_gene_evidence(self) -> None:
        bus_df = pd.DataFrame([{"barcode": "BC1", "umi": "U1", "ec": 0, "count": 4}])
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
        assert corrected[0, 1] == 0.0
        assert result.unique_viral.sum() == 1.0
        assert result.audit.ignored_read_multiplicity == 3

    def test_mixed_duplicate_host_virus_ec_preserves_selected_mass_and_upper_bound(self) -> None:
        bus_df = pd.DataFrame([{"barcode": "BC1", "umi": "U1", "ec": 0, "count": 6}])
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
        assert result.corrected[0, 0] == 0.5
        assert result.corrected[0, 1] == 0.5
        assert result.host_viral_ambiguous[0, 1] == 1.0
        assert result.host_viral_selected[0, 1] == 0.5
        assert result.viral_ambiguous_upper[0, 1] == 1.0

    def test_multiple_ecs_intersect_and_disjoint_collision_is_unresolved(self) -> None:
        bus = pd.DataFrame(
            [
                {"barcode": "BC1", "umi": "U1", "ec": 0, "count": 8},
                {"barcode": "BC1", "umi": "U1", "ec": 1, "count": 3},
                {"barcode": "BC1", "umi": "U2", "ec": 2, "count": 1},
                {"barcode": "BC1", "umi": "U2", "ec": 3, "count": 1},
            ]
        )
        molecules, audit = resolve_cb_umi_molecules(
            bus, {"BC1": 0}, {0: [0, 1], 1: [1, 2], 2: [0], 3: [2]}
        )
        assert molecules == [(0, (1,))]
        assert audit.input_molecules == 2
        assert audit.unique_molecules == 1
        assert audit.unresolved_molecules == 1
        assert audit.ignored_read_multiplicity == 9

    def test_streamed_bus_is_order_and_buffer_size_invariant(self, tmp_path) -> None:
        bus_df, barcode_to_idx, ec_map, viral_gene_indices, unique_counts = _toy_inputs()
        path = tmp_path / "output.bus.txt"
        ordered = bus_df.sort_values(["barcode", "umi", "ec"], kind="stable")
        path.write_text(
            "".join(
                f"{row.barcode}\t{row.umi}\t{row.ec}\t{row.count}\n"
                for row in ordered.itertuples(index=False)
            )
        )
        kwargs = dict(
            barcode_to_idx=barcode_to_idx,
            ec_map=ec_map,
            n_cells=2,
            n_genes=3,
            viral_gene_indices=viral_gene_indices,
            original_counts=unique_counts,
            method="equal",
        )
        expected = build_multimap_layers(bus_df.sample(frac=1, random_state=7), **kwargs)
        small = build_multimap_layers(path, bus_buffer_size=2, **kwargs)
        large = build_multimap_layers(path, bus_buffer_size=8192, **kwargs)
        np.testing.assert_allclose(small.corrected.toarray(), expected.corrected.toarray())
        np.testing.assert_allclose(large.corrected.toarray(), expected.corrected.toarray())
        assert small.audit == large.audit == expected.audit

    def test_streamed_bus_rejects_unsorted_molecule_keys(self, tmp_path) -> None:
        path = tmp_path / "unsorted.bus.txt"
        path.write_text("BC1\tU2\t0\t1\nBC1\tU1\t0\t1\n")
        with pytest.raises(ValueError, match="sorted"):
            build_multimap_layers(
                path,
                {"BC1": 0},
                {0: [0, 1]},
                n_cells=1,
                n_genes=2,
                viral_gene_indices=set(),
                original_counts=sparse.csr_matrix((1, 2)),
                method="equal",
            )


class TestMultimapEvidenceSummary:
    def test_evidence_columns_are_stable(self) -> None:
        empty = summarize_multimap_evidence(
            adata=None,
            group_by_virus={},
            config=RunConfig(multimap_method="equal", detection_threshold=1),
        )
        assert list(empty.columns) == MULTIMAP_EVIDENCE_COLUMNS

    def test_molecule_evidence_tiers_do_not_claim_validation(self) -> None:
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
        tiers = dict(zip(result["virus_name"], result["evidence_tier"]))
        assert tiers["StrongVirus"] == "candidate_unique"
        assert tiers["AmbiguousVirus"] == "candidate_virus_ambiguous"
        assert tiers["LowVirus"] == "candidate_host_virus_ambiguous"
        assert tiers["NoVirus"] == "not_detected"
        assert not {"probable", "strong"}.intersection(tiers.values())

    def test_host_virus_only_signal_is_ambiguity_labeled(self) -> None:
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
        assert result.loc[0, "evidence_tier"] == "candidate_host_virus_ambiguous"


class TestDetectionMatrixSelection:
    def test_v3_primary_call_uses_selected_method_x(self) -> None:
        import anndata as ad

        adata = ad.AnnData(X=sparse.csr_matrix([[0.0, 2.0]]))
        adata.layers["counts_unique_viral"] = sparse.csr_matrix([[0.0, 0.0]])
        selected = select_detection_matrix(
            adata, RunConfig(multimapping=True, multimap_primary_call="selected-method")
        )
        assert selected is adata.X

    def test_internal_config_cannot_switch_away_from_complete_x(self) -> None:
        import anndata as ad

        adata = ad.AnnData(X=sparse.csr_matrix([[0.0, 2.0]]))
        adata.layers["counts_unique_viral"] = sparse.csr_matrix([[0.0, 0.0]])
        selected = select_detection_matrix(
            adata, RunConfig(multimapping=True, multimap_primary_call="unique-only")
        )
        assert selected is adata.X

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

    def test_em_abundances_vectorised_matches_reference_loop(self) -> None:
        # The vectorised (sparse mat-vec) EM must produce numerically identical
        # abundances to the original per-EC Python loop, for arbitrary inputs
        # including degenerate (zero-theta) ECs. Reference implementation inlined.
        def _reference_em(ec_counts, unique_per_gene, pseudocount, max_iter, tol):
            unique_per_gene = np.asarray(unique_per_gene, dtype=float).reshape(-1)
            theta = unique_per_gene + float(pseudocount)
            items = [(np.asarray(g, dtype=int), float(c)) for g, c in ec_counts.items()]
            for _ in range(int(max_iter)):
                new = unique_per_gene.copy()
                for genes, count in items:
                    w = theta[genes]
                    s = float(w.sum())
                    if s <= 1e-12:
                        new[genes] += count / len(genes)
                    else:
                        new[genes] += count * w / s
                denom = float(theta.sum()) or 1.0
                if float(np.abs(new - theta).sum()) / denom < tol:
                    theta = new
                    break
                theta = new
            return theta

        rng = np.random.default_rng(0)
        n_genes = 40
        for _trial in range(5):
            unique = rng.integers(0, 50, size=n_genes).astype(float)
            # include a degenerate EC over genes with zero unique support + zero pseudocount path
            ec_counts = {}
            for _ in range(60):
                k = int(rng.integers(2, 5))
                genes = tuple(sorted(set(int(x) for x in rng.integers(0, n_genes, size=k))))
                if len(genes) >= 2:
                    ec_counts[genes] = ec_counts.get(genes, 0.0) + float(rng.integers(1, 20))
            kwargs = dict(unique_per_gene=unique, pseudocount=1.0, max_iter=100, tol=1e-10)
            ref = _reference_em(ec_counts=ec_counts, **kwargs)
            got = em_gene_abundances(ec_counts=ec_counts, **kwargs)
            np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-9)

    def test_em_abundances_empty_returns_unique_plus_pseudocount(self) -> None:
        theta = em_gene_abundances(
            {}, np.array([3.0, 7.0]), pseudocount=1.0, max_iter=100, tol=1e-9
        )
        np.testing.assert_allclose(theta, np.array([4.0, 8.0]))

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
            method="em-global",
        )
        # multi-gene ECs in the toy set: ec2 (count 4+2) + ec3 (count 3) = 9 total
        assert result.corrected.sum() == pytest.approx(9.0)
        assert result.method_diagnostics["global_model"]["converged"] is True
        assert result.method_diagnostics["global_model"]["iterations"] >= 1

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
        em = build_multimap_layers(bus_df, method="em-global", **kwargs).corrected.toarray()
        equal = build_multimap_layers(bus_df, method="equal", **kwargs).corrected.toarray()
        # column 2 = gene with zero unique support
        assert em[:, 2].sum() < equal[:, 2].sum()

    def test_em_cell_uses_local_evidence_with_global_shrinkage(self) -> None:
        theta = em_cell_abundances(
            {(0, 1): 2.0},
            unique_per_gene=np.array([10.0, 0.0]),
            global_theta=np.array([1.0, 10.0]),
            prior_strength=1.0,
            max_iter=100,
            tol=1e-9,
        )
        assert theta[0] > theta[1]
        assert theta.sum() == pytest.approx(13.0)

    def test_em_cell_without_local_compatible_support_remains_ambiguous(self) -> None:
        diagnostics = {}
        theta = em_cell_abundances(
            {(0, 1): 2.0},
            unique_per_gene=np.array([0.0, 0.0, 10.0]),
            global_theta=np.array([10.0, 1.0, 1.0]),
            prior_strength=1.0,
            max_iter=100,
            tol=1e-9,
            diagnostics=diagnostics,
        )
        np.testing.assert_array_equal(theta, np.zeros(3))
        assert diagnostics == {
            "iterations": 0,
            "converged": True,
            "fallback_reason": "no_local_compatible_unique_support",
        }

    def test_em_cell_unsupported_fallback_conserves_equal_ambiguity(self) -> None:
        bus = pd.DataFrame(
            [{"barcode": "BC1", "umi": f"U{i}", "ec": 0, "count": 1} for i in range(10)]
            + [{"barcode": "BC2", "umi": "U0", "ec": 1, "count": 1}]
        )
        result = build_multimap_layers(
            bus,
            {"BC1": 0, "BC2": 1},
            {0: [0], 1: [0, 1]},
            n_cells=2,
            n_genes=2,
            viral_gene_indices={1},
            original_counts=sparse.csr_matrix((2, 2)),
            method="em-cell",
        )
        np.testing.assert_allclose(result.corrected.toarray()[1], [0.5, 0.5])
        assert result.corrected.sum() == pytest.approx(result.audit.ambiguous_molecules)
        assert result.method_diagnostics["cell_models"]["fallback_reasons"] == [
            "no_local_compatible_unique_support"
        ]

    def test_em_reports_deterministic_max_iteration_fallback(self) -> None:
        diagnostics = {}
        em_gene_abundances(
            {(0, 1): 10.0},
            np.array([1.0, 0.0]),
            pseudocount=1.0,
            max_iter=1,
            tol=0.0,
            diagnostics=diagnostics,
        )
        assert diagnostics == {
            "iterations": 1,
            "converged": False,
            "fallback_reason": "max_iterations_reached",
        }
