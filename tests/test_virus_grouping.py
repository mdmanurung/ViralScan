"""Tests for the unified, boundary-aware viral gene-to-virus grouping.

Replaces the two drifted, buggy rules previously inlined in detection.py
(substring) and umap.py (prefix). Guards the specific real-data divergence
cases that motivated the fix.
"""

from __future__ import annotations

from viralscan.anellovirus import merged_name_map
from viralscan.constants import VIRUS_NAME_MAP
from viralscan.virus_grouping import group_genes_by_virus, virus_name_for_gene

# A synthetic map giving full control over boundary semantics.
SYN = {"AICHI": "Aichi virus", "TTV": "Torque teno virus", "HHV6": "HHV-6", "HHV6B": "HHV-6B"}


class TestBoundaryMatch:
    def test_exact_key_matches(self) -> None:
        assert virus_name_for_gene("AICHI", SYN) == "Aichi virus"

    def test_underscore_boundary_matches(self) -> None:
        assert virus_name_for_gene("AICHI_gp1", SYN) == "Aichi virus"

    def test_digit_boundary_matches(self) -> None:
        # TTV7_gp2 must group as TTV — the old prefix rule (key + "_") missed this.
        assert virus_name_for_gene("TTV7_gp2", SYN) == "Torque teno virus"

    def test_alpha_continuation_does_not_match(self) -> None:
        # 'AICHIX...' is NOT Aichi: the char after the key is a letter, not a boundary.
        assert virus_name_for_gene("AICHIX_gp1", SYN) == "AICHIX_gp1"

    def test_longest_key_wins(self) -> None:
        # HHV6B must beat HHV6 (HHV6 followed by 'B' is an alpha continuation anyway).
        assert virus_name_for_gene("HHV6B_U2", SYN) == "HHV-6B"
        assert virus_name_for_gene("HHV6_U2", SYN) == "HHV-6"

    def test_unmatched_returns_raw_gene_id(self) -> None:
        assert virus_name_for_gene("ZZZ_999", SYN) == "ZZZ_999"


class TestRealDataDivergenceCases:
    """The 151-gene divergence that motivated the unification (real VIRUS_NAME_MAP)."""

    def test_substring_overmatch_fixed_borf1(self) -> None:
        # 'EPSTEIN_HHV4_BORF1' must NOT be mislabeled Orf virus (key 'ORF' in 'BORF1').
        assert virus_name_for_gene("EPSTEIN_HHV4_BORF1") == VIRUS_NAME_MAP["EPSTEIN"]

    def test_substring_overmatch_fixed_bunyamw(self) -> None:
        # 'BUNYAMW_...' must be Bunyamwera, not Bunyavirus La Crosse (key 'BUNYA').
        assert virus_name_for_gene("BUNYAMW_BUNVsLgp1") == VIRUS_NAME_MAP["BUNYAMW"]

    def test_prefix_undermatch_fixed_ttv(self) -> None:
        # 'TTV7_gp2' must group as TTV (digit boundary) instead of falling through.
        assert virus_name_for_gene("TTV7_gp2") == VIRUS_NAME_MAP["TTV"]

    def test_nested_key_longest_wins(self) -> None:
        assert virus_name_for_gene("HUM_HERP6B_U63") == VIRUS_NAME_MAP["HUM_HERP6B"]


class TestGroupGenesByVirus:
    def test_groups_and_detected_set(self) -> None:
        genes = ["AICHI_gp1", "AICHI_gp2", "TTV7_gp1", "ZZZ_999"]
        groups, detected = group_genes_by_virus(genes, SYN)
        assert groups["Aichi virus"] == ["AICHI_gp1", "AICHI_gp2"]
        assert groups["Torque teno virus"] == ["TTV7_gp1"]
        assert groups["ZZZ_999"] == ["ZZZ_999"]  # raw fallback
        assert detected == {"Aichi virus", "Torque teno virus"}  # raw fallback excluded

    def test_input_order_preserved(self) -> None:
        groups, _ = group_genes_by_virus(["AICHI_b", "AICHI_a"], SYN)
        assert groups["Aichi virus"] == ["AICHI_b", "AICHI_a"]

    def test_empty_input(self) -> None:
        groups, detected = group_genes_by_virus([], SYN)
        assert groups == {}
        assert detected == set()


class TestMergedNameMapAnellovirus:
    """D.4 — merged_name_map() correctly resolves accession-keyed anellovirus gene IDs.

    Accession keys (e.g. ``"AB026929.1"``) come from the packaged TSV and are added
    by :func:`viralscan.anellovirus.anello_name_map`.  The boundary-aware prefix
    rule in :func:`~viralscan.virus_grouping.virus_name_for_gene` extends them to
    ``{acc}_geneN`` variants without any extra map entries.
    """

    # AB026929.1 is a clareaulab Betatorquevirus in the packaged TSV.
    # NC_002076.2 is a viralscan-refseq Alphatorquevirus.
    # NC_038345.1 is a clareaulab Betatorquevirus in the packaged TSV.

    def test_bare_genbank_accession_resolves(self) -> None:
        merged = merged_name_map()
        assert virus_name_for_gene("AB026929.1", merged) == "Betatorquevirus"

    def test_genbank_accession_gene_suffix_resolves(self) -> None:
        """``{acc}_gene1`` must resolve via boundary-aware prefix (no extra map entry)."""
        merged = merged_name_map()
        assert virus_name_for_gene("AB026929.1_gene1", merged) == "Betatorquevirus"

    def test_refseq_accession_resolves(self) -> None:
        """RefSeq NC_ accessions from the viralscan-refseq source still resolve."""
        merged = merged_name_map()
        assert virus_name_for_gene("NC_002076.2", merged) == "Alphatorquevirus"

    def test_refseq_accession_gene_suffix_resolves(self) -> None:
        merged = merged_name_map()
        assert virus_name_for_gene("NC_002076.2_gene1", merged) == "Alphatorquevirus"

    def test_clareaulab_betatorquevirus_refseq_accession(self) -> None:
        """NC_038345.1 appears in the packaged TSV as Betatorquevirus (clareaulab)."""
        merged = merged_name_map()
        assert virus_name_for_gene("NC_038345.1", merged) == "Betatorquevirus"

    def test_legacy_ttv_prefix_still_works(self) -> None:
        """Legacy TTV-prefix gene_ids from the bundled GTFs still resolve via VIRUS_NAME_MAP."""
        merged = merged_name_map()
        assert virus_name_for_gene("TTV7_gp2", merged) == VIRUS_NAME_MAP["TTV"]

    def test_unmapped_gene_id_returns_itself(self) -> None:
        """Gene IDs with no matching accession or prefix fall back to the raw ID."""
        merged = merged_name_map()
        assert virus_name_for_gene("COMPLETELY_UNKNOWN_XYZ", merged) == "COMPLETELY_UNKNOWN_XYZ"

    def test_merged_map_is_superset_of_virus_name_map(self) -> None:
        """merged_name_map() extends VIRUS_NAME_MAP — all legacy keys survive."""
        merged = merged_name_map()
        for key, value in VIRUS_NAME_MAP.items():
            assert merged[key] == value, f"Legacy key {key!r} missing or overwritten in merged map"

    def test_merged_map_does_not_mutate_virus_name_map(self) -> None:
        """merged_name_map() must return a new dict, not modify VIRUS_NAME_MAP."""
        original_len = len(VIRUS_NAME_MAP)
        _ = merged_name_map()
        assert len(VIRUS_NAME_MAP) == original_len
