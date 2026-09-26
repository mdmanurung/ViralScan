"""Tests for the unified, boundary-aware viral gene-to-virus grouping.

Replaces the two drifted, buggy rules previously inlined in detection.py
(substring) and umap.py (prefix). Guards the specific real-data divergence
cases that motivated the fix.
"""

from __future__ import annotations

import glob
import re

import pytest

from viralscan.anellovirus import merged_name_map
from viralscan.constants import VIRUS_GENE_ID_ALIASES, VIRUS_NAME_MAP
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


class TestAliasTier:
    """The alias tier exists for panel schemes that omit the separator.

    Measured motivation: against the panel the covid PBMC runs actually used
    (``references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf``,
    4,650 gene IDs), the strict boundary rule alone left 2,290 gene IDs (49.2 %)
    resolving to no virus name. Each became its own "virus" row in
    ``viral_summary.tsv`` — the covid run published 9 of 17 rows as bare gene IDs
    (``HHV1gp00p39``, ``CeHV2gUL24``, ``MPXV_gp132``), so per-virus aggregation,
    ``accession_breadth``, sibling cross-mapping and ``eve_risk`` were all
    computed per gene instead of per virus.
    """

    def test_concatenated_scheme_resolves(self) -> None:
        assert virus_name_for_gene("HHV1gp00p39") == "Human herpesvirus 1"
        assert virus_name_for_gene("HHV2p01") == "Human herpesvirus 2"
        assert virus_name_for_gene("MOCVgp001") == "Molluscum contagiosum virus"
        assert virus_name_for_gene("Ydvgp129") == "Yaba-like disease virus"

    def test_embedded_strain_token_resolves(self) -> None:
        # 'HHV5wt' = HHV-5 wild type; 'HHV8GK18' = HHV-8 Kaposi strain GK18.
        assert virus_name_for_gene("HHV5wtgp045") == "Human cytomegalovirus"
        assert virus_name_for_gene("HHV8GK18_gp56") == "Human herpesvirus 8"

    def test_bundled_ttv_gp_ids_resolve(self) -> None:
        """Regression: the 20 bundled RefSeq TTV GTFs use ``TTVgp1/2/3``.

        The boundary rule accepts a digit or ``_`` after the key but not a
        letter, so ``TTVgp1`` fell through to the raw-ID fallback while its
        sibling ``TTV_TTVgp1`` resolved — one genome reported two ways.
        """
        assert virus_name_for_gene("TTVgp1") == VIRUS_NAME_MAP["TTV"]
        assert virus_name_for_gene("TTVgp2") == VIRUS_NAME_MAP["TTV"]
        assert virus_name_for_gene("TTV7_gp2") == VIRUS_NAME_MAP["TTV"]

    def test_non_anellovirus_accession_keyed_gene_id_resolves(self) -> None:
        """``ncbi_fetch`` emits ``{acc}_geneN``; SARS-CoV-2 has no anello entry."""
        assert virus_name_for_gene("NC_045512.2_gene1") == "SARS coronavirus 2"

    def test_longest_alias_wins_over_shorter_sibling(self) -> None:
        assert virus_name_for_gene("HHV6B_something") == "Human herpesvirus 6b"
        assert virus_name_for_gene("HHV6_something") == "Human herpesvirus 6"

    def test_alias_never_overrides_a_strict_match(self) -> None:
        """The alias tier runs only after tier 1 misses, so it is purely additive.

        ``HUM_HERP6B`` is a tier-1 key and must keep winning over the ``HHV6B``
        alias even though the alias is what the concatenated panel needs.
        """
        assert virus_name_for_gene("HUM_HERP6B_U63") == VIRUS_NAME_MAP["HUM_HERP6B"]
        assert virus_name_for_gene("CERC_HERP_CeHV2gUL24") == "Cercopithecine herpesvirus"

    def test_strict_boundary_rule_is_not_weakened(self) -> None:
        """The over-match guards must survive the alias tier.

        The alias mechanism exists precisely so the boundary rule can stay
        strict; loosening it would re-open both original bugs.
        """
        assert virus_name_for_gene("AICHIX_gp1", SYN, aliases={}) == "AICHIX_gp1"
        assert virus_name_for_gene("EPSTEIN_HHV4_BORF1") == VIRUS_NAME_MAP["EPSTEIN"]
        assert virus_name_for_gene("BUNYAMW_BUNVsLgp1") == VIRUS_NAME_MAP["BUNYAMW"]

    def test_ambiguous_panel_tokens_keep_raw_fallback(self) -> None:
        """A token with no confident virus assignment must not be guessed.

        ``QKL08`` is a bare GenBank prefix in the Serratus panel. Inventing a
        name would be worse than the raw ID, which stays visible in the output.
        """
        assert virus_name_for_gene("QKL08_gp1") == "QKL08_gp1"

    def test_alias_tier_can_be_disabled(self) -> None:
        """``aliases={}`` restores pre-alias behaviour exactly (test harness use)."""
        assert virus_name_for_gene("TTVgp1", merged_name_map(), aliases={}) == "TTVgp1"
        assert virus_name_for_gene("TTVgp1", merged_name_map()) == "Torque teno virus"

    def test_aliases_are_not_empty_and_values_are_known_names(self) -> None:
        assert VIRUS_GENE_ID_ALIASES
        for prefix, name in VIRUS_GENE_ID_ALIASES.items():
            assert prefix and name
            assert not prefix.endswith("_"), (
                f"{prefix!r} has a trailing separator; use VIRUS_NAME_MAP"
            )


def _panel_gene_ids(gtf: str) -> set[str]:
    with open(gtf) as handle:
        return set(re.findall(r'gene_id "([^"]+)"', handle.read()))


class TestRealPanelResolution:
    """Measure name resolution against the panels the repo actually ships.

    Guards the class of defect rather than one token: a new panel that reverts
    to raw-gene-ID "viruses" must fail here, not in a published
    ``viral_summary.tsv``.
    """

    def _assert_mostly_resolved(self, gene_ids: set[str], floor: float) -> None:
        merged = merged_name_map()
        unresolved = [g for g in gene_ids if virus_name_for_gene(g, merged) == g]
        rate = len(unresolved) / len(gene_ids)
        assert rate <= floor, (
            f"{len(unresolved)}/{len(gene_ids)} gene IDs ({rate:.1%}) resolve to no virus "
            f"name, above the {floor:.1%} budget; e.g. {sorted(unresolved)[:5]}. Either add "
            "the panel's token to VIRUS_NAME_MAP (underscore-bounded) or "
            "VIRUS_GENE_ID_ALIASES (concatenated), or raise the floor deliberately."
        )

    def test_bundled_panel(self) -> None:
        gtfs = glob.glob("src/viralscan/data/*.gtf")
        if not gtfs:
            pytest.skip("bundled viral GTFs are gitignored and absent from this checkout")
        gene_ids: set[str] = set()
        for gtf in gtfs:
            gene_ids |= _panel_gene_ids(gtf)
        assert len(gene_ids) > 2000, "bundled GTFs not found; test would pass vacuously"
        self._assert_mostly_resolved(gene_ids, floor=0.02)

    def test_starsolo_packaged_panel(self) -> None:
        gtf = "references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf"
        try:
            gene_ids = _panel_gene_ids(gtf)
        except OSError:
            pytest.skip(f"{gtf} not present (local build artifact, not shipped)")
        assert len(gene_ids) > 4000, "panel looks truncated; test would pass vacuously"
        # Was 49.2% before the alias tier. The residual 13.5% is panel tokens
        # with no confident virus assignment (QKL08, HRCV, KPV, G128, SCV12).
        self._assert_mostly_resolved(gene_ids, floor=0.15)
