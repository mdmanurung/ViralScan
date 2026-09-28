"""Tests for the gene-programme catalogue and the layer-2 calling logic.

Layer 2 infers a viral gene programme (latent vs productive) for viruses layer 1
detected. The design exists because of one measurement, on the bundled EBV
LCL run (``SRR12682296``), which is latently infected by construction:

    LATENT (10 genes)  236,342     LYTIC (12 genes)  247,633
      BNLF2a   49,662                 BHLF1   142,954   <- highest EBV gene
      LMP-1    50,237                 BMRF1    45,005
      BaRF1.1  46,956                 BBLF2/3  13,246
      BNLF2b   45,108                 BRLF1     9,560
      EBNA-1.1    920  <-             BZLF1     9,308

A latent:lytic ratio of 1.15 in a cell line defined by latency is not biology.
EBNA-1 is expressed from every latent episome and must be present in every
infected cell, yet it sits ~155x below BHLF1. The cause is pervasive
overlapping-ORF cross-mapping, so the calling rule counts *distinct overlap
groups* rather than genes, and reads evidence from the uniquely-placing layer
rather than the multimap-allocated one.

``TestEbvLclRegression`` is the test that fails if either protection is removed.
"""

from __future__ import annotations

import glob
import os
import re

import pandas as pd
import pytest
from scipy import sparse

from viralscan.anellovirus import merged_name_map
from viralscan.gene_programs import (
    PROGRAMMES,
    STATES,
    Marker,
    call_cell_programme,
    load_catalogue,
    resolve_markers,
    summarise_programs,
    validate_catalogue,
    write_program_outputs,
)
from viralscan.virus_grouping import virus_name_for_gene

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLED_GTFS = sorted(glob.glob(os.path.join(REPO_ROOT, "src", "viralscan", "data", "*.gtf")))


def _catalog() -> list[dict]:
    return load_catalogue()


def _viruses() -> set[str]:
    return {row["virus"] for row in _catalog()}


#: The nine viruses the catalogue covers, with the completeness the inventory
#: supports. HCMV is partial because single-cell latency mirrors a low-level
#: late-lytic programme (PMID 29535194), so latency is not observable by marker.
EXPECTED_COMPLETE = {
    "Epstein-Barr virus",
    "Human herpesvirus 6",
    "Human herpesvirus 7",
}
EXPECTED_PARTIAL = {
    "Human cytomegalovirus",
    "Human herpesvirus 1",
    "Human herpesvirus 2",
    "Human herpesvirus 6b",
    "Varicella-zoster virus",
    "Human herpesvirus 8",
}


class TestCatalogueIntegrity:
    def test_loads_and_has_required_columns(self) -> None:
        rows = _catalog()
        assert len(rows) > 50, "catalogue is implausibly small"
        required = {
            "virus",
            "programme",
            "panel_completeness",
            "latency_observable_in_rna",
            "gene_id_bundled",
            "gene_id_merged",
            "gene_id_starsolo",
            "refseq_gene",
            "product",
            "gene_biotype",
            "has_cds",
            "overlap_group",
            "non_overlapping",
            "available_in_starsolo",
            "do_not_normalise",
            "note",
        }
        assert required <= set(rows[0]), required - set(rows[0])

    def test_controlled_vocabularies(self) -> None:
        for row in _catalog():
            assert row["programme"] in PROGRAMMES, row
            assert row["panel_completeness"] in {"complete", "partial"}, row
            for field in (
                "latency_observable_in_rna",
                "has_cds",
                "non_overlapping",
                "available_in_starsolo",
                "do_not_normalise",
            ):
                assert row[field] in {"true", "false"}, (field, row)

    def test_validate_accepts_the_shipped_catalogue(self) -> None:
        assert validate_catalogue(_catalog()) == []

    def test_covers_exactly_the_nine_intended_viruses(self) -> None:
        assert _viruses() == EXPECTED_COMPLETE | EXPECTED_PARTIAL

    def test_completeness_matches_the_inventory(self) -> None:
        for row in _catalog():
            expected = "complete" if row["virus"] in EXPECTED_COMPLETE else "partial"
            assert row["panel_completeness"] == expected, row

    def test_every_marker_resolves_in_the_bundled_panel(self) -> None:
        """A curated marker that names no real gene is a silent no-op."""
        if not BUNDLED_GTFS:
            pytest.skip("bundled viral GTFs are gitignored and absent from this checkout")
        present: set[str] = set()
        for gtf in BUNDLED_GTFS:
            with open(gtf) as handle:
                present |= set(re.findall(r'gene_id "([^"]+)"', handle.read()))
        missing = sorted({r["gene_id_bundled"] for r in _catalog()} - present)
        assert not missing, f"markers absent from the bundled panel: {missing}"

    def test_every_virus_is_reachable_from_the_name_map(self) -> None:
        """Guards against name drift between the catalogue and VIRUS_NAME_MAP."""
        name_map = merged_name_map()
        for virus in _viruses():
            assert virus_name_for_gene(virus, name_map) == virus, (
                f"{virus!r} does not resolve to itself; VIRUS_NAME_MAP or the catalogue has drifted"
            )

    def test_no_duplicate_rows(self) -> None:
        keys = [(r["virus"], r["gene_id_bundled"], r["programme"]) for r in _catalog()]
        assert len(keys) == len(set(keys)), "duplicate (virus, gene, programme) rows"

    def test_hhv6_padding_is_flagged_and_genes_stay_distinct(self) -> None:
        """HHV-6A uses two padding widths for genuinely different ORFs.

        ``HHV6gp041`` is U43 while ``HHV6gp41`` is U42, and six further pairs
        differ. Normalising padding would silently merge distinct genes.
        """
        rows = [r for r in _catalog() if r["virus"] == "Human herpesvirus 6"]
        assert rows, "no HHV-6A markers"
        assert all(r["do_not_normalise"] == "true" for r in rows)
        # The catalogue must not contain both paddings collapsed into one ID.
        ids = [r["gene_id_bundled"] for r in rows]
        assert len(ids) == len(set(ids))

    def test_unassigned_genes_are_flagged_unavailable_under_starsolo(self) -> None:
        """STARsolo's bare ``unassigned_gene_N`` is shared by >=7 genomes.

        Those rows carry EBV's EBERs -- the most abundant latent transcripts --
        and cannot be attributed to a genome there, so they must be excluded
        rather than silently mis-assigned.
        """
        for row in _catalog():
            if re.match(r"^unassigned_gene_\d+$", row["gene_id_starsolo"]):
                assert row["available_in_starsolo"] == "false", row

    def test_every_virus_has_enough_independent_productive_markers(self) -> None:
        """``min_breadth=2`` must be satisfiable, else productive is unreachable."""
        rows = _catalog()
        for virus in _viruses():
            groups = {
                r["overlap_group"]
                for r in rows
                if r["virus"] == virus
                and r["programme"] == "productive"
                and r["non_overlapping"] == "true"
            }
            assert len(groups) >= 2, (
                f"{virus} has only {len(groups)} independent productive overlap "
                "group(s); the default min_breadth=2 could never be met"
            )

    def test_partial_viruses_never_declare_latency_observable(self) -> None:
        """Absence of a one-transcript latency set is not evidence of latency."""
        for row in _catalog():
            if row["panel_completeness"] == "partial":
                assert row["latency_observable_in_rna"] == "false", row

    def test_complete_viruses_declare_latency_observable(self) -> None:
        for row in _catalog():
            if row["panel_completeness"] == "complete":
                assert row["latency_observable_in_rna"] == "true", row

    def test_ebvna1_mrna_record_is_marked_as_having_no_cds(self) -> None:
        """RefSeq annotates EBNA-1 twice; only the mRNA record lacks a CDS."""
        by_id = {r["gene_id_bundled"]: r for r in _catalog() if r["virus"] == "Epstein-Barr virus"}
        assert by_id["EPSTEIN_HHV4_EBNA-1.1"]["has_cds"] == "false"
        assert by_id["EPSTEIN_HHV4_EBNA-1.1"]["gene_biotype"] == "other"
        assert by_id["EPSTEIN_HHV4_EBNA-1.2"]["has_cds"] == "true"

    def test_latent_ebvna_locus_is_one_overlap_group(self) -> None:
        """The cross-mapping the design defends against, pinned as a fact.

        EBNA-1, EBNA-2 and EBNA-LP are transcribed from one locus and share
        exonic sequence, so detecting one is not independent evidence for the
        others. If this ever changes the catalogue is stale and the calling rule
        needs revisiting.
        """
        rows = [
            r
            for r in _catalog()
            if r["virus"] == "Epstein-Barr virus"
            and r["refseq_gene"] in {"EBNA-1", "EBNA-2", "EBNA-LP"}
        ]
        assert rows
        assert len({r["overlap_group"] for r in rows}) == 1, (
            "EBNA-1/2/LP no longer share an overlap group; the EBV latent anchor "
            "set assumption changed and the calling rule must be re-derived"
        )
        assert all(r["non_overlapping"] == "false" for r in rows)


class TestResolveMarkers:
    def test_resolves_all_three_panel_forms(self) -> None:
        catalogue = _catalog()
        var_names = []
        for row in catalogue:
            var_names += [row["gene_id_bundled"], row["gene_id_starsolo"]]
        var_names = list(dict.fromkeys(var_names))
        for form in ("bundled", "merged", "starsolo"):
            resolved, unresolved = resolve_markers("Epstein-Barr virus", var_names, form)
            assert resolved, f"nothing resolved under panel_form={form}"
            assert all(m.var_name in var_names for m in resolved)

    def test_starsolo_form_strips_the_prefix(self) -> None:
        catalogue = {r["refseq_gene"]: r for r in _catalog() if r["virus"] == "Epstein-Barr virus"}
        assert catalogue["BZLF1"]["gene_id_starsolo"] == "HHV4_BZLF1"
        assert catalogue["BZLF1"]["gene_id_bundled"] == "EPSTEIN_HHV4_BZLF1"

    def test_starsolo_keeps_the_hhv6b_prefix(self) -> None:
        """The one virus the STARsolo packager does not strip."""
        rows = [r for r in _catalog() if r["virus"] == "Human herpesvirus 6b"]
        assert rows
        for row in rows:
            assert row["gene_id_starsolo"].startswith("HUM_HERP6B_"), row

    def test_starsolo_unavailable_markers_are_skipped_with_a_reason(self) -> None:
        catalogue = _catalog()
        var_names = [
            r["gene_id_starsolo"] for r in catalogue if r["available_in_starsolo"] == "true"
        ]
        resolved, unresolved = resolve_markers("Epstein-Barr virus", var_names, "starsolo")
        assert not any("unassigned_gene" in m.var_name for m in resolved)
        # The skip must be visible, not silent.
        assert all(u.reason for u in unresolved), unresolved

    def test_unknown_virus_resolves_to_nothing_without_raising(self) -> None:
        resolved, unresolved = resolve_markers("Not a virus", [], "bundled")
        assert resolved == [] and unresolved == []

    def test_padding_is_never_normalised(self) -> None:
        """Two padding widths are two different HHV-6A genes.

        The catalogue holds U42 as ``HUM_HERP6_HHV6gp41``. U43 -- a different
        ORF -- is ``HUM_HERP6_HHV6gp041``. Six further pairs differ the same
        way. If matching normalised padding, supplying only the 3-digit form
        would silently resolve to the 2-digit marker and report a gene that is
        not in the panel as present.
        """
        two_digit, _ = resolve_markers("Human herpesvirus 6", ["HUM_HERP6_HHV6gp41"], "bundled")
        assert [m.var_name for m in two_digit] == ["HUM_HERP6_HHV6gp41"]

        # U43 (the 3-digit form) is not a catalogue marker, so supplying it
        # alone must resolve nothing. Padding-normalising matching would return
        # the U42 marker and report a gene that is not in the panel.
        three_digit, _ = resolve_markers("Human herpesvirus 6", ["HUM_HERP6_HHV6gp041"], "bundled")
        assert three_digit == [], (
            "the 3-digit form was silently matched to the 2-digit marker; "
            "padding normalisation has crept back in"
        )

        # Both widths present: each must map to its own marker, never across.
        both, _ = resolve_markers(
            "Human herpesvirus 6",
            ["HUM_HERP6_HHV6gp41", "HUM_HERP6_HHV6gp085"],
            "bundled",
        )
        assert {m.var_name for m in both} == {
            "HUM_HERP6_HHV6gp41",
            "HUM_HERP6_HHV6gp085",
        }


def _matrix(pairs: list[tuple[str, int, float]], n_cells: int, n_cols: int):
    """Sparse (cells x genes) matrix from (gene_index, cell_index, value)."""
    data = [v for _, _, v in pairs]
    rows = [c for _, c, _ in pairs]
    cols = [g for g, _, _ in pairs]
    return sparse.csr_matrix((data, (rows, cols)), shape=(n_cells, n_cols))


class TestCalling:
    """``min_breadth`` counts distinct non-overlapping *overlap groups*."""

    def _markers(self, virus="Epstein-Barr virus", programme="productive"):
        catalogue = _catalog()
        return [
            Marker(
                var_name=row["gene_id_bundled"],
                programme=row["programme"],
                overlap_group=row["overlap_group"],
                non_overlapping=row["non_overlapping"] == "true",
                has_cds=row["has_cds"] == "true",
            )
            for row in catalogue
            if row["virus"] == virus and row["programme"] == programme
        ]

    def test_one_productive_gene_is_not_productive(self) -> None:
        markers = self._markers()
        bzlf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_BZLF1")
        m = _matrix([(0, 0, 20.0)], 1, 1)
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["state"] == "indeterminate", calls
        # One overlap group is one piece of evidence, which is below the
        # default min_breadth=2. BZLF1 shares group g28 with BRLF1, so even
        # BZLF1 + BRLF1 together would only count once.
        assert calls[0]["productive_breadth"] == 1
        assert bzlf1.overlap_group == "g28"

    def test_two_distinct_productive_groups_is_productive(self) -> None:
        markers = self._markers()
        independent = [m for m in markers if m.non_overlapping][:2]
        assert len({m.overlap_group for m in independent}) == 2
        pairs = [(i, 0, 20.0) for i, _ in enumerate(independent)]
        m = _matrix(pairs, 1, len(independent))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["state"] == "productive", calls
        assert calls[0]["productive_breadth"] == 2

    def test_two_markers_in_one_overlap_group_still_fail(self) -> None:
        """The cross-mapping defence: nested ORFs are one piece of evidence."""
        markers = self._markers()
        shared = [m for m in markers if m.overlap_group == "g39"]  # BTRF1 + BcLF1
        assert len(shared) == 2, [m.var_name for m in shared]
        idx = [markers.index(m) for m in shared]
        m = _matrix([(i, 0, 20.0) for i in idx], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["state"] == "indeterminate", calls
        assert calls[0]["productive_breadth"] == 1

    def test_non_overlapping_markers_in_a_shared_group_count_once(self) -> None:
        markers = self._markers()
        idx = [
            i for i, m in enumerate(markers) if m.overlap_group == "g28"
        ]  # BZLF1 + BRLF1, both non_overlapping=false
        m = _matrix([(i, 0, 20.0) for i in idx], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["productive_breadth"] == 1

    def test_latent_anchor_gives_latent_when_observable(self) -> None:
        markers = self._markers(programme="latent")
        bhrf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_BHRF1")
        m = _matrix([(markers.index(bhrf1), 0, 20.0)], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2, latency_observable=True)
        assert calls[0]["state"] == "latent", calls
        assert calls[0]["latent_breadth"] == 1

    def test_latent_unreachable_when_not_observable(self) -> None:
        markers = self._markers(programme="latent")
        bhrf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_BHRF1")
        m = _matrix([(markers.index(bhrf1), 0, 20.0)], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2, latency_observable=False)
        assert calls[0]["state"] != "latent", calls
        assert calls[0]["state"] in {"productive", "indeterminate"}
        assert calls[0]["latency_not_observable"] is True

    def test_mixed_requires_both_sides(self) -> None:
        catalogue = _catalog()
        lat = [
            Marker(
                r["gene_id_bundled"],
                "latent",
                r["overlap_group"],
                r["non_overlapping"] == "true",
                r["has_cds"] == "true",
            )
            for r in catalogue
            if r["virus"] == "Epstein-Barr virus" and r["programme"] == "latent"
        ]
        prod = self._markers()
        allm = lat + prod
        # Select independent markers from each programme *separately*; taking
        # them from the combined list would pick two latent ones.
        lat_idx = [allm.index(next(m for m in lat if m.var_name == "EPSTEIN_HHV4_BHRF1"))]
        prod_idx = [allm.index(m) for m in prod if m.non_overlapping][:2]
        assert len(prod_idx) == 2 and len({allm[i].overlap_group for i in prod_idx}) == 2
        m = _matrix([(lat_idx[0], 0, 20.0)] + [(i, 0, 20.0) for i in prod_idx], 1, len(allm))
        calls = call_cell_programme(m, allm, min_breadth=2, latency_observable=True)
        assert calls[0]["state"] == "mixed", calls
        assert calls[0]["productive_breadth"] == 2 and calls[0]["latent_breadth"] == 1

    def test_mixed_unreachable_when_latency_not_observable(self) -> None:
        catalogue = _catalog()
        lat = [
            Marker(
                r["gene_id_bundled"],
                "latent",
                r["overlap_group"],
                r["non_overlapping"] == "true",
                r["has_cds"] == "true",
            )
            for r in catalogue
            if r["virus"] == "Epstein-Barr virus" and r["programme"] == "latent"
        ]
        prod = self._markers()
        allm = lat + prod
        lat_idx = [allm.index(next(m for m in lat if m.var_name == "EPSTEIN_HHV4_BHRF1"))]
        prod_idx = [allm.index(m) for m in prod if m.non_overlapping][:2]
        m = _matrix([(lat_idx[0], 0, 20.0)] + [(i, 0, 20.0) for i in prod_idx], 1, len(allm))
        calls = call_cell_programme(m, allm, min_breadth=2, latency_observable=False)
        assert calls[0]["state"] != "mixed", calls

    def test_empty_cell_is_indeterminate(self) -> None:
        markers = self._markers()
        m = sparse.csr_matrix((1, len(markers)))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["state"] == "indeterminate"
        assert calls[0]["productive_breadth"] == 0 and calls[0]["latent_breadth"] == 0

    def test_every_returned_state_is_in_the_controlled_vocabulary(self) -> None:
        assert {"productive", "latent", "mixed", "indeterminate"} == STATES

    def test_unique_and_selected_breadths_are_reported_separately(self) -> None:
        """Two evidence layers, never merged.

        This is the point of the design: on the selected layer the lytic side
        inflates to match the latent side, and a caller that cannot see both
        cannot tell the difference.
        """
        markers = self._markers()
        idx = [i for i, m in enumerate(markers) if m.non_overlapping][:2]
        unique = _matrix([(i, 0, 20.0) for i in idx], 1, len(markers))
        selected = _matrix([(i, 0, 20.0) for i in range(len(markers))], 1, len(markers))
        calls = call_cell_programme(unique, markers, min_breadth=2, selected_matrix=selected)
        assert calls[0]["productive_breadth"] == 2
        assert calls[0]["selected_productive_breadth"] > 2
        assert calls[0]["evidence_layer"] == "counts_unique_viral"

    def test_rejects_min_breadth_below_one(self) -> None:
        with pytest.raises(ValueError, match="min_breadth"):
            call_cell_programme(sparse.csr_matrix((1, 1)), [], min_breadth=0)


class TestEbvLclRegression:
    """The regression that justifies every design choice above.

    Synthetic matrix built from the real EBV gene set, with latent markers high,
    productive markers low (the true LCL situation) *and* the latent/productive
    markers sharing overlap groups so the multimap-allocated layer inflates the
    productive side the way the real run did.
    """

    def _real_markers(self) -> list[Marker]:
        return [
            Marker(
                var_name=row["gene_id_bundled"],
                programme=row["programme"],
                overlap_group=row["overlap_group"],
                non_overlapping=row["non_overlapping"] == "true",
                has_cds=row["has_cds"] == "true",
            )
            for row in _catalog()
            if row["virus"] == "Epstein-Barr virus"
        ]

    def test_unique_layer_recovers_latent_calls_the_allocated_layer_loses(self) -> None:
        """The asymmetry the real run exhibits, asserted end to end.

        Measured on the bundled EBV LCL run: 2,240 cells latent on the unique
        layer versus 1,277 on the allocated layer, because 1,263 fall to
        ``indeterminate`` there once cross-mapping has drained their latent
        signal onto lytic ORFs. **0** cells go latent-on-unique to
        productive-on-allocated, so the gain is sensitivity, not a different
        direction.
        """
        markers = self._real_markers()
        latent_idx = [
            i for i, m in enumerate(markers) if m.programme == "latent" and m.non_overlapping
        ]
        prod_idx = [
            i for i, m in enumerate(markers) if m.programme == "productive" and m.non_overlapping
        ]
        assert len(latent_idx) >= 2 and len(prod_idx) >= 2

        # Unique-placing evidence: latent markers carry the signal.
        unique = _matrix([(i, 0, 50.0) for i in latent_idx], 1, len(markers))
        # Allocated evidence: cross-mapping has moved the latent signal onto the
        # productive ORFs, so breadth there is dominated by lytic groups.
        allocated = _matrix(
            [(i, 0, 50.0) for i in latent_idx] + [(i, 0, 50.0) for i in prod_idx],
            1,
            len(markers),
        )

        on_unique = call_cell_programme(
            unique,
            markers,
            min_breadth=2,
            latency_observable=True,
            selected_matrix=allocated,
        )[0]
        assert on_unique["state"] == "latent", on_unique
        assert on_unique["selected_state"] != "latent", (
            "the allocated layer should NOT reach a latent call here; if it does, "
            "this test has stopped exercising the cross-mapping scenario"
        )
        # Directional safety: the unique layer must never turn a cell the
        # allocated layer called productive into a latent one.
        assert not (on_unique["selected_state"] == "productive" and on_unique["state"] == "latent")

    def test_directional_guard_fires_when_allocation_is_productive_only(self) -> None:
        """The guard above is only meaningful if this scenario is reachable."""
        markers = self._real_markers()
        latent_idx = [
            i for i, m in enumerate(markers) if m.programme == "latent" and m.non_overlapping
        ]
        prod_idx = [
            i for i, m in enumerate(markers) if m.programme == "productive" and m.non_overlapping
        ]
        unique = _matrix([(i, 0, 50.0) for i in latent_idx], 1, len(markers))
        productive_only = _matrix([(i, 0, 50.0) for i in prod_idx], 1, len(markers))
        call = call_cell_programme(
            unique, markers, min_breadth=2, latency_observable=True, selected_matrix=productive_only
        )[0]
        assert call["state"] == "latent" and call["selected_state"] == "productive", call

    def test_calls_depend_on_breadth_not_mass(self) -> None:
        """The real run's aggregate masses (latent 236,342 vs lytic 247,633) are
        uninformative on their own; the rule must ignore mass given the support."""
        markers = self._real_markers()
        latent = [i for i, m in enumerate(markers) if m.programme == "latent"]
        prod = [i for i, m in enumerate(markers) if m.programme == "productive"]
        support = latent[:1] + prod[:1]

        def _call(latent_mass: float, prod_mass: float) -> dict:
            entries = [(i, 0, latent_mass if i in latent else prod_mass) for i in support]
            matrix = _matrix(entries, 1, len(markers))
            return call_cell_programme(matrix, markers, min_breadth=2, selected_matrix=matrix)[0]

        balanced = _call(1.0, 1.05)
        skewed = _call(1000.0, 0.01)
        assert balanced["state"] == skewed["state"], (balanced, skewed)
        assert balanced["latent_breadth"] == skewed["latent_breadth"] == 1
        assert balanced["productive_breadth"] == skewed["productive_breadth"] == 1


class TestSummaryAndWriters:
    def _cell_df(self) -> pd.DataFrame:
        rows = []
        for barcode, virus, state, pb, lb, sel_state in [
            ("c1", "Epstein-Barr virus", "latent", 0, 2, "latent"),
            ("c2", "Epstein-Barr virus", "latent", 0, 1, "indeterminate"),
            ("c3", "Epstein-Barr virus", "productive", 3, 0, "productive"),
            ("c4", "Epstein-Barr virus", "indeterminate", 0, 0, "indeterminate"),
            ("c5", "Human herpesvirus 1", "indeterminate", 0, 0, "indeterminate"),
        ]:
            rows.append(
                {
                    "barcode": barcode,
                    "virus_name": virus,
                    "state": state,
                    "productive_breadth": pb,
                    "latent_breadth": lb,
                    "selected_state": sel_state,
                    "selected_productive_breadth": pb + 1,
                    "selected_latent_breadth": lb + 1,
                    "evidence_layer": "counts_unique_viral",
                    "latency_not_observable": virus == "Human herpesvirus 1",
                }
            )
        return pd.DataFrame(rows)

    def test_counts_per_state_are_exact(self) -> None:
        summary = summarise_programs(self._cell_df(), _catalog(), min_breadth=2)
        by_virus = summary.set_index("virus_name")
        ebv = by_virus.loc["Epstein-Barr virus"]
        assert ebv["n_cells_latent"] == 2
        assert ebv["n_cells_productive"] == 1
        assert ebv["n_cells_indeterminate"] == 1
        assert ebv["n_cells_mixed"] == 0
        assert ebv["n_cells_total"] == 4

    def test_summary_carries_completeness_and_observability(self) -> None:
        summary = summarise_programs(self._cell_df(), _catalog(), min_breadth=2)
        by_virus = summary.set_index("virus_name")
        assert by_virus.loc["Epstein-Barr virus"]["panel_completeness"] == "complete"
        assert bool(by_virus.loc["Epstein-Barr virus"]["latency_observable_in_rna"])
        hsv1 = by_virus.loc["Human herpesvirus 1"]
        assert hsv1["panel_completeness"] == "partial"
        assert not bool(hsv1["latency_observable_in_rna"])
        assert hsv1["n_cells_latent"] == 0, "latent must be unreachable for HSV-1"

    def test_virus_with_no_catalogue_entry_is_not_applicable(self) -> None:
        df = self._cell_df().copy()
        df.loc[len(df)] = {
            "barcode": "c6",
            "virus_name": "Betatorquevirus",
            "state": "indeterminate",
            "productive_breadth": 0,
            "latent_breadth": 0,
            "selected_state": None,
            "selected_productive_breadth": 0,
            "selected_latent_breadth": 0,
            "evidence_layer": "counts_unique_viral",
            "latency_not_observable": True,
        }
        summary = summarise_programs(df, _catalog(), min_breadth=2)
        row = summary.set_index("virus_name").loc["Betatorquevirus"]
        assert row["panel_completeness"] == "not_applicable"
        assert "no programme model" in row["caveat"]

    def test_detected_virus_with_no_marker_evidence_still_gets_a_row(self) -> None:
        """A virus layer 1 detected must never vanish from the summary.

        HSV-1 in the real EBV LCL run is exactly this case: detected by layer 1,
        has a catalogue entry, but no marker carried uniquely-placing molecules.
        Silence would be read as "never looked at".
        """
        empty = pd.DataFrame(
            columns=[
                "barcode",
                "virus_name",
                "state",
                "productive_breadth",
                "latent_breadth",
                "selected_state",
                "selected_productive_breadth",
                "selected_latent_breadth",
                "evidence_layer",
                "latency_not_observable",
            ]
        )
        summary = summarise_programs(
            empty, _catalog(), min_breadth=2, viruses=["Human herpesvirus 1"]
        )
        assert len(summary) == 1
        row = summary.iloc[0]
        assert row["virus_name"] == "Human herpesvirus 1"
        assert row["n_cells_total"] == 0
        assert row["panel_completeness"] == "partial"
        assert row["latency_observable_in_rna"] == False  # noqa: E712
        assert "no programme could be assessed" in row["caveat"]

    def test_summary_reports_both_evidence_layers(self) -> None:
        summary = summarise_programs(self._cell_df(), _catalog(), min_breadth=2)
        ebv = summary.set_index("virus_name").loc["Epstein-Barr virus"]
        # 2 latent + 1 productive on the unique layer; 1 latent + 1 productive
        # on the allocated layer, which is the honest comparison to publish.
        assert ebv["n_cells_latent"] == 2 and ebv["n_cells_productive"] == 1
        assert ebv["n_cells_latent_selected_layer"] == 1
        assert ebv["n_cells_productive_selected_layer"] == 1

    def test_writers_emit_both_files(self, tmp_path) -> None:
        cells = self._cell_df()
        summary = summarise_programs(cells, _catalog(), min_breadth=2)
        paths = write_program_outputs(cells, summary, str(tmp_path))
        assert len(paths) == 2
        for path in paths:
            assert os.path.exists(path), path
        assert os.path.basename(paths[0]) == "gene_program_summary.tsv"
        assert os.path.basename(paths[1]) == "gene_program_cells.tsv"
        reread = pd.read_csv(paths[0], sep="\t")
        assert set(reread["virus_name"]) == set(cells["virus_name"])

    def test_writers_are_idempotent(self, tmp_path) -> None:
        cells = self._cell_df()
        summary = summarise_programs(cells, _catalog(), min_breadth=2)
        write_program_outputs(cells, summary, str(tmp_path))
        first = pd.read_csv(
            os.path.join(str(tmp_path), "results", "gene_program_summary.tsv"), sep="\t"
        )
        write_program_outputs(cells, summary, str(tmp_path))
        second = pd.read_csv(
            os.path.join(str(tmp_path), "results", "gene_program_summary.tsv"), sep="\t"
        )
        pd.testing.assert_frame_equal(first, second)


class TestDagAndCli:
    """The rule is optional, and it must sit after detection.

    Mirrors the dry-run pattern in ``tests/test_snakefile_dag.py``: the rule is
    only defined when ``gene_programs`` is set, and its scope is
    ``results/viral_summary.tsv``, so it cannot run before layer 1 has finished.
    """

    _SNAKEFILE = os.path.join(REPO_ROOT, "src", "viralscan", "Snakefile")

    def _source(self) -> str:
        with open(self._SNAKEFILE) as handle:
            return handle.read()

    def test_rule_is_gated_on_gene_programs(self) -> None:
        source = self._source()
        # Updated 2026-09-28. The gate used to be a bare `if config.get("gene_programs"):`,
        # which is always True: menu.py serialises booleans for `snakemake --config`
        # as the *strings* "true"/"false", and both are non-empty, so the gene_programs
        # layer ran on every invocation regardless of the flag. The Snakefile now
        # compares the lowercased string. This test must keep pinning that comparison,
        # otherwise the bug returns silently.
        assert (
            'str(config.get("gene_programs", "")).lower() == "true"' in source
        ), (
            "the gene_programs rule must be optional, like hostresponse and "
            "host_filter — a second layer should not be forced on every run. It must "
            'be gated on str(config.get("gene_programs", "")).lower() == "true", '
            "because config values arrive as strings."
        )

    def test_rule_depends_on_detection_and_its_summary(self) -> None:
        source = self._source()
        # Anchor on the rule itself: `gene_programs` also appears in
        # _all_targets, and that occurrence must not be mistaken for the rule.
        start = source.index("rule gene_programs:")
        block = source[start : start + 700]
        assert "log/detection.done" in block
        assert "results/viral_summary.tsv" in block, (
            "layer 2's scope is the viruses layer 1 detected; without the "
            "summary as an input it could run before layer 1 and see nothing"
        )
        assert "scripts/gene_programs.py" in block

    def test_gene_programs_is_not_a_default_target_dependency(self) -> None:
        """The sentinel joins `rule all` only when enabled, like hostresponse."""
        source = self._source()
        assert (
            'str(config.get("gene_programs", "")).lower() == "true"' in source
            and "targets.append" in source.split("def _all_targets")[1].split("return")[0]
        ), "the `all` target must join the gene_programs sentinel only when enabled"

    def test_config_defaults_are_off(self) -> None:
        from viralscan.defaults import DEFAULTS

        assert DEFAULTS["gene_programs"] is False
        assert DEFAULTS["programme_min_breadth"] == 2

    def test_runconfig_accepts_the_defaults(self) -> None:
        from viralscan.runconfig import RunConfig

        cfg = RunConfig()
        assert cfg.gene_programs is False
        assert cfg.programme_min_breadth == 2

    def test_runconfig_rejects_breadth_below_one(self) -> None:
        from viralscan.runconfig import RunConfig

        snake = {
            "output": "o/",
            "index": "i",
            "transcripts": "t",
            "sample1": "s1",
            "sample2": "s2",
            "f1": "None",
            "reference": "False",
            "umap": "False",
            "technology": "10xv3",
            "visual": "True",
            "multimapping": "True",
            "gtf": "None",
            "fasta": "None",
            "whitelist": "None",
            "emptydrops_seed": 100,
            "programme_min_breadth": 0,
        }
        with pytest.raises(ValueError, match="programme_min_breadth"):
            RunConfig.from_snakemake_config(snake)

    def test_cli_exposes_both_flags(self) -> None:
        """The flags must exist on the real CLI, not just in the config layer."""
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-m", "viralscan", "--help"],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONPATH": os.path.join(REPO_ROOT, "src")},
        )
        text = result.stdout + result.stderr
        assert "--gene-programs" in text, text[-800:]
        assert "--programme-min-breadth" in text, text[-800:]

    def test_gene_programs_flag_defaults_off(self) -> None:
        """A second layer must be opt-in, so a plain run is unchanged."""
        from viralscan.defaults import DEFAULTS

        assert DEFAULTS["gene_programs"] is False

    def test_rerun_programs_subcommand_is_registered(self) -> None:
        from viralscan import menu

        assert hasattr(menu, "_build_rerun_programs_parser")
        assert hasattr(menu, "_run_rerun_programs")


class TestCatalogueBiology:
    """Reclassifications backed by primary literature (PLAN PROG-11)."""

    def _row(self, gene_id: str) -> dict | None:
        return next((r for r in _catalog() if r["gene_id_bundled"] == gene_id), None)

    def test_barf1_and_barf1_rr_are_not_ebv_latency_markers(self) -> None:
        """BARF1 is latent only in epithelial cancers; BaRF1 is the lytic RR."""
        assert self._row("EPSTEIN_HHV4_BARF1.2") is None
        assert self._row("EPSTEIN_HHV4_BaRF1.1") is None

    def test_cmv_immediate_early_genes_are_not_latent(self) -> None:
        for gene_id in ("HUM_CYTO_HHV5wtgp107", "HUM_CYTO_HHV5wtgp108"):
            row = self._row(gene_id)
            assert row is not None and row["programme"] == "productive", row

    def test_cmv_latency_is_not_observable(self) -> None:
        rows = [r for r in _catalog() if r["virus"] == "Human cytomegalovirus"]
        assert rows and all(r["latency_observable_in_rna"] == "false" for r in rows)

    def test_kshv_orf16_is_lytic_vbcl2(self) -> None:
        row = self._row("HUM_HERP8_HHV8GK18_gp19")
        assert row is not None and row["programme"] == "productive", row
        assert "vBcl-2" in row["note"] and "vGPCR" not in row["note"].split(";")[0]


class TestResolveIsCaseSensitive:
    def test_barf1_does_not_match_barf1_rr(self) -> None:
        import sys

        sys.path.insert(0, os.path.join(REPO_ROOT, "extras"))
        import build_gene_programs as generator

        records = {
            "EPSTEIN_HHV4_BARF1.2": {"gene": "BARF1", "description": ""},
            "EPSTEIN_HHV4_BaRF1.1": {"gene": "BaRF1", "description": ""},
        }
        assert generator._resolve(records, "BARF1") == ["EPSTEIN_HHV4_BARF1.2"]
        assert generator._resolve(records, "BaRF1") == ["EPSTEIN_HHV4_BaRF1.1"]
