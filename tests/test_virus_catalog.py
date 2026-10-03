"""Tests for src/viralscan/virus_catalog.py (PLAN CAT-03).

Two failures this catalogue exists to prevent, both of which show up as wrong
rows in ``viral_summary.tsv``:

* a genome-scoped gene ID resolving to no virus, so a bare accession is
  reported as though it were a virus name;
* a segmented virus fragmenting, so one influenza strain is reported as eight
  separate "viruses" and N strains as 8 x N.
"""

from __future__ import annotations

import csv

import pytest

from viralscan import virus_catalog
from viralscan.virus_grouping import group_genes_by_virus, virus_name_for_gene

CATALOGUE_PRESENT = virus_catalog.catalogue_path().exists()
requires_catalogue = pytest.mark.skipif(
    not CATALOGUE_PRESENT, reason="virus_catalog.tsv is not packaged in this checkout"
)

_ROWS = [
    # Two segments of one influenza strain: different accessions, one species.
    {
        "accession": "NC_002016",
        "accession_version": "NC_002016.1",
        "species": "Influenza A virus",
        "segment": "7",
    },
    {
        "accession": "NC_002017",
        "accession_version": "NC_002017.1",
        "species": "Influenza A virus",
        "segment": "4",
    },
    # A different species, to prove grouping does not over-merge.
    {
        "accession": "NC_045512",
        "accession_version": "NC_045512.2",
        "species": "Severe acute respiratory syndrome coronavirus 2",
        "segment": "",
    },
    # A row with no species name must be skipped, not mapped to "".
    {"accession": "XX_000000", "accession_version": "XX_000000.1", "species": "", "segment": ""},
]


@pytest.fixture
def catalogue(tmp_path, monkeypatch):
    path = tmp_path / "virus_catalog.tsv"
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "accession_version", "species", "segment"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(_ROWS)
    monkeypatch.setattr(virus_catalog, "_default_tsv_path", lambda: path)
    for cached in (
        virus_catalog.catalog_name_map,
        virus_catalog.accession_metadata,
        virus_catalog.merged_name_map,
    ):
        cached.cache_clear()
    yield path
    for cached in (
        virus_catalog.catalog_name_map,
        virus_catalog.accession_metadata,
        virus_catalog.merged_name_map,
    ):
        cached.cache_clear()


class TestNameMap:
    def test_maps_versioned_and_bare_accessions(self, catalogue) -> None:
        name_map = virus_catalog.catalog_name_map()
        assert name_map["NC_002016.1"] == "Influenza A virus"
        assert name_map["NC_002016"] == "Influenza A virus"

    def test_row_without_a_species_is_skipped(self, catalogue) -> None:
        """An empty species must not become an empty virus name."""
        name_map = virus_catalog.catalog_name_map()
        assert "XX_000000.1" not in name_map
        assert "" not in name_map.values()

    def test_missing_catalogue_degrades_to_empty(self, tmp_path, monkeypatch) -> None:
        """A checkout without the catalogue keeps working on the legacy map."""
        monkeypatch.setattr(virus_catalog, "_default_tsv_path", lambda: tmp_path / "absent.tsv")
        virus_catalog.catalog_name_map.cache_clear()
        virus_catalog.merged_name_map.cache_clear()
        assert virus_catalog.load_catalogue() == []
        assert virus_catalog.catalog_name_map() == {}
        assert virus_catalog.merged_name_map(), "legacy entries must still be present"
        virus_catalog.catalog_name_map.cache_clear()
        virus_catalog.merged_name_map.cache_clear()


class TestSegmentGrouping:
    """The reason this module exists."""

    def test_segments_of_one_strain_group_into_one_virus(self, catalogue) -> None:
        genes = ["NC_002016.1_gene1", "NC_002017.1_gene1"]
        groups, detected = group_genes_by_virus(genes, virus_catalog.merged_name_map())
        assert list(groups) == ["Influenza A virus"]
        assert len(groups["Influenza A virus"]) == 2
        assert detected == {"Influenza A virus"}

    def test_different_species_do_not_merge(self, catalogue) -> None:
        genes = ["NC_002016.1_gene1", "NC_045512.2_gene1"]
        groups, _ = group_genes_by_virus(genes, virus_catalog.merged_name_map())
        assert len(groups) == 2

    def test_real_gene_tokens_resolve_too(self, catalogue) -> None:
        """Not just the `_gene1` placeholder: CAT-01 emits real gene tokens."""
        name_map = virus_catalog.merged_name_map()
        assert virus_name_for_gene("NC_045512.2_orf1ab", name_map) == (
            "Severe acute respiratory syndrome coronavirus 2"
        )

    def test_segments_for_species_orders_by_segment_number(self, catalogue) -> None:
        assert virus_catalog.segments_for_species("Influenza A virus") == [
            "NC_002017.1",  # segment 4
            "NC_002016.1",  # segment 7
        ]


class TestPrecedence:
    def test_anellovirus_genus_wins_over_the_catalogue_species(self, catalogue) -> None:
        """ANDET-05 wants one genus label per anellovirus genome, not a species."""
        name_map = virus_catalog.merged_name_map()
        assert name_map.get("NC_002076.2") == "Alphatorquevirus"


@requires_catalogue
class TestPackagedCatalogue:
    """Assertions about the catalogue actually shipped in this checkout."""

    def test_every_row_has_an_accession_and_species(self) -> None:
        rows = virus_catalog.load_catalogue()
        assert rows, "the packaged catalogue must not be empty"
        assert all(r["accession_version"] for r in rows)
        assert all(r["species"] for r in rows), "a row without a species names nothing"

    def test_accession_versions_are_unique(self) -> None:
        versions = [r["accession_version"] for r in virus_catalog.load_catalogue()]
        assert len(versions) == len(set(versions))

    def test_influenza_a_is_one_virus_across_all_eight_segments(self) -> None:
        """One isolate's eight segments group into one virus.

        Since the MECH-A merge the catalogue holds several Influenza A isolates,
        so the eight segments are those of the RefSeq A/Puerto Rico/8/1934 set
        (NC_002016-NC_002023), selected by their isolate-level taxid. The
        strain text is not a usable key: two PR8 segments spell it
        "A/Puerto Rico/8/1934(H1N1)".
        """
        rows = virus_catalog.load_catalogue()
        pr8 = next(r for r in rows if r["accession_version"] == "NC_002016.1")
        segments = [r["accession_version"] for r in rows if r["taxid"] == pr8["taxid"]]
        assert len(segments) == 8, segments
        assert {r["segment"] for r in rows if r["taxid"] == pr8["taxid"]} == {
            str(n) for n in range(1, 9)
        }
        groups, _ = group_genes_by_virus(
            [f"{a}_gene1" for a in segments], virus_catalog.merged_name_map()
        )
        assert list(groups) == ["Influenza A virus"]

    def test_panel_is_a_known_scope(self) -> None:
        panels = {r["panel"] for r in virus_catalog.load_catalogue()}
        assert panels <= {"shipped", "max", "broad", "legacy"}, panels

    def test_broad_discovery_list_matches_the_catalogue_exactly(self) -> None:
        """CAT-09: the frozen broad-discovery list is the whole catalogue, one row each."""
        path = virus_catalog.catalogue_path().with_name("broad_discovery_accessions.tsv")
        with open(path, newline="") as handle:
            listed = list(csv.DictReader(handle, delimiter="\t"))
        got = [(r["accession_version"], r["panel"]) for r in listed]
        assert len(got) == len({a for a, _ in got}), "duplicate accession in the broad list"
        want = [(r["accession_version"], r["panel"]) for r in virus_catalog.load_catalogue()]
        assert sorted(got) == sorted(want)
