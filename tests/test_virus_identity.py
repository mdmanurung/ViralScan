"""Tests for the per-Run Virus Identity table (PLAN MECH-A step 2)."""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from viralscan import virus_catalog
from viralscan.virus_identity import (
    CATALOGUED,
    COMBINED,
    HOST,
    LEGACY_PREFIX,
    UNCATALOGUED,
    VIRUS_ONLY,
    VirusIdentityTable,
    _Catalogue,
    _default_anello_genus,
    _segment_label,
    assign_virus_keys,
    build_identity_table,
    normalise_strain,
    read_t2g,
)


def _row(accession_version: str, **fields: str) -> dict[str, str]:
    row = {
        "accession": accession_version.split(".")[0],
        "accession_version": accession_version,
        "species": "",
        "genus": "",
        "family": "",
        "segment": "",
        "isolate": "",
        "refseq": "true",
        "taxid": "",
        "organism": "",
        "strain": "",
        "common_name": "",
        "sibling_group": "",
        "risk_class": "",
        "role": "",
    }
    row.update(fields)
    return row


EBV1 = _row(
    "NC_007605.1",
    species="Human gammaherpesvirus 4",
    organism="Human gammaherpesvirus 4",
    family="Orthoherpesviridae",
    taxid="10376",
    common_name="Epstein-Barr virus",
    sibling_group="HHV-4",
)
EBV2 = _row(
    "NC_009334.1",
    species="Human herpesvirus 4 type 2",
    organism="Human herpesvirus 4 type 2",
    family="Orthoherpesviridae",
    taxid="12509",
    sibling_group="HHV-4",
)
HPV45_A = _row(
    "EF202156.1",
    species="Alphapapillomavirus 7",
    organism="Human papillomavirus type 45",
    taxid="10593",
    refseq="false",
)
HPV45_B = _row(
    "NC_075269.1",
    species="Alphapapillomavirus 7",
    organism="human papillomavirus type 45",
    taxid="10593",
)
# One isolate, own taxid, strain spelled two ways across its segments.
GOOSE = [
    _row(
        "NC_007357.1",
        species="Influenza A virus",
        organism="Influenza A virus (A/goose/Guangdong/1/1996(H5N1))",
        taxid="93838",
        segment="1",
        strain="A/Goose/Guangdong/1/96(H5N1)",
    ),
    _row(
        "NC_007358.1",
        species="Influenza A virus",
        organism="Influenza A virus (A/goose/Guangdong/1/1996(H5N1))",
        taxid="93838",
        segment="2",
        strain="A/goose/Guangdong/1/1996",
    ),
]
# Two isolates sharing the species-level taxid 11320.
SPECIES_LEVEL = [
    _row(
        "PQ000001.1",
        species="Influenza A virus",
        organism="Influenza A virus",
        taxid="11320",
        segment="1",
        strain="A/California/147/2024(H5N1)",
        refseq="false",
    ),
    _row(
        "PQ000002.1",
        species="Influenza A virus",
        organism="Influenza A virus",
        taxid="11320",
        segment="2",
        strain="A/California/147/2024",
        refseq="false",
    ),
    _row(
        "PQ000003.1",
        species="Influenza A virus",
        organism="Influenza A virus",
        taxid="11320",
        segment="1",
        strain="A/Washington/239/2024",
        refseq="false",
    ),
]
ANELLO_BETA = _row(
    "PQ438050.1",
    species="Anelloviridae sp.",
    family="Anelloviridae",
    taxid="2055263",
    risk_class="eve",
)
ANELLO_ALPHA = _row(
    "PQ438051.1",
    species="Anelloviridae sp.",
    family="Anelloviridae",
    taxid="2055263",
    risk_class="eve",
)
ANELLO_GENUS = {"PQ438050.1": "Betatorquevirus", "PQ438051.1": "Alphatorquevirus"}

CATALOGUE = [EBV1, EBV2, HPV45_A, HPV45_B, *GOOSE, *SPECIES_LEVEL, ANELLO_BETA, ANELLO_ALPHA]


def _keys(rows):
    return dict(zip((r["accession_version"] for r in rows), assign_virus_keys(rows, ANELLO_GENUS)))


class TestVirusKeys:
    def test_one_taxid_is_one_virus_and_sibling_taxids_stay_apart(self):
        keys = _keys(CATALOGUE)
        assert keys["EF202156.1"] == keys["NC_075269.1"] == "taxid:10593"
        assert keys["NC_007605.1"] == "taxid:10376"
        assert keys["NC_009334.1"] == "taxid:12509"

    def test_isolate_taxid_groups_segments_whatever_the_strain_spelling(self):
        keys = _keys(CATALOGUE)
        assert keys["NC_007357.1"] == keys["NC_007358.1"] == "taxid:93838"

    def test_species_level_taxid_splits_into_isolates(self):
        keys = _keys(CATALOGUE)
        assert keys["PQ000001.1"] == keys["PQ000002.1"] == "taxid:11320|a/california/147/2024"
        assert keys["PQ000003.1"] == "taxid:11320|a/washington/239/2024"

    def test_anelloviruses_key_by_genus_not_by_catch_all_taxid(self):
        keys = _keys(CATALOGUE)
        assert keys["PQ438050.1"] == "genus:Betatorquevirus"
        assert keys["PQ438051.1"] == "genus:Alphatorquevirus"

    def test_anellovirus_outside_the_table_falls_back_to_catalogue_genus_then_family(self):
        rows = [
            _row("MZ000001.1", family="Anelloviridae", genus="Gammatorquevirus", taxid="1"),
            _row("MZ000002.1", family="Anelloviridae", taxid="1"),
        ]
        assert assign_virus_keys(rows, {}) == ["genus:Gammatorquevirus", "genus:Anelloviridae"]

    def test_a_lineage_token_that_is_not_a_genus_falls_back_to_the_family(self):
        rows = [_row("MZ000003.1", family="Anelloviridae", genus="Small anellovirus", taxid="1")]
        assert assign_virus_keys(rows, {}) == ["genus:Anelloviridae"]

    def test_row_without_taxid_keys_by_accession(self):
        assert assign_virus_keys([_row("AB000001.1")], {}) == ["accession:AB000001"]

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("A/Puerto Rico/8/1934(H1N1)", "a/puerto rico/8/1934"),
            ("A/Puerto Rico/8/1934", "a/puerto rico/8/1934"),
            ("  A/Perth/16/2009 (H3N2) ", "a/perth/16/2009"),
            ("B/Lee/1940", "b/lee/1940"),
        ],
    )
    def test_normalise_strain(self, text, expected):
        assert normalise_strain(text) == expected

    def test_segment_labels_compare_across_spellings(self):
        assert _segment_label({"segment": "L RNA"}) == _segment_label({"segment": "L"}) == "l"
        assert _segment_label({"segment": "segment 4"}) == "4"


class TestDisplayNames:
    def test_curated_name_wins_then_organism_then_isolate(self):
        cat = _Catalogue(CATALOGUE, ANELLO_GENUS)
        assert cat.names["taxid:10376"] == "Epstein-Barr virus"
        assert cat.names["taxid:12509"] == "Human herpesvirus 4 type 2"
        # The RefSeq row names a key with two records.
        assert cat.names["taxid:10593"] == "human papillomavirus type 45"
        assert cat.names["taxid:11320|a/washington/239/2024"] == (
            "Influenza A virus (A/Washington/239/2024)"
        )
        assert cat.names["genus:Betatorquevirus"] == "Betatorquevirus"

    def test_anelloviruses_without_a_genus_are_named_as_unassigned(self):
        rows = [_row("MZ000002.1", family="Anelloviridae", taxid="1")]
        assert _Catalogue(rows, {}).names == {
            "genus:Anelloviridae": "Anelloviridae (genus unassigned)"
        }

    def test_a_shared_name_is_disambiguated_by_key(self):
        rows = [
            _row("AB000001.1", organism="Same virus", taxid="1"),
            _row("AB000002.1", organism="Same virus", taxid="2"),
        ]
        names = _Catalogue(rows, {}).names
        assert names == {"taxid:1": "Same virus [taxid:1]", "taxid:2": "Same virus [taxid:2]"}


def _write_t2g(path: Path, lines: list[str]) -> Path:
    path.write_text("".join(line + "\n" for line in lines))
    return path


def _combined_t2g(tmp_path: Path) -> Path:
    return _write_t2g(
        tmp_path / "t2g.txt",
        [
            "ENST1.1\tENSG1.1\t\t\tENST1.1\t1\t400\t+",
            "ENST2.1\tENSG2.1\t\t\tENST1.1\t1\t300\t+",
            "EBV_tx1\tNC_007605.1_gene1\tBZLF1\t\tNC_007605.1\t10\t900\t+",
            "EBV_tx2\tNC_007605.1_gene2\tBRLF1\t\tNC_007605.1\t900\t1900\t+",
            "EBV2_tx1\tNC_009334_gene1\t\t\tNC_009334\t10\t900\t+",
            "MY_tx1\tMY_VIRUS_g1\t\t\tMY000001.1\t1\t500\t+",
            "ZZ_tx1\tMYSTERY_g1\t\t\tZZ999999.1\t1\t500\t+",
        ],
    )


class TestReadT2g:
    def test_strict_tab_split_keeps_column_five_on_host_rows(self, tmp_path):
        t2g = read_t2g(_combined_t2g(tmp_path))
        assert not t2g.legacy
        assert t2g.accession["ENSG1.1"] == "ENST1.1"
        assert t2g.accession["NC_007605.1_gene1"] == "NC_007605.1"
        assert t2g.structural_host == {"ENSG1.1", "ENSG2.1"}

    def test_a_short_row_marks_the_t2g_legacy(self, tmp_path):
        t2g = read_t2g(_write_t2g(tmp_path / "t2g.txt", ["tx1 GENE_A GENE_A", "tx2 HOST HOST"]))
        assert t2g.legacy
        assert t2g.accession == {"GENE_A": "", "HOST": ""}


class TestBuildIdentityTable:
    def _table(self, tmp_path, gtf_genes=("MY_VIRUS_g1",)):
        return build_identity_table(
            _combined_t2g(tmp_path),
            gtf_genes,
            catalogue_rows=CATALOGUE,
            anello_genus=ANELLO_GENUS,
        )

    def test_every_status(self, tmp_path):
        genes = self._table(tmp_path).by_gene()
        assert genes["NC_007605.1_gene1"].status == CATALOGUED
        assert genes["NC_009334_gene1"].status == CATALOGUED  # bare accession
        assert genes["MY_VIRUS_g1"].status == UNCATALOGUED
        assert genes["MYSTERY_g1"].status == HOST
        assert genes["ENSG1.1"].status == HOST
        assert [g for g, row in genes.items() if row.viral] == [
            "NC_007605.1_gene1",
            "NC_007605.1_gene2",
            "NC_009334_gene1",
            "MY_VIRUS_g1",
        ]

    def test_catalogued_gene_carries_its_catalogue_decisions(self, tmp_path):
        ebv = self._table(tmp_path).by_gene()["NC_007605.1_gene1"]
        assert ebv.virus_key == "taxid:10376"
        assert ebv.virus_name == "Epstein-Barr virus"
        assert ebv.taxid == "10376"
        assert ebv.sibling_group == "HHV-4"
        assert ebv.family == "Orthoherpesviridae"

    def test_ebv_types_are_two_viruses_in_one_sibling_group(self, tmp_path):
        table = self._table(tmp_path)
        assert table.groups() == {
            "taxid:10376": ["NC_007605.1_gene1", "NC_007605.1_gene2"],
            "taxid:12509": ["NC_009334_gene1"],
            "accession:MY000001.1": ["MY_VIRUS_g1"],
        }
        genes = table.by_gene()
        assert genes["NC_009334_gene1"].sibling_group == genes["NC_007605.1_gene1"].sibling_group

    def test_uncatalogued_gene_is_named_by_accession_and_warned(self, tmp_path, caplog):
        with caplog.at_level(logging.WARNING, logger="viralscan.virus_identity"):
            gene = self._table(tmp_path).by_gene()["MY_VIRUS_g1"]
        assert (gene.virus_key, gene.virus_name) == ("accession:MY000001.1", "MY000001.1")
        assert "1 viral gene(s) on 1 genome(s) are not in the virus catalogue" in caplog.text

    def test_host_gene_in_a_combined_gtf_stays_host(self, tmp_path):
        genes = self._table(tmp_path, gtf_genes=("MY_VIRUS_g1", "ENSG1.1", "ENSG2.1")).by_gene()
        assert genes["ENSG1.1"].status == HOST
        assert genes["ENSG2.1"].status == HOST

    def test_index_kind(self, tmp_path):
        assert self._table(tmp_path).index_kind == COMBINED
        virus_only = _write_t2g(
            tmp_path / "virus_t2g.txt",
            ["EBV_tx1\tNC_007605.1_gene1\t\t\tNC_007605.1\t10\t900\t+"],
        )
        table = build_identity_table(virus_only, (), catalogue_rows=CATALOGUE, anello_genus={})
        assert table.index_kind == VIRUS_ONLY

    def test_no_viral_gene_is_an_error(self, tmp_path):
        host_only = _write_t2g(tmp_path / "t2g.txt", ["ENST1.1\tENSG1.1\t\t\tENST1.1\t1\t400\t+"])
        with pytest.raises(ValueError, match="No viral gene") as exc:
            build_identity_table(host_only, ["ENSG1.1"], catalogue_rows=CATALOGUE, anello_genus={})
        assert "1 GTF gene(s) were treated as host" in str(exc.value)

    def test_legacy_t2g_uses_the_gtf_gene_set_and_prefix_names(self, tmp_path):
        legacy = _write_t2g(
            tmp_path / "t2g.txt",
            ["tx1\tEPSTEIN_HHV4_BZLF1\tBZLF1", "tx2\tHOST_GENE\tHOST_GENE"],
        )
        genes = build_identity_table(legacy, ["EPSTEIN_HHV4_BZLF1"]).by_gene()
        assert genes["EPSTEIN_HHV4_BZLF1"].status == LEGACY_PREFIX
        assert genes["EPSTEIN_HHV4_BZLF1"].virus_name == "Epstein-Barr virus"
        assert genes["HOST_GENE"].status == HOST

    def test_tsv_round_trip(self, tmp_path):
        table = self._table(tmp_path)
        path = table.write_tsv(tmp_path / "results" / "virus_identity.tsv")
        assert VirusIdentityTable.read_tsv(path) == table

    def test_read_tsv_rejects_other_tables(self, tmp_path):
        other = tmp_path / "other.tsv"
        other.write_text("gene_id\tviral\nX\ttrue\n")
        with pytest.raises(ValueError, match="not a Virus Identity table"):
            VirusIdentityTable.read_tsv(other)


@pytest.fixture(scope="module")
def catalogue():
    rows = virus_catalog.load_catalogue()
    return rows, _Catalogue(rows, _default_anello_genus())


class TestPackagedCatalogue:
    def test_no_virus_holds_the_same_segment_twice(self, catalogue):
        rows, cat = catalogue
        segments: dict[str, Counter[str]] = defaultdict(Counter)
        for row in rows:
            _, key = cat.lookup(row["accession_version"])
            if _segment_label(row) and not key.startswith("genus:"):
                segments[key][_segment_label(row)] += 1
        assert {k: c for k, c in segments.items() if max(c.values()) > 1} == {}

    def test_pr8_is_one_virus_of_eight_segments(self, catalogue):
        rows, cat = catalogue
        keys = {cat.lookup(f"NC_0020{n}.1")[1] for n in range(16, 24)}
        assert keys == {"taxid:211044"}

    def test_no_display_name_needs_a_key_suffix(self, catalogue):
        _, cat = catalogue
        assert [name for name in cat.names.values() if " [" in name] == []

    def test_legacy_display_names_are_kept(self, catalogue):
        _, cat = catalogue
        assert cat.names[cat.lookup("NC_007605.1")[1]] == "Epstein-Barr virus"
        assert cat.names[cat.lookup("NC_001806")[1]] == "Human herpesvirus 1"


# Golden tests over the real indexes. They live outside the repo (≈470k rows
# each) and are skipped where they are absent.
_WORK = Path("/exports/para-lipg-hpc/mdmanurung/ViralScan")
_GOLDEN = {
    "v2": _WORK / "ebv_latest_ref_2026-09-27/full_panel/t2g_v2.txt",
    "final": Path(
        "/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_final/build/panel.t2g"
    ),
    "max": _WORK / "viral_panel_max_2026-09-28/t2g_max.txt",
}


@pytest.fixture(scope="module")
def golden_tables():
    tables = {
        name: build_identity_table(path, ()) for name, path in _GOLDEN.items() if path.exists()
    }
    if not tables:
        pytest.skip("real t2g files not available on this machine")
    return tables


class TestGoldenIndexes:
    def test_every_viral_gene_is_catalogued_and_every_host_gene_is_ensembl(self, golden_tables):
        for name, table in golden_tables.items():
            assert set(table.status_counts()) == {CATALOGUED, HOST}, name
            assert table.index_kind == COMBINED
            for gene in table.genes:
                assert gene.gene_id.startswith("ENSG") == (gene.status == HOST), (name, gene)

    def test_ebv_types_are_separate_viruses(self, golden_tables):
        for name in ("final", "max"):
            if name not in golden_tables:
                continue
            groups = golden_tables[name].groups()
            assert "taxid:10376" in groups and "taxid:12509" in groups, name
            genes = golden_tables[name].by_gene()
            assert {genes[g].sibling_group for g in groups["taxid:12509"]} == {"HHV-4"}

    def test_max_panel_groupings(self, golden_tables):
        if "max" not in golden_tables:
            pytest.skip("max panel t2g not available")
        table = golden_tables["max"]
        accessions: dict[str, set[str]] = defaultdict(set)
        segments: dict[str, Counter[str]] = defaultdict(Counter)
        for gene in table.genes:
            if gene.viral:
                accessions[gene.virus_key].add(gene.genome_accession)
        for gene in {g.genome_accession: g for g in table.genes if g.viral and g.segment}.values():
            segments[gene.virus_key][_segment_label({"segment": gene.segment})] += 1
        assert accessions["taxid:10593"] == {"EF202156.1", "NC_075269.1"}  # HPV45
        assert len(accessions["taxid:2697049"]) == 18  # SARS-CoV-2
        assert sum(segments["taxid:211044"].values()) == 8  # PR8
        assert {k: c for k, c in segments.items() if max(c.values()) > 1} == {}
        assert len(set(table.virus_names().values())) == len(table.groups())
