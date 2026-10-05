"""Tests for the named HPV ORF catalogue.

The catalogue exists because the reference ViralScan ships today cannot answer
the question it is being asked. Its HPV gene IDs are RefSeq ``locus_tag`` values
(``HpV16gp1`` ... ``HpV16gp8``) across four accessions, so nothing in the
kallisto ``t2g`` says which ORF is E6 and which is L1, and an
oncogene-versus-capsid contrast in oropharyngeal tissue is not expressible.

What these tests defend
-----------------------
* **The names are the record's, not a table's.** Every row states which
  qualifier the name came from (``/gene`` or ``/product``). A coordinate-derived
  naming scheme would have nothing to state, and would be wrong: papillomavirus
  genomes are linearised circles cut at the submitter's choice of point, so the
  same E6 ORF sits at 7125-7601 in HPV16 and 105-581 in HPV18.
* **The identity of the oncogenes is checkable.** E6 and E7 must be present for
  every high-risk genotype, must be the expected size relative to L1, and must
  not be satisfiable by an isoform such as E6* that does not encode the
  oncoprotein.
* **Coordinates are usable.** Within the genome, consistent in strand, with
  ``cds_length_nt`` and ``n_exons`` agreeing with the exon blocks.
* **The catalogue cannot reintroduce the defect it fixes.** ``gene_id`` is unique
  across the merged reference; the bare symbol ``E6`` occurring 16 times is
  exactly what ``ncbi_fetch._genbank_to_gtf`` produces today.

Nothing here needs the network. ``TestAgainstNcbiCache`` checks the catalogue
against the on-disk NCBI cache and skips per accession; ``TestCatalogueHasNotDrifted``
re-fetches and is gated behind ``@pytest.mark.network``.
"""

from __future__ import annotations

import csv
import importlib.util
import os
import re
from pathlib import Path

import pytest

from viralscan import hpv_genes
from viralscan.hpv_genes import (
    CAPSID_GENES,
    HIGH_RISK_GENOTYPES,
    REQUIRED_COLUMNS,
    HpvOrf,
    capsid_gene_ids,
    early_gene_ids,
    gene_ids_for,
    genotypes,
    load_catalogue,
    oncogene_gene_ids,
    oncogene_locus_gene_ids,
    orfs,
    split_oncogene_vs_capsid,
    validate_catalogue,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "extras" / "build_hpv_reference.py"
SPEC = importlib.util.spec_from_file_location("build_hpv_reference", SCRIPT_PATH)
assert SPEC is not None
build_hpv_reference = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(build_hpv_reference)

CANONICAL = build_hpv_reference.CANONICAL_GENES
EXPECTED_GENOTYPES = {
    "1",
    "2",
    "16",
    "18",
    "31",
    "33",
    "35",
    "39",
    "45",
    "51",
    "52",
    "56",
    "58",
    "59",
    "66",
    "68",
}

#: External ground truth: the IARC Group 1 / clinical 14-type high-risk HPV
#: list. Kept as a literal so the test fails if the catalogue's own constant
#: drifts from the clinical list (it previously carried 69 in place of 68).
IARC_HIGH_RISK_14 = {
    "16",
    "18",
    "31",
    "33",
    "35",
    "39",
    "45",
    "51",
    "52",
    "56",
    "58",
    "59",
    "66",
    "68",
}


def _rows() -> list[dict]:
    return load_catalogue()


def _by_genotype() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for row in _rows():
        out.setdefault(row["genotype"], []).append(row)
    return out


def _names(genotype: str) -> set[str]:
    return {row["canonical_gene_name"] for row in _by_genotype().get(genotype, [])}


class TestCatalogueIntegrity:
    def test_loads_and_has_required_columns(self) -> None:
        rows = _rows()
        assert len(rows) > 100, f"catalogue is implausibly small: {len(rows)} rows"
        missing = set(REQUIRED_COLUMNS) - set(rows[0])
        assert not missing, missing

    def test_validate_catalogue_accepts_the_shipped_catalogue(self) -> None:
        assert validate_catalogue(_rows()) == []

    def test_every_canonical_name_is_in_the_build_s_allow_list(self) -> None:
        """A row whose identity is unknown would be counted as if it meant
        something, so the build has no generic-ORF fallback and neither may the
        shipped table."""
        for row in _rows():
            assert row["canonical_gene_name"] in CANONICAL, row

    def test_controlled_vocabularies(self) -> None:
        for row in _rows():
            assert row["gene_class"] in hpv_genes.GENE_CLASSES, row
            assert row["name_source"] in {"gene", "product"}, row
            assert row["source"] in {"refseq", "insdc"}, row
            assert row["topology"] in {"circular", "linear"}, row
            assert row["strand"] in {"+", "-"}, row
            for field in ("high_risk", "spans_origin"):
                assert row[field] in {"true", "false"}, (field, row)

    def test_every_accession_has_at_least_one_gene(self) -> None:
        for genotype, rows in _by_genotype().items():
            assert rows, genotype

    def test_covers_the_expected_genotypes(self) -> None:
        assert set(_by_genotype()) == EXPECTED_GENOTYPES

    def test_all_fourteen_high_risk_genotypes_are_present(self) -> None:
        present = set(genotypes(_rows(), high_risk_only=True))
        assert present == set(HIGH_RISK_GENOTYPES), (
            f"missing high-risk genotypes: {sorted(set(HIGH_RISK_GENOTYPES) - present)}"
        )
        assert len(HIGH_RISK_GENOTYPES) == 14

    def test_high_risk_set_matches_the_external_iarc_list(self) -> None:
        # Guard against the catalogue's own constant drifting from the clinical
        # 14-type list: both the loader and the build script must agree with
        # the hardcoded IARC Group 1 set above (68, not 69).
        assert set(HIGH_RISK_GENOTYPES) == IARC_HIGH_RISK_14
        assert set(build_hpv_reference.HIGH_RISK_GENOTYPES) == IARC_HIGH_RISK_14

    def test_high_risk_column_matches_the_declared_set(self) -> None:
        for row in _rows():
            expected = row["genotype"] in HIGH_RISK_GENOTYPES
            assert (row["high_risk"] == "true") is expected, row

    def test_no_duplicate_accession_and_gene_pairs(self) -> None:
        keys = [(row["accession_version"], row["canonical_gene_name"]) for row in _rows()]
        duplicates = [k for k in set(keys) if keys.count(k) > 1]
        assert not duplicates, f"duplicate (accession, gene) rows: {duplicates}"

    def test_gene_id_is_unique_across_the_merged_reference(self) -> None:
        """The defect being fixed. ``ncbi_fetch._genbank_to_gtf`` takes
        ``gene_id`` from ``/gene=``, so merging two HPV genomes yields 16 ORFs on
        9 gene IDs with every genome's E6 collapsed onto one row."""
        ids = [row["gene_id"] for row in _rows()]
        assert len(ids) == len(set(ids)), "gene_id is not unique"
        assert len(ids) > len(set(row["accession_version"] for row in _rows())) * 5

    def test_gene_id_is_namespaced_by_accession(self) -> None:
        for row in _rows():
            assert row["gene_id"].startswith(row["accession_version"] + "_"), row

    def test_gene_id_carries_no_bare_orf_symbol(self) -> None:
        """``*`` and ``^`` are spelled out: a kallisto ``t2g`` is an unquoted
        two-column TSV and bustools has to agree with kallisto on every byte."""
        for row in _rows():
            prefix = row["accession_version"] + "_"
            symbol = row["gene_id"][len(prefix) :]
            assert re.fullmatch(r"[A-Za-z0-9_]+", symbol), row

    def test_gene_ids_are_not_the_opaque_locus_tags(self) -> None:
        """The current index's ``HpV16gp1``-style IDs are what this catalogue
        exists to replace; none may reappear."""
        for row in _rows():
            assert not re.search(r"gp\d+$", row["gene_id"]), row


class TestCoordinates:
    def test_coordinates_lie_within_the_genome(self) -> None:
        for row in _rows():
            start = int(row["genome_start"])
            end = int(row["genome_end"])
            length = int(row["genome_length"])
            assert 1 <= start <= end <= length, row

    def test_cds_length_matches_the_exon_blocks(self) -> None:
        for row in _rows():
            blocks = [tuple(int(v) for v in b.split("..")) for b in row["exon_blocks"].split(",")]
            assert int(row["n_exons"]) == len(blocks), row
            assert int(row["cds_length_nt"]) == sum(e - s + 1 for s, e in blocks), row
            assert int(row["genome_start"]) == min(s for s, _ in blocks), row
            assert int(row["genome_end"]) == max(e for _, e in blocks), row

    def test_one_strand_per_genome(self) -> None:
        """Papillomavirus transcribes its early and late ORFs from a single
        origin, so a genome whose CDS features disagree on strand has been
        misparsed."""
        for version, rows in _group_by_version().items():
            strands = {row["strand"] for row in rows}
            assert len(strands) == 1, (version, strands)

    def test_shipped_catalogue_is_all_on_the_plus_strand(self) -> None:
        """Recorded as a fact about this build. A future reverse-complemented
        submission would legitimately fail this, and the failure is the point:
        read it before quoting coordinates."""
        assert {row["strand"] for row in _rows()} == {"+"}

    def test_genome_length_is_consistent_within_a_genome(self) -> None:
        for version, rows in _group_by_version().items():
            assert len({row["genome_length"] for row in rows}) == 1, version
            assert len({row["topology"] for row in rows}) == 1, version

    def test_papillomavirus_genomes_are_about_eight_kilobases(self) -> None:
        for row in _rows():
            assert 7500 <= int(row["genome_length"]) <= 8200, row

    def test_origin_spanning_orfs_are_flagged_and_look_the_way_one_expects(self) -> None:
        """HPV16's ``E1^E4`` is ``join(1..16, 2494..2756)`` in a circular record,
        so its bounding box is 1-2756 and would read as a 2756-nt gene without the
        flag."""
        spanning = [row for row in _rows() if row["spans_origin"] == "true"]
        assert spanning, "no origin-spanning ORF recorded; NC_001526.4 E1^E4 should be"
        for row in spanning:
            assert row["topology"] == "circular", row
            assert int(row["n_exons"]) > 1, row
            blocks = row["exon_blocks"].split(",")
            assert blocks[0].startswith("1.."), row
            span = int(row["genome_end"]) - int(row["genome_start"]) + 1
            assert span > int(row["cds_length_nt"]), (
                f"{row['gene_id']}: bounding box {span} nt is far larger than the "
                f"CDS {row['cds_length_nt']} nt, so the flag is load-bearing"
            )


def _group_by_version() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for row in _rows():
        out.setdefault(row["accession_version"], []).append(row)
    return out


class TestOncogeneVersusCapsid:
    """The contrast the catalogue exists to make expressible."""

    def test_oncogene_and_capsid_vocabularies_are_exactly_e6_e7_and_l1_l2(self) -> None:
        """Pinned rather than derived, because a widening of either set would
        silently change what a positive call means."""
        oncogenes = {"E6", "E6*", "E7", "E7*"}
        capsids = {"L1", "L2"}
        isoforms = {"E6*", "E7*"}
        assert set(hpv_genes.ONCOGENE_GENES) == oncogenes
        assert set(hpv_genes.CAPSID_GENES) == capsids
        assert set(hpv_genes.ONCOGENE_LOCUS_GENES) == isoforms
        assert not oncogenes & capsids
        assert isoforms < oncogenes

    def test_e6_and_e7_exist_for_every_high_risk_genotype(self) -> None:
        for genotype in sorted(HIGH_RISK_GENOTYPES):
            names = _names(genotype)
            assert "E6" in names, f"HPV-{genotype} has no E6"
            assert "E7" in names, f"HPV-{genotype} has no E7"

    def test_capsid_genes_exist_for_every_genotype(self) -> None:
        """Without L1 and L2 the contrast is not computable for that genotype."""
        for genotype in sorted(EXPECTED_GENOTYPES):
            names = _names(genotype)
            missing = CAPSID_GENES - names
            assert names >= CAPSID_GENES, f"HPV-{genotype} missing {missing}"

    def test_oncogene_and_capsid_gene_ids_are_disjoint(self) -> None:
        onco = set(oncogene_gene_ids(_rows()))
        capsid = set(capsid_gene_ids(_rows()))
        assert onco and capsid
        assert not onco & capsid

    def test_two_oncogenes_per_genotype(self) -> None:
        assert len(oncogene_gene_ids(_rows())) == 2 * len(EXPECTED_GENOTYPES)

    def test_oncogene_isoforms_are_not_reported_as_the_oncogene(self) -> None:
        """E6* lacks the PDZ-binding motif and E7* lacks the LXCXE pRb-binding
        motif, so neither encodes an oncoprotein. Their reads are
        indistinguishable from E6/E7, so letting them into the oncogene class
        would let E6* alone produce a confident positive call."""
        for row in _rows():
            if row["canonical_gene_name"] in {"E6*", "E7*"}:
                assert row["gene_class"] == "oncogene_locus", row
            if row["gene_class"] == "oncogene":
                assert row["canonical_gene_name"] in {"E6", "E7"}, row

    def test_isoforms_reach_the_locus_helper_but_not_the_oncogene_helper(self) -> None:
        onco = set(oncogene_gene_ids(_rows()))
        locus = set(oncogene_locus_gene_ids(_rows()))
        assert locus > onco
        assert not locus & set(capsid_gene_ids(_rows()))

    def test_l1_is_far_larger_than_either_oncogene(self) -> None:
        """An E6/E7/L1 mix-up would break this. Note the invariant that *looks*
        obvious and is false: L1 is not the largest ORF, because E1 is."""
        for genotype, rows in _by_genotype().items():
            spans = {row["canonical_gene_name"]: int(row["cds_length_nt"]) for row in rows}
            for oncogene in ("E6", "E7"):
                if oncogene in spans:
                    assert spans["L1"] >= 3 * spans[oncogene], (genotype, spans)

    def test_e1_is_longer_than_l1_which_is_why_l1_is_not_the_longest_orf(self) -> None:
        """Pinned because getting it backwards fails every record. E1 is the
        replication helicase at ~650 aa; L1 runs ~505-570 aa."""
        for genotype, rows in _by_genotype().items():
            if "E1" not in {row["canonical_gene_name"] for row in rows}:
                continue
            spans = {row["canonical_gene_name"]: int(row["cds_length_nt"]) for row in rows}
            assert spans["E1"] > spans["L1"], (genotype, spans)

    def test_every_genotype_has_a_usable_oncogene_to_capsid_contrast(self) -> None:
        for genotype in sorted(HIGH_RISK_GENOTYPES):
            onco = oncogene_gene_ids(_rows(), [genotype])
            capsid = capsid_gene_ids(_rows(), [genotype])
            assert onco, genotype
            assert capsid, genotype

    def test_capsid_genes_are_larger_than_the_oncogenes_they_are_contrasted_against(self) -> None:
        rows = _rows()
        mean_onco = sum(
            int(r["cds_length_nt"]) for r in rows if r["gene_class"] == "oncogene"
        ) / max(1, len(oncogene_gene_ids(rows)))
        mean_capsid = sum(
            int(r["cds_length_nt"]) for r in rows if r["gene_class"] == "late_capsid"
        ) / max(1, len(capsid_gene_ids(rows)))
        assert mean_capsid > 2 * mean_onco, (mean_onco, mean_capsid)


class TestNameProvenance:
    def test_every_row_states_which_qualifier_the_name_came_from(self) -> None:
        for row in _rows():
            assert row["name_source"] in {"gene", "product"}, row
            assert row["product_as_in_ncbi"] or row["name_source"] == "gene", row

    def test_gene_sourced_names_are_the_verbatim_gene_qualifier(self) -> None:
        for row in _rows():
            if row["name_source"] == "gene":
                assert row["product_as_in_ncbi"] == row["canonical_gene_name"], row

    def test_product_sourced_names_appear_in_the_product_qualifier(self) -> None:
        """HPV45's record carries no ``/gene`` at all, so its names come from
        ``/product``. Token matching must not be substring matching: ``E6`` must
        not be read out of ``E6*``."""
        for row in _rows():
            if row["name_source"] != "product":
                continue
            symbol = row["canonical_gene_name"]
            assert re.search(
                rf"(?:^|[^A-Za-z0-9]){re.escape(symbol)}(?:[^A-Za-z0-9]|$)",
                row["product_as_in_ncbi"],
            ), row

    def test_refseq_and_insdc_are_both_represented(self) -> None:
        sources = {row["source"] for row in _rows()}
        assert sources == {"refseq", "insdc"}

    def test_refseq_is_used_where_it_exists(self) -> None:
        """RefSeq has complete genomes for only 4 of the 14 high-risk types, so
        the other 10 must be INSDC. This is the documented gap, not a choice."""
        refseq = {row["genotype"] for row in _rows() if row["source"] == "refseq"}
        assert refseq == {"1", "2", "16", "18", "31", "33"}, refseq


class TestAgainstNcbiCache:
    """Coordinates against the sequences ViralScan already has on disk.

    Skipped per accession when the FASTA is not cached, so the suite stays fast
    and offline. The cached FASTAs are the ones ``ncbi_fetch`` wrote with a
    SHA-256 sidecar, i.e. the sequences an actual reference build would use.
    """

    def _cached_fasta_length(self, accession: str) -> int | None:
        from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR

        path = Path(DEFAULT_CACHE_DIR) / accession / f"{accession}.fasta"
        if not path.exists() or not path.stat().st_size:
            return None
        return sum(
            len(line.strip()) for line in path.read_text().splitlines() if not line.startswith(">")
        )

    def test_coordinates_lie_within_the_cached_fasta(self) -> None:
        checked = 0
        for version, rows in _group_by_version().items():
            length = self._cached_fasta_length(version)
            if length is None:
                continue
            checked += 1
            assert int(rows[0]["genome_length"]) == length, (
                f"{version}: catalogue says {rows[0]['genome_length']} nt, cached "
                f"FASTA is {length} nt"
            )
            for row in rows:
                assert int(row["genome_end"]) <= length, (version, row)
        if checked == 0:
            pytest.skip(
                "no HPV FASTA in the ncbi_fetch cache; run "
                "extras/build_hpv_reference.py --emit-reference to populate it"
            )


class TestParser:
    """Unit tests for the GenBank reader, on synthetic records.

    These pin the two parsing decisions that a coordinate table would have
    sidestepped and that are easy to get wrong: a record's topology, and a CDS
    that wraps the linearisation point of a circular record.
    """

    HEADER = (
        "LOCUS       TESTHPV                  1000 bp    DNA     circular SYN 01-JAN-2020\n"
        "DEFINITION  synthetic.\n"
        "VERSION     TESTHPV.1\n"
    )

    def _record(self, features: str) -> str:
        return self.HEADER + "FEATURES             Location/Qualifiers\n" + features + ("ORIGIN\n")

    def _cds(self, location: str, qualifiers: str) -> str:
        return f"     CDS             {location}\n" + "".join(
            f"                     {q}\n" for q in qualifiers.split("|")
        )

    def test_reads_locus_length_and_topology(self) -> None:
        header = build_hpv_reference.parse_genbank(self._record(self._cds("1..300", '/gene="E6"')))
        assert header["locus"] == "TESTHPV"
        assert header["accession"] == "TESTHPV"
        assert header["accession_version"] == "TESTHPV.1"
        assert header["length"] == 1000
        assert header["topology"] == "circular"

    def test_linear_topology_is_distinguished_from_circular(self) -> None:
        text = self._record(self._cds("1..300", '/gene="E6"')).replace(
            "DNA     circular", "DNA     linear  "
        )
        assert build_hpv_reference.parse_genbank(text)["topology"] == "linear"

    def test_a_record_with_no_topology_is_an_error(self) -> None:
        text = self._record(self._cds("1..300", '/gene="E6"')).replace(
            "DNA     circular", "DNA          "
        )
        with pytest.raises(build_hpv_reference.BuildError, match="topology"):
            build_hpv_reference.parse_genbank(text)

    def test_spliced_cds_keeps_its_blocks(self) -> None:
        header = build_hpv_reference.parse_genbank(
            self._record(self._cds("join(1..16,2494..2756)", '/gene="E1^E4"'))
        )
        cds = header["cds"][0]
        assert cds["exon_blocks"] == [(1, 16), (2494, 2756)]
        assert cds["strand"] == "+"

    def test_complement_location_is_read_as_minus_strand(self) -> None:
        header = build_hpv_reference.parse_genbank(
            self._record(self._cds("complement(100..300)", '/gene="E7"'))
        )
        assert header["cds"][0]["strand"] == "-"
        assert header["cds"][0]["exon_blocks"] == [(100, 300)]

    def test_wrapped_qualifier_values_are_not_truncated(self) -> None:
        """A truncated ``/translation`` would silently drop the motif checks, and
        the first 60 characters of an E6 or E7 product carry no motif at all."""
        first, rest = "MVIIIIIAAAA", "BBBCCCCDDDDEEEEFFFFGGGGHHHHIIIIJJJJKKKK"
        header = build_hpv_reference.parse_genbank(
            self._record(self._cds("1..300", f'/gene="E6"|/translation="{first}{rest}"'))
        )
        assert header["cds"][0]["qualifiers"]["translation"] == first + rest

    def test_wrapped_location_is_rejoined(self) -> None:
        text = self._record(
            "     CDS             join(1..16,\n"
            "                     2494..2756)\n"
            '                     /gene="E1^E4"\n'
        )
        cds = build_hpv_reference.parse_genbank(text)["cds"][0]
        assert cds["exon_blocks"] == [(1, 16), (2494, 2756)]

    def test_multiple_cds_do_not_share_qualifiers(self) -> None:
        header = build_hpv_reference.parse_genbank(
            self._record(self._cds("1..300", '/gene="E6"') + self._cds("400..600", '/gene="E7"'))
        )
        names = [c["qualifiers"]["gene"] for c in header["cds"]]
        assert names == ["E6", "E7"]
        blocks = [c["exon_blocks"] for c in header["cds"]]
        assert blocks == [[(1, 300)], [(400, 600)]]

    def test_a_record_with_no_cds_is_an_error(self) -> None:
        with pytest.raises(build_hpv_reference.BuildError, match="no CDS"):
            build_hpv_reference.parse_genbank(self._record(""))


class TestNameResolution:
    def test_gene_qualifier_wins_when_it_is_the_only_source(self) -> None:
        assert build_hpv_reference.resolve_gene_name(
            {"gene": "E6", "product": "transforming protein E6"}
        ) == ("E6", "gene", "E6")

    def test_product_is_used_when_gene_is_absent(self) -> None:
        name, source, verbatim = build_hpv_reference.resolve_gene_name({"product": "E6 protein"})
        assert (name, source, verbatim) == ("E6", "product", "E6 protein")

    def test_product_token_matching_ignores_surrounding_words(self) -> None:
        for product in (
            "putative transforming protein E6",
            "early protein E6",
            "major capsid protein L1",
            "late protein L1",
            "putative E4 protein",
            "cell cycle modulating protein E1^E4",
        ):
            name, source, _ = build_hpv_reference.resolve_gene_name({"product": product})
            assert source == "product", product
            assert name in CANONICAL, product

    def test_e6_star_is_not_misread_as_e6(self) -> None:
        """Substring matching would turn ``protein E6*`` into ``E6`` and report a
        non-transforming isoform as the oncoprotein."""
        name, _, _ = build_hpv_reference.resolve_gene_name({"product": "protein E6*"})
        assert name == "E6*"

    def test_an_isoform_in_product_overrides_its_parent_in_gene(self) -> None:
        """``NC_001526.4`` annotates the truncated E6* CDS with ``/gene="E6"`` and
        ``/product="protein E6*"``. Taking /gene at face value would report the
        isoform as the oncoprotein."""
        name, source, _ = build_hpv_reference.resolve_gene_name(
            {"gene": "E6", "product": "protein E6*"}
        )
        assert (name, source) == ("E6*", "product")

    def test_e1_caret_e4_is_not_misread_as_e1(self) -> None:
        name, _, _ = build_hpv_reference.resolve_gene_name(
            {"gene": "E1^E4", "product": "cell cycle modulating protein E1^E4"}
        )
        assert name == "E1^E4"

    def test_gene_and_product_disagreeing_unrelatedly_is_an_error(self) -> None:
        with pytest.raises(build_hpv_reference.BuildError, match="disagrees"):
            build_hpv_reference.resolve_gene_name({"gene": "E6", "product": "L1 protein"})

    def test_a_product_naming_two_symbols_is_an_error(self) -> None:
        with pytest.raises(build_hpv_reference.BuildError, match="ambiguous"):
            build_hpv_reference.resolve_gene_name({"product": "E6 protein and L1 protein"})

    def test_a_generic_product_with_no_symbol_is_an_error(self) -> None:
        with pytest.raises(build_hpv_reference.BuildError, match="no recognisable"):
            build_hpv_reference.resolve_gene_name(
                {"product": "hypothetical protein", "protein_id": "YP_1.1"}
            )

    def test_an_unannotated_cds_is_an_error_not_a_guess(self) -> None:
        with pytest.raises(build_hpv_reference.BuildError, match="no recognisable"):
            build_hpv_reference.resolve_gene_name({"protein_id": "BAA12345.1"})


class TestSafeGeneIds:
    def test_exotic_symbols_are_spelled_out(self) -> None:
        assert build_hpv_reference.safe_gene_id("NC_001526.4", "E6") == "NC_001526.4_E6"
        assert build_hpv_reference.safe_gene_id("NC_001526.4", "E6*") == "NC_001526.4_E6_star"
        assert build_hpv_reference.safe_gene_id("NC_001526.4", "E1^E4") == "NC_001526.4_E1_E4"

    def test_the_same_orf_in_two_genomes_gets_two_ids(self) -> None:
        """The whole point. A bare ``E6`` from two genomes is one kallisto gene."""
        a = build_hpv_reference.safe_gene_id("NC_001526.4", "E6")
        b = build_hpv_reference.safe_gene_id("NC_001357.1", "E6")
        assert a != b

    def test_no_symbol_is_silently_mangled(self) -> None:
        for symbol in CANONICAL:
            gene_id = build_hpv_reference.safe_gene_id("ACC.1", symbol)
            assert re.fullmatch(r"ACC\.1_[A-Za-z0-9_]+", gene_id), symbol


class TestOrfOrderIsRotationTolerant:
    """The linearisation hazard, as a test rather than a comment."""

    def test_a_record_cut_inside_e1_still_passes(self) -> None:
        """``NC_001526.4`` reports its ORFs as E1, E2, E5, L2, L1, E6, E7 while
        ``NC_001357.1`` reports E6, E7, E1, ... Both are the canonical order."""
        rotated = ["E1", "E2", "E5", "L2", "L1", "E6", "E7"]
        assert build_hpv_reference._check_orf_order("16", "NC_001526.4", rotated) == []

    def test_the_canonical_order_passes(self) -> None:
        canonical = ["E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"]
        assert build_hpv_reference._check_orf_order("18", "NC_001357.1", canonical) == []

    def test_a_reversed_order_is_rejected(self) -> None:
        reversed_order = ["E7", "E6", "E1", "E2", "E4", "E5", "L2", "L1"]
        problems = build_hpv_reference._check_orf_order("16", "X.1", reversed_order)
        assert problems and "not a rotation" in problems[0]

    def test_a_l1_l2_swap_is_rejected(self) -> None:
        problems = build_hpv_reference._check_orf_order(
            "16", "X.1", ["E6", "E7", "E1", "E2", "E4", "E5", "L1", "L2"]
        )
        assert problems

    def test_the_shipped_catalogue_is_consistent_with_the_builds_order_check(self) -> None:
        for version, rows in _group_by_version().items():
            names = [
                row["canonical_gene_name"]
                for row in sorted(rows, key=lambda r: int(r["genome_start"]))
            ]
            assert (
                build_hpv_reference._check_orf_order(rows[0]["genotype"], version, names) == []
            ), version


class TestOpaquePanelIdsCannotBeRescued:
    """The gap this catalogue does *not* close, pinned so it is not forgotten."""

    def test_an_opaque_gene_id_does_not_match_a_catalogue_gene_id(self) -> None:
        catalogue_ids = {row["gene_id"] for row in _rows()}
        for opaque in ("HpV16gp1", "HpV16gp2", "HpV1agp1", "Hpv1gp01"):
            assert opaque not in catalogue_ids

    def test_the_catalogue_names_the_orfs_the_index_only_numbered(self) -> None:
        """``HpV16gp1`` is E6 and ``HpV16gp2`` is E7 in the shipped index. The
        catalogue can say so only because the *record* says so, which is why the
        index has to be rebuilt rather than relabelled after the fact."""
        rows = [r for r in _rows() if r["genotype"] == "16"]
        by_name = {r["canonical_gene_name"]: r for r in rows}
        assert by_name["E6"]["genome_start"] == "7125"
        assert by_name["E7"]["genome_start"] == "7604"
        assert by_name["L1"]["canonical_gene_name"] == "L1"
        assert by_name["E6"]["gene_class"] == "oncogene"
        assert by_name["L1"]["gene_class"] == "late_capsid"


class TestModuleApi:
    def test_typed_orfs_round_trip(self) -> None:
        records = orfs(_rows())
        assert len(records) == len(_rows())
        assert all(isinstance(r, HpvOrf) for r in records)
        sample = next(r for r in records if r.canonical_gene_name == "E6")
        assert sample.is_oncogene and not sample.is_capsid
        assert sample.is_from_oncogene_locus

    def test_capsid_orfs_report_as_capsid(self) -> None:
        record = next(r for r in orfs(_rows()) if r.canonical_gene_name == "L1")
        assert record.is_capsid and not record.is_oncogene
        assert not record.is_from_oncogene_locus

    def test_gene_ids_for_filters_by_genotype(self) -> None:
        rows = _rows()
        everything = gene_ids_for(rows, "late_capsid")
        only_16 = gene_ids_for(rows, "late_capsid", ["16"])
        assert set(only_16) < set(everything)
        assert len(only_16) == 2

    def test_early_gene_ids_exclude_the_oncogene_locus(self) -> None:
        early = set(early_gene_ids(_rows()))
        assert early
        assert not early & set(oncogene_locus_gene_ids(_rows()))
        assert not early & set(capsid_gene_ids(_rows()))

    def test_genotypes_are_sorted_numerically(self) -> None:
        assert genotypes(_rows())[0] == "1"
        assert genotypes(_rows())[-1] == "68"

    def test_split_reports_both_sides(self) -> None:
        rows = _rows()
        counts = {row["gene_id"]: 10.0 for row in rows if row["gene_class"] == "oncogene"}
        result = split_oncogene_vs_capsid(counts, rows, ["16"])
        assert result["oncogene_locus_total"] == 20.0
        assert result["capsid_total"] == 0.0
        assert result["oncogene_to_capsid_ratio"] is None
        assert result["genotypes"] == ["16"]

    def test_split_reports_a_ratio_only_when_both_sides_are_non_zero(self) -> None:
        rows = _rows()
        counts = {
            row["gene_id"]: 5.0 for row in rows if row["gene_class"] in {"oncogene", "late_capsid"}
        }
        result = split_oncogene_vs_capsid(counts, rows)
        assert result["oncogene_locus_total"] > 0
        assert result["capsid_total"] > 0
        assert result["oncogene_to_capsid_ratio"] == pytest.approx(1.0)

    def test_capsid_only_is_not_reported_as_a_zero_ratio(self) -> None:
        """Capsid-but-no-oncogene is the ordinary state of a productive infection,
        not a ratio of zero."""
        rows = _rows()
        counts = {row["gene_id"]: 3.0 for row in rows if row["gene_class"] == "late_capsid"}
        result = split_oncogene_vs_capsid(counts, rows)
        assert result["oncogene_locus_total"] == 0.0
        assert result["capsid_total"] > 0
        assert result["oncogene_to_capsid_ratio"] is None

    def test_unknown_genotype_yields_nothing_rather_than_everything(self) -> None:
        rows = _rows()
        assert split_oncogene_vs_capsid({}, rows, ["999"])["oncogene_to_capsid_ratio"] is None
        assert split_oncogene_vs_capsid({}, rows, ["999"])["genotypes"] == ["999"]

    def test_validate_reports_a_broken_row_rather_than_raising(self) -> None:
        rows = [dict(_rows()[0])]
        rows[0]["gene_class"] = "made_up"
        rows[0]["high_risk"] = "yes"
        rows[0]["genome_end"] = "99999999"
        problems = validate_catalogue(rows)
        assert any("gene_class" in p for p in problems)
        assert any("high_risk" in p for p in problems)
        assert any("exceeds genome_length" in p for p in problems)

    def test_validate_reports_a_duplicate_pair(self) -> None:
        row = dict(_rows()[0])
        problems = validate_catalogue([row, dict(row)])
        assert any("duplicate" in p for p in problems)

    def test_validate_reports_an_empty_catalogue(self) -> None:
        assert validate_catalogue([]) == ["HPV gene catalogue is empty"]


class TestCatalogueHasNotDrifted:
    """Re-derives one record from NCBI and compares it to the shipped row.

    Gated behind ``@pytest.mark.network`` and skipped without ``NCBI_EMAIL``, so
    it never runs in the default suite.
    """

    @pytest.mark.network
    def test_hpv16_e6_row_still_matches_ncbi(self) -> None:
        email = os.environ.get("NCBI_EMAIL")
        if not email:
            pytest.skip("NCBI_EMAIL is not set")
        from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR

        text = build_hpv_reference.fetch_genbank(
            "NC_001526.4", email, os.environ.get("NCBI_API_KEY"), Path(DEFAULT_CACHE_DIR)
        )
        header = build_hpv_reference.parse_genbank(text)
        shipped = next(
            r
            for r in _rows()
            if r["accession_version"] == "NC_001526.4" and r["canonical_gene_name"] == "E6"
        )
        cds = next(
            c
            for c in header["cds"]
            if build_hpv_reference.resolve_gene_name(c["qualifiers"])[0] == "E6"
        )
        blocks = cds["exon_blocks"]
        assert str(header["length"]) == shipped["genome_length"]
        assert str(min(s for s, _ in blocks)) == shipped["genome_start"]
        assert str(max(e for _, e in blocks)) == shipped["genome_end"]

    @pytest.mark.network
    def test_every_selected_accession_is_still_resolvable_at_ncbi(self) -> None:
        email = os.environ.get("NCBI_EMAIL")
        if not email:
            pytest.skip("NCBI_EMAIL is not set")
        from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR

        shipped = {row["accession_version"] for row in _rows()}
        selected = {accession for _, accession, _, _ in build_hpv_reference.SELECTED_RECORDS}
        assert shipped == selected, shipped.symmetric_difference(selected)
        for _, accession, _, _ in build_hpv_reference.SELECTED_RECORDS:
            text = build_hpv_reference.fetch_genbank(
                accession,
                email,
                os.environ.get("NCBI_API_KEY"),
                Path(DEFAULT_CACHE_DIR),
            )
            names = {
                build_hpv_reference.resolve_gene_name(c["qualifiers"])[0]
                for c in build_hpv_reference.parse_genbank(text)["cds"]
            }
            assert {"E6", "E7", "L1", "L2"} <= names, (accession, sorted(names))


def test_tsv_is_tab_separated_with_a_single_header() -> None:
    path = REPO_ROOT / "src" / "viralscan" / "data" / "hpv_genes.tsv"
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    assert header == list(build_hpv_reference.TSV_COLUMNS)
    assert rows
    for row in rows:
        assert len(row) == len(header), row
