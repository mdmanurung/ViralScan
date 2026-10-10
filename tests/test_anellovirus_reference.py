"""Anelloviridae panel gene-structure regression tests.

The bug these lock down
-----------------------
Every one of the panel's 2,042 genomes used to reach a built reference as a
single placeholder gene ``{accession}_gene1`` spanning the whole genome, because
``scripts/build_bundled_panel_ref.py`` and
``build_anellovirus_reference`` both discarded the GenBank-derived GTF that
``ncbi_fetch._fetch_one`` had already written and regenerated a placeholder from
the FASTA.  NCBI does carry real CDS features.

A whole-genome transcript shares sequence with every other genome in the panel,
so reads cross-map in proportion to conservation and the panel's entire viral
signal is absorbed by whichever genome is most conserved.  Measured on a COVID
scRNA-seq run: 1,169,272 anellovirus UMI in total, 1,167,103 of them (99.8 %)
in ``MW455439.1_gene1``, reported as "Alphatorquevirus" against 1,241 UMI in
Betatorquevirus.  That ratio is a conservation rank, not a measurement.

Everything in the default (non-``network``) selection here reads the packaged
``anellovirus_genes.tsv`` and committed GenBank fixtures, so the suite is fast
and needs no internet.  ``@pytest.mark.network`` tests hit NCBI live.
"""

from __future__ import annotations

import csv
import re
import textwrap
from pathlib import Path

import pytest

from viralscan.anellovirus import (
    annotated_accessions,
    candidate_gene_ids,
    gtf_text_for,
    load_accession_table,
    load_gene_table,
    merged_name_map,
)
from viralscan.scripts.ncbi_fetch import (
    NCBIFetchError,
    _genbank_to_gtf,
    fetch_genbank,
)
from viralscan.virus_grouping import virus_name_for_gene

PANEL_TSV = Path("src/viralscan/data/anellovirus_accessions.tsv")
GENE_TSV = Path("src/viralscan/data/anellovirus_genes.tsv")
GENUSES = (
    "Alphatorquevirus",
    "Anelloviridae",
    "Betatorquevirus",
    "Gammatorquevirus",
    "Gyrovirus",
    "Hetorquevirus",
    "Memtorquevirus",
    "Samektorquevirus",
)

# A representative accession per genus, drawn from the panel.
REPRESENTATIVE = {
    "Alphatorquevirus": "NC_002076.2",
    "Anelloviridae": "AB303555.1",
    "Betatorquevirus": "AY823988.1",
    "Gammatorquevirus": "AB303552.1",
    "Gyrovirus": "KF294862.1",
    "Hetorquevirus": "MN769649.1",
    "Memtorquevirus": "MN774938.1",
    "Samektorquevirus": "PP728781.1",
}

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
ACCESSION_PREFIX_RE = re.compile(r"^[A-Z]{1,3}_?\d+\.\d+_")

#: Accessions whose CDS count is asserted to be exact.  Kept as a dict so a
#: re-annotation that changes the count shows up as a named diff rather than an
#: opaque total.
EXPECTED_CDS = {
    "NC_002076.2": 3,
    "AB303552.1": 4,
    "KF294862.1": 3,
    "PP728781.1": 2,
}


def _panel_rows() -> list[dict[str, str]]:
    with open(PANEL_TSV, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _gene_rows() -> list[dict[str, str]]:
    with open(GENE_TSV, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _exon_blocks(exons: str) -> list[tuple[int, int]]:
    return [tuple(int(v) for v in piece.split(":")) for piece in exons.split(",") if piece]


# ── Panel integrity ──────────────────────────────────────────────────────────


class TestPanelLoads:
    def test_panel_has_2041_accessions(self) -> None:
        # 2042 before 2026-09-28, when CAT-05 dropped AB303562.1: it is
        # byte-identical to NC_038359.1 (both 3,187 bp, taxid 2065053, md5
        # 5e78d2d32ea911d61962cd8ad6255a4f) and a duplicate sequence has broken
        # `kallisto index` before. The RefSeq copy is retained because it carries
        # 4 curated gene features the GenBank record lacks.
        assert len(_panel_rows()) == 2041

    def test_no_duplicate_accessions(self) -> None:
        accessions = [row["accession"] for row in _panel_rows()]
        assert len(set(accessions)) == len(accessions)

    def test_every_row_declares_one_of_the_eight_genera(self) -> None:
        genera = {row["genus"] for row in _panel_rows()}
        assert genera == set(GENUSES)

    def test_load_accession_table_matches_the_file(self) -> None:
        assert len(load_accession_table()) == 2041

    def test_documented_panel_counts_match_the_packaged_tables(self) -> None:
        """Counts quoted in docs and help text (2,041 / 2,021 / 2,511 / 1,994) must not drift."""
        import csv
        import importlib.resources
        from pathlib import Path

        data = Path(str(importlib.resources.files("viralscan.data")))
        panel = {row["accession"] for row in _panel_rows()}
        bundled = set()
        for gtf in data.glob("*.gtf"):
            for line in gtf.read_text().splitlines():
                seqname = line.split("\t", 1)[0]
                if seqname in panel:
                    bundled.add(seqname)
        # Genomes with no bundled GTF: the ones `analysis.py` must add from the table.
        assert len(panel) - len(bundled) == 2021
        with (data / "anellovirus_genes.tsv").open() as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        assert len(genes) == 2511
        assert len({row["accession"] for row in genes}) == 1994
        # The 4 stale gene rows of the CAT-05 duplicate AB303562.1 are gone (2,515 / 1,995
        # before 2026-10-10); every gene-table accession is a panel accession.
        assert {row["accession"] for row in genes} - panel == set()
        assert "AB303562.1" not in {row["accession"] for row in genes}


# ── The central regression: real genes, not placeholders ──────────────────────


class TestRealGeneStructure:
    def test_catalogue_ships_and_is_non_trivial(self) -> None:
        rows = _gene_rows()
        assert len(rows) > 2000
        assert len(annotated_accessions()) > 1900

    def test_no_placeholder_gene_ids_survive(self) -> None:
        """The regression itself: ``{accession}_gene1`` must not be a gene ID."""
        placeholders = [row["gene_id"] for row in _gene_rows() if row["gene_id"].endswith("_gene1")]
        assert placeholders == []

    def test_every_gene_id_is_accession_scoped(self) -> None:
        for row in _gene_rows():
            assert row["gene_id"].startswith(f"{row['accession']}_"), row
            assert ACCESSION_PREFIX_RE.match(row["gene_id"]), row

    def test_gene_ids_are_unique_across_the_whole_panel(self) -> None:
        ids = [row["gene_id"] for row in _gene_rows()]
        assert len(set(ids)) == len(ids)

    def test_bare_ncbi_symbol_is_never_a_gene_id(self) -> None:
        """``ORF1`` is 150 different genomes' /product; it must not be a gene ID."""
        for row in _gene_rows():
            assert row["gene_id"] != row["gene"].lower()
            assert "/" not in row["gene_id"]

    def test_representative_subset_yields_real_genes_in_every_genus(self) -> None:
        rows = {row["accession"]: row for row in _gene_rows()}
        for genus, accession in REPRESENTATIVE.items():
            assert accession in rows, f"{genus} representative {accession} not annotated"
            row = rows[accession]
            assert row["gene_id"] != f"{accession}_gene1"
            assert int(row["n_exons"]) >= 1
            assert _exon_blocks(row["exons"]), row

    def test_spliced_genes_have_more_than_one_exon(self) -> None:
        spliced = [row for row in _gene_rows() if int(row["n_exons"]) > 1]
        assert spliced, "no multi-exon gene in the panel — the CDS path did not run"
        for row in spliced:
            blocks = _exon_blocks(row["exons"])
            assert len(blocks) == int(row["n_exons"]) > 1, row
            assert all(end > start for start, end in blocks), row

    def test_torque_teno_virus_has_its_three_known_cds(self) -> None:
        """NC_002076.2 is the reference TTV record: orf1 (Rep), orf2/4 and orf2/5.

        Its spliced VP2 and VP3 are the canonical demonstration that the CDS path
        is running: the whole-genome placeholder can never express a two-exon
        gene.
        """
        rows = [row for row in _gene_rows() if row["accession"] == "NC_002076.2"]
        assert len(rows) == EXPECTED_CDS["NC_002076.2"]
        by_symbol = {row["gene_symbol"]: row for row in rows}
        assert set(by_symbol) == {"orf1", "orf2/4", "orf2/5"}
        assert by_symbol["orf2/4"]["n_exons"] == "2"
        assert by_symbol["orf2/5"]["n_exons"] == "2"
        assert by_symbol["orf1"]["n_exons"] == "1"

    def test_exact_cds_counts_for_known_records(self) -> None:
        counts: dict[str, int] = {}
        for row in _gene_rows():
            counts[row["accession"]] = counts.get(row["accession"], 0) + 1
        for accession, expected in EXPECTED_CDS.items():
            assert counts.get(accession) == expected, accession

    def test_panel_wide_gene_count_is_plausible(self) -> None:
        """Anelloviridae are small: most genomes are 1–5 ORFs, none has dozens."""
        counts: dict[str, int] = {}
        for row in _gene_rows():
            counts[row["accession"]] = counts.get(row["accession"], 0) + 1
        assert max(counts.values()) <= 10, max(counts.items(), key=lambda kv: kv[1])
        assert sum(1 for n in counts.values() if n == 1) / len(counts) > 0.6


# ── ORF1 / Rep presence ──────────────────────────────────────────────────────


class TestRepPresence:
    @staticmethod
    def _is_rep(row: dict[str, str]) -> bool:
        return (
            row["gene_symbol"].strip().lower().startswith("orf1")
            or row["product"].strip().lower().startswith("orf1")
            or "rep" in row["product"].lower()
        )

    def test_orf1_or_rep_is_present_where_ncbi_annotates_it(self) -> None:
        rows = _gene_rows()
        marked = [row for row in rows if self._is_rep(row)]
        assert len(marked) > 1500, "ORF1/Rep should be the dominant annotation"
        assert len(marked) / len(rows) > 0.6

    def test_orf1_is_the_dominant_single_product(self) -> None:
        counts: dict[str, int] = {}
        for row in _gene_rows():
            counts[row["product"]] = counts.get(row["product"], 0) + 1
        assert counts.get("ORF1", 0) > 1500

    def test_orf1_is_a_single_exon_in_the_reference_ttv(self) -> None:
        row = next(
            row
            for row in _gene_rows()
            if row["accession"] == "NC_002076.2" and row["gene_symbol"] == "orf1"
        )
        assert row["n_exons"] == "1"
        assert row["exons"] == "589:2901"
        assert row["gene_id"] == "NC_002076.2_TTVgp3"
        assert row["product"] == "VP1"

    def test_orf1_present_in_every_genus_that_carries_one(self) -> None:
        """Anelloviridae (unclassified) and Gyrovirus annotate differently.

        The two genera absent from the hit list below are the two whose records
        use neither a Rep/ORF1 product nor an ORF1 symbol, so requiring all eight
        would be asserting an NCBI annotation that does not exist.
        """
        genus_of = {row["accession"]: row["genus"] for row in _panel_rows()}
        seen: set[str] = set()
        for row in _gene_rows():
            if self._is_rep(row):
                seen.add(genus_of.get(row["accession"], "?"))
        assert "Alphatorquevirus" in seen
        assert "Betatorquevirus" in seen
        assert "Gammatorquevirus" in seen
        assert "Hetorquevirus" in seen
        assert "Samektorquevirus" in seen
        assert "Memtorquevirus" in seen


# ── Coordinate and topology validity ─────────────────────────────────────────


class TestCoordinates:
    def test_every_exon_lies_within_the_genome(self) -> None:
        for row in _gene_rows():
            length = int(row["genome_length"])
            assert length > 0, row
            for start, end in _exon_blocks(row["exons"]):
                assert 1 <= start <= end <= length, row

    def test_strand_is_plus_or_minus(self) -> None:
        for row in _gene_rows():
            assert row["strand"] in ("+", "-"), row

    def test_all_exons_of_a_gene_share_one_strand(self) -> None:
        gtf = gtf_text_for(sorted(annotated_accessions())[:200])
        by_gene: dict[str, set[str]] = {}
        for line in gtf.strip().splitlines():
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            gene_id = GENE_ID_RE.search(cols[8])
            assert gene_id
            by_gene.setdefault(gene_id.group(1), set()).add(cols[6])
        for gene_id, strands in by_gene.items():
            assert len(strands) == 1, (gene_id, strands)

    def test_the_panels_one_origin_spanning_gene_keeps_transcript_order(self) -> None:
        """The one real wrap in the panel, verified against its GenBank record.

        ``KU243129.1`` (2,824 bp, ``ss-DNA``, ``circular``) annotates a CDS as
        ``join(2677..2824,1..80)``: the first exon sits at the end of the linear
        representation and the second back at the origin. NCBI writes the
        intervals in transcript order, so a coordinate sort would emit
        ``1..80`` before ``2677..2824`` and silently transpose the gene's two
        exons. Exactly one catalogued gene wraps; the assertion is on that gene,
        not on a synthetic fixture, so a future re-annotation that adds a second
        wrap fails loudly instead of passing.
        """
        rows = _gene_rows()
        spanning = [row for row in rows if row["origin_spanning"] == "true"]
        assert len(spanning) == 1, [(r["accession"], r["gene_id"]) for r in spanning]
        row = spanning[0]
        assert row["accession"] == "KU243129.1"
        assert row["topology"] == "circular"
        assert row["exons"] == "2677:2824,1:80"
        assert row["genome_length"] == "2824"
        assert row["n_exons"] == "2"
        gtf = gtf_text_for(["KU243129.1"])
        blocks = [
            (line.split("\t")[3], line.split("\t")[4])
            for line in gtf.strip().splitlines()
            if GENE_ID_RE.search(line.split("\t")[8])
            and GENE_ID_RE.search(line.split("\t")[8]).group(1) == row["gene_id"]
        ]
        assert blocks == [("2677", "2824"), ("1", "80")]

    def test_internal_joins_are_not_mistaken_for_wraps(self) -> None:
        """``NC_002076.2`` VP2/VP3 are internal splices, not origin wraps."""
        for symbol in ("orf2/4", "orf2/5"):
            row = next(
                r
                for r in _gene_rows()
                if r["accession"] == "NC_002076.2" and r["gene_symbol"] == symbol
            )
            assert row["origin_spanning"] == "false"
            assert int(row["n_exons"]) == 2

    def test_topology_is_recorded_as_declared(self) -> None:
        rows = _gene_rows()
        topologies = {row["topology"] for row in rows}
        assert topologies <= {"linear", "circular"}
        assert "circular" in topologies, "expected some records to declare circular"
        assert "linear" in topologies, (
            "most Anelloviridae records declare linear despite ssDNA biology, so the "
            "declaration must be recorded rather than inferred"
        )


# ── Genogroup honesty ────────────────────────────────────────────────────────


class TestGenogroup:
    def test_source_genotype_is_never_invented(self) -> None:
        """NCBI carries no ``/genogroup`` qualifier for Anelloviridae, full stop.

        The full 2,042-accession panel carries no ``/genogroup`` qualifier at
        all. Exactly two accessions carry NCBI's own ``/genotype`` qualifier,
        copied verbatim into ``source_genotype``: ``NC_014081.1`` ("6") and
        ``NC_014094.1`` ("28") -- both contradict their own ``/organism``
        species number, which is the concrete evidence that genogroup must
        never be inferred from organism or back-filled from free text. Any
        non-empty value outside this exact pair is a fabricated retrieval
        fact.
        """
        non_empty = {
            row["accession"]: row["source_genotype"]
            for row in _gene_rows()
            if row["source_genotype"]
        }
        assert non_empty == {"NC_014081.1": "6", "NC_014094.1": "28"}

    def test_retrieved_source_qualifiers_are_carried_verbatim(self) -> None:
        rows = _gene_rows()
        assert any(row["source_isolate"] for row in rows)
        assert any(row["source_strain"] for row in rows)


# ── GTF emission and the downstream naming chain ─────────────────────────────


class TestGtfEmission:
    def test_gtf_text_for_emits_real_genes_not_placeholders(self) -> None:
        gtf = gtf_text_for(["NC_002076.2"])
        ids = set(GENE_ID_RE.findall(gtf))
        assert "NC_002076.2_gene1" not in ids
        assert ids == {
            "NC_002076.2_TTVgp1",
            "NC_002076.2_TTVgp2",
            "NC_002076.2_TTVgp3",
        }

    def test_spliced_gene_appears_once_per_exon(self) -> None:
        gtf = gtf_text_for(["NC_002076.2"])
        rows = [line.split("\t") for line in gtf.strip().splitlines()]
        vp3 = [row for row in rows if GENE_ID_RE.search(row[8]) and "VP3" in row[8]]
        assert len(vp3) == 2
        assert (vp3[0][3], vp3[0][4]) == ("353", "711")
        assert (vp3[1][3], vp3[1][4]) == ("2564", "3077")

    def test_gtf_seqname_is_the_accession(self) -> None:
        gtf = gtf_text_for(["NC_002076.2"])
        for line in gtf.strip().splitlines():
            assert line.split("\t")[0] == "NC_002076.2"

    def test_uncovered_accession_falls_back_to_a_placeholder(self) -> None:
        gtf = gtf_text_for(
            ["NOT_IN_CATALOGUE.1"],
            fasta_texts={"NOT_IN_CATALOGUE.1": ">NOT_IN_CATALOGUE.1\nACGTACGTAC\n"},
        )
        assert set(GENE_ID_RE.findall(gtf)) == {"NOT_IN_CATALOGUE.1_gene1"}
        assert "whole_genome" in gtf

    def test_uncovered_accession_without_fasta_text_raises(self) -> None:
        """SW-18: it used to return '' and the genome silently lost its GTF rows."""
        with pytest.raises(ValueError, match="NOT_IN_CATALOGUE.1"):
            gtf_text_for(["NC_002076.2", "NOT_IN_CATALOGUE.1"])

    def test_a_generator_of_accessions_still_gets_placeholders(self) -> None:
        """Both loops iterate `accessions`; a generator used to empty the second."""
        gtf = gtf_text_for(
            (a for a in ["NOT_IN_CATALOGUE.1"]),
            fasta_texts={"NOT_IN_CATALOGUE.1": ">NOT_IN_CATALOGUE.1\nACGTACGTAC\n"},
        )
        assert set(GENE_ID_RE.findall(gtf)) == {"NOT_IN_CATALOGUE.1_gene1"}

    def test_real_genes_resolve_to_their_genus(self) -> None:
        name_map = merged_name_map()
        genus_of = {row["accession"]: row["genus"] for row in _panel_rows()}
        for accession in list(REPRESENTATIVE.values())[:4]:
            expected = genus_of[accession]
            for gene_id in GENE_ID_RE.findall(gtf_text_for([accession])):
                assert virus_name_for_gene(gene_id, name_map) == expected

    def test_candidate_gene_ids_includes_the_real_gene_names(self) -> None:
        candidates = candidate_gene_ids()
        assert "NC_002076.2_TTVgp3" in candidates
        assert "NC_002076.2_gene1" in candidates

    def test_build_reference_uses_the_catalogue(self, tmp_path) -> None:
        from viralscan.scripts.build_reference import _anellovirus_gtf

        fasta = tmp_path / "anello.fa"
        fasta.write_text(">NC_002076.2 Torque teno virus\n" + "ACGT" * 20 + "\n")
        gtf = tmp_path / "anello.gtf"
        annotated, placeholder = _anellovirus_gtf(fasta, gtf)
        assert (annotated, placeholder) == (1, 0)
        text = gtf.read_text()
        assert "NC_002076.2_TTVgp3" in text
        assert "_gene1" not in text


# ── GenBank converter, on committed fixtures (offline) ──────────────────────

TTC_RECORD = textwrap.dedent(
    """\
    LOCUS       NC_002076               3852 bp ss-DNA     circular VRL 13-AUG-2018
    VERSION     NC_002076.2
    FEATURES             Location/Qualifiers
         source          1..3852
                         /organism="Torque teno virus 1"
                         /mol_type="genomic DNA"
                         /isolate="VT416"
                         /db_xref="taxon:687340"
                         /note="genotype 1;
                         vector:pTZ18U"
         CDS             join(353..711,2564..3077)
                         /gene="orf2/5"
                         /locus_tag="TTVgp1"
                         /product="VP3"
                         /protein_id="NP_817120.1"
         CDS             join(353..711,2374..2875)
                         /gene="orf2/4"
                         /locus_tag="TTVgp2"
                         /product="VP2"
                         /protein_id="NP_817121.1"
         CDS             589..2901
                         /gene="orf1"
                         /locus_tag="TTVgp3"
                         /product="VP1"
                         /protein_id="NP_817122.1"
    ORIGIN
    //
    """
)

CIRCULAR_RECORD = textwrap.dedent(
    """\
    LOCUS       KU243129               2824 bp ss-DNA     circular VRL 08-SEP-2017
    VERSION     KU243129.1
    FEATURES             Location/Qualifiers
         source          1..2824
                         /organism="Anellovirus"
                         /isolate="MDJHem2"
         CDS             join(2677..2824,1..80)
                         /product="hypothetical protein"
                         /protein_id="APR63556.1"
    ORIGIN
    //
    """
)

NO_CDS_RECORD = textwrap.dedent(
    """\
    LOCUS       KP343852                2359 bp    DNA     linear   VRL 01-JAN-2020
    VERSION     KP343852.1
    FEATURES             Location/Qualifiers
         source          1..2359
                         /organism="Torque teno virus"
                         /isolate="S85"
    ORIGIN
    //
    """
)


class TestGenbankConverterFixtures:
    def test_reference_ttv_gene_ids_match_the_bundled_panel_convention(self) -> None:
        """``/locus_tag`` wins, so the IDs are ``<accession>_TTVgpN``.

        The bundled RefSeq GTF for the same genome writes ``TTV_TTVgp1``. The only
        difference is the virus token: accession rather than ``TTV``, because a
        2,042-genome panel cannot share one token.
        """
        gtf = _genbank_to_gtf(TTC_RECORD, "NC_002076.2")
        assert list(dict.fromkeys(GENE_ID_RE.findall(gtf))) == [
            "NC_002076.2_TTVgp1",
            "NC_002076.2_TTVgp2",
            "NC_002076.2_TTVgp3",
        ]

    def test_spliced_vp2_and_vp3_keep_both_exons(self) -> None:
        gtf = _genbank_to_gtf(TTC_RECORD, "NC_002076.2")
        rows = [line.split("\t") for line in gtf.strip().splitlines()]
        exons = [(int(row[3]), int(row[4])) for row in rows]
        assert exons == [
            (353, 711),
            (2564, 3077),
            (353, 711),
            (2374, 2875),
            (589, 2901),
        ]

    def test_origin_spanning_join_is_emitted_in_transcript_order(self) -> None:
        """A circular genome's wrapping join must not be coordinate-sorted.

        ``KU243129.1`` is the panel's one real wrap: exon 1 sits at the end of the
        linear representation and exon 2 back at the origin. NCBI writes the
        intervals in transcript order, so sorting ascending would emit
        ``1..80`` before ``2677..2824`` and silently transpose the gene.
        """
        gtf = _genbank_to_gtf(CIRCULAR_RECORD, "KU243129.1")
        rows = [line.split("\t") for line in gtf.strip().splitlines()]
        assert [(row[3], row[4]) for row in rows] == [("2677", "2824"), ("1", "80")]
        assert 'origin_spanning "true"' in gtf

    def test_record_without_cds_still_raises_for_the_converter(self) -> None:
        with pytest.raises(NCBIFetchError):
            _genbank_to_gtf(NO_CDS_RECORD, "KP343852.1")

    def test_catalogue_rows_carry_the_source_qualifiers(self) -> None:
        from viralscan.scripts.ncbi_fetch import catalogue_rows

        rows = catalogue_rows("NC_002076.2", TTC_RECORD)
        assert {row["source_isolate"] for row in rows} == {"VT416"}
        assert {row["source_genotype"] for row in rows} == {""}
        assert {row["topology"] for row in rows} == {"circular"}
        assert {row["genome_length"] for row in rows} == {3852}


# ── Generator behaviour (offline, on a synthetic cache) ───────────────────────


class TestGenerator:
    def _cache(self, tmp_path: Path, accession: str, record: str) -> Path:
        import hashlib

        acc_dir = tmp_path / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        path = acc_dir / f"{accession}.gb"
        path.write_text(record)
        path.with_suffix(".gb.sha256").write_text(hashlib.sha256(record.encode()).hexdigest())
        return tmp_path

    def test_audit_rejects_an_out_of_range_exon(self) -> None:
        import sys
        from pathlib import Path as _Path

        sys.path.insert(0, str(_Path("extras").resolve()))
        from build_anellovirus_genes import audit

        row = {
            "accession": "ACC.1",
            "gene_id": "ACC.1_orf1",
            "n_exons": 1,
            "exons": "1:9999",
            "genome_length": 3000,
            "strand": "+",
        }
        assert audit([row]) != []

    def test_audit_rejects_a_duplicate_gene_id(self) -> None:
        import sys
        from pathlib import Path as _Path

        sys.path.insert(0, str(_Path("extras").resolve()))
        from build_anellovirus_genes import audit

        row = {
            "accession": "ACC.1",
            "gene_id": "ACC.1_orf1",
            "n_exons": 1,
            "exons": "1:900",
            "genome_length": 3000,
            "strand": "+",
        }
        assert audit([row, dict(row, accession="ACC.2")]) != []

    def test_generator_runs_offline_from_a_seeded_cache(self, tmp_path) -> None:
        import sys
        from pathlib import Path as _Path

        sys.path.insert(0, str(_Path("extras").resolve()))
        from build_anellovirus_genes import main

        cache = self._cache(tmp_path / "cache", "NC_002076.2", TTC_RECORD)
        out = tmp_path / "subset.tsv"
        code = main(
            [
                "--email",
                "test@example.org",
                "--accessions",
                "NC_002076.2",
                "--cache-dir",
                str(cache),
                "--out",
                str(out),
            ]
        )
        assert code == 0
        with open(out, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 3
        assert {row["gene_id"] for row in rows} == {
            "NC_002076.2_TTVgp1",
            "NC_002076.2_TTVgp2",
            "NC_002076.2_TTVgp3",
        }

    def test_generator_refuses_to_write_when_coverage_is_too_low(self, tmp_path) -> None:
        import sys
        from pathlib import Path as _Path

        sys.path.insert(0, str(_Path("extras").resolve()))
        from build_anellovirus_genes import main

        cache = self._cache(tmp_path / "cache", "NC_002076.2", TTC_RECORD)
        out = tmp_path / "subset.tsv"
        code = main(
            [
                "--email",
                "test@example.org",
                "--accessions",
                "NC_002076.2",
                "KP343852.1",
                "--cache-dir",
                str(cache),
                "--out",
                str(out),
                "--min-coverage",
                "0.9",
            ]
        )
        assert code == 1
        assert not out.exists()


# ── Live NCBI (opt-in) ───────────────────────────────────────────────────────


@pytest.mark.network
class TestLiveNcbi:
    """Opt-in: ``PYTHONPATH=src pytest -m network tests/test_anellovirus_reference.py``."""

    def test_fetch_genbank_returns_a_cacheable_flatfile(self, tmp_path) -> None:
        path, text = fetch_genbank(
            "NC_002076.2",
            email="viralscan-test@example.org",
            cache_dir=tmp_path,
        )
        assert "FEATURES" in text
        assert path.with_suffix(".gb.sha256").exists()
        again_path, again_text = fetch_genbank(
            "NC_002076.2",
            email="viralscan-test@example.org",
            cache_dir=tmp_path,
        )
        assert again_path == path
        assert again_text == text

    def test_live_record_still_yields_three_real_genes(self, tmp_path) -> None:
        _path, text = fetch_genbank(
            "NC_002076.2",
            email="viralscan-test@example.org",
            cache_dir=tmp_path,
        )
        gtf = _genbank_to_gtf(text, "NC_002076.2")
        ids = GENE_ID_RE.findall(gtf)
        assert "NC_002076.2_orf1" in ids
        assert "NC_002076.2_gene1" not in ids
        assert gtf.count("\texon\t") >= 4, "spliced VP2/VP3 exons are missing"

    def test_live_catalogue_matches_the_shipped_rows(self, tmp_path) -> None:
        from viralscan.scripts.ncbi_fetch import catalogue_rows

        _path, text = fetch_genbank(
            "NC_002076.2",
            email="viralscan-test@example.org",
            cache_dir=tmp_path,
        )
        live = catalogue_rows("NC_002076.2", text)
        shipped = {
            row["gene_id"]: row for row in load_gene_table() if row["accession"] == "NC_002076.2"
        }
        assert shipped, "NC_002076.2 is missing from the shipped catalogue"
        for row in live:
            assert row["gene_id"] in shipped, row["gene_id"]
            assert int(shipped[str(row["gene_id"])]["n_exons"]) == row["n_exons"]
