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
EBNA-1 is expressed from the latent episome in latency I-III, so a latency I-III cell
line should show it well above a lytic gene, yet it sits ~155x below BHLF1. The cause is pervasive
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
    KINETIC_CLASSES,
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
#: HHV-6A and HHV-7 are partial: their U90/U86 are immediate-early (PMID
#: 12706083, 10573164), no HHV-7 latency transcript is established, and the
#: remaining HHV-6A latent genes do not separate latency from productive
#: infection (U94 is immediate-early in HHV-6B, PMID 33627386).
EXPECTED_COMPLETE = {
    "Epstein-Barr virus",
    "Human herpesvirus 8",  # PROG-08, 2026-10-03
}
EXPECTED_PARTIAL = {
    "Human cytomegalovirus",
    "Human herpesvirus 1",
    "Human herpesvirus 2",
    "Human herpesvirus 6A",
    "Human herpesvirus 6b",
    "Human herpesvirus 7",
    "Varicella-zoster virus",
    "Human herpesvirus 8",
}


#: Words that name a protein's function. Used by :func:`_function_contradictions`.
_FUNCTION_WORDS = (
    "polymerase",
    "helicase",
    "primase",
    "kinase",
    "thymidine",
    "capsid",
    "protease",
    "tegument",
    "glycoprotein",
    "processivity",
    "catalytic",
    "transactivator",
)


def _function_contradictions(rows: list[dict]) -> list[tuple[str, str]]:
    """``(gene, word)`` for each function word in a note that RefSeq's product text lacks.

    Rows whose product is just the gene name, or names no function ("protein Zta"), have
    nothing to contradict and are skipped; clauses that
    start with "not" ("not a polymerase subunit") are negations and are ignored.
    """
    out = []
    for row in rows:
        product = row["product"].lower()
        if row["product"] == row["refseq_gene"] or not any(w in product for w in _FUNCTION_WORDS):
            continue
        note = re.sub(r"\bnot [^;,)]*", "", row["note"].lower())
        out += [(row["refseq_gene"], w) for w in _FUNCTION_WORDS if w in note and w not in product]
    return out


def _generator():
    import sys

    sys.path.insert(0, os.path.join(REPO_ROOT, "extras"))
    import build_gene_programs as generator

    return generator


def _genes(spans: dict[str, tuple[int, int]], seqname: str = "seq") -> dict[str, dict]:
    return {
        g: {"blocks": {span}, "span_start": span[0], "span_end": span[1], "seqname": seqname}
        for g, span in spans.items()
    }


class TestKmerGrouping:
    """AUD-01: groups are cliques of the k-mer link graph, not connected components."""

    def test_a_chain_is_not_one_group(self) -> None:
        """A~B and B~C but A!~C (the EBV BMRF1-...-EBNA-1 chain): never one group."""
        g = _generator()
        records = _genes({"A": (1, 100), "B": (50, 150), "C": (120, 220)})
        sharing = {("A", "B"): 0.9, ("B", "C"): 0.9}
        groups = g.build_overlap_groups(records, sharing, 0.5)
        assert groups["A"] == groups["B"] != groups["C"]
        # every pair inside a group is linked
        by_group: dict[str, list[str]] = {}
        for gene, label in groups.items():
            by_group.setdefault(label, []).append(gene)
        for members in by_group.values():
            for i, a in enumerate(members):
                for b in members[i + 1 :]:
                    assert sharing.get((a, b), 0) >= 0.5

    def test_threshold_and_determinism(self) -> None:
        g = _generator()
        records = _genes({"A": (1, 100), "B": (200, 300)})
        sharing = {("A", "B"): 0.4}
        assert (
            g.build_overlap_groups(records, sharing, 0.5)["A"]
            != (g.build_overlap_groups(records, sharing, 0.5)["B"])
        )
        assert g.build_overlap_groups(records, sharing, 0.25) == {"A": "g1", "B": "g1"}
        shuffled = dict(reversed(list(records.items())))
        assert g.build_overlap_groups(shuffled, sharing, 0.25) == {"A": "g1", "B": "g1"}

    def test_sequence_identical_copies_group_without_shared_coordinates(self, tmp_path) -> None:
        """A-02: a copy at other coordinates, on the other strand, shares every canonical 31-mer."""
        g = _generator()
        import random

        rng = random.Random(3)
        x = "".join(rng.choice("ACGT") for _ in range(120))
        y = "".join(rng.choice("ACGT") for _ in range(120))
        rc = x.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        genome = x + y + rc  # copy of x, reverse-complemented, 240 bp downstream
        cache = tmp_path / "cache" / "NC_FAKE.1"
        cache.mkdir(parents=True)
        (cache / "NC_FAKE.1.fasta").write_text(f">NC_FAKE.1\n{genome}\n")
        gtf = tmp_path / "x.gtf"

        def line(feature, start, end, gene, strand):
            return (
                f'NC_FAKE.1\tt\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; '
                f'transcript_id "{gene}";'
            )

        gtf.write_text(
            "\n".join(
                [
                    line("gene", 1, 120, "X", "+"),
                    line("exon", 1, 120, "X", "+"),
                    line("gene", 121, 240, "Y", "+"),
                    line("exon", 121, 240, "Y", "+"),
                    line("gene", 241, 360, "XCOPY", "-"),
                    line("exon", 241, 360, "XCOPY", "-"),
                ]
            )
            + "\n"
        )
        records = g.parse_gtf(str(gtf))
        kmers = g.gene_kmers(str(gtf), str(tmp_path / "cache"))
        assert len(kmers["X"]) == 90 and kmers["X"] == kmers["XCOPY"]
        sharing = g.kmer_sharing(kmers)
        assert sharing == {("X", "XCOPY"): 1.0}
        groups = g.build_overlap_groups(records, sharing)
        assert groups["X"] == groups["XCOPY"] != groups["Y"]
        # exon geometry alone cannot see it
        assert g.build_overlap_groups(records)["X"] != g.build_overlap_groups(records)["XCOPY"]

    def test_audit_flags_cross_programme_sharing_in_different_groups(self) -> None:
        g = _generator()
        programmes = {"L": "latent", "P": "productive", "Q": "productive"}
        groups = {"L": "g1", "P": "g2", "Q": "g3"}
        sharing = {("L", "P"): 0.6, ("P", "Q"): 0.95}
        report = g.audit_markers(programmes, groups, sharing, 0.5)
        assert report["cross_programme_split"] == [("L", "P")]
        assert report["same_programme_split"] == [("P", "Q")]
        assert report["mixed_group"] == []
        groups["P"] = "g1"
        assert g.audit_markers(programmes, groups, sharing, 0.5)["mixed_group"] == [("L", "P")]

    def test_non_overlapping_is_false_when_linked_to_a_marker_in_another_group(self) -> None:
        g = _generator()
        rows = [
            {"virus": "V", "gene_id_bundled": "A", "overlap_group": "g1"},
            {"virus": "V", "gene_id_bundled": "B", "overlap_group": "g2"},
            {"virus": "V", "gene_id_bundled": "C", "overlap_group": "g3"},
        ]
        g._fill_non_overlapping(rows, {"V": {frozenset({"A", "B"})}})
        assert [r["non_overlapping"] for r in rows] == ["false", "false", "true"]

    def test_parse_gtf_reads_product_from_the_cds_row(self, tmp_path) -> None:
        """A-04 root cause: RefSeq puts ``product`` on CDS rows, never on the gene row."""
        g = _generator()
        gtf = tmp_path / "x.gtf"
        gtf.write_text(
            'S\tt\tgene\t11\t310\t.\t+\t.\tgene_id "G"; gene "BGLF4";\n'
            'S\tt\tCDS\t11\t310\t.\t+\t0\tgene_id "G"; transcript_id "T"; '
            'product "tegument serine/threonine protein kinase";\n'
            'S\tt\tgene\t411\t610\t.\t+\t.\tgene_id "H"; product "from the gene row";\n'
            'S\tt\tCDS\t411\t610\t.\t+\t0\tgene_id "H"; transcript_id "U"; product "from the CDS";\n'
        )
        records = g.parse_gtf(str(gtf))
        assert records["G"]["product"] == "tegument serine/threonine protein kinase"
        assert records["H"]["product"] == "from the gene row"  # the gene row wins


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
            "kinetic_class",
            "kinetic_pmid",
            "note",
        }
        assert required <= set(rows[0]), required - set(rows[0])

    def test_controlled_vocabularies(self) -> None:
        for row in _catalog():
            assert row["programme"] in PROGRAMMES, row
            assert row["kinetic_class"] in KINETIC_CLASSES, row
            assert row["panel_completeness"] in {"complete", "partial"}, row
            for field in (
                "latency_observable_in_rna",
                "has_cds",
                "non_overlapping",
                "available_in_starsolo",
                "do_not_normalise",
            ):
                assert row[field] in {"true", "false"}, (field, row)

    def test_classified_rows_cite_a_pmid_and_unclassified_do_not(self) -> None:
        """A kinetic class without a primary source is a guess; ``unclassified`` is not."""
        for row in _catalog():
            cited = bool(row["kinetic_pmid"].strip())
            assert cited == (row["kinetic_class"] != "unclassified"), row
            for pmid in filter(None, row["kinetic_pmid"].split(";")):
                assert pmid.isdigit(), row

    def test_validate_rejects_an_unknown_kinetic_class(self) -> None:
        rows = [dict(r) for r in _catalog()]
        rows[0]["kinetic_class"] = "very_late"
        assert any("kinetic_class" in e for e in validate_catalogue(rows))

    def test_cited_biology_corrections_2026_10_02(self) -> None:
        """Each case was checked against retrieved primary literature (see PLAN PROG-11)."""
        by = {(r["virus"], r["refseq_gene"]): r for r in _catalog()}
        ebv, h6a, h6b, h7 = (
            "Epstein-Barr virus",
            "Human herpesvirus 6A",
            "Human herpesvirus 6b",
            "Human herpesvirus 7",
        )
        # EBV BHRF1 and BNLF2a/b are early lytic (PMID 29864140), never latent.
        # They are removed rather than moved: BNLF2a/b lie inside LMP-1's overlap
        # group and BHRF1 reads may come from latent EBNA-LP transcripts.
        assert not {"BHRF1", "BNLF2a", "BNLF2b"} & {g for v, g in by if v == ebv}
        # BcLF1 is the late major capsid protein, not a polymerase/lytic switch.
        assert by[(ebv, "BcLF1")]["kinetic_class"] == "late"
        assert "capsid" in by[(ebv, "BcLF1")]["note"]
        # Roseolovirus IE1/IE2 are immediate-early, not latent markers.
        for key in ((h6a, "U90"), (h6a, "U86"), (h6b, "U95"), (h7, "U90")):
            assert by[key]["programme"] == "productive", key
            assert by[key]["kinetic_class"] == "immediate_early", key
        # HHV-7 and HHV-6A have no latent marker that separates the states.
        for row in _catalog():
            if row["virus"] in {h6a, h7}:
                assert row["panel_completeness"] == "partial", row
                assert row["latency_observable_in_rna"] == "false", row
        # VZV ORF4 is IE4 (ICP27 homologue); IE62 is ORF62. KSHV ORF17 is the
        # protease; MTA is ORF57.
        zoster = by[("Varicella-zoster virus", "ORF4")]["note"]
        assert "IE4" in zoster and "not IE62" in zoster
        orf17 = by[("Human herpesvirus 8", "ORF17")]
        assert "protease" in orf17["note"] and orf17["kinetic_class"] == "late"
        # No latent and productive marker may share an overlap group: a latent
        # read there would count as productive breadth.
        groups: dict[tuple[str, str], set[str]] = {}
        for row in _catalog():
            groups.setdefault((row["virus"], row["overlap_group"]), set()).add(row["programme"])
        assert not [k for k, v in groups.items() if len(v) > 1], groups

    def test_curated_notes_do_not_contradict_the_refseq_product_2026_10(self) -> None:
        """AUD-03 (A-04): notes written from memory said BALF5 is the processivity subunit,
        BGLF4 the thymidine kinase, BBLF4 the polymerase and so on; RefSeq's CDS product says
        otherwise. A function word in a note must also be in the product text, unless the note
        negates it. Programme and kinetic class are separate columns and are not checked here."""
        assert _function_contradictions(_catalog()) == []
        # the check has teeth: the pre-fix notes are flagged
        assert _function_contradictions(
            [
                {
                    "refseq_gene": "BGLF4",
                    "product": "tegument serine/threonine protein kinase",
                    "note": "viral thymidine kinase; lytic",
                },
                {
                    "refseq_gene": "BALF5",
                    "product": "DNA polymerase catalytic subunit",
                    "note": "lytic polymerase processivity subunit",
                },
            ]
        ) == [("BGLF4", "thymidine"), ("BALF5", "processivity")]

    def test_curated_note_wording_fixes_2026_10(self) -> None:
        """A-05 / A-16: U90 cites PMID 12388818 (not 12367753) for IE1; LMP-1 and EBNA-1 wording."""
        by = {(r["virus"], r["refseq_gene"]): r for r in _catalog()}
        ebv = "Epstein-Barr virus"
        u90 = by[("Human herpesvirus 6A", "U90")]
        assert "12388818" in u90["note"] and "12367753" not in u90["note"]
        assert "12367753" in u90["kinetic_pmid"]  # still the class source (Biomni-verified)
        assert "independent" not in by[(ebv, "LMP-1")]["note"]
        assert "every infected cell" not in by[(ebv, "EBNA-1")]["note"]
        # HHV-6B U94 is latency-associated and immediate-early: not a contradiction, unchanged
        h6b = by[("Human herpesvirus 6b", "U94")]
        assert (h6b["programme"], h6b["kinetic_class"]) == ("latent", "immediate_early")

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
        rows = [r for r in _catalog() if r["virus"] == "Human herpesvirus 6A"]
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

    def test_latent_ebna_groups_follow_shared_sequence(self) -> None:
        """The cross-mapping the design defends against, pinned as a fact (AUD-01).

        Groups are cliques of the 31-mer link graph (>= 50 % of the smaller gene). The two
        EBNA-1 records of one locus share 96 % and EBNA-2 / EBNA-LP share 91 %, so each pair is
        one group; EBNA-1 shares only 24 % with EBNA-2 and 57 % with EBNA-LP, which a clique
        cannot join, so EBNA-1 is a group of its own. None is independent. BMRF1, which the old
        coordinate closure chained to all four, shares no group with any latent marker.
        """
        ebv = [r for r in _catalog() if r["virus"] == "Epstein-Barr virus"]
        group = {r["gene_id_bundled"]: r["overlap_group"] for r in ebv}
        ebna1 = {group[g] for g in group if g.startswith("EPSTEIN_HHV4_EBNA-1.")}
        assert len(ebna1) == 1 and len(ebna1 | {group["EPSTEIN_HHV4_EBNA-2"]}) == 2
        assert group["EPSTEIN_HHV4_EBNA-2"] == group["EPSTEIN_HHV4_EBNA-LP"]
        ebna_rows = [r for r in ebv if r["refseq_gene"] in {"EBNA-1", "EBNA-2", "EBNA-LP"}]
        assert all(r["non_overlapping"] == "false" for r in ebna_rows)
        latent_groups = {r["overlap_group"] for r in ebv if r["programme"] == "latent"}
        assert group["EPSTEIN_HHV4_BMRF1"] not in latent_groups


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

    def test_bare_starsolo_id_is_never_a_fallback_for_unavailable_rows(self) -> None:
        """A-10: under panel_form='bundled' the bare STARsolo ID (shared by >=7 genomes) is skipped."""
        row = dict(next(r for r in _catalog() if r["virus"] == "Epstein-Barr virus"))
        row.update(
            gene_id_bundled="EPSTEIN_unassigned_gene_3",
            gene_id_merged="EPSTEIN_unassigned_gene_3",
            gene_id_starsolo="unassigned_gene_3",
            available_in_starsolo="false",
        )
        for form in ("bundled", "merged", "starsolo"):
            resolved, unresolved = resolve_markers(
                "Epstein-Barr virus", {"unassigned_gene_3"}, form, catalogue=[row]
            )
            assert resolved == [], form
            assert len(unresolved) == 1, form
        row["available_in_starsolo"] = "true"  # the fallback still works for available rows
        resolved, _ = resolve_markers(
            "Epstein-Barr virus", {"unassigned_gene_3"}, "bundled", catalogue=[row]
        )
        assert [m.var_name for m in resolved] == ["unassigned_gene_3"]

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
        two_digit, _ = resolve_markers("Human herpesvirus 6A", ["HUM_HERP6_HHV6gp41"], "bundled")
        assert [m.var_name for m in two_digit] == ["HUM_HERP6_HHV6gp41"]

        # U43 (the 3-digit form) is not a catalogue marker, so supplying it
        # alone must resolve nothing. Padding-normalising matching would return
        # the U42 marker and report a gene that is not in the panel.
        three_digit, _ = resolve_markers("Human herpesvirus 6A", ["HUM_HERP6_HHV6gp041"], "bundled")
        assert three_digit == [], (
            "the 3-digit form was silently matched to the 2-digit marker; "
            "padding normalisation has crept back in"
        )

        # Both widths present: each must map to its own marker, never across.
        both, _ = resolve_markers(
            "Human herpesvirus 6A",
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


class TestKshvLatency:
    """PROG-08: the cited KSHV latency set makes a KSHV latent call reachable."""

    def _latent(self):
        return {
            r["refseq_gene"]: r
            for r in _catalog()
            if r["virus"] == "Human herpesvirus 8" and r["programme"] == "latent"
        }

    def test_latent_set_is_the_cited_cluster_plus_k12_and_lana2(self) -> None:
        assert set(self._latent()) == {"ORF73", "ORF72", "ORF71", "K12", "vIRF-3"}

    def test_k1_is_not_a_latency_anchor(self) -> None:
        assert "K1" not in self._latent()

    def test_lana_vcyclin_vflip_are_one_breadth_unit(self) -> None:
        """One promoter, 3'-coterminal mRNAs (PMID 9733875): never independent evidence."""
        lat = self._latent()
        assert len({lat[g]["overlap_group"] for g in ("ORF73", "ORF72", "ORF71")}) == 1
        assert lat["K12"]["overlap_group"] != lat["ORF73"]["overlap_group"]
        assert lat["vIRF-3"]["overlap_group"] not in {
            lat["ORF73"]["overlap_group"],
            lat["K12"]["overlap_group"],
        }

    def test_every_kshv_latent_row_cites_a_pmid(self) -> None:
        assert all(r["kinetic_pmid"] for r in self._latent().values())


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
        # default min_breadth=2. BZLF1 shares its group with BRLF1, so even
        # BZLF1 + BRLF1 together would only count once.
        assert calls[0]["productive_breadth"] == 1
        brlf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_BRLF1")
        assert bzlf1.overlap_group == brlf1.overlap_group

    def test_two_distinct_productive_groups_is_productive(self) -> None:
        markers = self._markers()
        independent = [m for m in markers if m.non_overlapping][:2]
        assert len({m.overlap_group for m in independent}) == 2
        pairs = [(markers.index(m), 0, 20.0) for m in independent]
        m = _matrix(pairs, 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["state"] == "productive", calls
        assert calls[0]["productive_breadth"] == 2

    def test_two_markers_in_one_overlap_group_still_fail(self) -> None:
        """The cross-mapping defence: nested ORFs are one piece of evidence."""
        markers = self._markers()
        bglf4 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_BGLF4")
        shared = [m for m in markers if m.overlap_group == bglf4.overlap_group]  # + BBLF1
        assert len(shared) == 2, [m.var_name for m in shared]
        idx = [markers.index(m) for m in shared]
        m = _matrix([(i, 0, 20.0) for i in idx], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["state"] == "indeterminate", calls
        assert calls[0]["productive_breadth"] == 1

    def test_non_overlapping_markers_in_a_shared_group_count_once(self) -> None:
        markers = self._markers()
        bzlf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_BZLF1")
        idx = [
            i for i, m in enumerate(markers) if m.overlap_group == bzlf1.overlap_group
        ]  # BZLF1 + BRLF1, both non_overlapping=false
        assert len(idx) == 2 and not any(markers[i].non_overlapping for i in idx)
        m = _matrix([(i, 0, 20.0) for i in idx], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2)
        assert calls[0]["productive_breadth"] == 1

    def test_latent_anchor_gives_latent_when_observable(self) -> None:
        markers = self._markers(programme="latent")
        bhrf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_LMP-1")
        m = _matrix([(markers.index(bhrf1), 0, 20.0)], 1, len(markers))
        calls = call_cell_programme(m, markers, min_breadth=2, latency_observable=True)
        assert calls[0]["state"] == "latent", calls
        assert calls[0]["latent_breadth"] == 1

    def test_kshv_cluster_counts_once_and_k12_adds_breadth(self) -> None:
        markers = self._markers(virus="Human herpesvirus 8", programme="latent")
        by_gene = {
            r["gene_id_bundled"]: r["refseq_gene"]
            for r in _catalog()
            if r["virus"] == "Human herpesvirus 8"
        }
        idx = {by_gene[m.var_name]: i for i, m in enumerate(markers)}
        cluster = _matrix([(idx[g], 0, 5.0) for g in ("ORF73", "ORF72", "ORF71")], 1, len(markers))
        calls = call_cell_programme(cluster, markers, min_breadth=2, latency_observable=True)
        assert calls[0]["state"] == "latent" and calls[0]["latent_breadth"] == 1, calls
        two = _matrix([(idx["ORF73"], 0, 5.0), (idx["K12"], 0, 5.0)], 1, len(markers))
        calls = call_cell_programme(two, markers, min_breadth=2, latency_observable=True)
        assert calls[0]["latent_breadth"] == 2, calls

    def test_latent_unreachable_when_not_observable(self) -> None:
        markers = self._markers(programme="latent")
        bhrf1 = next(m for m in markers if m.var_name == "EPSTEIN_HHV4_LMP-1")
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
        lat_idx = [allm.index(next(m for m in lat if m.var_name == "EPSTEIN_HHV4_LMP-1"))]
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
        lat_idx = [allm.index(next(m for m in lat if m.var_name == "EPSTEIN_HHV4_LMP-1"))]
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
        assert {"productive", "latent", "reactivating", "mixed", "indeterminate"} == STATES

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
    only defined when ``gene_programs`` is set, and it waits on ``log/detection.done``
    (which stands for ``results/viral_summary.tsv``), so it cannot run before layer 1.
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
        assert 'str(config.get("gene_programs", "")).lower() == "true"' in source, (
            "the gene_programs rule must be optional, like hostresponse and "
            "host_filter — a second layer should not be forced on every run. It must "
            'be gated on str(config.get("gene_programs", "")).lower() == "true", '
            "because config values arrive as strings."
        )

    def test_rule_depends_on_detection_not_on_a_phantom_summary(self) -> None:
        source = self._source()
        # Anchor on the rule itself: `gene_programs` also appears in
        # _all_targets, and that occurrence must not be mistaken for the rule.
        block = source[source.index("rule gene_programs:") :].split("script:")[0]
        inputs = "\n".join(
            line
            for line in block.split("input:")[1].split("output:")[0].splitlines()
            if not line.strip().startswith("#")
        )
        assert "log/detection.done" in inputs, (
            "layer 2's scope is the viruses layer 1 detected; detection.done is what "
            "orders it after layer 1 (it stands for results/viral_summary.tsv)"
        )
        assert "viral_summary.tsv" not in inputs, (
            "no rule produces viral_summary.tsv, so listing it as an input is a "
            "MissingInputException on every fresh run (AUD-08)"
        )

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
        # PROG-12: the defaults are the pre-PROG-12 rule (any nonzero, one latent group).
        assert DEFAULTS["programme_latent_min_breadth"] == 1
        assert DEFAULTS["programme_min_umi"] == 0.0

    def test_runconfig_accepts_the_defaults(self) -> None:
        from viralscan.runconfig import RunConfig

        cfg = RunConfig()
        assert cfg.gene_programs is False
        assert cfg.programme_min_breadth == 2
        assert cfg.programme_latent_min_breadth == 1
        assert cfg.programme_min_umi == 0.0

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
        assert "--programme-latent-min-breadth" in text, text[-800:]
        assert "--programme-min-umi" in text, text[-800:]

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


def test_no_model_placeholder_row_is_not_counted_as_a_cell() -> None:
    """PROG-19: run_one's barcode="" marker row must not read as one cell."""
    placeholder = pd.DataFrame(
        [
            {
                "barcode": "",
                "virus_name": "Monkeypox virus",
                "state": "not_applicable",
                "productive_breadth": 0,
                "latent_breadth": 0,
                "selected_state": None,
                "selected_productive_breadth": 0,
                "selected_latent_breadth": 0,
                "evidence_layer": "counts_unique_viral",
                "latency_not_observable": True,
            }
        ]
    )
    for viruses in (None, ["Monkeypox virus"]):
        row = summarise_programs(placeholder, _catalog(), viruses=viruses).iloc[0]
        assert row["virus_name"] == "Monkeypox virus"
        assert row["n_cells_total"] == 0
        assert row["panel_completeness"] == "not_applicable"
        assert row["caveat"] == "no programme model defined for this virus"


def _legacy_decide(prod: int, lat: int, min_breadth: int, latency_observable: bool) -> str:
    """``_decide`` exactly as it stood before PROG-12, frozen for the equivalence test."""
    if prod >= min_breadth and lat >= 1 and latency_observable:
        return "mixed"
    if prod >= min_breadth:
        return "productive"
    if lat >= 1 and prod < min_breadth and latency_observable:
        return "latent"
    return "indeterminate"


def _m(name: str, programme: str, group: str, kinetic: str = "unclassified") -> Marker:
    return Marker(name, programme, group, True, True, kinetic)


class TestReactivatingAndParameters:
    """PLAN PROG-12: the ``reactivating`` state, latent breadth and the UMI floor."""

    IE1 = _m("ie1", "productive", "g1", "immediate_early")
    IE2 = _m("ie2", "productive", "g2", "immediate_early")
    EARLY = _m("early", "productive", "g3", "early")
    LAT = _m("lat", "latent", "g4", "latent")
    LAT2 = _m("lat2", "latent", "g5", "latent")
    MARKERS = [IE1, IE2, EARLY, LAT, LAT2]

    def _call(self, values: dict[str, float], **kwargs):
        pairs = [
            (i, 0, values[m.var_name]) for i, m in enumerate(self.MARKERS) if m.var_name in values
        ]
        matrix = _matrix(pairs, 1, len(self.MARKERS))
        kwargs.setdefault("selected_matrix", matrix)
        return call_cell_programme(matrix, self.MARKERS, **kwargs)[0]

    def test_state_is_declared(self) -> None:
        assert "reactivating" in STATES

    def test_undeclared_state_raises_value_error_not_assert(self, monkeypatch) -> None:
        """A-14: the STATES guard must survive ``python -O``, so it is a ValueError."""
        from viralscan import gene_programs

        monkeypatch.setattr(gene_programs, "_decide", lambda *a, **k: "bogus")
        with pytest.raises(ValueError, match="bogus"):
            self._call({"ie1": 3.0})

    def test_ie_only_single_group_is_reactivating(self) -> None:
        call = self._call({"ie1": 3.0})
        assert call["state"] == "reactivating", call
        assert call["selected_state"] == "reactivating"

    def test_ie_plus_later_evidence_is_not_reactivating(self) -> None:
        assert self._call({"ie1": 3.0, "early": 1.0})["state"] == "productive"
        assert self._call({"ie1": 3.0}, min_breadth=3)["state"] == "reactivating"
        # an early-class group alongside the IE one: not IE-only
        assert self._call({"ie1": 3.0, "early": 1.0}, min_breadth=3)["state"] == "indeterminate"

    def test_two_ie_groups_are_productive_not_reactivating(self) -> None:
        """HHV-7 has two IE groups; at min_breadth=2 that stays ``productive``."""
        assert self._call({"ie1": 3.0, "ie2": 3.0})["state"] == "productive"

    def test_latent_evidence_wins_over_ie_only(self) -> None:
        assert self._call({"ie1": 3.0, "lat": 2.0})["state"] == "latent"

    def test_latent_marker_kinetic_class_is_ignored(self) -> None:
        """Only productive markers define IE-only evidence (VZV ORF4, HHV-6B U94)."""
        ie_latent = _m("ie_lat", "latent", "g9", "immediate_early")
        markers = [self.IE1, ie_latent]
        matrix = _matrix([(0, 0, 1.0), (1, 0, 1.0)], 1, 2)
        # AUD-19 (D-01): latency not observable => never reactivating; the latent
        # marker's IE class must not turn this into a call either.
        assert call_cell_programme(matrix, markers, latency_observable=False)[0]["state"] == (
            "indeterminate"
        )
        assert call_cell_programme(_matrix([(0, 0, 1.0)], 1, 2), markers)[0]["state"] == (
            "reactivating"
        )
        only_latent = call_cell_programme(_matrix([(1, 0, 1.0)], 1, 2), markers)[0]
        assert only_latent["state"] == "latent"

    def test_not_applicable_virus_never_reactivates(self) -> None:
        call = self._call({"ie1": 3.0}, reactivating_applicable=False)
        assert call["state"] == "indeterminate", call

    def test_hcmv_is_not_applicable_in_the_summary(self) -> None:
        cells = pd.DataFrame(
            [
                {
                    "barcode": "c1",
                    "virus_name": v,
                    "state": "reactivating" if v == "Epstein-Barr virus" else "indeterminate",
                    "productive_breadth": 1,
                    "latent_breadth": 0,
                    "selected_state": "indeterminate",
                    "selected_productive_breadth": 0,
                    "selected_latent_breadth": 0,
                    "evidence_layer": "counts_unique_viral",
                    "latency_not_observable": False,
                }
                for v in ("Epstein-Barr virus", "Human cytomegalovirus")
            ]
        )
        summary = summarise_programs(cells, _catalog(), latent_min_breadth=2, min_umi=3.0)
        by_virus = summary.set_index("virus_name")
        assert by_virus.loc["Epstein-Barr virus", "n_cells_reactivating"] == 1
        assert by_virus.loc["Human cytomegalovirus", "n_cells_reactivating"] == "not_applicable"
        assert by_virus["latent_min_breadth"].eq(2).all()
        assert by_virus["marker_min_umi"].eq(3.0).all()

    def test_run_one_never_emits_reactivating_for_hcmv(self) -> None:
        """Real UL122 (HCMV) stays ``indeterminate``; real BZLF1 (EBV) is reactivating."""
        import anndata as ad
        import numpy as np

        from viralscan.scripts.gene_programs import run_one

        catalogue = _catalog()

        def run(virus: str, gene_id: str) -> str:
            row = next(r for r in catalogue if r["gene_id_bundled"] == gene_id)
            assert row["virus"] == virus and row["kinetic_class"] == "immediate_early"
            # Full index (AUD-19 D-02: a one-gene index is a reduced marker set and
            # can no longer call reactivating); only ``gene_id`` carries molecules.
            ids = [r["gene_id_bundled"] for r in catalogue if r["virus"] == virus]
            dense = np.zeros((1, len(ids)), dtype=np.float32)
            dense[0, ids.index(gene_id)] = 5.0
            x = sparse.csr_matrix(dense)
            adata = ad.AnnData(x)
            adata.var_names = ids
            adata.obs_names = ["c1"]
            adata.layers["counts_unique_viral"] = x.copy()
            df = run_one(adata, [virus], catalogue, min_breadth=2, form="bundled")
            return str(df["state"].iloc[0])

        assert run("Human cytomegalovirus", "HUM_CYTO_HHV5wtgp107") == "indeterminate"
        assert run("Epstein-Barr virus", "EPSTEIN_HHV4_BZLF1") == "reactivating"

    def test_run_one_keeps_cells_with_selected_layer_only_evidence(self) -> None:
        """AUD-02 (A-03): unique 0 / allocated >0 is the BZLF1 comparison; it must not be dropped."""
        import anndata as ad
        import numpy as np

        from viralscan.scripts.gene_programs import run_one

        virus = "Epstein-Barr virus"
        catalogue = _catalog()
        names = [
            r["gene_id_bundled"]
            for r in catalogue
            if r["virus"] == virus
            and r["programme"] == "productive"
            and r["non_overlapping"] == "true"
        ][:3]
        assert len(names) == 3
        # c0: unique 0, allocated on all three (selected-only); c1: unique and allocated;
        # c2: nothing on either layer (the only cell that may be skipped)
        unique = sparse.csr_matrix(np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]], dtype=np.float32))
        alloc = sparse.csr_matrix(np.array([[2, 2, 2], [1, 1, 1], [0, 0, 0]], dtype=np.float32))
        adata = ad.AnnData(alloc)
        adata.var_names = names
        adata.obs_names = ["c0", "c1", "c2"]
        adata.layers["counts_unique_viral"] = unique
        df = run_one(adata, [virus], catalogue, min_breadth=2, form="bundled")
        assert list(df["barcode"]) == ["c0", "c1"]
        c0 = df.set_index("barcode").loc["c0"]
        assert c0["state"] == "indeterminate" and c0["productive_breadth"] == 0
        assert c0["selected_state"] == "productive" and c0["selected_productive_breadth"] == 3
        summary = summarise_programs(
            df, catalogue, viruses=[virus], marker_resolution=df.attrs["marker_resolution"]
        ).iloc[0]
        assert summary["n_cells_productive"] == 1
        assert summary["n_cells_productive_selected_layer"] == 2

    def test_marker_carries_the_catalogue_kinetic_class(self) -> None:
        catalogue = _catalog()
        resolved, _ = resolve_markers(
            "Epstein-Barr virus", ["EPSTEIN_HHV4_BZLF1"], catalogue=catalogue
        )
        assert [m.kinetic_class for m in resolved] == ["immediate_early"]

    # -- parameters -------------------------------------------------------

    def test_latent_min_breadth_raises_the_latent_bar(self) -> None:
        one = {"lat": 1.0}
        two = {"lat": 1.0, "lat2": 1.0}
        assert self._call(one)["state"] == "latent"
        assert self._call(one, latent_min_breadth=2)["state"] == "indeterminate"
        assert self._call(two, latent_min_breadth=2)["state"] == "latent"

    def test_min_umi_floor_drops_thin_evidence(self) -> None:
        assert self._call({"lat": 1.0})["state"] == "latent"
        assert self._call({"lat": 1.0}, min_umi=2.0)["state"] == "indeterminate"
        assert self._call({"lat": 2.0}, min_umi=2.0)["state"] == "latent"

    def test_default_floor_still_counts_fractional_allocated_values(self) -> None:
        """The allocated layer holds fractions; the default must count a 0.3."""
        call = self._call({"lat": 0.3})
        assert call["state"] == "latent" and call["selected_state"] == "latent", call

    @pytest.mark.parametrize("kwargs", [{"latent_min_breadth": 0}, {"min_umi": -1.0}])
    def test_rejects_out_of_range_parameters(self, kwargs) -> None:
        with pytest.raises(ValueError):
            call_cell_programme(sparse.csr_matrix((1, 1)), [], **kwargs)

    # -- defaults reproduce today -----------------------------------------

    def test_default_decision_table_matches_the_pre_prog12_rule(self) -> None:
        from viralscan.gene_programs import _decide

        for prod in range(5):
            for lat in range(4):
                for min_breadth in (1, 2, 3):
                    for observable in (True, False):
                        assert _decide(prod, lat, min_breadth, observable) == _legacy_decide(
                            prod, lat, min_breadth, observable
                        )

    def test_reactivating_only_ever_replaces_indeterminate(self) -> None:
        """Random matrices over every virus's real markers, defaults throughout."""
        import numpy as np

        catalogue = _catalog()
        for virus in sorted(_viruses()):
            resolved, _ = resolve_markers(
                virus, [r["gene_id_bundled"] for r in catalogue], catalogue=catalogue
            )
            observable = next(
                r["latency_observable_in_rna"] == "true" for r in catalogue if r["virus"] == virus
            )
            rng = np.random.default_rng(7)
            shape = (300, len(resolved))
            dense = (rng.random(shape) < 0.15) * rng.integers(1, 5, shape)
            matrix = sparse.csr_matrix(dense.astype(float))
            calls = call_cell_programme(matrix, resolved, latency_observable=observable)
            for i, call in enumerate(calls):
                prod, lat = call["productive_breadth"], call["latent_breadth"]
                old = _legacy_decide(prod, lat, 2, observable)
                if call["state"] == "reactivating":
                    assert old == "indeterminate", (virus, i, call)
                else:
                    assert call["state"] == old, (virus, i, call, old)
