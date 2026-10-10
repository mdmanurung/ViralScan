"""REF-13: the build-time viral GTF normaliser, checked through ngs_tools' real cDNA path."""

from __future__ import annotations

import random
import re
import sys
from pathlib import Path

import pytest

from viralscan import gtf_normalise
from viralscan.gtf_normalise import (
    cds_joins_from_genbank,
    normalise_gtf_file,
    normalise_viral_gtf,
)

ROOT = Path(__file__).resolve().parents[1]
ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
SEQ = "NC_TEST.1"
GENOME = "".join(random.Random(7).choice("ACGT") for _ in range(900))
COMP = str.maketrans("ACGT", "TGCA")


def row(feature, start, end, gene, strand="+", seq=SEQ, **attrs):
    extra = "".join(f' {k} "{v}";' for k, v in attrs.items())
    return f'{seq}\tfixture\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}";{extra}'


def cds_gene(
    gene, gs, ge, blocks, stop=None, strand="+", protein="P1", tx="unassigned_transcript_1"
):
    rows = [row("gene", gs, ge, gene, strand)]
    for s, e in blocks:
        rows.append(row("CDS", s, e, gene, strand, transcript_id=tx, protein_id=protein))
    if stop:
        rows.append(row("stop_codon", *stop, gene, strand, transcript_id=tx, protein_id=protein))
    return rows


def cdna(lines, tmp_path):
    """``{gene: {transcript: (header, sequence)}}`` from the installed ngs_tools, as kb ref does."""
    ngs = pytest.importorskip("ngs_tools")
    gtf, fasta, out = tmp_path / "x.gtf", tmp_path / "g.fa", tmp_path / "cdna.fa"
    gtf.write_text("\n".join(lines) + "\n")
    fasta.write_text(f">{SEQ}\n{GENOME}\n")
    genes, transcripts = ngs.gtf.genes_and_transcripts_from_gtf(str(gtf), use_version=True)
    ngs.fasta.split_genomic_fasta_to_cdna(str(fasta), str(out), genes, transcripts)
    result: dict[str, dict[str, tuple[str, str]]] = {}
    with ngs.fasta.Fasta(str(out), "r") as handle:
        for entry in handle:
            result.setdefault(entry.attributes["gene_id"], {})[entry.name] = (
                entry.header,
                entry.sequence,
            )
    return result


def exons(lines, gene):
    """Exon blocks of a gene in genomic order (minus-strand rows are written 5' to 3')."""
    return sorted(
        (int(c[3]), int(c[4]))
        for c in (line.split("\t") for line in lines)
        if len(c) == 9 and c[2] == "exon" and dict(ATTR.findall(c[8]))["gene_id"] == gene
    )


def revcomp(text):
    return text.translate(COMP)[::-1]


def test_gap_free_gene_cdna_is_byte_identical(tmp_path):
    raw = cds_gene("G", 21, 200, [(41, 180)], stop=(181, 183))
    new = normalise_viral_gtf(raw)
    assert cdna(new, tmp_path) == cdna(raw, tmp_path)
    assert exons(new, "G") == [(21, 200)]
    assert set(cdna(new, tmp_path)["G"]) == {"G"}


def test_exon_is_clipped_to_the_gene_row_when_the_stop_runs_past_it(tmp_path):
    raw = cds_gene("G", 21, 180, [(41, 180)], stop=(181, 183))  # HAV: stop 3 nt past the gene row
    new = normalise_viral_gtf(raw)
    assert exons(new, "G") == [(21, 180)]
    assert cdna(new, tmp_path) == cdna(raw, tmp_path)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_single_protein_intron_is_removed(strand, tmp_path):
    raw = cds_gene("G", 11, 300, [(31, 100), (201, 260)], stop=(261, 263), strand=strand)
    new = normalise_viral_gtf(raw)
    assert exons(new, "G") == [(11, 100), (201, 300)]
    expected = GENOME[10:100] + GENOME[200:300]
    assert cdna(new, tmp_path)["G"]["G"][1] == (expected if strand == "+" else revcomp(expected))
    assert len(cdna(raw, tmp_path)["G"]["G"][1]) == 290  # today: the intron stays in


def test_repeat_copies_never_join_and_use_flatfile_joins(tmp_path):
    # HHV-6B: two terminal-repeat copies, the spliced CDS collapsed to one row per copy.
    raw = [
        row("gene", 11, 300, "G"),
        row("CDS", 31, 260, "G", product="X"),
        row("gene", 411, 700, "G"),
        row("CDS", 431, 660, "G", product="X"),
    ]
    joins = {SEQ: [[(31, 100), (201, 260)], [(431, 500), (601, 660)]]}
    new = normalise_viral_gtf(raw, joins)
    assert exons(new, "G") == [(11, 100), (201, 300), (411, 500), (601, 700)]
    transcripts = cdna(new, tmp_path)["G"]
    assert {k: v[1] for k, v in transcripts.items()} == {
        "G": GENOME[10:100] + GENOME[200:300],
        "G-c2": GENOME[410:500] + GENOME[600:700],
    }
    # without the flatfile the copies still do not join: one span transcript per copy
    spans = cdna(normalise_viral_gtf(raw), tmp_path)["G"]
    assert {k: v[1] for k, v in spans.items()} == {"G": GENOME[10:300], "G-c2": GENOME[410:700]}
    # each CDS row follows its own copy
    retagged = [line for line in new if "\tCDS\t" in line]
    assert [dict(ATTR.findall(c.split("\t")[8]))["transcript_id"] for c in retagged] == [
        "G",
        "G-c2",
    ]


def test_origin_wrapping_gene_keeps_the_last_gene_row_and_is_recorded(tmp_path):
    # HBV P/S: ngs_tools concatenates exons in genomic order, so a wrap cannot be represented.
    raw = [
        row("gene", 301, 400, "G", part="1"),
        row("gene", 1, 100, "G", part="2"),
        row("CDS", 301, 400, "G", transcript_id="t", protein_id="P", part="1"),
        row("CDS", 1, 97, "G", transcript_id="t", protein_id="P", part="2"),
        row("stop_codon", 98, 100, "G", transcript_id="t", protein_id="P"),
    ]
    new = normalise_viral_gtf(raw)
    assert exons(new, "G") == [(1, 100)]
    assert 'viralscan_norm "wrap_last_row+span"' in "\n".join(new)
    assert cdna(new, tmp_path) == cdna(raw, tmp_path)


def test_multi_protein_gene_keeps_its_span(tmp_path):
    # adenovirus 3'-coterminal family: the gap is another protein's territory, not an intron
    raw = [
        row("gene", 11, 300, "G"),
        row("CDS", 31, 100, "G", transcript_id="t1", protein_id="P1"),
        row("CDS", 201, 260, "G", transcript_id="t2", protein_id="P2"),
    ]
    new = normalise_viral_gtf(raw)
    assert exons(new, "G") == [(11, 300)]
    assert cdna(new, tmp_path) == cdna(raw, tmp_path)


def test_frameshift_gap_is_not_an_intron():
    raw = cds_gene("G", 11, 300, [(31, 100), (102, 260)])
    assert exons(normalise_viral_gtf(raw), "G") == [(11, 300)]


def test_retained_intron_isoform_stays_indexable(tmp_path, monkeypatch):
    monkeypatch.setattr(gtf_normalise, "RETAINED_INTRON_GENES", frozenset({"G"}))
    raw = cds_gene("G", 11, 300, [(31, 100), (201, 260)], stop=(261, 263))
    transcripts = cdna(normalise_viral_gtf(raw), tmp_path)["G"]
    assert {k: v[1] for k, v in transcripts.items()} == {
        "G": GENOME[10:100] + GENOME[200:300],
        "G-span": GENOME[10:300],
    }
    assert transcripts["G-span"][1] == cdna(raw, tmp_path)["G"]["G"][1]  # today's transcript


def test_transcript_ids_are_gene_scoped_and_never_shared(tmp_path):
    raw = cds_gene("G1", 11, 100, [(21, 90)]) + cds_gene("G2", 201, 300, [(211, 290)])
    assert {dict(ATTR.findall(c.split("\t")[8]))["transcript_id"] for c in raw if "CDS" in c} == {
        "unassigned_transcript_1"
    }
    new = normalise_viral_gtf(raw)
    owners: dict[str, set[str]] = {}
    for line in new:
        attrs = dict(ATTR.findall(line.split("\t")[8]))
        if "transcript_id" in attrs:
            owners.setdefault(attrs["transcript_id"], set()).add(attrs["gene_id"])
    assert owners == {"G1": {"G1"}, "G2": {"G2"}}
    assert set(cdna(new, tmp_path)) == {"G1", "G2"}


def test_genes_with_exons_and_gene_only_genes(tmp_path):
    spliced = [
        row("gene", 11, 300, "E"),
        row("transcript", 11, 300, "E", transcript_id="E.1"),
        row("exon", 11, 100, "E", transcript_id="E.1"),
        row("exon", 201, 300, "E", transcript_id="E.1"),
    ]
    only = [row("gene", 401, 500, "O")]
    new = normalise_viral_gtf(spliced + only)
    assert new[: len(spliced)] == spliced  # a gene that already has exons is untouched
    assert exons(new, "O") == [(401, 500)]  # a gene-only gene keeps its span, now as an exon


def test_orphan_cds_joins_the_gene_row_at_its_locus(tmp_path):
    # B19V: CDS rows say unassigned_gene_1, the gene row at the same span is gp4 with no CDS
    raw = [row("gene", 51, 150, "REAL")] + [
        row("CDS", 51, 147, "ORPHAN", transcript_id="t", protein_id="P"),
        row("stop_codon", 148, 150, "ORPHAN", transcript_id="t", protein_id="P"),
    ]
    new = normalise_viral_gtf(raw)
    assert all(dict(ATTR.findall(c.split("\t")[8]))["gene_id"] == "REAL" for c in new)
    assert set(cdna(new, tmp_path)) == {"REAL"}
    assert cdna(new, tmp_path) == cdna(raw[:1], tmp_path)


def test_cds_without_any_gene_row_gets_one():
    raw = [row("CDS", 51, 147, "ORPHAN", transcript_id="t", protein_id="P")]
    new = normalise_viral_gtf(raw)
    assert [c.split("\t")[2] for c in new].count("gene") == 1
    assert exons(new, "ORPHAN") == [(51, 147)]


def test_normalising_twice_changes_nothing():
    raw = (
        cds_gene("A", 11, 300, [(31, 100), (201, 260)], stop=(261, 263))
        + cds_gene("B", 401, 500, [(411, 490)])
        + [row("gene", 601, 700, "C"), row("CDS", 611, 690, "ORPHAN", transcript_id="t")]
        + [row("gene", 11, 50, "W", part="1"), row("gene", 801, 850, "W", part="2")]
        + [row("CDS", 11, 50, "W", transcript_id="t", protein_id="P")]
        + ["# comment", ""]
    )
    once = normalise_viral_gtf(raw)
    assert normalise_viral_gtf(once) == once


GENBANK = """LOCUS       KX000001                 900 bp    DNA     linear   VRL 01-JAN-2000
VERSION     KX000001.1
FEATURES             Location/Qualifiers
     source          1..900
     CDS             complement(join(31..100,
                     201..260))
                     /product="X"
     CDS             401..490
                     /product="Y"
ORIGIN
//
"""


def test_flatfile_joins_are_read_from_genbank_text_and_the_cache(tmp_path):
    assert cds_joins_from_genbank(GENBANK) == [[(31, 100), (201, 260)]]
    acc = "KX000001.1"
    (tmp_path / acc).mkdir()
    (tmp_path / acc / f"{acc}.gb").write_text(GENBANK)
    gtf = tmp_path / "x.gtf"
    gtf.write_text(
        "\n".join(
            [
                row("gene", 11, 300, "G", "-", seq=acc),
                row("CDS", 31, 260, "G", "-", seq=acc, product="X"),
            ]
        )
        + "\n"
    )
    assert exons(normalise_gtf_file(gtf, tmp_path), "G") == [(11, 100), (201, 300)]
    assert exons(normalise_gtf_file(gtf, tmp_path / "empty"), "G") == [(11, 300)]


def test_overlap_groups_use_the_indexed_exons(tmp_path):
    sys.path.insert(0, str(ROOT / "extras"))
    import build_gene_programs as generator

    gtf = tmp_path / "x.gtf"
    gtf.write_text(
        "\n".join(
            cds_gene("A", 11, 300, [(31, 100), (201, 260)], stop=(261, 263))
            + cds_gene("B", 120, 180, [(125, 175)])
        )
        + "\n"
    )
    records = generator.parse_gtf(str(gtf))
    assert records["A"]["blocks"] == {(11, 100), (201, 300)}
    groups = generator.build_overlap_groups(records)
    assert groups["A"] != groups["B"]  # B sits in A's intron: no shared read sequence


REAL = ROOT / "src" / "viralscan" / "data" / "Human_herpesvirus6B.gtf"


@pytest.mark.skipif(
    not REAL.is_file() or not (Path.home() / ".cache/viralscan/ncbi/AF157706.1").is_dir(),
    reason="gitignored HHV-6B GTF or its cached flatfile is absent",
)
def test_hhv6b_collapsed_joins_come_back_from_the_flatfile():
    new = normalise_gtf_file(REAL)
    assert exons(new, "HUM_HERP6B_U66") == [(97675, 99704), (102912, 103784)]
    assert len(exons(new, "HUM_HERP6B_DR1")) == 4  # two copies, each spliced
