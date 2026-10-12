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


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(("gap", "spliced"), [(10, False), (11, True)])
def test_intron_boundary_is_eleven_bp(strand, gap, spliced):
    # MIN_INTRON_BP: a 10-bp gap is a frameshift (span kept), an 11-bp gap is an intron.
    second = (100 + gap + 1, 100 + gap + 60)
    raw = cds_gene("G", 11, 300, [(31, 100), second], strand=strand)
    expected = [(11, 100), (second[0], 300)] if spliced else [(11, 300)]
    assert exons(normalise_viral_gtf(raw), "G") == expected


@pytest.mark.parametrize("strand", ["+", "-"])
def test_multi_exon_rows_follow_the_transcript_on_both_strands(strand):
    raw = cds_gene("G", 11, 300, [(31, 100), (201, 260)], stop=(261, 263), strand=strand)
    new = normalise_viral_gtf(raw)
    numbered = [
        (int(c[3]), dict(ATTR.findall(c[8]))["exon_number"])
        for c in (line.split("\t") for line in new)
        if c[2] == "exon"
    ]
    # exon_number counts 5' to 3': descending coordinates on the minus strand
    assert numbered == ([(11, "1"), (201, "2")] if strand == "+" else [(201, "1"), (11, "2")])
    assert {line.split("\t")[6] for line in new} == {strand}
    retagged = [line for line in new if "\tCDS\t" in line or "\tstop_codon\t" in line]
    assert len(retagged) == 3
    assert all('transcript_id "G"' in line for line in retagged)


@pytest.mark.parametrize(
    "bad",
    [
        [row("gene", 300, 100, "G"), row("CDS", 120, 200, "G", protein_id="P")],  # start > end
        [row("gene", 0, 100, "G")],  # GTF is 1-based
    ],
)
def test_invalid_coordinates_raise(bad):
    with pytest.raises(ValueError, match="invalid coordinates"):
        normalise_viral_gtf(bad)


def test_gene_on_two_strands_is_skipped_with_a_warning(caplog):
    raw = [row("gene", 1, 100, "G", "+"), row("gene", 1, 100, "G", "-"), row("CDS", 5, 50, "G")]
    with caplog.at_level("WARNING", logger="viralscan"):
        new = normalise_viral_gtf(raw)
    assert new == raw  # no transcript is invented for it
    assert "gene G has rows on two sequences/strands" in caplog.text


def test_native_exon_wrap_gene_is_tagged_but_not_rewritten():
    # EBV LMP-2A: native transcript + exon rows, gene rows and transcripts with part 1/2
    wrap = [
        row("gene", 301, 400, "W", part="1"),
        row("gene", 1, 100, "W", part="2"),
        row("transcript", 301, 400, "W", transcript_id="t", part="1"),
        row("transcript", 1, 100, "W", transcript_id="t", part="2"),
        row("exon", 301, 400, "W", transcript_id="t", exon_number="1"),
        row("exon", 11, 100, "W", transcript_id="t", exon_number="2"),
    ]
    plain = [row("gene", 501, 600, "P"), row("exon", 501, 600, "P", transcript_id="p")]
    new = normalise_viral_gtf(wrap + plain)
    untagged = [line for line in new if "\ttranscript\t" not in line]
    assert untagged == [line for line in wrap + plain if "\ttranscript\t" not in line]
    tagged = [line for line in new if 'viralscan_norm "wrap_native"' in line]
    assert [line.split("\t")[2] for line in tagged] == ["transcript", "transcript"]
    assert all(line.endswith(";") for line in tagged)
    assert normalise_viral_gtf(new) == new  # idempotent: no double tag


def test_missing_flatfile_is_warned_about(tmp_path, caplog):
    gtf = tmp_path / "x.gtf"
    gtf.write_text(row("gene", 11, 300, "G", seq="KX000001.1") + "\n")
    with caplog.at_level("WARNING", logger="viralscan"):
        normalise_gtf_file(gtf, tmp_path / "empty")
    assert "no cached NCBI flatfile for KX000001.1" in caplog.text
    caplog.clear()
    (tmp_path / "KX000001.1").mkdir()
    (tmp_path / "KX000001.1" / "KX000001.1.gb").write_text(GENBANK)
    with caplog.at_level("WARNING", logger="viralscan"):
        normalise_gtf_file(gtf, tmp_path)
    assert caplog.text == ""


def test_overlap_groups_do_not_cross_contigs(tmp_path):
    sys.path.insert(0, str(ROOT / "extras"))
    import build_gene_programs as generator

    lines = [
        row("gene", 11, 300, "A", seq="contig1"),
        row("exon", 31, 260, "A", seq="contig1"),
        row("gene", 11, 300, "B", seq="contig2"),
        row("exon", 31, 260, "B", seq="contig2"),
    ]
    gtf = tmp_path / "x.gtf"
    gtf.write_text("\n".join(lines) + "\n")
    groups = generator.build_overlap_groups(generator.parse_gtf(str(gtf)))
    assert groups["A"] != groups["B"]  # different sequences cannot share reads


def test_docstring_figures_are_the_ones_measured_for_aud13_a12():
    # A-12: the docstring said 2,517 genes and ~102 kb over ~102 genes; the audit measured these.
    doc = gtf_normalise.__doc__ or ""
    for stale in ("2,517", "~102 kb", "~102 genes"):
        assert stale not in doc
    for measured in ("2,522", "101,346 bp", "105 spliced", "103 genes"):
        assert measured in doc


@pytest.mark.research
def test_docstring_figures_match_the_local_gtf_corpus():
    data = ROOT / "src" / "viralscan" / "data"
    cache = Path.home() / ".cache/viralscan/ncbi"
    gtfs = sorted(data.glob("*.gtf"))
    if len(gtfs) != 393 or not cache.is_dir():
        pytest.skip("the 393-GTF local corpus and the NCBI flatfile cache are required")
    no_exon = spliced = intron = 0
    spliced_genes = set()
    for path in gtfs:
        raw = [ln.split("\t") for ln in path.read_text().splitlines() if ln and ln[0] != "#"]
        genes = {(c[0], dict(ATTR.findall(c[8])).get("gene_id")) for c in raw if c[2] == "gene"}
        with_exon = {(c[0], dict(ATTR.findall(c[8])).get("gene_id")) for c in raw if c[2] == "exon"}
        no_exon += len(genes - with_exon)
        exons: dict[str, list[tuple[int, int]]] = {}
        spliced_tx = {}
        for line in normalise_gtf_file(path, cache):
            c = line.split("\t")
            if len(c) != 9 or c[2] not in ("transcript", "exon"):
                continue
            attrs = dict(ATTR.findall(c[8]))
            if c[2] == "exon":
                exons.setdefault(attrs["transcript_id"], []).append((int(c[3]), int(c[4])))
            elif "spliced" in attrs.get("viralscan_norm", ""):
                spliced_tx[attrs["transcript_id"]] = attrs["gene_id"]
        for tx, gene in spliced_tx.items():
            blocks = sorted(exons[tx])
            spliced += 1
            spliced_genes.add((path.name, gene))
            intron += sum(b[0] - a[1] - 1 for a, b in zip(blocks, blocks[1:]))
    assert (no_exon, spliced, len(spliced_genes), intron) == (2522, 105, 103, 101346)
