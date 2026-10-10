"""Pure-function tests for scripts/fix_thin_gtfs.py and scripts/model_hpv_cds_gtfs.py (PLAN CAT-38/39)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fix_thin_gtfs as fix  # noqa: E402
import model_hpv_cds_gtfs as model  # noqa: E402


def row(feature, start, end, gene, tx="", gbkey="CDS"):
    attrs = f'gene_id "{gene}"; transcript_id "{tx}"; gbkey "{gbkey}";'
    return "\t".join(["S.1", "RefSeq", feature, str(start), str(end), ".", "+", ".", attrs])


def test_spanning_transcript_is_dropped_and_gene_is_respanned():
    lines = [
        row("gene", 1, 1000, "G"),
        row("transcript", 1, 1000, "G", "T1", "misc_RNA"),
        row("exon", 1, 1000, "G", "T1", "misc_RNA"),
        row("transcript", 100, 400, "G", "T2", "CDS"),
        row("exon", 100, 400, "G", "T2", "CDS"),
        row("CDS", 100, 400, "G", "T2", "CDS"),
    ]
    out = fix.drop_spanning(lines, genome_length=1000)
    assert [c.split("\t")[2] for c in out] == ["gene", "transcript", "exon", "CDS"]
    gene = out[0].split("\t")
    assert (gene[3], gene[4]) == ("100", "400")
    assert fix.drop_spanning(out, 1000) == out  # idempotent


def test_gene_with_only_a_spanning_transcript_disappears():
    lines = [
        row("gene", 1, 1000, "HAVgs1"),
        row("transcript", 1, 1000, "HAVgs1", "T1", "misc_RNA"),
        row("exon", 1, 1000, "HAVgs1", "T1", "misc_RNA"),
        row("gene", 50, 500, "real"),
        row("transcript", 50, 500, "real", "T2", "CDS"),
        row("exon", 50, 500, "real", "T2", "CDS"),
    ]
    out = fix.drop_spanning(lines, 1000)
    assert all('gene_id "HAVgs1"' not in line for line in out)
    assert len(out) == 3


def test_long_cds_is_never_dropped():
    """Only misc_RNA / prim_transcript are removed: a genome-length polyprotein CDS stays."""
    lines = [
        row("gene", 10, 990, "poly"),
        row("transcript", 10, 990, "poly", "T1", "CDS"),
        row("exon", 10, 990, "poly", "T1", "CDS"),
    ]
    assert fix.drop_spanning(lines, 1000) == lines


def test_translate_and_orf_window():
    genome = "TAA" + "CC" + "ATG" + "AAA" + "GGG" + "TGA" + "AAAAAA"
    assert model.translate("ATGAAATGA") == "MK*"
    start, end = model.orf_around(genome, hit_start=9, hit_end=14)
    assert (start, end) == (6, 17)  # ATG at 6, stop codon TGA ends at 17
    assert model.translate(genome[start - 1 : end]) == "MKG*"


def test_orf_returns_none_without_stop_or_atg():
    assert model.orf_around("CC" + "AAA" * 10, 3, 30) is None  # runs off the end, no stop
    assert model.orf_around("TAA" + "AAA" * 4 + "TGA", 4, 15) is None  # window has no ATG


def test_model_gene_flags_length_mismatch_and_minus_strand():
    genome = "TAA" + "CC" + "ATG" + "AAA" * 30 + "TGA" + "GG"
    hsp = ["G", "6", str(5 + 3 * 31), "70.0", "31", "1", "31", "1e-30", "99"]
    ok = model.model_gene("G", "M" + "K" * 30, hsp, genome)
    assert ok["status"] == "modelled" and ok["aa"] == 31
    short = model.model_gene("G", "M" + "K" * 300, hsp, genome)
    assert short["status"] == "length_mismatch"
    minus = model.model_gene(
        "G", "MK", ["G", "50", "10", "60", "3", "1", "3", "1e-5", "20"], genome
    )
    assert minus["status"] == "minus_strand_hit"
    assert model.model_gene("G", "MK", None, genome)["status"] == "no_hit"


def test_sibling_proteins_skip_spliced_products(tmp_path):
    genome = "ATG" + "AAA" * 4 + "TGA" + "ATG" + "CCC" + "TAA" + "ATG" + "GGG" + "TAA"
    gtf = tmp_path / "s.gtf"
    gtf.write_text(
        "\n".join(
            [
                # E6: one CDS row -> kept; E7: exon-only (older RefSeq style) -> kept
                'S\tR\tCDS\t1\t18\t.\t+\t0\tgene "E6"; transcript_id "t1";',
                'S\tR\texon\t19\t27\t.\t+\t.\tgene "E7"; transcript_id "t2";',
                # E1^E4-like: two CDS rows in one transcript -> skipped
                'S\tR\tCDS\t28\t30\t.\t+\t0\tgene "E1"; transcript_id "t3";',
                'S\tR\tCDS\t31\t33\t.\t+\t0\tgene "E1"; transcript_id "t3";',
            ]
        )
        + "\n"
    )
    proteins = model.sibling_proteins(gtf, genome)
    assert proteins == {"E6": "MKKKK", "E7": "MP"}
