"""The one FASTA reader and the one GTF reader (PLAN SEQIO-01)."""

import gzip

import pytest

from viralscan.fasta_reader import FastaError, iter_fasta, read_fasta
from viralscan.gtf_reader import (
    gene_id_of,
    gene_ids,
    gtf_gene_ids,
    iter_gtf_fields,
    parse_attributes,
)

ROW = "c1\ts\tgene\t1\t10\t.\t+\t.\t{}\n"


def _write(tmp_path, name, text):
    path = tmp_path / name
    if name.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("name", ["x.fa", "x.fa.gz"])
def test_fasta_plain_and_gz_read_alike(tmp_path, name):
    path = _write(tmp_path, name, ">a first desc\nACGT\nacgt\n>b\nGG\n")
    assert list(iter_fasta(path)) == [("a", "ACGTacgt"), ("b", "GG")]


def test_fasta_blank_lines_crlf_and_missing_final_newline(tmp_path):
    path = tmp_path / "x.fa"
    path.write_bytes(b">a\r\n\r\nAC\r\nGT\r\n\r\n>b\tdesc\r\nTT")
    assert read_fasta(path) == {"a": "ACGT", "b": "TT"}


def test_fasta_upper_and_empty_record(tmp_path):
    path = _write(tmp_path, "x.fa", ">a\nacgtn\n>b\n>c\nA\n")
    assert read_fasta(path, upper=True) == {"a": "ACGTN", "b": "", "c": "A"}
    assert read_fasta(path)["a"] == "acgtn"


def test_fasta_text_before_first_header_is_ignored(tmp_path):
    assert read_fasta(_write(tmp_path, "x.fa", "AAAA\n>a\nAC\n")) == {"a": "AC"}


def test_fasta_empty_file_has_no_records(tmp_path):
    assert read_fasta(_write(tmp_path, "x.fa", "")) == {}


def test_fasta_repeated_id_is_the_callers_rule(tmp_path):
    path = _write(tmp_path, "x.fa", ">a\nAC\n>a\nGGG\n>b\nT\n")
    assert list(iter_fasta(path)) == [("a", "AC"), ("a", "GGG"), ("b", "T")]
    assert read_fasta(path) == {"a": "GGG", "b": "T"}
    with pytest.raises(FastaError, match="Duplicate FASTA ID 'a'"):
        read_fasta(path, unique=True)


@pytest.mark.parametrize("header", [">", ">   "])
def test_fasta_header_without_identifier_raises(tmp_path, header):
    with pytest.raises(FastaError, match="without an identifier"):
        read_fasta(_write(tmp_path, "x.fa", f"{header}\nAC\n"))


def test_gtf_fields_skip_comments_blank_and_short_rows():
    lines = [
        "#comment\n",
        "\n",
        'c1\ts\tgene\t1\t10\t.\t+\tgene_id "G1";\n',  # eight columns
        ROW.format('gene_id "G1";'),
        ROW.format('gene_id "G2";').rstrip("\n") + "\t\n",  # a trailing tab: ten columns
    ]
    assert [f[8] for f in iter_gtf_fields(lines)] == ['gene_id "G1";', 'gene_id "G2";']
    assert [f[8] for f in iter_gtf_fields(lines, exact=True)] == ['gene_id "G1";']


def test_gtf_attributes_quoted_unquoted_and_missing():
    assert parse_attributes('gene_id "G1"; note "a b"; n 3; bare x;') == {
        "gene_id": "G1",
        "note": "a b",
    }
    assert parse_attributes("gene_id G1;") == {}
    assert parse_attributes("") == {}


def test_gtf_attribute_rules_differ_on_repeats_and_prefixes():
    """Two rules, both kept: the dict parse and the ``gene_id`` search disagree on these inputs."""
    text = 'ref_gene_id "R"; gene_id "G1"; gene_id "G2";'
    assert parse_attributes(text)["gene_id"] == "G2"
    assert gene_id_of(text) == "R"
    assert gene_id_of('gene_id "";') is None
    assert parse_attributes('gene_id "";') == {"gene_id": ""}
    assert gene_id_of('gene_id  "G1";') is None
    assert parse_attributes('gene_id  "G1";') == {"gene_id": "G1"}


def test_gtf_gene_ids_plain_gz_and_lines(tmp_path):
    text = (
        "#c\n"
        + ROW.format('gene_id "A"; gene_name "x";')
        + ROW.format('gene_name "y";')
        + ROW.format('gene_id "B";')
        + ROW.format('gene_id "A";')
    )
    assert gene_ids(text.splitlines(keepends=True)) == {"A", "B"}
    assert gtf_gene_ids(_write(tmp_path, "x.gtf", text)) == {"A", "B"}
    assert gtf_gene_ids(_write(tmp_path, "x.gtf.gz", text)) == {"A", "B"}
