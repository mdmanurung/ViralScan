"""Tests for viral read-evidence tracing and extraction (pure logic)."""

from __future__ import annotations

import pytest

from viralscan.evidence import (
    _parse_blast_output,
    _parse_coverage_output,
    _run,
    cb_umi_geometry,
    extract_viral_reads,
    have_tools,
    viral_assigned_keys,
    viral_equivalence_classes,
)
from tests._fastq import write_fastq


class TestRunSurfacesErrors:
    def test_failing_command_raises_with_stderr(self) -> None:
        # The fail-fast contract: a non-zero exit raises RuntimeError that
        # INCLUDES the tool's stderr, never a bare swallowed CalledProcessError.
        with pytest.raises(RuntimeError) as exc:
            _run(["sh", "-c", "echo boom-message 1>&2; exit 3"])
        assert "exit 3" in str(exc.value)
        assert "boom-message" in str(exc.value)

    def test_capture_returns_stdout(self) -> None:
        assert _run(["printf", "hello"], capture=True) == b"hello"

    def test_have_tools_reports_missing(self) -> None:
        assert have_tools(["definitely-not-a-real-binary-xyz"]) == [
            "definitely-not-a-real-binary-xyz"
        ]
        assert have_tools(["sh"]) == []


class TestGeometry:
    def test_named_chemistries(self) -> None:
        assert cb_umi_geometry("10xv2") == (16, 10)
        assert cb_umi_geometry("10xv3") == (16, 12)
        # Drop-seq is case-insensitive and returns the correct (12, 8) geometry.
        # Regression guard: the old _cb_umi_lengths() silently returned (16,12)
        # for any non-10x chemistry.
        assert cb_umi_geometry("DROPSEQ") == (12, 8)
        assert cb_umi_geometry("dropseq") == (12, 8)
        assert cb_umi_geometry("DropSeq") == (12, 8)

    def test_explicit_geometry_string(self) -> None:
        # Explicit kallisto "bc:umi:seq" triplets bypass the named-chemistry table.
        assert cb_umi_geometry("0,0,16:0,16,28:1,0,0") == (16, 12)  # 10xv3 layout
        assert cb_umi_geometry("0,0,12:0,12,20:1,0,0") == (12, 8)   # Drop-seq layout

    def test_unknown_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown technology"):
            cb_umi_geometry("nanopore-bonkers")


class TestViralEcs:
    def test_selects_ecs_touching_viral_genes(self) -> None:
        ec_map = {0: [0], 1: [1], 2: [0, 1], 3: [3]}
        viral = {1}
        assert viral_equivalence_classes(ec_map, viral) == {1, 2}


class TestViralKeys:
    def test_collects_barcode_umi_for_viral_ecs(self) -> None:
        lines = [
            "AAAA\tUUUU1\t2\t1",  # viral EC -> kept
            "AAAA\tUUUU2\t0\t1",  # host EC  -> skipped
            "CCCC\tUUUU3\t2\t3",  # viral EC -> kept
            "garbage line",
        ]
        keys = viral_assigned_keys(lines, viral_ecs={2})
        assert keys == {("AAAA", "UUUU1"), ("CCCC", "UUUU3")}


class TestParseCoverageOutput:
    """_parse_coverage_output parses samtools coverage TSV without the binary."""

    _HEADER = "#rname\tstart\tend\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
    _ROW_HIT = "virus_A\t0\t29903\t412\t27100\t90.6\t4.1\t37\t60"
    _ROW_ZERO = "virus_B\t0\t2000\t0\t0\t0.0\t0.0\t0\t0"

    def test_parses_covered_references(self) -> None:
        text = "\n".join([self._HEADER, self._ROW_HIT, self._ROW_ZERO])
        rows = _parse_coverage_output(text)
        assert len(rows) == 1
        assert rows[0]["rname"] == "virus_A"
        assert rows[0]["numreads"] == "412"

    def test_excludes_zero_read_references(self) -> None:
        text = "\n".join([self._HEADER, self._ROW_ZERO])
        assert _parse_coverage_output(text) == []

    def test_empty_output(self) -> None:
        assert _parse_coverage_output("") == []

    def test_header_hash_stripped(self) -> None:
        text = "\n".join([self._HEADER, self._ROW_HIT])
        rows = _parse_coverage_output(text)
        # Column key must be "rname" not "#rname"
        assert "rname" in rows[0]
        assert "#rname" not in rows[0]


class TestParseBlastOutput:
    """_parse_blast_output parses blastn -outfmt 6 text without the binary."""

    _BLAST_TEXT = (
        "read1\tNC_045512.2\t98.5\t150\t1e-80\n"
        "read1\tNC_002549.1\t72.0\t130\t1e-20\n"  # second hit for read1 — ignored
        "read2\tNC_045512.2\t95.0\t148\t1e-75\n"
        "read3\tNC_002549.1\t60.0\t100\t1e-10\n"
    )

    def test_best_hit_per_read(self) -> None:
        rows = _parse_blast_output(self._BLAST_TEXT)
        # read1 appears twice; only the first (best) hit must be returned
        read_ids = [r["read"] for r in rows]
        assert read_ids.count("read1") == 1

    def test_all_reads_present(self) -> None:
        rows = _parse_blast_output(self._BLAST_TEXT)
        assert {r["read"] for r in rows} == {"read1", "read2", "read3"}

    def test_fields_populated(self) -> None:
        rows = _parse_blast_output(self._BLAST_TEXT)
        r1 = next(r for r in rows if r["read"] == "read1")
        assert r1["subject"] == "NC_045512.2"
        assert r1["pident"] == "98.5"
        assert r1["length"] == "150"

    def test_empty_output(self) -> None:
        assert _parse_blast_output("") == []

    def test_short_lines_skipped(self) -> None:
        # Lines with fewer than 4 columns are silently ignored
        text = "only\ttwo\n" + "read1\tNC_045512.2\t98.5\t150\t1e-80\n"
        rows = _parse_blast_output(text)
        assert len(rows) == 1
        assert rows[0]["read"] == "read1"


class TestExtractReads:
    def test_extracts_only_viral_assigned_cdna_reads(self, tmp_path) -> None:
        # 10xv2 geometry: CB=16, UMI=10 -> first 26 bases of R1
        cb = "A" * 16
        umi_hit = "C" * 10
        umi_miss = "G" * 10
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(
            r1,
            [("r1", cb + umi_hit), ("r2", cb + umi_miss), ("r3", cb + umi_hit)],
        )
        write_fastq(
            r2,
            [("r1", "ACGTACGTAC"), ("r2", "TTTTTTTTTT"), ("r3", "GGGGCCCCGG")],
        )
        out = tmp_path / "evidence.fasta"
        stats = extract_viral_reads(
            str(r1), str(r2), keys={(cb, umi_hit)}, cb_len=16, umi_len=10, out_fasta=str(out)
        )
        assert stats.total_reads == 3
        assert stats.viral_reads == 2  # r1 and r3 (umi_hit), not r2 (umi_miss)
        body = out.read_text()
        assert "ACGTACGTAC" in body and "GGGGCCCCGG" in body
        assert "TTTTTTTTTT" not in body

    def test_handles_gzipped_input(self, tmp_path) -> None:
        cb, umi = "A" * 16, "C" * 12  # 10xv3
        r1 = tmp_path / "R1.fastq.gz"
        r2 = tmp_path / "R2.fastq.gz"
        write_fastq(r1, [("r1", cb + umi)])
        write_fastq(r2, [("r1", "ACGTACGTACGT")])
        out = tmp_path / "ev.fasta"
        stats = extract_viral_reads(
            str(r1), str(r2), keys={(cb, umi)}, cb_len=16, umi_len=12, out_fasta=str(out)
        )
        assert stats.viral_reads == 1
        assert "ACGTACGTACGT" in out.read_text()
