"""Tests for viral read-evidence tracing and extraction (pure logic)."""

from __future__ import annotations

import gzip

import pytest

from viralscan.evidence import (
    _run,
    cb_umi_geometry,
    extract_viral_reads,
    have_tools,
    viral_assigned_keys,
    viral_equivalence_classes,
)


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
        assert cb_umi_geometry("DROPSEQ") == (12, 8)  # case-insensitive, non-10x supported

    def test_explicit_geometry_string(self) -> None:
        assert cb_umi_geometry("0,0,16:0,16,28:1,0,0") == (16, 12)

    def test_unknown_raises(self) -> None:
        with pytest.raises(ValueError):
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


def _write_fastq(path, records):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


class TestExtractReads:
    def test_extracts_only_viral_assigned_cdna_reads(self, tmp_path) -> None:
        # 10xv2 geometry: CB=16, UMI=10 -> first 26 bases of R1
        cb = "A" * 16
        umi_hit = "C" * 10
        umi_miss = "G" * 10
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        _write_fastq(
            r1,
            [("r1", cb + umi_hit), ("r2", cb + umi_miss), ("r3", cb + umi_hit)],
        )
        _write_fastq(
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
        _write_fastq(r1, [("r1", cb + umi)])
        _write_fastq(r2, [("r1", "ACGTACGTACGT")])
        out = tmp_path / "ev.fasta"
        stats = extract_viral_reads(
            str(r1), str(r2), keys={(cb, umi)}, cb_len=16, umi_len=12, out_fasta=str(out)
        )
        assert stats.viral_reads == 1
        assert "ACGTACGTACGT" in out.read_text()
