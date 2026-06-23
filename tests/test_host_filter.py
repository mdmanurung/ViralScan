"""Unit tests for host_filter CB/UMI geometry and FASTQ-pair filtering.

Tests run without Snakemake, minimap2, STARsolo, or kallisto — all tool
I/O is simulated with synthetic temp FASTQ files or passed as in-memory data.

These are regression guards for PLAN S1/S6:
- host_filter._TECH_PARAMS was 10x-only and fell back to (16,12) for unknown
  techs, silently mis-slicing Drop-seq and other non-10x barcodes.
- cb_umi_geometry() now raises ValueError for unknown techs (fails loudly).
- filter_fastq_pairs() takes all params explicitly — testable without Snakemake.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from viralscan.scripts.host_filter import filter_fastq_pairs
from viralscan.evidence import cb_umi_geometry
from tests._fastq import write_fastq


# cb_umi_geometry() tests (named chemistries, explicit triplets, unknown-tech error)
# live in tests/test_evidence.py::TestGeometry.  The FASTQ-filter tests below
# exercise the actual geometry use-path through filter_fastq_pairs().


class TestFilterFastqPairs:
    """filter_fastq_pairs keeps pairs whose (CB, UMI) is NOT in host_mapped."""

    def test_keeps_unmapped_drops_mapped(self, tmp_path: Path) -> None:
        # Drop-seq geometry: CB=12 bases, UMI=8 bases
        cb = "A" * 12
        umi_host = "C" * 8   # mapped to host — should be dropped
        umi_virus = "G" * 8  # not host-mapped — should be kept

        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        out_r1 = str(tmp_path / "out_R1.fastq.gz")
        out_r2 = str(tmp_path / "out_R2.fastq.gz")

        write_fastq(
            r1,
            [
                ("pair1", cb + umi_host),  # host-mapped → drop
                ("pair2", cb + umi_virus),  # not host-mapped → keep
                ("pair3", cb + umi_host),  # host-mapped → drop
            ],
        )
        write_fastq(
            r2,
            [
                ("pair1", "ACGT" * 5),
                ("pair2", "TTTT" * 5),
                ("pair3", "GGGG" * 5),
            ],
        )

        cb_len, umi_len = cb_umi_geometry("dropseq")
        host_mapped = {(cb, umi_host)}
        kept, total = filter_fastq_pairs(
            str(r1), str(r2), out_r1, out_r2, cb_len, umi_len, host_mapped
        )

        assert total == 3
        assert kept == 1

        # Verify only pair2's cDNA sequence is in the output
        with gzip.open(out_r2, "rt") as fh:
            content = fh.read()
        assert "TTTT" * 5 in content
        assert "ACGT" * 5 not in content
        assert "GGGG" * 5 not in content

    def test_empty_host_set_keeps_all(self, tmp_path: Path) -> None:
        cb = "A" * 16
        umi = "C" * 12
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        out_r1 = str(tmp_path / "out_R1.fastq.gz")
        out_r2 = str(tmp_path / "out_R2.fastq.gz")

        write_fastq(r1, [("r1", cb + umi), ("r2", cb + umi)])
        write_fastq(r2, [("r1", "ACGT"), ("r2", "TGCA")])

        kept, total = filter_fastq_pairs(
            str(r1), str(r2), out_r1, out_r2, 16, 12, set()
        )
        assert kept == 2
        assert total == 2

    def test_all_host_mapped_drops_all(self, tmp_path: Path) -> None:
        cb = "G" * 16
        umi = "T" * 12
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        out_r1 = str(tmp_path / "out_R1.fastq.gz")
        out_r2 = str(tmp_path / "out_R2.fastq.gz")

        write_fastq(r1, [("r1", cb + umi)])
        write_fastq(r2, [("r1", "AAAA")])

        kept, total = filter_fastq_pairs(
            str(r1), str(r2), out_r1, out_r2, 16, 12, {(cb, umi)}
        )
        assert kept == 0
        assert total == 1

    def test_gzipped_input(self, tmp_path: Path) -> None:
        cb = "C" * 12
        umi_keep = "A" * 8
        umi_drop = "T" * 8

        r1 = tmp_path / "R1.fastq.gz"
        r2 = tmp_path / "R2.fastq.gz"
        out_r1 = str(tmp_path / "out_R1.fastq.gz")
        out_r2 = str(tmp_path / "out_R2.fastq.gz")

        write_fastq(r1, [("a", cb + umi_keep), ("b", cb + umi_drop)])
        write_fastq(r2, [("a", "CCCCCCCC"), ("b", "GGGGGGGG")])

        kept, total = filter_fastq_pairs(
            str(r1), str(r2), out_r1, out_r2, 12, 8, {(cb, umi_drop)}
        )
        assert kept == 1
        assert total == 2

        with gzip.open(out_r2, "rt") as fh:
            out_content = fh.read()
        assert "CCCCCCCC" in out_content
        assert "GGGGGGGG" not in out_content

    def test_returns_counts_tuple(self, tmp_path: Path) -> None:
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(r1, [("x", "A" * 20)])
        write_fastq(r2, [("x", "T" * 20)])
        result = filter_fastq_pairs(
            str(r1), str(r2),
            str(tmp_path / "o1.fastq.gz"),
            str(tmp_path / "o2.fastq.gz"),
            16, 12, set()
        )
        kept, total = result
        assert isinstance(kept, int)
        assert isinstance(total, int)

    def test_truncated_r1_raises(self, tmp_path: Path) -> None:
        # A file that ends after the @header line (mid-record) must raise, not
        # silently emit a malformed record with empty sequence/quality lines.
        # R2 gets 2 complete records so the R2 EOF check doesn't fire first.
        r1 = tmp_path / "R1.fastq"
        r2 = tmp_path / "R2.fastq"
        write_fastq(r1, [("good", "A" * 20)])
        write_fastq(r2, [("good", "T" * 20), ("extra", "G" * 20)])
        # Append a truncated record (header only, no seq/qual) to R1
        with open(r1, "a") as fh:
            fh.write("@truncated\n")
        with pytest.raises(ValueError, match="Truncated FASTQ"):
            filter_fastq_pairs(
                str(r1), str(r2),
                str(tmp_path / "o1.fastq.gz"),
                str(tmp_path / "o2.fastq.gz"),
                16, 12, set(),
            )
