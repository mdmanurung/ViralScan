import gzip
from pathlib import Path

import pytest

from scripts.audit_fastq_pair import FastqAuditError, audit_pair


def test_audit_fastq_pair_hashes_and_validates(tmp_path: Path) -> None:
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    read1.write_bytes(b"@spot/1\nAC\n+\nII\n@next/1\nGT\n+\nJJ\n")
    read2.write_bytes(b"@spot/2\nTG\n+\nII\n@next/2\nCA\n+\nJJ\n")
    result = audit_pair(read1, read2, "sample")
    assert result["status"] == "valid"
    assert result["records"] == 2
    assert result["read1"]["bytes"] == read1.stat().st_size
    assert len(result["read2"]["sha256"]) == 64
    assert len(result["read2"]["md5"]) == 32


def test_audit_fastq_pair_rejects_mismatched_ids(tmp_path: Path) -> None:
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    read1.write_bytes(b"@spot/1\nAC\n+\nII\n")
    read2.write_bytes(b"@other/2\nTG\n+\nII\n")
    with pytest.raises(FastqAuditError, match="identifiers differ"):
        audit_pair(read1, read2, "sample")


def test_audit_fastq_pair_reads_gzip_and_validates_10xv2_geometry(
    tmp_path: Path,
) -> None:
    read1 = tmp_path / "r1.fastq.gz"
    read2 = tmp_path / "r2.fastq.gz"
    with gzip.open(read1, "wb") as handle:
        handle.write(b"@spot/1\n" + b"A" * 26 + b"\n+\n" + b"I" * 26 + b"\n")
    with gzip.open(read2, "wb") as handle:
        handle.write(b"@spot/2\nACGT\n+\nIIII\n")

    result = audit_pair(read1, read2, "sample", chemistry="10xv2")

    assert result["status"] == "valid"
    assert result["records"] == 1
    assert result["chemistry"] == "10xv2"
    assert result["read1"]["compression"] == "gzip"
    assert result["read1"]["min_sequence_bases"] == 26
    assert result["read1"]["stored_bytes"] == read1.stat().st_size


def test_audit_fastq_pair_rejects_short_barcode_umi_read(tmp_path: Path) -> None:
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    read1.write_bytes(b"@spot/1\n" + b"A" * 25 + b"\n+\n" + b"I" * 25 + b"\n")
    read2.write_bytes(b"@spot/2\nACGT\n+\nIIII\n")

    with pytest.raises(FastqAuditError, match="10xv2 requires at least 26"):
        audit_pair(read1, read2, "sample", chemistry="10xv2")
