#!/usr/bin/env python3
"""Stream, hash, and structurally validate one paired FASTQ input."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import BinaryIO


class FastqAuditError(RuntimeError):
    """Raised when a FASTQ pair violates the frozen input contract."""


REQUIRED_R1_BASES = {"10xv2": 26, "10xv3": 28}


def _record(handle: BinaryIO, path: Path, record_number: int) -> tuple[bytes, ...] | None:
    lines = tuple(handle.readline() for _ in range(4))
    if not any(lines):
        return None
    if any(not line for line in lines):
        raise FastqAuditError(f"{path}: truncated record {record_number}")
    header, sequence, plus, quality = lines
    if not header.startswith(b"@"):
        raise FastqAuditError(f"{path}: invalid header at record {record_number}")
    if not plus.startswith(b"+"):
        raise FastqAuditError(f"{path}: invalid plus line at record {record_number}")
    if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
        raise FastqAuditError(f"{path}: sequence/quality length mismatch at record {record_number}")
    return lines


def _pair_id(header: bytes) -> bytes:
    identifier = header[1:].split(maxsplit=1)[0]
    if identifier.endswith((b"/1", b"/2")):
        identifier = identifier[:-2]
    return identifier


def _open_fastq(path: Path) -> BinaryIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rb")
    return path.open("rb")


def audit_pair(
    read1: Path,
    read2: Path,
    sample_id: str,
    *,
    chemistry: str | None = None,
) -> dict[str, object]:
    """Read a pair once while computing hashes and validating record pairing."""

    if chemistry is not None and chemistry not in REQUIRED_R1_BASES:
        raise FastqAuditError(f"unsupported chemistry: {chemistry}")

    hashes = {
        "read1": {"sha256": hashlib.sha256(), "md5": hashlib.md5()},
        "read2": {"sha256": hashlib.sha256(), "md5": hashlib.md5()},
    }
    byte_counts = {"read1": 0, "read2": 0}
    sequence_min: dict[str, int | None] = {"read1": None, "read2": None}
    sequence_max = {"read1": 0, "read2": 0}
    record_number = 0
    with _open_fastq(read1) as left, _open_fastq(read2) as right:
        while True:
            left_record = _record(left, read1, record_number + 1)
            right_record = _record(right, read2, record_number + 1)
            if left_record is None and right_record is None:
                break
            if left_record is None or right_record is None:
                raise FastqAuditError("paired FASTQs have different record counts")
            record_number += 1
            if _pair_id(left_record[0]) != _pair_id(right_record[0]):
                raise FastqAuditError(f"paired identifiers differ at record {record_number}")
            left_sequence_bases = len(left_record[1].rstrip(b"\r\n"))
            if chemistry is not None and left_sequence_bases < REQUIRED_R1_BASES[chemistry]:
                raise FastqAuditError(
                    f"{chemistry} requires at least {REQUIRED_R1_BASES[chemistry]} "
                    f"R1 bases; found {left_sequence_bases} at record {record_number}"
                )
            for key, record in (("read1", left_record), ("read2", right_record)):
                sequence_bases = len(record[1].rstrip(b"\r\n"))
                if sequence_min[key] is None:
                    sequence_min[key] = sequence_bases
                else:
                    sequence_min[key] = min(sequence_min[key], sequence_bases)
                sequence_max[key] = max(sequence_max[key], sequence_bases)
                for line in record:
                    byte_counts[key] += len(line)
                    hashes[key]["sha256"].update(line)
                    hashes[key]["md5"].update(line)

    def read_result(key: str, path: Path) -> dict[str, object]:
        return {
            "bytes": byte_counts[key],
            "stored_bytes": path.stat().st_size,
            "compression": "gzip" if path.name.endswith(".gz") else "none",
            "min_sequence_bases": sequence_min[key] or 0,
            "max_sequence_bases": sequence_max[key],
            "sha256": hashes[key]["sha256"].hexdigest(),
            "md5": hashes[key]["md5"].hexdigest(),
        }

    return {
        "schema_version": "1.0.0",
        "sample_id": sample_id,
        "chemistry": chemistry,
        "required_r1_bases": REQUIRED_R1_BASES.get(chemistry),
        "status": "valid",
        "records": record_number,
        "read1": read_result("read1", read1),
        "read2": read_result("read2", read2),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--chemistry", choices=("10xv2", "10xv3"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit_pair(
        args.read1,
        args.read2,
        args.sample_id,
        chemistry=args.chemistry,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
