#!/usr/bin/env python3
"""Freeze a validated prefix of an existing paired FASTQ for a bounded pilot."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
from pathlib import Path

from audit_fastq_pair import FastqAuditError, _open_fastq, _pair_id, _record, audit_pair

from viralscan.run_safety import sha256_file


def subset(read1: Path, read2: Path, output: Path, pairs: int) -> dict:
    if pairs <= 0 or output.exists():
        raise ValueError("Pair limit must be positive and output must be a new directory")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        destinations = [stage / "R1.fastq.gz", stage / "R2.fastq.gz"]
        with _open_fastq(read1) as left, _open_fastq(read2) as right:
            with (
                gzip.open(destinations[0], "wb", compresslevel=1) as out1,
                gzip.open(destinations[1], "wb", compresslevel=1) as out2,
            ):
                for number in range(1, pairs + 1):
                    records = [_record(left, read1, number), _record(right, read2, number)]
                    if any(record is None for record in records):
                        raise FastqAuditError(f"Input ended before requested {pairs} pairs")
                    if _pair_id(records[0][0]) != _pair_id(records[1][0]):
                        raise FastqAuditError(f"Paired identifiers differ at record {number}")
                    out1.writelines(records[0])
                    out2.writelines(records[1])
        audit = audit_pair(*destinations, sample_id=output.name)
        receipt = {
            "label": "explore",
            "sampling": "first intact read pairs, not random or virus-selected",
            "pairs": pairs,
            "script_sha256": sha256_file(Path(__file__)),
            "sources": [
                {
                    "path": str(p.resolve()),
                    "stored_bytes": p.stat().st_size,
                    "mtime_ns": p.stat().st_mtime_ns,
                }
                for p in (read1, read2)
            ],
            "subset_sha256": {p.name: sha256_file(p) for p in destinations},
            "audit": audit,
            "limitation": "Whole-source gzip hashes are not recomputed; copied pairs are validated and hashed",
        }
        (stage / "subset.json").write_text(json.dumps(receipt, indent=2) + "\n")
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pairs", type=int, default=1_000_000)
    args = parser.parse_args()
    print(json.dumps(subset(args.read1, args.read2, args.output, args.pairs), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
