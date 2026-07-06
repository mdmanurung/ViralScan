#!/usr/bin/env python3
"""Validate and summarize the reference-strategy benchmark TSV."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from viralscan.reference_strategy import (
    BenchmarkContractError,
    validate_results,
    write_parsed_benchmark_results,
    write_run_status_from_results,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, help="Benchmark run directory used as provenance input")
    parser.add_argument("--out", type=Path, help="Output benchmark TSV path")
    parser.add_argument("--validate", type=Path, required=True, help="results/reference_strategy_benchmark.tsv")
    parser.add_argument("--run-status-out", type=Path)
    args = parser.parse_args(argv)

    try:
        if args.run_dir and args.out and args.out != args.validate:
            raise BenchmarkContractError("--out and --validate must point to the same TSV in this implementation")
        if args.run_dir:
            write_parsed_benchmark_results(args.run_dir, args.validate)
        rows = validate_results(args.validate)
        if args.run_status_out:
            write_run_status_from_results(args.validate, args.run_status_out)
    except (OSError, BenchmarkContractError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    complete = sum(1 for row in rows if row["status"] == "complete")
    failed = [row for row in rows if row["status"] != "complete"]
    print(f"Validated {len(rows)} rows ({complete} complete, {len(failed)} incomplete/failed)")
    for row in failed:
        print(
            f"{row['dataset']}\t{row['method']}\t{row['reference_strategy']}\t{row['status']}\t{row['failure_reason']}"
        )
    return 0 if complete == 12 else 1


if __name__ == "__main__":
    raise SystemExit(main())
