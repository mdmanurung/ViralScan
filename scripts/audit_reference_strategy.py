#!/usr/bin/env python3
"""Audit the matched-reference benchmark manifest and write reference_audit.tsv."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from viralscan.reference_strategy import (
    BenchmarkContractError,
    PANEL_ID,
    audit_manifest,
    audit_fastqs,
    default_manifest,
    load_manifest,
    validate_manifest,
    validate_audit_rows,
    validate_fastq_audit_rows,
    write_fastq_audit,
    write_json,
    write_reference_audit,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--expected-panel", default=PANEL_ID)
    parser.add_argument("--fail-on-single-virus", action="store_true")
    parser.add_argument("--no-sha256", action="store_true", help="Skip file hashing for a fast structural audit")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Write audit evidence but do not fail on missing paths/provenance/target checks.",
    )
    parser.add_argument("--init-template", action="store_true", help="Create reference_manifest.json if missing")
    args = parser.parse_args(argv)

    run_dir = args.run_dir
    manifest_path = args.manifest or run_dir / "reference_manifest.json"
    if args.init_template and not manifest_path.exists():
        write_json(manifest_path, default_manifest())

    try:
        manifest = load_manifest(manifest_path)
        validate_manifest(
            manifest,
            expected_panel=args.expected_panel,
            fail_on_single_virus=args.fail_on_single_virus,
            require_provenance=not args.allow_incomplete,
        )
        rows = audit_manifest(
            manifest,
            compute_sha256=not args.no_sha256,
            base_dir=manifest_path.parent,
        )
        write_reference_audit(run_dir / "reference_audit.tsv", rows)
        fastq_rows = audit_fastqs(
            manifest,
            base_dir=manifest_path.parent,
            compute_sha256=not args.no_sha256,
        )
        write_fastq_audit(run_dir / "fastq_audit.tsv", fastq_rows)
        manifest["fastq_audit"] = [row.__dict__ for row in fastq_rows]
        write_json(run_dir / "reference_manifest.json", manifest)
        if not args.allow_incomplete:
            validate_audit_rows(
                rows,
                require_hashes=not args.no_sha256,
                require_target_presence=True,
            )
            validate_fastq_audit_rows(
                fastq_rows,
                require_hashes=not args.no_sha256,
                require_source_urls=True,
            )
    except (OSError, BenchmarkContractError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {run_dir / 'reference_audit.tsv'}")
    print(f"Wrote {run_dir / 'fastq_audit.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
