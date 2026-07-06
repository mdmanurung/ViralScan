#!/usr/bin/env python3
"""Prepare command and SLURM manifests for the 12-row reference-strategy benchmark."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from viralscan.reference_strategy import (
    BenchmarkContractError,
    audit_manifest,
    audit_fastqs,
    default_manifest,
    load_manifest,
    validate_audit_rows,
    validate_commands,
    validate_fastq_audit_rows,
    validate_manifest,
    write_commands,
    write_fastq_audit,
    write_reference_audit,
    write_json,
    write_slurm_array,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--init-template", action="store_true")
    args = parser.parse_args(argv)

    manifest_path = args.manifest or args.run_dir / "reference_manifest.json"
    try:
        args.run_dir.mkdir(parents=True, exist_ok=True)
        if args.init_template and not manifest_path.exists():
            write_json(manifest_path, default_manifest())
        manifest = load_manifest(manifest_path)
        validate_manifest(manifest, fail_on_single_virus=True, require_provenance=True)
        audit_rows = audit_manifest(manifest, compute_sha256=True, base_dir=manifest_path.parent)
        write_reference_audit(args.run_dir / "reference_audit.tsv", audit_rows)
        validate_audit_rows(audit_rows, require_hashes=True, require_target_presence=True)
        fastq_rows = audit_fastqs(manifest, base_dir=manifest_path.parent)
        write_fastq_audit(args.run_dir / "fastq_audit.tsv", fastq_rows)
        validate_fastq_audit_rows(fastq_rows, require_source_urls=True)
        manifest["fastq_audit"] = [row.__dict__ for row in fastq_rows]
        write_json(args.run_dir / "reference_manifest.json", manifest)
        write_commands(args.run_dir, manifest)
        slurm = write_slurm_array(args.run_dir)
        validate_commands(args.run_dir / "commands.jsonl")
    except (OSError, KeyError, BenchmarkContractError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {args.run_dir / 'commands.jsonl'}")
    print(f"Wrote {slurm}")
    print(f"Submit with: sbatch --array=0-11 {slurm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
