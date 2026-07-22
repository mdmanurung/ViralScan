#!/usr/bin/env python3
"""Fail CI if restricted biological outputs or institutional paths enter ship scope."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

RESTRICTED_SUFFIXES = (
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".h5ad",
    ".bam",
    ".cram",
    ".vcf",
    ".vcf.gz",
)
SHIP_TEXT = {"Dockerfile", "Singularity.def", "pyproject.toml", "MANIFEST.in"}
INSTITUTIONAL_PREFIXES = ("/exports/", "/lustre/", "/gpfs/", "/scratch/")
SYNTHETIC_FIXTURE_SHA256 = {
    "tests/data/evidence_tiny/R1.fastq": "b3e4b9d8fbe251e1430bde40b50fe3e770da100407cdd543c5eff896a6f84d7f",
    "tests/data/evidence_tiny/R2.fastq": "05d6b421d286dce208c9d9f6fb6782b65d0d5fd20415d4064cb7da2036b2bf5e",
}


def _is_pinned_synthetic_fixture(path: Path) -> bool:
    expected = SYNTHETIC_FIXTURE_SHA256.get(path.as_posix())
    if expected is None:
        return False
    try:
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return False
    return observed == expected


def tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"], check=True, stdout=subprocess.PIPE
    ).stdout.decode()
    return [Path(value) for value in result.split("\0") if value]


def violations(files: list[Path]) -> list[str]:
    problems: list[str] = []
    for path in files:
        lower = path.name.lower()
        if lower.endswith(RESTRICTED_SUFFIXES) and not _is_pinned_synthetic_fixture(path):
            problems.append(f"restricted biological artifact is tracked: {path}")
        if not (str(path).startswith("src/") or str(path) in SHIP_TEXT):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for prefix in INSTITUTIONAL_PREFIXES:
            if prefix in text:
                problems.append(f"institutional absolute path {prefix!r} in ship scope: {path}")
    return problems


def main() -> int:
    problems = violations(tracked_files())
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("data-governance check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
