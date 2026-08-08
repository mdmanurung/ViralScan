#!/usr/bin/env python3
"""Fail CI if restricted biological outputs or institutional paths enter ship scope."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

try:
    from scripts.governance_utils import INSTITUTIONAL_PREFIXES, load_json
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from governance_utils import INSTITUTIONAL_PREFIXES, load_json  # type: ignore[no-redef]

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
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCOPE = REPO_ROOT / "config" / "public_ship_scope.json"
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


def violations(files: list[Path], text_files: list[Path]) -> list[str]:
    problems: list[str] = []
    for path in files:
        lower = path.name.lower()
        if lower.endswith(RESTRICTED_SUFFIXES) and not _is_pinned_synthetic_fixture(path):
            problems.append(f"restricted biological artifact is tracked: {path}")

    for path in text_files:
        if not path.is_file():
            problems.append(f"allowlisted text file is missing: {path}")
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            problems.append(f"allowlisted text file is not UTF-8: {path}")
            continue
        except OSError as exc:
            problems.append(f"allowlisted text file is unreadable: {path}: {exc}")
            continue
        for prefix in INSTITUTIONAL_PREFIXES:
            if prefix in text:
                problems.append(f"institutional absolute path {prefix!r} in ship scope: {path}")
    return problems


def scope_text_files(scope_path: Path = DEFAULT_SCOPE) -> tuple[list[Path], list[str]]:
    if not scope_path.is_file():
        return [], [f"ship-scope config not found: {scope_path}"]
    try:
        scope = load_json(scope_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return [], [f"cannot load ship-scope config {scope_path}: {exc}"]
    values = scope.get("text_files") if isinstance(scope, dict) else None
    if not isinstance(values, list) or not all(isinstance(value, str) for value in values):
        return [], ["ship-scope text_files must be a list of relative paths"]
    return [REPO_ROOT / value for value in values], []


def main() -> int:
    text_files, problems = scope_text_files()
    problems.extend(violations(tracked_files(), text_files))
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("data-governance check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
