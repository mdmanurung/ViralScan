#!/usr/bin/env python3
"""Fail CI if restricted biological outputs or institutional paths enter ship scope."""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

try:
    from scripts.check_ship_scope import load_ship_scope, ship_scope_source_paths
    from scripts.governance_utils import INSTITUTIONAL_PREFIXES
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from check_ship_scope import (  # type: ignore[no-redef]
        load_ship_scope,
        ship_scope_source_paths,
    )
    from governance_utils import INSTITUTIONAL_PREFIXES  # type: ignore[no-redef]

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
SYNTHETIC_FIXTURE_SHA256 = {
    "tests/data/evidence_tiny/R1.fastq": "b3e4b9d8fbe251e1430bde40b50fe3e770da100407cdd543c5eff896a6f84d7f",
    "tests/data/evidence_tiny/R2.fastq": "05d6b421d286dce208c9d9f6fb6782b65d0d5fd20415d4064cb7da2036b2bf5e",
}


def _is_pinned_synthetic_fixture(path: Path, repo_root: Path) -> bool:
    try:
        relative = path if not path.is_absolute() else path.relative_to(repo_root)
    except ValueError:
        return False
    expected = SYNTHETIC_FIXTURE_SHA256.get(relative.as_posix())
    if expected is None:
        return False
    try:
        observed = hashlib.sha256((repo_root / relative).read_bytes()).hexdigest()
    except OSError:
        return False
    return observed == expected


def tracked_files(git_root: Path = REPO_ROOT) -> list[Path]:
    result = subprocess.run(
        ["git", "-C", str(git_root), "ls-files", "-z"],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout.decode("utf-8")
    return [Path(value) for value in result.split("\0") if value]


def violations(
    files: list[Path], text_files: list[Path], *, repo_root: Optional[Path] = None
) -> list[str]:
    repo_root = (repo_root or Path.cwd()).resolve()
    problems: list[str] = []
    for path in files:
        lower = path.name.lower()
        if lower.endswith(RESTRICTED_SUFFIXES) and not _is_pinned_synthetic_fixture(
            path, repo_root
        ):
            problems.append(f"restricted biological artifact is tracked: {path}")

    for path in text_files:
        full_path = path if path.is_absolute() else repo_root / path
        if not full_path.is_file():
            problems.append(f"allowlisted text file is missing: {path}")
            continue
        try:
            text = full_path.read_text(encoding="utf-8")
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


def scope_text_files(
    scope_path: Optional[Path] = None, *, repo_root: Path = REPO_ROOT
) -> tuple[list[Path], list[str]]:
    scope_path = scope_path or repo_root / "config/public_ship_scope.json"
    scope, errors = load_ship_scope(scope_path, repo_root=repo_root)
    if not scope:
        return [], errors
    sources, mapping_errors = ship_scope_source_paths(scope)
    return [Path(value) for value in sorted(sources)], sorted(set(errors + mapping_errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--git-root", type=Path)
    parser.add_argument("--scope", type=Path)
    args = parser.parse_args()

    git_root = args.git_root
    if git_root is None:
        configured_git_root = os.environ.get("VIRALSCAN_GIT_ROOT")
        git_root = Path(configured_git_root) if configured_git_root else args.repo_root
    text_files, problems = scope_text_files(args.scope, repo_root=args.repo_root)
    try:
        files = tracked_files(git_root)
    except (OSError, UnicodeDecodeError, subprocess.CalledProcessError) as exc:
        problems.append(f"cannot enumerate tracked files from git root {git_root}: {exc}")
        files = []
    problems.extend(violations(files, text_files, repo_root=args.repo_root))
    if problems:
        print("\n".join(sorted(set(problems))), file=sys.stderr)
        return 1
    print("data-governance check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
