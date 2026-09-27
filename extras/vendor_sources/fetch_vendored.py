#!/usr/bin/env python3
"""Re-fetch the third-party reference sequences vendored for the WP4H panel review.

These files are deliberately untracked (see `.gitignore`): the provenance that
matters is the pinned upstream commit recorded in `VENDOR_MANIFEST.tsv`, not a
blob in this repository.  This script reproduces the working tree from those
pins and verifies every file against the recorded sha256, so a partial or
tampered vendor directory is detected rather than silently used.

Usage
-----
    python extras/vendor_sources/fetch_vendored.py            # fetch + verify
    python extras/vendor_sources/fetch_vendored.py --verify   # verify only
    python extras/vendor_sources/fetch_vendored.py --list     # what is pinned

Nothing here is committed to git and nothing is uploaded anywhere; it only reads
from raw.githubusercontent.com at fixed commits.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "VENDOR_MANIFEST.tsv"
RAW = "https://raw.githubusercontent.com/{repo}/{commit}/{path}"
TIMEOUT = 120


def load_manifest() -> list[dict[str, str]]:
    if not MANIFEST.is_file():
        sys.exit(f"missing manifest: {MANIFEST}")
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        sys.exit(f"manifest is empty: {MANIFEST}")
    return rows


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def column(rows: list[dict[str, str]], name: str) -> str:
    for row in rows:
        if name in row:
            return name
    keys = ", ".join(sorted(rows[0]))
    sys.exit(f"manifest has no column {name!r}; it has: {keys}")


def target_for(row: dict[str, str], local_column: str) -> Path:
    """`local_path` in the manifest is relative to this directory."""
    return HERE / row[local_column]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="check only, do not download")
    parser.add_argument("--list", action="store_true", help="print the pinned inventory")
    args = parser.parse_args()

    rows = load_manifest()
    repo_col = column(rows, "source_repo")
    commit_col = column(rows, "pinned_commit")
    path_col = column(rows, "path_in_repo")
    local_col = column(rows, "local_path")
    sha_col = column(rows, "sha256")
    size_col = column(rows, "bytes")

    if args.list:
        for row in sorted(rows, key=lambda r: (r[repo_col], r[path_col])):
            print(f"{row[repo_col]}@{row[commit_col][:8]}  {row[path_col]}")
        print(f"\n{len(rows)} pinned files")
        return 0

    ok = missing = mismatch = 0
    problems: list[str] = []
    for row in rows:
        repo, commit = row[repo_col], row[commit_col]
        dest = target_for(row, local_col)
        url = RAW.format(repo=repo, commit=commit, path=row[path_col])
        if not args.verify and (not dest.is_file() or sha256_of(dest) != row[sha_col]):
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
                    payload = response.read()
            except (urllib.error.URLError, TimeoutError) as exc:
                problems.append(f"FETCH FAILED {url}: {exc}")
                missing += 1
                continue
            dest.write_bytes(payload)
        if not dest.is_file():
            problems.append(f"MISSING   {dest}")
            missing += 1
            continue
        observed = sha256_of(dest)
        if observed != row[sha_col]:
            problems.append(
                f"SHA MISMATCH {dest}: expected {row[sha_col][:12]} observed {observed[:12]}"
            )
            mismatch += 1
            continue
        if int(row[size_col] or 0) != dest.stat().st_size:
            problems.append(f"SIZE MISMATCH {dest}")
            mismatch += 1
            continue
        ok += 1

    for line in problems:
        print(line, file=sys.stderr)
    print(f"verified {ok}/{len(rows)}  missing={missing}  mismatch={mismatch}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
