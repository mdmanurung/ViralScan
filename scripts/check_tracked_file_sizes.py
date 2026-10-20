#!/usr/bin/env python3
"""Fail if a git-tracked or staged file under analysis/ or results/ exceeds 1 MB.

Sizes are the indexed blob sizes (so staged content counts). Exceptions go in
config/large_file_allowlist.txt, one repo-relative path per line.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

DEFAULT_LIMIT = 1024 * 1024
DEFAULT_ROOTS = ("analysis", "results")
DEFAULT_ALLOWLIST = "config/large_file_allowlist.txt"


def _git(repo: Path, args: list[str], stdin: bytes | None = None) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repo), *args], input=stdin, check=True, capture_output=True
    ).stdout


def load_allowlist(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    lines = (ln.strip() for ln in path.read_text().splitlines())
    return {ln for ln in lines if ln and not ln.startswith("#")}


def oversized(
    repo: Path, roots: tuple[str, ...], limit: int, allow: set[str]
) -> list[tuple[str, int]]:
    """Return (path, bytes) of indexed files under roots larger than limit and not allowlisted."""
    out = _git(repo, ["ls-files", "--stage", "-z", "--", *roots])
    entries = [e.split(b"\t", 1) for e in out.split(b"\0") if e]
    shas = b"".join(meta.split()[1] + b"\n" for meta, _ in entries)
    sizes = _git(repo, ["cat-file", "--batch-check=%(objectsize)"], shas).split()
    bad = []
    for (_, name), size in zip(entries, sizes, strict=True):
        p = name.decode()
        if int(size) > limit and p not in allow:
            bad.append((p, int(size)))
    return bad


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--max-bytes", type=int, default=DEFAULT_LIMIT)
    ap.add_argument("--allowlist", type=Path, default=None)
    ap.add_argument("roots", nargs="*", default=list(DEFAULT_ROOTS))
    args = ap.parse_args(argv)
    allow = load_allowlist(args.allowlist or args.repo / DEFAULT_ALLOWLIST)
    try:
        bad = oversized(args.repo, tuple(args.roots), args.max_bytes, allow)
    except (subprocess.CalledProcessError, OSError) as exc:
        print(f"git failed: {exc}", file=sys.stderr)
        return 2
    for p, size in bad:
        print(f"{p}: {size / 1e6:.2f} MB > {args.max_bytes / 1e6:.2f} MB", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
