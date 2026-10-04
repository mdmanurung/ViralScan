#!/usr/bin/env python3
"""CAT-13: fail if any built wheel/sdist exceeds the size limit (default 60 MB)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_LIMIT_MB = 60.0  # PLAN CAT-13: PyPI's 60 MB project limit


def oversized(dist_dir: Path, limit_mb: float) -> list[tuple[Path, float]]:
    files = [*dist_dir.glob("*.whl"), *dist_dir.glob("*.tar.gz")]
    return [(f, f.stat().st_size / 1e6) for f in files if f.stat().st_size > limit_mb * 1e6]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dist_dir", nargs="?", default="dist", type=Path)
    ap.add_argument("--max-mb", type=float, default=DEFAULT_LIMIT_MB)
    args = ap.parse_args(argv)
    if not any(args.dist_dir.glob("*.whl")) and not any(args.dist_dir.glob("*.tar.gz")):
        print(f"no wheel/sdist in {args.dist_dir}", file=sys.stderr)
        return 2
    bad = oversized(args.dist_dir, args.max_mb)
    for f, mb in bad:
        print(f"{f}: {mb:.1f} MB > {args.max_mb:g} MB", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
