#!/usr/bin/env python3
"""Write or check a per-GTF sha256 manifest for the bundled viral references (PLAN REF-13, partial)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from viralscan.reference_strategy import sha256_file

HEADER = "file\tsha256\tsize_bytes\tsource_date"
DEFAULT_DIR = Path(__file__).resolve().parents[1] / "src" / "viralscan" / "data"


def build_rows(data_dir: Path) -> list[str]:
    # source_date is always "unknown": provenance dates are not recoverable (REF-13 partial).
    return [
        f"{p.name}\t{sha256_file(p)}\t{p.stat().st_size}\tunknown"
        for p in sorted(data_dir.glob("*.gtf"))
    ]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DIR)
    ap.add_argument("-o", "--output", type=Path, help="write manifest here (default: stdout)")
    ap.add_argument("--check", type=Path, metavar="MANIFEST", help="verify data dir against MANIFEST")
    args = ap.parse_args(argv)
    if not args.data_dir.is_dir():
        print(f"data dir not found: {args.data_dir}", file=sys.stderr)
        return 2
    rows = build_rows(args.data_dir)
    if args.check:
        # Compare file/sha256/size only; source_date is not re-derivable.
        want_d = {r[0]: r for r in (l.split("\t")[:3] for l in args.check.read_text().splitlines()[1:] if l)}
        got_d = {r[0]: r for r in (l.split("\t")[:3] for l in rows)}
        bad = [f"missing: {k}" for k in want_d.keys() - got_d.keys()]
        bad += [f"extra: {k}" for k in got_d.keys() - want_d.keys()]
        bad += [f"mismatch: {k}" for k in want_d.keys() & got_d.keys() if want_d[k] != got_d[k]]
        for b in sorted(bad):
            print(b, file=sys.stderr)
        return 1 if bad else 0
    text = "\n".join([HEADER, *rows]) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
