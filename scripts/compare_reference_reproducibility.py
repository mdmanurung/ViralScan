#!/usr/bin/env python3
"""Compare two ``reference_reproducibility.json`` files (PLAN REF-04).

Exit 0 when the deterministic content is identical (``content_sha256`` and every
per-file digest), 1 with one line per differing field otherwise. The binary
kallisto ``.idx`` is deliberately not part of this audit; see
``docs/reference_maintenance.md``.

    python scripts/compare_reference_reproducibility.py rebuild_a.json rebuild_b.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any


def compare_reproducibility(first: dict[str, Any], second: dict[str, Any]) -> list[str]:
    """One ``field: a != b`` line per difference; empty when identical."""
    diffs = []
    if first.get("content_sha256") != second.get("content_sha256"):
        diffs.append(
            f"content_sha256: {first.get('content_sha256')} != {second.get('content_sha256')}"
        )
    files_a, files_b = first.get("content_files", {}), second.get("content_files", {})
    for name in sorted(files_a.keys() | files_b.keys()):
        if files_a.get(name) != files_b.get(name):
            diffs.append(f"content_files.{name}: {files_a.get(name)} != {files_b.get(name)}")
    return diffs


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    args = parser.parse_args(argv)
    diffs = compare_reproducibility(
        json.loads(args.first.read_text()), json.loads(args.second.read_text())
    )
    for line in diffs:
        print(line)
    if not diffs:
        print("identical")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
