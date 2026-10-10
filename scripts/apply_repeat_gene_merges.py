#!/usr/bin/env python3
"""Merge 31-mer-identical repeat-copy genes in an existing kb index's t2g (PLAN PROG-13).

HSV-1 LAT / RL2 (ICP0) / RS1 (ICP4) and VZV ORF62 / ORF63 each exist twice in the genome and
RefSeq annotates each copy as a gene, so reads from them are multi-gene and never reach the
uniquely-placing layer. This rewrites the gene column of the index's ``t2g.txt`` (or
``panel.t2g``) in place, per ``src/viralscan/data/repeat_gene_merges.tsv``. The index
(``index.idx``) is untouched. New builds already do this
(``scripts/build_bundled_panel_ref.py``); use this for an index built before.

Idempotent: a second run changes nothing. An index whose t2g names none of the listed
genes (for example one built without HSV-1 or VZV) is left alone.

Usage: python scripts/apply_repeat_gene_merges.py --index-dir DIR [--merges TSV]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from viralscan.repeat_merges import (  # noqa: E402
    apply_repeat_merges,
    find_t2g,
    load_repeat_merges,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--index-dir", required=True, type=Path, help="kb index directory")
    parser.add_argument("--merges", type=Path, default=None, help="override the packaged table")
    args = parser.parse_args(argv)
    try:
        t2g = find_t2g(args.index_dir)
    except FileNotFoundError as exc:
        sys.exit(f"ERROR: {exc}")
    changed = apply_repeat_merges(t2g, load_repeat_merges(args.merges))
    print(
        f"{t2g}: {changed} row(s) rewritten"
        if changed
        else f"{t2g}: already merged / nothing to do"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
