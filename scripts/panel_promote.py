#!/usr/bin/env python3
"""PANEL-01 (WP2): apply an accepted promotions table to the virus catalogue's overlay columns.

    python scripts/panel_promote.py analysis/panel_expansion/promotions_round1.tsv
        [--catalogue src/viralscan/data/virus_catalog.tsv] [--check]

`promotions_round1.tsv` has one row per accession: `accession`, `curation_id` (the row of the curation sheet that
accepted it), `panel`, `sibling_group` (blank leaves the catalogue value), `decided_by`. Only `panel` and
`sibling_group` are written; every other column and the row order stay as they are. An accession that is not in the
catalogue is an error: add it first with `extras/build_virus_catalog.py`. The table is the record of who accepted which
accession (as `index_exclusions.tsv` is for omissions). `--check` reports what would change and writes nothing.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

DEFAULT_CATALOGUE = (
    Path(__file__).resolve().parents[1] / "src" / "viralscan" / "data" / "virus_catalog.tsv"
)
OVERLAYS = ("panel", "sibling_group")


def apply(
    rows: list[dict[str, str]], promotions: list[dict[str, str]]
) -> list[tuple[str, str, str, str]]:
    """Edit `rows` in place; return (accession, column, old, new) for every change."""
    by_acc = {r["accession"]: r for r in rows}
    missing = sorted({p["accession"] for p in promotions} - by_acc.keys())
    if missing:
        sys.exit(
            f"not in the catalogue (add with extras/build_virus_catalog.py): {', '.join(missing)}"
        )
    changes = []
    for p in promotions:
        row = by_acc[p["accession"]]
        for col in OVERLAYS:
            if p.get(col) and row[col] != p[col]:
                changes.append((p["accession"], col, row[col], p[col]))
                row[col] = p[col]
    return changes


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("promotions", type=Path)
    ap.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    with args.catalogue.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    with args.promotions.open(newline="") as handle:
        promotions = list(csv.DictReader(handle, delimiter="\t"))
    changes = apply(rows, promotions)
    print(f"{len(changes)} changes in {len({c[0] for c in changes})} accessions", file=sys.stderr)
    if not args.check:
        with args.catalogue.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
