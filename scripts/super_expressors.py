#!/usr/bin/env python3
"""Per-virus super-expressor summary over a frozen cell set (PLAN SE-01).

Follows the HHV-6B CAR-T study (Lareau et al. 2023, Nature): a super-expressor is a cell with at least
``--min-umi`` (default 10) viral UMIs. For every virus in ``per_cell_viral.tsv`` it reports the number and
percentage of super-expressor cells, the share of that virus's UMIs that sit in them, and the maximum viral UMIs
per 10,000 UMIs. The denominator is the cell set (``--cells``), so cells with no viral UMI count as zero; without
``--cells`` the ``is_called_cell`` column is used. Exploratory: no detection claim is made here.

    python scripts/super_expressors.py --per-cell-viral results/per_cell_viral.tsv --cells cells.tsv \
        --dataset ebv --sample GSM4796271 --out super_expressors.tsv
"""

import argparse
import contextlib
import csv
import sys
from collections import defaultdict
from pathlib import Path

SENSITIVITY = (3, 5, 10, 50)


def read_cells(path):
    with open(path, newline="") as fh:
        rows = [line.strip().split("\t")[0] for line in fh if line.strip()]
    return {r for r in rows if r != "barcode"}


def summarise(per_cell_viral, cells=None, min_umi=10):
    """Return one dict per virus. ``cells`` is the denominator set (None = is_called_cell)."""
    by_virus = defaultdict(dict)
    called = set()
    with open(per_cell_viral, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            bc = row["barcode"]
            if row.get("is_called_cell", "").lower() == "true":
                called.add(bc)
            total = float(row["molecules_total_est"] or 0)
            n = float(row["viral_molecules_total_est"] or 0)
            by_virus[row["virus_name"]][bc] = (n, total)
    universe = set(cells) if cells is not None else called
    out = []
    for virus, per in sorted(by_virus.items()):
        vals = {bc: v for bc, v in per.items() if bc in universe}
        counts = [v[0] for v in vals.values()]
        mol = sum(counts)
        se = [v for v in vals.values() if v[0] >= min_umi]
        mol_se = sum(v[0] for v in se)
        row = {
            "virus": virus,
            "n_cells": len(universe),
            "n_cells_viral_umi": sum(1 for c in counts if c >= 1),
            "min_umi": min_umi,
            "n_super_expressors": len(se),
            "pct_cells_super": round(100 * len(se) / len(universe), 4) if universe else 0.0,
            "viral_molecules": round(mol, 3),
            "viral_molecules_in_super": round(mol_se, 3),
            "share_molecules_in_super": round(mol_se / mol, 4) if mol else 0.0,
            "max_viral_per_10k": round(
                max((1e4 * n / t for n, t in vals.values() if t), default=0.0), 1
            ),
        }
        for t in SENSITIVITY:
            row[f"n_ge_{t}"] = sum(1 for c in counts if c >= t)
        out.append(row)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--per-cell-viral", required=True, type=Path)
    ap.add_argument(
        "--cells", type=Path, help="frozen cell list (one barcode per line); default is_called_cell"
    )
    ap.add_argument("--min-umi", type=int, default=10)
    ap.add_argument("--dataset", default="")
    ap.add_argument("--sample", default="")
    ap.add_argument("--out", type=Path, help="TSV path (default stdout)")
    args = ap.parse_args(argv)
    rows = summarise(
        args.per_cell_viral, read_cells(args.cells) if args.cells else None, args.min_umi
    )
    for r in rows:
        r.update(dataset=args.dataset, sample=args.sample)
    fields = ["dataset", "sample"] + [
        k for k in (rows[0] if rows else {}) if k not in ("dataset", "sample")
    ]
    with contextlib.ExitStack() as stack:
        fh = stack.enter_context(open(args.out, "w", newline="")) if args.out else sys.stdout
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
