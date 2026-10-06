#!/usr/bin/env python3
"""DSR-16: score every arm of a sample over the SAME called cells.

Arms call cells on different matrices (kb for combined, STARsolo host for twostep), so their
own `pct_infected_called` use different denominators. The reference set of a sample is the
`combined_off` emptyDrops set (full kb matrix, no read filter, so it depends on no arm). Each
arm's `per_cell_viral.tsv` is then counted over exactly those barcodes.

    python scripts/dsr_common_cells.py <round_dir> > common_cells_summary.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

REFERENCE_ARM = "combined_off"
COLS = ["dataset", "sample", "arm", "virus", "n_reference_cells", "infected_in_reference",
        "pct_infected_reference", "molecules_in_reference"]


def reference_cells(sample_dir: Path) -> set[str]:
    with (sample_dir / "results" / "called_cells.tsv").open() as fh:
        next(fh)
        return {line.strip() for line in fh if line.strip()}


def score_arm(sample_dir: Path, cells: set[str]):
    infected: dict[str, int] = defaultdict(int)
    molecules: dict[str, float] = defaultdict(float)
    with (sample_dir / "results" / "per_cell_viral.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n = float(row["viral_molecules_total_est"] or 0)
            if row["barcode"] in cells and n > 0:
                infected[row["virus_name"]] += 1
                molecules[row["virus_name"]] += n
    return infected, molecules


def rows(round_dir: Path):
    for ref_summary in sorted(round_dir.glob(f"runs/*/{REFERENCE_ARM}/*/*/results/called_cells.tsv")):
        ref_dir = ref_summary.parent.parent
        dataset, sample = ref_dir.parent.parent.parent.name, ref_dir.name
        cells = reference_cells(ref_dir)
        for arm_dir in sorted((round_dir / "runs" / dataset).iterdir()):
            sdir = arm_dir / sample / sample
            if arm_dir.name.endswith(".meta") or not (sdir / "results" / "per_cell_viral.tsv").is_file():
                continue
            infected, molecules = score_arm(sdir, cells)
            for virus in sorted(infected, key=lambda v: -molecules[v]):
                yield {
                    "dataset": dataset, "sample": sample, "arm": arm_dir.name, "virus": virus,
                    "n_reference_cells": len(cells), "infected_in_reference": infected[virus],
                    "pct_infected_reference": f"{100 * infected[virus] / len(cells):.4f}",
                    "molecules_in_reference": f"{molecules[virus]:.1f}",
                }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    args = ap.parse_args(argv)
    out = csv.DictWriter(sys.stdout, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    out.writeheader()
    out.writerows(rows(args.round_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
