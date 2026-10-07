#!/usr/bin/env python3
"""DSR-16: score every arm of a sample over the SAME called cells.

Arms call cells on different matrices (kb for combined, STARsolo host for twostep), so their
own `pct_infected_called` use different denominators. The reference set of a sample is the
`combined_off` emptyDrops set (full kb matrix, no read filter, so it depends on no arm). Each
arm's `per_cell_viral.tsv` is then counted over exactly those barcodes.

    python scripts/dsr_common_cells.py <round_dir> > common_cells_summary.tsv

A sample with a vendor cell call (cellranger filtered barcodes) uses that instead:
    --reference-cells sfl_tonsil/x223=<...>/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

REFERENCE_ARM = "combined_off"
COLS = ["dataset", "sample", "arm", "virus", "n_reference_cells", "infected_in_reference",
        "pct_infected_reference", "molecules_in_reference", "evidence_status"]


def reference_cells(sample_dir: Path) -> set[str]:
    with (sample_dir / "results" / "called_cells.tsv").open() as fh:
        next(fh)
        return {line.strip() for line in fh if line.strip()}


def vendor_cells(path: Path) -> set[str]:
    """Barcodes from a cellranger barcodes.tsv[.gz] or sample_filtered_barcodes.csv (`-1` stripped)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return {re.sub(r"-\d+$", "", re.split(r"[,\t]", line.strip())[-1]) for line in fh if line.strip()}


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


def reference_sets(round_dir: Path, overrides: dict[str, Path]):
    seen = set()
    for ref_summary in sorted(round_dir.glob(f"runs/*/{REFERENCE_ARM}/*/*/results/called_cells.tsv")):
        ref_dir = ref_summary.parent.parent
        dataset, sample = ref_dir.parent.parent.parent.name, ref_dir.name
        key = f"{dataset}/{sample}"
        seen.add(key)
        yield dataset, sample, vendor_cells(overrides[key]) if key in overrides else reference_cells(ref_dir)
    for key in sorted(set(overrides) - seen):  # vendor-called sample whose combined_off has not finished
        dataset, sample = key.split("/")
        yield dataset, sample, vendor_cells(overrides[key])


def evidence_status(dataset, sample, arm, virus, verdicts):
    """F-028: a twostep call is reported only with a viral_best read verdict.

    The verdict comes from this arm's evidence directory, else the pre-v2 twostep one, else combined_off's.
    A combined arm merely calling the virus is not enough: HPV77 is called by both on host CAG reads.
    """
    if not arm.startswith("twostep"):
        return "not_gated"
    slug = re.sub(r"[^A-Za-z0-9]+", "_", virus).strip("_")
    for a in (arm, arm.removesuffix("_v2"), REFERENCE_ARM):
        verdict = verdicts.get(f"{dataset}__{sample}__{a}__{slug}")
        if verdict is not None:
            return "verified" if verdict == "viral_best" else f"rejected_{verdict}"
    return "needs_evidence"


def rows(round_dir: Path, overrides: dict[str, Path] | None = None, verdicts: dict[str, str] | None = None):
    for dataset, sample, cells in reference_sets(round_dir, overrides or {}):
        for arm_dir in sorted((round_dir / "runs" / dataset).iterdir()):
            sdir = arm_dir / sample / sample
            if arm_dir.name.endswith(".meta") or not (sdir / "results" / "per_cell_viral.tsv").is_file():
                continue
            # a twostep_v2 copy inherits v1's results until its detection finishes, which writes the marker
            v2_done = (round_dir / "runs" / dataset / "twostep_v2" / sample / "run_complete.json").is_file()
            if (arm_dir.name == "twostep_v2" and not v2_done) or (arm_dir.name == "twostep" and v2_done):
                continue
            infected, molecules = score_arm(sdir, cells)
            for virus in sorted(infected, key=lambda v: -molecules[v]):
                yield {
                    "dataset": dataset, "sample": sample, "arm": arm_dir.name, "virus": virus,
                    "n_reference_cells": len(cells), "infected_in_reference": infected[virus],
                    "pct_infected_reference": f"{100 * infected[virus] / len(cells):.4f}",
                    "molecules_in_reference": f"{molecules[virus]:.1f}",
                    "evidence_status": evidence_status(dataset, sample, arm_dir.name, virus, verdicts or {}),
                }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument("--reference-cells", action="append", default=[], metavar="DATASET/SAMPLE=PATH")
    ap.add_argument("--verdicts", type=Path, help="evidence/verdicts.tsv from dsr02_verdicts.py")
    args = ap.parse_args(argv)
    verdicts = {}
    if args.verdicts:
        with args.verdicts.open() as fh:
            verdicts = {r["evidence_dir"]: r["verdict"] for r in csv.DictReader(fh, delimiter="\t")}
    overrides = {k: Path(v) for k, v in (a.split("=", 1) for a in args.reference_cells)}
    out = csv.DictWriter(sys.stdout, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    out.writeheader()
    out.writerows(rows(args.round_dir, overrides, verdicts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
