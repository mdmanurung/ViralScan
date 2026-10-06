#!/usr/bin/env python3
"""DSR-02: list every (run, virus) call that needs a read-level check.

Walks <round>/runs/<dataset>/<arm>/<sample>/<sample>/results/viral_summary.tsv and keeps each
virus with >= MIN_MOLECULES (run-plan rule). Anellovirus calls are always kept, whatever the
count, because no anellovirus call has ever been read-validated (F-022). `expected` marks the
virus the dataset is meant to contain; it only orders the table, it never skips a check.

    python scripts/dsr02_enumerate_calls.py <round_dir> [--min-molecules 3] > calls.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

MIN_MOLECULES = 3
ANELLO = re.compile(r"anello|torque|alphatorque|betatorque|gammatorque|\bTTV\b|\bTTMV\b|\bTTMDV\b", re.I)
EXPECTED = {
    "ebv": r"epstein|gammaherpesvirus 4",
    "hhv6b": r"betaherpesvirus 6b|herpesvirus 6b",
    "hsv1": r"alphaherpesvirus 1|herpes simplex virus 1",
    "hpv16": r"papillomavirus type 16|hpv16|alphapapillomavirus 9",
    "kshv_gse190558": r"gammaherpesvirus 8|kaposi",
    "kshv_ebv_gse154900": r"gammaherpesvirus 8|kaposi|epstein|gammaherpesvirus 4",
    "hnscc_gse164690": r"papillomavirus",
    "covid": r"sars|coronavirus 2",
}


def calls(round_dir: Path, min_molecules: int = MIN_MOLECULES):
    for summary in sorted(round_dir.glob("runs/*/*/*/*/results/viral_summary.tsv")):
        sample_dir = summary.parent.parent
        dataset, arm = summary.relative_to(round_dir / "runs").parts[:2]
        with summary.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                n = float(row["viral_molecules_total_est"] or 0)
                anello = bool(ANELLO.search(row["virus_name"]))
                if n < min_molecules and not (anello and n > 0):
                    continue
                pat = EXPECTED.get(dataset)
                yield {
                    "dataset": dataset,
                    "arm": arm,
                    "sample": sample_dir.name,
                    "run_dir": str(sample_dir),
                    "virus": row["virus_name"],
                    "molecules": f"{n:.1f}",
                    "anellovirus": str(anello).lower(),
                    "expected": str(bool(pat and re.search(pat, row["virus_name"], re.I))).lower(),
                }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument("--min-molecules", type=int, default=MIN_MOLECULES)
    args = ap.parse_args(argv)
    rows = list(calls(args.round_dir, args.min_molecules))
    cols = ["dataset", "arm", "sample", "run_dir", "virus", "molecules", "anellovirus", "expected"]
    out = csv.DictWriter(sys.stdout, fieldnames=cols, delimiter="\t", lineterminator="\n")
    out.writeheader()
    out.writerows(rows)
    print(f"{len(rows)} calls", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
