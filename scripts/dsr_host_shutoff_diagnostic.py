#!/usr/bin/env python3
"""Q17: do host+virus emptyDrops calls lose infected cells that the host-only anchor keeps?

    python scripts/dsr_host_shutoff_diagnostic.py <round_dir> --anchors DIR --out-dir DIR

Strong host shutoff (HSV-1, lytic KSHV, EBV lytic, HHV-6B) shrinks a cell's host UMI count, so a call on the
host+virus matrix can drop it while the host-only anchor (`dsr06_reference_cells.py --source host_called`) keeps it.
Per library, with the round's `combined_off` call (`results/called_cells.tsv`) and anchor:

  combined_only = called on host+virus, not in the anchor (the cells an anchor-based analysis drops)
  anchor_only   = in the anchor, not called on host+virus

Viral load is `viral_molecules_total_est` of the expected virus in `results/per_cell_viral.tsv` (a barcode with no
row has 0). Also reports the barcode format of both lists (length, `-N` suffix): they must match.
Writes `host_shutoff_diagnostic.tsv` (one row per library x virus) and `host_shutoff_top_dropped.tsv`.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# dataset -> expected virus names (exact `virus_name`); other libraries fall back to the top virus by molecules.
EXPECTED = {
    "ebv": ["Epstein-Barr virus"],
    "hhv6b": ["Human herpesvirus 6b"],
    "hsv1": ["Human herpesvirus 1"],
    "hpv16": ["Human papillomavirus 16"],
    "hnscc_gse164690": ["Human papillomavirus 16"],
    "kshv_gse190558": ["Human herpesvirus 8"],
    "kshv_ebv_gse154900": ["Human herpesvirus 8", "Epstein-Barr virus"],
}
COLS = [
    "dataset", "sample", "virus", "expected", "n_anchor", "n_combined_called", "n_both", "n_combined_only",
    "n_anchor_only", "combined_only_with_virus", "combined_only_ge10", "combined_only_ge100",
    "combined_only_max", "combined_only_molecules", "combined_called_molecules",
    "pct_combined_called_molecules_dropped", "anchor_only_with_virus", "anchor_only_ge10", "anchor_only_ge100",
    "barcode_format", "flag",
]  # fmt: skip
TOP_COLS = ["dataset", "sample", "virus", "rank", "barcode", "viral_molecules_total_est", "total_molecules"]


def read_barcodes(path: Path) -> set[str]:
    with path.open() as fh:
        next(fh)
        return {line.strip() for line in fh if line.strip()}


def barcode_format(cells: set[str]) -> str:
    lens = sorted({len(b) for b in cells})
    suffix = any(re.search(r"-\d+$", b) for b in cells)
    return f"len={','.join(map(str, lens))};suffix={'yes' if suffix else 'no'}"


def inner_dir(root: Path, sample: str) -> Path | None:
    if (root / sample).is_dir():
        return root / sample
    subs = [p for p in root.iterdir() if (p / "results").is_dir()] if root.is_dir() else []
    return subs[0] if len(subs) == 1 else None


def load_viral(path: Path) -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    """virus -> barcode -> molecules, and each barcode's total molecules."""
    load: dict[str, dict[str, float]] = defaultdict(dict)
    total: dict[str, float] = {}
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            load[r["virus_name"]][r["barcode"]] = float(r["viral_molecules_total_est"] or 0)
            total[r["barcode"]] = float(r["molecules_total_est"] or 0)
    return load, total


def diagnose(round_dir: Path, anchors: Path, expected: dict[str, list[str]] | None = None):
    """Yield (summary row, top-dropped rows) per library x virus."""
    expected = EXPECTED if expected is None else expected
    for anchor in sorted(anchors.glob("*__*.tsv")):
        if anchor.name == "index.tsv":
            continue
        dataset, sample = anchor.stem.split("__", 1)
        inner = inner_dir(round_dir / "runs" / dataset / "combined_off" / sample, sample)
        if inner is None or not (inner / "results" / "per_cell_viral.tsv").is_file():
            continue
        anc, called = read_barcodes(anchor), read_barcodes(inner / "results" / "called_cells.tsv")
        load, total = load_viral(inner / "results" / "per_cell_viral.tsv")
        viruses = expected.get(dataset) or sorted(load, key=lambda v: -sum(load[v].values()))[:1]
        only_c, only_a = called - anc, anc - called
        fmt_a, fmt_c = barcode_format(anc), barcode_format(called)
        fmt = fmt_a if fmt_a == fmt_c else f"MISMATCH anchor[{fmt_a}] combined[{fmt_c}]"
        for virus in viruses:
            v = load.get(virus, {})
            dropped = sorted(only_c, key=lambda b: -v.get(b, 0))
            d = [v.get(b, 0) for b in only_c]
            a = [v.get(b, 0) for b in only_a]
            called_mol = sum(v.get(b, 0) for b in called)
            ge10 = sum(x >= 10 for x in d)
            row = {
                "dataset": dataset, "sample": sample, "virus": virus, "expected": str(dataset in expected),
                "n_anchor": len(anc), "n_combined_called": len(called), "n_both": len(anc & called),
                "n_combined_only": len(only_c), "n_anchor_only": len(only_a),
                "combined_only_with_virus": sum(x > 0 for x in d), "combined_only_ge10": ge10,
                "combined_only_ge100": sum(x >= 100 for x in d), "combined_only_max": f"{max(d, default=0):.1f}",
                "combined_only_molecules": f"{sum(d):.1f}", "combined_called_molecules": f"{called_mol:.1f}",
                "pct_combined_called_molecules_dropped": f"{100 * sum(d) / called_mol:.2f}" if called_mol else "NA",
                "anchor_only_with_virus": sum(x > 0 for x in a), "anchor_only_ge10": sum(x >= 10 for x in a),
                "anchor_only_ge100": sum(x >= 100 for x in a), "barcode_format": fmt,
                "flag": "HIGH_LOAD_CELLS_LOST" if ge10 else "",
            }  # fmt: skip
            top = [
                {"dataset": dataset, "sample": sample, "virus": virus, "rank": i + 1, "barcode": b,
                 "viral_molecules_total_est": v.get(b, 0), "total_molecules": total.get(b, 0)}
                for i, b in enumerate(dropped[:20]) if v.get(b, 0) > 0
            ]  # fmt: skip
            yield row, top


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument("--anchors", type=Path, required=True, help="dir of <ds>__<sample>.tsv anchors")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, tops = [], []
    for row, top in diagnose(args.round_dir, args.anchors):
        rows.append(row)
        tops += top
    for name, cols, data in (("host_shutoff_diagnostic", COLS, rows), ("host_shutoff_top_dropped", TOP_COLS, tops)):
        with (args.out_dir / f"{name}.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
    print(f"{len(rows)} library x virus rows, {sum(r['flag'] != '' for r in rows)} flagged", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
