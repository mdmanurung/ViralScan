#!/usr/bin/env python3
"""Per-library virus/protein UMI summary of one arm-P ``kb count --aa`` output (read-only).

Generalises the AA-02 EBV/HSV-1 summary. Needs ``bustools`` on PATH (BUS records give the
protein-level and ambiguity counts; the H5AD gives per-virus UMIs).

Expected virus: ``--expected-taxa 10376[;37296]`` or ``--expected-table expected_virus.tsv``
(columns ``sample`` and ``taxon_ids``, as written by ``make_configs.py``) with ``--sample``.
``--cells FILE`` restricts everything to a frozen called-cell list (one barcode per line;
a ``-1``-style suffix is ignored); without it all barcodes in the unfiltered output count.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_barcodes import normalize_barcode, read_barcodes  # noqa: E402


def taxon_ids(spec: str) -> list[str]:
    return [
        t if t.startswith("taxon_") else f"taxon_{t}"
        for t in spec.replace(",", ";").split(";")
        if t.strip()
    ]


def read_proteins(path: Path) -> dict[str, tuple[str, str, str, str]]:
    with open(path, newline="") as handle:
        return {
            r["target"]: (r["virus_id"], r["species"], r["product"], r["gene"])
            for r in csv.DictReader(handle, delimiter="\t")
        }


def summarize(
    adata,
    transcripts: list[str],
    ec_lines: Iterable[str],
    bus_rows: Iterable[list[str]],
    proteins: dict,
    expected: list[str],
    cells: set[str] | None = None,
    top: int = 10,
) -> dict:
    """``bus_rows``: ``bustools text`` rows (barcode, umi, ec, count); ``cells``: normalised barcodes."""
    if cells is not None:
        adata = adata[[normalize_barcode(b) in cells for b in adata.obs_names]]
    ec = {
        int(a): [transcripts[int(i)] for i in b.split(",")]
        for a, b in (line.rstrip("\n").split("\t")[:2] for line in ec_lines)
    }
    X = adata.X.tocsc() if hasattr(adata.X, "tocsc") else adata.X
    per_virus = {}
    for i, v in enumerate(adata.var_names):
        col = X[:, i]
        total = int(col.sum())
        if total > 0:
            per_virus[v] = (total, int((col > 0).sum()), int((col >= 10).sum()))
    prot_umi: dict = defaultdict(Counter)
    prot_cells: dict = defaultdict(lambda: defaultdict(set))
    ambiguous: Counter = Counter()
    n_records = 0
    for barcode, _umi, e, _count in (r[:4] for r in bus_rows):
        if cells is not None and normalize_barcode(barcode) not in cells:
            continue
        n_records += 1
        targets = ec[int(e)]
        viruses = {proteins[t][0] for t in targets}
        if len(viruses) == 1:
            v = next(iter(viruses))
            for t in targets:
                prot_umi[v][t] += 1
                prot_cells[v][t].add(barcode)
        else:
            ambiguous[tuple(sorted(viruses))] += 1
    label = {v: next((p[1] for p in proteins.values() if p[0] == v), "-") for v in per_virus}

    def top_proteins(v, n=4):
        return [
            [t.split("__")[-1], proteins[t][2] or proteins[t][3], c, len(prot_cells[v][t])]
            for t, c in prot_umi[v].most_common(n)
        ]

    result = {
        "restricted_to_cells": cells is not None,
        "barcodes_in_matrix": int(adata.n_obs),
        "expected": {
            v: {
                "label": label.get(v, "-"),
                "umis": per_virus.get(v, (0, 0, 0))[0],
                "barcodes": per_virus.get(v, (0, 0, 0))[1],
                "barcodes_ge10_umis": per_virus.get(v, (0, 0, 0))[2],
                "distinct_proteins": len(prot_umi[v]),
                "top_proteins": top_proteins(v),
            }
            for v in expected
        },
        "total_virus_umis": int(X.sum()),
        "barcodes_any_virus": int((X.sum(axis=1) > 0).sum()),
        "viruses_nonzero": len(per_virus),
        "bus_records": n_records,
        "ambiguous_multi_virus_records": sum(ambiguous.values()),
        "top_ambiguous": [[list(k), n] for k, n in ambiguous.most_common(5)],
        "top_viruses": [
            [v, label[v], *c, top_proteins(v, 2)]
            for v, c in sorted(per_virus.items(), key=lambda kv: -kv[1][0])[:top]
        ],
    }
    return result


def bus_text_rows(bus: Path) -> Iterable[list[str]]:
    with subprocess.Popen(
        ["bustools", "text", "-p", str(bus)], stdout=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line.rstrip("\n").split("\t")
    if proc.returncode:
        raise RuntimeError(f"bustools text failed for {bus} (exit {proc.returncode})")


def expected_from_table(path: Path, sample: str) -> list[str]:
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] == sample:
                return taxon_ids(row["taxon_ids"])
    raise SystemExit(f"{sample!r} not in {path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--counts-dir", type=Path, required=True, help="kb count output directory")
    parser.add_argument(
        "--reference-dir", type=Path, required=True, help="directory holding proteins.tsv"
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sample", required=True, help="library name; used for output file names")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--expected-taxa", help="e.g. 10376 or '37296;10376'")
    group.add_argument("--expected-table", type=Path, help="TSV with sample and taxon_ids columns")
    parser.add_argument("--cells", type=Path, help="frozen called-cell barcode list")
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()
    d = args.counts_dir
    expected = (
        taxon_ids(args.expected_taxa)
        if args.expected_taxa
        else expected_from_table(args.expected_table, args.sample)
        if args.expected_table
        else []
    )
    cells = {normalize_barcode(b) for b in read_barcodes(args.cells)} if args.cells else None
    result = summarize(
        ad.read_h5ad(d / "counts_unfiltered" / "adata.h5ad"),
        (d / "transcripts.txt").read_text().split(),
        (d / "matrix.ec").open(),
        bus_text_rows(d / "output.unfiltered.bus"),
        read_proteins(args.reference_dir / "proteins.tsv"),
        expected,
        cells,
        args.top,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / f"{args.sample}.summary.json"
    out.write_text(json.dumps({args.sample: result}, indent=1) + "\n")
    print(json.dumps({args.sample: result}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
