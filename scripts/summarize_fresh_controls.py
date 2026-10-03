#!/usr/bin/env python3
"""Sanitized LVC-13 tables from a ``compare_legacy_v2_v3.py fresh-vs-archive`` report.

Writes ``fresh_control_stack_comparison.tsv``: per sample and virus, the fresh v2
total, the fresh v3 host-conservative total and both archived values, for every
virus non-zero in any of the four. Descriptive only; no new metric.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

COLUMNS = ["sample_id", "virus", "fresh_v2", "fresh_v3_host_conservative", "archived_v2", "archived_v3_host_conservative"]


def stack_rows(report: dict) -> list[dict]:
    cells: dict[tuple[str, str], dict[str, float | None]] = defaultdict(dict)
    for rec in report["records"]:
        if rec["record_type"] != "virus":
            continue
        key = (rec["sample_id"], rec["identifier"])
        fresh, archived = ("fresh_v2", "archived_v2") if rec["stack"] == "v2" else ("fresh_v3_host_conservative", "archived_v3_host_conservative")
        cells[key][fresh] = rec["fresh_value"]
        cells[key][archived] = rec["archived_value"]
    rows = []
    for (sample, virus), vals in sorted(cells.items()):
        row = {"sample_id": sample, "virus": virus, **{c: vals.get(c) for c in COLUMNS[2:]}}
        if any(v for v in vals.values()):
            rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = stack_rows(json.loads(args.report.read_text(encoding="utf-8")))
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
