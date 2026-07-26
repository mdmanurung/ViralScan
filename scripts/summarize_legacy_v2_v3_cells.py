#!/usr/bin/env python3
"""Summarize barcode and nonzero-cell overlap for validated v2/v3 rows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


class CellSummaryError(RuntimeError):
    """Raised when a planned row cannot be summarized safely."""


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def summarize_row(plan: dict[str, str], row_dir: Path) -> dict[str, Any]:
    """Return barcode and nonzero-cell overlap for one successful row."""

    import anndata as ad
    import numpy as np
    import pandas as pd

    run_id = plan["run_id"]
    status_path = row_dir / "status.json"
    if not status_path.is_file():
        raise CellSummaryError(f"{run_id}: missing status.json")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "success":
        raise CellSummaryError(f"{run_id}: row status is not success")

    legacy_path = Path(plan["adata_multimap_path"])
    per_cell_path = row_dir / "v3" / "per_cell.tsv"
    if not legacy_path.is_file() or not per_cell_path.is_file():
        raise CellSummaryError(f"{run_id}: required cell input is missing")

    legacy = ad.read_h5ad(legacy_path, backed="r")
    try:
        if "is_viral" not in legacy.var or "counts_combined" not in legacy.layers:
            raise CellSummaryError(f"{run_id}: legacy viral layer contract is missing")
        legacy_barcodes = set(map(str, legacy.obs_names))
        viral_columns = np.flatnonzero(
            np.asarray(legacy.var["is_viral"], dtype=bool)
        )
        legacy_mass = np.asarray(
            legacy.layers["counts_combined"][:, viral_columns].sum(axis=1)
        ).ravel()
        legacy_nonzero = {
            str(barcode)
            for barcode, value in zip(legacy.obs_names, legacy_mass, strict=True)
            if float(value) > 0
        }
    finally:
        legacy.file.close()

    v3 = pd.read_csv(per_cell_path, sep="\t")
    required = {
        "barcode",
        "v3_unique_viral",
        "v3_equal_viral",
        "v3_host_conservative_viral",
    }
    if not required.issubset(v3.columns):
        raise CellSummaryError(f"{run_id}: v3 per-cell columns are incomplete")
    if v3["barcode"].duplicated().any():
        raise CellSummaryError(f"{run_id}: duplicate v3 barcodes")
    v3_barcodes = set(v3["barcode"].astype(str))
    endpoint_sets = {
        "unique": set(
            v3.loc[v3["v3_unique_viral"] > 0, "barcode"].astype(str)
        ),
        "equal": set(v3.loc[v3["v3_equal_viral"] > 0, "barcode"].astype(str)),
        "host_conservative": set(
            v3.loc[v3["v3_host_conservative_viral"] > 0, "barcode"].astype(str)
        ),
    }

    output: dict[str, Any] = {
        "run_id": run_id,
        "logical_id": plan.get("logical_id", ""),
        "sample_class": plan.get("sample_class", ""),
        "n_legacy_cells": len(legacy_barcodes),
        "n_v3_cells": len(v3_barcodes),
        "n_barcode_overlap": len(legacy_barcodes & v3_barcodes),
        "barcode_jaccard": _jaccard(legacy_barcodes, v3_barcodes),
        "legacy_nonzero_viral_cells": len(legacy_nonzero),
    }
    for endpoint, cells in endpoint_sets.items():
        output[f"v3_{endpoint}_nonzero_viral_cells"] = len(cells)
        output[f"legacy_v3_{endpoint}_nonzero_overlap"] = len(
            legacy_nonzero & cells
        )
        output[f"legacy_v3_{endpoint}_nonzero_jaccard"] = _jaccard(
            legacy_nonzero, cells
        )
    return output


def summarize_cells(raw_manifest: Path, run_root: Path, output: Path) -> Path:
    """Summarize every planned row and write one sanitized row per run."""

    plans = sorted(_read_tsv(raw_manifest), key=lambda row: row["run_id"])
    rows = [summarize_row(plan, run_root / plan["run_id"]) for plan in plans]
    if not rows:
        raise CellSummaryError("manifest has no planned rows")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-manifest", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(summarize_cells(args.raw_manifest, args.run_root, args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
