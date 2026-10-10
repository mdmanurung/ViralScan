#!/usr/bin/env python3
"""CELLS-01: freeze ONE emptyDrops call per sample and score every arm over it (user, 2026-10-07).

    python scripts/dsr06_reference_cells.py <round_dir> [--force]

The reference call of a sample is its `combined_off` emptyDrops result (full kb host+virus matrix, no read
filter, so it depends on no other arm; DSR-16). This script

  1. copies each `combined_off/.../results/called_cells.tsv` to `<round>/reference_cells/<ds>__<sample>.tsv`
     and records its provenance in `reference_cells/index.tsv` (source hash, emptyDrops parameters, barcodes
     tested, knee, inflection, code SHA). An existing frozen list is never overwritten unless its content is
     identical, or --force is given;
  2. writes `<round>/cell_level/<ds>__<sample>.tsv.gz`: every (arm, called cell, virus) with an estimated
     viral molecule count > 0, restricted to the frozen cells, for every operative arm. Cell-level virus
     analyses read this file, never an arm's own `is_called_cell`.

The arms' own cell calls (for example the STARsolo host matrix of twostep_v2) stay in their run directories as
diagnostics; they are not used for any rate.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ARMS = ("combined_off", "combined_artefact", "twostep_v2")
REFERENCE_ARM = "combined_off"
INDEX_COLS = [
    "dataset", "sample", "n_cells", "source", "source_sha256", "frozen_sha256", "method", "emptydrops_fdr",
    "emptydrops_lower", "emptydrops_niters", "emptydrops_seed", "barcodes_tested", "knee", "inflection",
    "code_sha", "frozen_utc",
]  # fmt: skip
CELL_COLS = [
    "arm",
    "barcode",
    "virus_name",
    "viral_molecules_total_est",
    "molecules_total_est",
    "viral_fraction",
]


def _common():
    spec = importlib.util.spec_from_file_location(
        "dsr_common_cells", REPO / "scripts" / "dsr_common_cells.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _config_value(config: Path, key: str) -> str:
    m = re.search(rf"^{key}:\s*(\S+)", config.read_text(), re.M) if config.is_file() else None
    return m.group(1) if m else ""


def provenance(round_dir: Path, root: Path, inner: Path) -> dict[str, object]:
    """emptyDrops settings and the knee / inflection / barcodes tested behind a frozen list."""
    manifest = root / "run_manifest.json"
    method = (
        json.loads(manifest.read_text()).get("cell_calling", {}).get("method", "")
        if manifest.is_file()
        else ""
    )
    tested, knee, inflection = 0, "", ""
    drops = inner / "kb-python" / "counts_unfiltered" / "emptydrops_cells.tsv"
    if drops.is_file():
        with drops.open(newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tested += 1
                knee, inflection = r.get("knee", ""), r.get("inflection", "")
    code_sha = (
        round_dir
        / "runs"
        / root.parent.parent.name
        / f"{REFERENCE_ARM}"
        / f"{root.name}.meta"
        / "code_sha.txt"
    )
    cfg = inner / "config.yaml"
    return {
        "method": method,
        "emptydrops_fdr": _config_value(cfg, "emptydrops_fdr"),
        "emptydrops_lower": _config_value(cfg, "emptydrops_lower"),
        "emptydrops_niters": _config_value(cfg, "emptydrops_niters"),
        "emptydrops_seed": _config_value(cfg, "emptydrops_seed"),
        "barcodes_tested": tested,
        "knee": knee,
        "inflection": inflection,
        "code_sha": code_sha.read_text().strip() if code_sha.is_file() else "",
    }


def freeze(round_dir: Path, force: bool = False) -> list[dict[str, object]]:
    """Write the frozen lists and return one index row per sample."""
    common = _common()
    out_dir = round_dir / "reference_cells"
    out_dir.mkdir(exist_ok=True)
    rows = []
    for called in sorted(round_dir.glob(f"runs/*/{REFERENCE_ARM}/*/*/results/called_cells.tsv")):
        inner = called.parent.parent
        root = inner.parent
        dataset, sample = root.parent.parent.name, root.name
        cells = sorted(common.reference_cells(inner))
        target = out_dir / f"{dataset}__{sample}.tsv"
        body = "barcode\n" + "".join(f"{b}\n" for b in cells)
        if target.is_file() and target.read_text() != body and not force:
            raise SystemExit(
                f"{target} exists with different content; refusing to overwrite (use --force)"
            )
        if not target.is_file() or force:
            target.write_text(body)
        rows.append(
            {
                "dataset": dataset,
                "sample": sample,
                "n_cells": len(cells),
                "source": str(called.relative_to(round_dir)),
                "source_sha256": sha256(called),
                "frozen_sha256": sha256(target),
                **provenance(round_dir, root, inner),
                "frozen_utc": datetime.fromtimestamp(target.stat().st_mtime, timezone.utc).strftime(
                    "%Y-%m-%dT%H:%M:%SZ"
                ),
            }
        )
    with (out_dir / "index.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=INDEX_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def export_cell_level(round_dir: Path) -> dict[tuple[str, str], int]:
    """One gz table per sample: viral molecule estimates of the frozen cells, every operative arm."""
    common = _common()
    out_dir = round_dir / "cell_level"
    out_dir.mkdir(exist_ok=True)
    written: dict[tuple[str, str], int] = {}
    for dataset, sample, cells in common.reference_sets(round_dir, {}):
        n = 0
        with gzip.open(out_dir / f"{dataset}__{sample}.tsv.gz", "wt", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(CELL_COLS)
            for arm in ARMS:
                inner = common.arm_inner(round_dir, dataset, arm, sample)
                if inner is None or not common.arm_is_final(round_dir, dataset, arm, sample):
                    continue
                with (inner / "results" / "per_cell_viral.tsv").open(newline="") as src:
                    for r in csv.DictReader(src, delimiter="\t"):
                        if r["barcode"] in cells and float(r["viral_molecules_total_est"] or 0) > 0:
                            w.writerow([arm, r["barcode"], r["virus_name"], r["viral_molecules_total_est"],
                                        r["molecules_total_est"], r["viral_fraction"]])  # fmt: skip
                            n += 1
        written[dataset, sample] = n
    return written


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("round_dir", type=Path)
    ap.add_argument(
        "--force", action="store_true", help="overwrite a frozen list whose content differs"
    )
    args = ap.parse_args(argv)
    rows = freeze(args.round_dir.resolve(), args.force)
    written = export_cell_level(args.round_dir.resolve())
    for r in rows:
        print(
            f"{r['dataset']}\t{r['sample']}\t{r['n_cells']}\t{written.get((r['dataset'], r['sample']), 0)}"
        )
    print(f"{len(rows)} samples frozen", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
