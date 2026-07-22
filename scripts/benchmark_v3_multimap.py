#!/usr/bin/env python3
"""Benchmark the v3 molecule allocator from retained BUS/reference intermediates."""

from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import anndata as ad

from viralscan.multimapping import build_multimap_layers
from viralscan.scripts.multimap import (
    load_barcodes,
    load_genes,
    load_transcripts,
    prepare_resolved_bus,
    read_ec,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-dir", type=Path, required=True)
    parser.add_argument("--t2g", type=Path, required=True)
    parser.add_argument("--method", default="host-conservative")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    kb = args.sample_dir / "kb-python"
    counts = kb / "counts_unfiltered"
    required = {
        "bus": kb / "output.bus",
        "ec": kb / "matrix.ec",
        "transcripts": kb / "transcripts.txt",
        "barcodes": counts / "cells_x_genes.barcodes.txt",
        "genes": counts / "cells_x_genes.genes.txt",
        "gene_names": counts / "cells_x_genes.genes.names.txt",
        "adata": counts / "adata.h5ad",
        "run_info": kb / "run_info.json",
        "viral_ids": args.sample_dir / "log" / "analysis.txt",
        "t2g": args.t2g,
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing benchmark inputs: " + ", ".join(missing))

    resolved_bus = args.output.parent / "output.resolved.sorted.bus"
    resolved_text = args.output.parent / "output.resolved.sorted.bus.txt"
    if not resolved_bus.is_file() or not resolved_text.is_file():
        prepare_resolved_bus(
            required["bus"],
            resolved_bus,
            resolved_text,
            whitelist=None,
            threads=args.threads,
        )

    started = time.monotonic()
    adata = ad.read_h5ad(required["adata"])
    barcode_to_idx, n_cells = load_barcodes(required["barcodes"])
    gene_ids, _ = load_genes(required["genes"], required["gene_names"], adata.n_vars)
    transcripts, t2g = load_transcripts(required["transcripts"], required["t2g"])
    ec_map = read_ec(required["ec"], transcripts, t2g, gene_ids)
    viral_ids = set(required["viral_ids"].read_text().splitlines())
    viral_indices = {i for i, gene in enumerate(gene_ids) if gene in viral_ids}
    setup_seconds = time.monotonic() - started

    allocation_started = time.monotonic()
    layers = build_multimap_layers(
        resolved_text,
        barcode_to_idx,
        ec_map,
        n_cells,
        adata.n_vars,
        viral_indices,
        adata.X,
        method=args.method,
    )
    allocation_seconds = time.monotonic() - allocation_started
    audit = layers.audit
    audit.validate(float(layers.corrected.sum()))
    if audit.resolved_molecules > audit.input_molecules:
        raise AssertionError("Resolved molecule count exceeds selected CB-UMI input.")

    bus_records = int(json.loads(required["run_info"].read_text())["n_pseudoaligned"])

    result = {
        "schema_version": "3.0.0",
        "sample": "SRR12682296",
        "bus_records": bus_records,
        "method": args.method,
        "setup_seconds": setup_seconds,
        "allocation_seconds": allocation_seconds,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "n_cells": n_cells,
        "n_genes": adata.n_vars,
        "n_ec": len(ec_map),
        "audit": vars(audit),
        "allocated_ambiguous_mass": float(layers.corrected.sum()),
        "unique_molecule_mass": float(layers.unique.sum()),
        "selected_matrix_mass": float((layers.unique + layers.corrected).sum()),
        "legacy_bustools_matrix_mass": float(adata.X.sum()),
        "method_diagnostics": layers.method_diagnostics,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    staging = args.output.with_suffix(args.output.suffix + ".tmp")
    staging.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    staging.replace(args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
