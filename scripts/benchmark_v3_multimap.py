#!/usr/bin/env python3
"""Benchmark the v3 molecule allocator from retained BUS/reference intermediates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import resource
import shutil
import sys
import tempfile
import time
import traceback
import uuid
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import viralscan.constants as virus_constants_module
import viralscan.multimapping as multimapping_module
import viralscan.scripts.multimap as multimap_worker_module
import viralscan.virus_grouping as virus_grouping_module
from viralscan.constants import VIRUS_NAME_MAP
from viralscan.multimapping import build_multimap_layers
from viralscan.scripts.multimap import (
    load_barcodes,
    load_genes,
    load_transcripts,
    prepare_resolved_bus,
    read_ec,
)
from viralscan.virus_grouping import virus_name_for_gene

FROZEN_ENDPOINT_METHODS = ("equal", "host-conservative")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse the retained-BUS benchmark command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sample-dir", type=Path, required=True)
    parser.add_argument("--whitelist", type=Path, required=True)
    parser.add_argument("--t2g", type=Path, required=True)
    parser.add_argument(
        "--method",
        choices=FROZEN_ENDPOINT_METHODS,
        default="host-conservative",
        help="Primary frozen comparison endpoint; all frozen endpoint metrics are emitted.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args(argv)


def endpoint_mass_metrics(layers: object) -> dict[str, dict[str, float | str]]:
    """Return mass accounting for each prespecified v3 comparison endpoint."""
    unique_mass = float(layers.unique.sum())
    endpoints: dict[str, tuple[str, float]] = {
        "v3-unique": ("unique-only", 0.0),
        "v3-equal": ("equal", float(layers.equal.sum())),
        "v3-host-conservative": (
            "host-conservative",
            float(layers.host_conservative.sum()),
        ),
    }
    return {
        endpoint: {
            "allocation_method": method,
            "unique_molecule_mass": unique_mass,
            "allocated_ambiguous_mass": ambiguous_mass,
            "selected_matrix_mass": unique_mass + ambiguous_mass,
        }
        for endpoint, (method, ambiguous_mass) in endpoints.items()
    }


def _atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.with_suffix(path.suffix + ".tmp")
    staging.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    staging.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact_record(path: Path) -> dict[str, int | str]:
    return {"size_bytes": path.stat().st_size, "sha256": _sha256(path)}


def _scientific_parameter_hash(args: argparse.Namespace) -> str:
    payload = {
        "schema_version": "3.0.0",
        "sample_id": args.sample_id,
        "sample_dir": str(args.sample_dir.resolve()),
        "whitelist": str(args.whitelist.resolve()),
        "t2g": str(args.t2g.resolve()),
        "method": args.method,
        "minimum_tmpdir_to_raw_bus_ratio": 3,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(canonical).hexdigest()


def _require_tmpdir_capacity(
    raw_bus: Path,
    tmpdir: Path,
    *,
    minimum_ratio: int = 3,
) -> dict[str, int | str]:
    required_bytes = raw_bus.stat().st_size * minimum_ratio
    available_bytes = shutil.disk_usage(tmpdir).free
    if available_bytes < required_bytes:
        raise OSError(
            "TMPDIR capacity preflight requires at least "
            f"{required_bytes} bytes ({minimum_ratio}x raw BUS); "
            f"{available_bytes} bytes available at {tmpdir}."
        )
    return {
        "tmpdir": str(tmpdir),
        "raw_bus_bytes": raw_bus.stat().st_size,
        "required_tmpdir_bytes": required_bytes,
        "available_tmpdir_bytes": available_bytes,
        "minimum_ratio": minimum_ratio,
    }


def _tree_bytes(root: Path) -> int:
    return sum(path.stat().st_size for path in root.rglob("*") if path.is_file())


def _loaded_module_path(name: str, module: Any) -> Path:
    """Resolve the on-disk artifact backing an imported module."""
    module_file = getattr(module, "__file__", None)
    if module_file is None:
        module_file = getattr(getattr(module, "__spec__", None), "origin", None)
    if not module_file:
        raise RuntimeError(f"Loaded module {name!r} has no on-disk origin to fingerprint.")
    path = Path(module_file).resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Loaded module {name!r} fingerprint artifact does not exist: {path}"
        )
    return path


def _code_fingerprints() -> dict[str, dict[str, int | str]]:
    paths = {
        "benchmark_v3_multimap": Path(__file__).resolve(),
        "multimapping": _loaded_module_path("viralscan.multimapping", multimapping_module),
        "multimap_worker": _loaded_module_path(
            "viralscan.scripts.multimap", multimap_worker_module
        ),
        "virus_grouping": _loaded_module_path(
            "viralscan.virus_grouping", virus_grouping_module
        ),
        "virus_constants": _loaded_module_path(
            "viralscan.constants", virus_constants_module
        ),
    }
    return {name: _artifact_record(path) for name, path in sorted(paths.items())}


def _virus_grouping_fingerprint() -> str:
    canonical = json.dumps(VIRUS_NAME_MAP, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(canonical).hexdigest()


def _matrix_values(matrix: Any) -> np.ndarray:
    return matrix.data if sparse.issparse(matrix) else np.asarray(matrix)


def _validate_matrix(name: str, matrix: Any, shape: tuple[int, int]) -> None:
    if matrix.shape != shape:
        raise ValueError(f"{name} shape {matrix.shape} does not match declared {shape}.")
    values = _matrix_values(matrix)
    if not np.isfinite(values).all():
        raise ValueError(f"{name} contains non-finite values.")
    if (values < 0).any():
        raise ValueError(f"{name} contains negative values.")


def _write_tsv(path: Path, frame: pd.DataFrame) -> None:
    staging = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(staging, sep="\t", index=False)
    staging.replace(path)


def _h5_matrix_shape(node: h5py.Dataset | h5py.Group) -> tuple[int, int]:
    raw_shape = node.shape if isinstance(node, h5py.Dataset) else node.attrs.get("shape")
    if raw_shape is None:
        raise ValueError("Written H5AD matrix node has no declared shape.")
    return tuple(int(value) for value in raw_shape)


def _h5_matrix_rows(
    node: h5py.Dataset | h5py.Group,
    start: int,
    stop: int,
) -> np.ndarray | sparse.csr_matrix:
    if isinstance(node, h5py.Dataset):
        return np.asarray(node[start:stop])
    encoding = node.attrs.get("encoding-type")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding != "csr_matrix":
        raise ValueError(f"Written H5AD matrix encoding {encoding!r} is not supported.")
    shape = _h5_matrix_shape(node)
    indptr = np.asarray(node["indptr"][start : stop + 1], dtype=np.int64)
    offset = int(indptr[0])
    end = int(indptr[-1])
    data = np.asarray(node["data"][offset:end])
    indices = np.asarray(node["indices"][offset:end])
    return sparse.csr_matrix(
        (data, indices, indptr - offset),
        shape=(stop - start, shape[1]),
    )


def _validate_written_h5ad(
    path: Path,
    shape: tuple[int, int],
    *,
    row_chunk_size: int = 4096,
) -> None:
    layer_names = (
        "counts_unique",
        "counts_ambiguous_allocated",
        "counts_multimap_equal",
        "counts_multimap_host_conservative",
    )
    with h5py.File(path, "r") as verified:
        nodes = {"X": verified["X"]}
        for layer in layer_names:
            if layer not in verified["layers"]:
                raise ValueError(f"Written H5AD layer {layer!r} is missing.")
            nodes[layer] = verified["layers"][layer]
        for name, node in nodes.items():
            if _h5_matrix_shape(node) != shape:
                raise ValueError(f"Written H5AD matrix {name!r} is misaligned.")
        for start in range(0, shape[0], row_chunk_size):
            stop = min(start + row_chunk_size, shape[0])
            chunks = {
                name: _h5_matrix_rows(node, start, stop)
                for name, node in nodes.items()
            }
            for name, matrix in chunks.items():
                _validate_matrix(f"written {name}", matrix, (stop - start, shape[1]))
            delta = (
                chunks["X"]
                - chunks["counts_unique"]
                - chunks["counts_ambiguous_allocated"]
            )
            values = _matrix_values(delta)
            if values.size and float(np.max(np.abs(values))) > 1e-9:
                raise ValueError("Written H5AD X does not equal unique plus allocated ambiguity.")


def _write_benchmark_artifacts(
    *,
    output_dir: Path,
    source_adata: ad.AnnData,
    barcode_to_idx: dict[str, int],
    gene_ids: list[str],
    gene_names: list[str],
    viral_indices: set[int],
    layers: object,
    method: str,
) -> dict[str, Path]:
    shape = (len(barcode_to_idx), len(gene_ids))
    expected_barcodes = [
        barcode for barcode, _index in sorted(barcode_to_idx.items(), key=lambda item: item[1])
    ]
    observed_barcodes = [str(barcode).removesuffix("-1") for barcode in source_adata.obs_names]
    if observed_barcodes != expected_barcodes:
        raise ValueError("Source AnnData barcode order does not match the declared BUS barcode universe.")
    if list(source_adata.var_names) != gene_ids:
        raise ValueError("Source AnnData feature order does not match the declared gene universe.")

    matrices = {
        "counts_unique": layers.unique,
        "counts_ambiguous_allocated": layers.corrected,
        "counts_multimap_equal": layers.equal,
        "counts_multimap_host_conservative": layers.host_conservative,
    }
    for name, matrix in matrices.items():
        _validate_matrix(name, matrix, shape)
    selected = layers.unique + layers.corrected
    _validate_matrix("X", selected, shape)

    exported = ad.AnnData(
        X=selected,
        obs=source_adata.obs.copy(),
        var=source_adata.var.copy(),
    )
    exported.obs_names = expected_barcodes
    exported.var["gene_id"] = gene_ids
    exported.var["gene_name"] = gene_names
    exported.var["is_viral"] = [index in viral_indices for index in range(len(gene_ids))]
    for name, matrix in matrices.items():
        exported.layers[name] = matrix
    exported.uns["count_schema_version"] = "3.0.0"
    exported.uns["quantification_unit"] = "bustools-resolved-cb-umi-molecule"
    exported.uns["multimap_method"] = method
    exported.uns["molecule_audit"] = {
        **vars(layers.audit),
        "allocated_ambiguous_mass": float(layers.corrected.sum()),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = output_dir / "adata_v3.h5ad"
    h5ad_staging = h5ad_path.with_suffix(".h5ad.tmp")
    exported.write_h5ad(h5ad_staging)
    h5ad_staging.replace(h5ad_path)
    _validate_written_h5ad(h5ad_path, shape)

    viral_columns = sorted(viral_indices)
    unique_viral = np.asarray(layers.unique[:, viral_columns].sum(axis=0)).reshape(-1)
    equal_viral = np.asarray(
        (layers.unique + layers.equal)[:, viral_columns].sum(axis=0)
    ).reshape(-1)
    conservative_viral = np.asarray(
        (layers.unique + layers.host_conservative)[:, viral_columns].sum(axis=0)
    ).reshape(-1)
    per_feature = pd.DataFrame(
        {
            "gene_id": [gene_ids[index] for index in viral_columns],
            "virus_name": [
                virus_name_for_gene(gene_ids[index]) for index in viral_columns
            ],
            "v3_unique": unique_viral,
            "v3_equal": equal_viral,
            "v3_host_conservative": conservative_viral,
        }
    )
    if per_feature.empty:
        per_virus = pd.DataFrame(
            columns=[
                "virus_name",
                "n_features",
                "v3_unique",
                "v3_equal",
                "v3_host_conservative",
            ]
        )
    else:
        per_virus = (
            per_feature.groupby("virus_name", sort=True, as_index=False)
            .agg(
                n_features=("gene_id", "size"),
                v3_unique=("v3_unique", "sum"),
                v3_equal=("v3_equal", "sum"),
                v3_host_conservative=("v3_host_conservative", "sum"),
            )
        )
    per_cell = pd.DataFrame(
        {
            "barcode": expected_barcodes,
            "v3_unique_viral": np.asarray(
                layers.unique[:, viral_columns].sum(axis=1)
            ).reshape(-1),
            "v3_equal_viral": np.asarray(
                (layers.unique + layers.equal)[:, viral_columns].sum(axis=1)
            ).reshape(-1),
            "v3_host_conservative_viral": np.asarray(
                (layers.unique + layers.host_conservative)[:, viral_columns].sum(axis=1)
            ).reshape(-1),
        }
    )
    per_virus_path = output_dir / "per_virus.tsv"
    per_cell_path = output_dir / "per_cell.tsv"
    _write_tsv(per_virus_path, per_virus)
    _write_tsv(per_cell_path, per_cell)
    count_audit_path = output_dir / "count_audit.json"
    _atomic_json(
        count_audit_path,
        {
            **vars(layers.audit),
            "allocated_ambiguous_mass": float(layers.corrected.sum()),
        },
    )
    return {
        "adata_v3.h5ad": h5ad_path,
        "per_virus.tsv": per_virus_path,
        "per_cell.tsv": per_cell_path,
        "count_audit.json": count_audit_path,
    }


def _write_failure(
    *,
    args: argparse.Namespace,
    command: list[str],
    attempt_id: str,
    scientific_parameter_hash: str,
    context: dict[str, Any],
    error: Exception,
    started: float,
) -> None:
    output_dir = args.output.parent
    stderr_path = output_dir / "stderr.log"
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_staging = stderr_path.with_suffix(".log.tmp")
    stderr_staging.write_text(traceback.format_exc())
    stderr_staging.replace(stderr_path)
    scratch_path = context.get("scratch_path")
    failure = {
        "schema_version": "3.0.0",
        "sample_id": args.sample_id,
        "status": "failed",
        "stage": str(context.get("stage", "initialization")),
        "command": command,
        "exit_code": int(getattr(error, "returncode", 1)),
        "stderr_path": str(stderr_path),
        "attempt_id": attempt_id,
        "scientific_parameter_hash": scientific_parameter_hash,
        "error_type": type(error).__name__,
        "message": str(error),
        "scratch_path": str(scratch_path) if scratch_path else None,
        "scratch_retained": bool(scratch_path and Path(scratch_path).is_dir()),
        "total_wall_seconds": time.monotonic() - started,
    }
    _atomic_json(output_dir / "failure.json", failure)
    _atomic_json(output_dir / "status.json", failure)


def _run(
    args: argparse.Namespace,
    command: list[str],
    context: dict[str, Any],
    attempt_id: str,
    scientific_parameter_hash: str,
    started: float,
) -> int:

    context["stage"] = "validate_inputs"
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
        "whitelist": args.whitelist,
        "t2g": args.t2g,
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing benchmark inputs: " + ", ".join(missing))

    context["stage"] = "fingerprint_inputs_before"
    input_fingerprints_before = {
        name: _artifact_record(path) for name, path in sorted(required.items())
    }
    code_fingerprints = _code_fingerprints()
    grouping_fingerprint = _virus_grouping_fingerprint()
    before_fingerprint_path = args.output.parent / "input_fingerprints_before.json"
    _atomic_json(
        before_fingerprint_path,
        {
            "schema_version": "3.0.0",
            "attempt_id": attempt_id,
            "scientific_parameter_hash": scientific_parameter_hash,
            "inputs": input_fingerprints_before,
            "code": code_fingerprints,
            "virus_grouping_sha256": grouping_fingerprint,
        },
    )

    context["stage"] = "preflight_scratch"
    resolved_bus = args.output.parent / "output.resolved.sorted.bus"
    tmpdir = (
        Path(os.environ["TMPDIR"])
        if os.environ.get("TMPDIR")
        else Path(tempfile.gettempdir())
    )
    preflight = _require_tmpdir_capacity(required["bus"], tmpdir)
    scratch_path = Path(
        tempfile.mkdtemp(
            prefix=f"viralscan-{args.sample_id}-",
            dir=tmpdir,
        )
    )
    context["scratch_path"] = str(scratch_path)
    corrected_bus = scratch_path / "output.corrected.bus"
    resolved_text = scratch_path / "output.resolved.sorted.bus.txt"
    context["stage"] = "prepare_resolved_bus"
    if not resolved_bus.is_file() or not resolved_text.is_file():
        prepare_resolved_bus(
            required["bus"],
            resolved_bus,
            resolved_text,
            whitelist=str(required["whitelist"]),
            threads=args.threads,
            corrected_bus=corrected_bus,
        )
    if not resolved_bus.is_file() or not resolved_text.is_file():
        raise FileNotFoundError("BUS preparation did not produce the resolved binary and text boundary.")

    context["stage"] = "load_count_inputs"
    setup_started = time.monotonic()
    adata = ad.read_h5ad(required["adata"])
    barcode_to_idx, n_cells = load_barcodes(required["barcodes"])
    gene_ids, gene_names = load_genes(required["genes"], required["gene_names"], adata.n_vars)
    transcripts, t2g = load_transcripts(required["transcripts"], required["t2g"])
    ec_map = read_ec(required["ec"], transcripts, t2g, gene_ids)
    viral_ids = set(required["viral_ids"].read_text().splitlines())
    viral_indices = {i for i, gene in enumerate(gene_ids) if gene in viral_ids}
    setup_seconds = time.monotonic() - setup_started

    context["stage"] = "allocate_molecules"
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
    endpoint_metrics = endpoint_mass_metrics(layers)
    for endpoint in ("v3-equal", "v3-host-conservative"):
        allocated = float(endpoint_metrics[endpoint]["allocated_ambiguous_mass"])
        if not math.isclose(
            allocated,
            float(audit.ambiguous_molecules),
            rel_tol=0.0,
            abs_tol=1e-6,
        ):
            raise AssertionError(
                f"{endpoint} allocated mass does not equal audited ambiguous molecules."
            )

    bus_records = int(json.loads(required["run_info"].read_text())["n_pseudoaligned"])
    context["stage"] = "write_scientific_outputs"
    artifact_paths = _write_benchmark_artifacts(
        output_dir=args.output.parent,
        source_adata=adata,
        barcode_to_idx=barcode_to_idx,
        gene_ids=gene_ids,
        gene_names=gene_names,
        viral_indices=viral_indices,
        layers=layers,
        method=args.method,
    )
    artifact_paths[before_fingerprint_path.name] = before_fingerprint_path
    scratch_bytes = _tree_bytes(scratch_path)

    result = {
        "schema_version": "3.0.0",
        "sample": args.sample_id,
        "sample_id": args.sample_id,
        "resolved_bus": str(resolved_bus),
        "bus_records": bus_records,
        "method": args.method,
        "primary_endpoint": f"v3-{args.method}",
        "comparison_endpoints": endpoint_metrics,
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
        "artifacts": sorted(artifact_paths),
        "attempt_id": attempt_id,
        "scientific_parameter_hash": scientific_parameter_hash,
        "resource_metrics": {
            **preflight,
            "scratch_bytes_before_cleanup": scratch_bytes,
            "total_wall_seconds_before_hashing": time.monotonic() - started,
        },
    }
    _atomic_json(args.output, result)
    artifact_paths[args.output.name] = args.output
    if resolved_bus.is_file():
        artifact_paths[resolved_bus.name] = resolved_bus
    context["stage"] = "fingerprint_inputs_after"
    input_fingerprints_after = {
        name: _artifact_record(path) for name, path in sorted(required.items())
    }
    after_fingerprint_path = args.output.parent / "input_fingerprints_after.json"
    _atomic_json(
        after_fingerprint_path,
        {
            "schema_version": "3.0.0",
            "attempt_id": attempt_id,
            "scientific_parameter_hash": scientific_parameter_hash,
            "inputs": input_fingerprints_after,
        },
    )
    artifact_paths[after_fingerprint_path.name] = after_fingerprint_path
    if input_fingerprints_after != input_fingerprints_before:
        raise RuntimeError("One or more benchmark inputs changed during execution.")
    context["stage"] = "hash_outputs"
    hashes = {
        "schema_version": "3.0.0",
        "inputs": input_fingerprints_before,
        "inputs_before": input_fingerprints_before,
        "inputs_after": input_fingerprints_after,
        "code": code_fingerprints,
        "virus_grouping_sha256": grouping_fingerprint,
        "scientific_parameter_hash": scientific_parameter_hash,
        "outputs": {
            name: _artifact_record(path) for name, path in sorted(artifact_paths.items())
        },
    }
    _atomic_json(args.output.parent / "hashes.json", hashes)
    _atomic_json(
        args.output.parent / "status.json",
        {
            "schema_version": "3.0.0",
            "sample_id": args.sample_id,
            "status": "success",
            "stage": "complete",
            "command": command,
            "exit_code": 0,
            "stderr_path": None,
            "attempt_id": attempt_id,
            "scientific_parameter_hash": scientific_parameter_hash,
            "allocator_passes": 1,
            "scratch_path": str(scratch_path),
            "scratch_retained": False,
            "scratch_bytes_before_cleanup": scratch_bytes,
            "total_wall_seconds": time.monotonic() - started,
            "artifacts": sorted([*artifact_paths, "hashes.json"]),
        },
    )
    (args.output.parent / "failure.json").unlink(missing_ok=True)
    print(json.dumps(result, indent=2, sort_keys=True))
    context["stage"] = "cleanup_scratch"
    shutil.rmtree(scratch_path)
    context["scratch_path"] = None
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    command = [sys.argv[0], *(list(argv) if argv is not None else sys.argv[1:])]
    attempt_id = uuid.uuid4().hex
    scientific_parameter_hash = _scientific_parameter_hash(args)
    context: dict[str, Any] = {"stage": "initialization", "scratch_path": None}
    started = time.monotonic()
    try:
        return _run(
            args,
            command,
            context,
            attempt_id,
            scientific_parameter_hash,
            started,
        )
    except Exception as error:
        _write_failure(
            args=args,
            command=command,
            attempt_id=attempt_id,
            scientific_parameter_hash=scientific_parameter_hash,
            context=context,
            error=error,
            started=started,
        )
        raise


if __name__ == "__main__":
    raise SystemExit(main())
