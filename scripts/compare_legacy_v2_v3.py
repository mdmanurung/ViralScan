#!/usr/bin/env python3
"""Read-only inventory and comparison harness for archived ViralScan 2.2 trees.

Tracked manifests produced by this module contain only paths relative to the
declared archive root. The raw manifest is intentionally written below the
gitignored benchmark run root and retains absolute paths for execution.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml


class InventoryContractError(ValueError):
    """Raised when the archive does not satisfy the frozen cohort contract."""


class LegacyReconstructionError(ValueError):
    """Raised when retained v2 artifacts cannot be reconstructed safely."""


class RunRowError(ValueError):
    """Raised when a manifest row cannot be executed without violating its contract."""


@dataclass(frozen=True)
class InventoryContract:
    """Expected inventory cardinalities for fail-closed cohort discovery."""

    technical_rows: int = 44
    logical_inputs: int = 42
    chemistry_counts: Mapping[str, int] = field(default_factory=lambda: {"10xv2": 6, "10xv3": 38})


@dataclass(frozen=True)
class LegacySummary:
    """Machine-readable values retained in a v2.2 ``summary.txt``."""

    gene_totals: Mapping[str, float]
    virus_totals: Mapping[str, float]
    total_viral_load: Optional[float]
    no_viral_genes: bool


@dataclass(frozen=True)
class LegacyReconstruction:
    """Reconstructed v2 calls and an explicit summary comparison audit."""

    called_gene_totals: Mapping[str, float]
    reconstructed_virus_totals: Mapping[str, float]
    reconstructed_total_viral_load: float
    summary: LegacySummary
    comparisons: tuple[dict[str, Any], ...]


SANITIZED_FIELDS = [
    "run_id",
    "logical_id",
    "technical_repeat_group",
    "chemistry",
    "sample_class",
    "expected_target",
    "sample1_name",
    "sample2_name",
    "config_path",
    "summary_path",
    "output_bus_path",
    "matrix_ec_path",
    "transcripts_path",
    "whitelist_path",
    "adata_multimap_path",
    "barcodes_path",
    "genes_path",
    "gene_names_path",
    "adata_path",
    "run_info_path",
    "analysis_path",
    "config_bytes",
    "config_sha256",
    "summary_bytes",
    "summary_sha256",
    "output_bus_bytes",
    "output_bus_sha256",
    "matrix_ec_bytes",
    "matrix_ec_sha256",
    "transcripts_bytes",
    "transcripts_sha256",
    "whitelist_bytes",
    "whitelist_sha256",
    "adata_multimap_bytes",
    "adata_multimap_sha256",
    "barcodes_bytes",
    "barcodes_sha256",
    "genes_bytes",
    "genes_sha256",
    "gene_names_bytes",
    "gene_names_sha256",
    "adata_bytes",
    "adata_sha256",
    "run_info_bytes",
    "run_info_sha256",
    "analysis_bytes",
    "analysis_sha256",
    "index_bytes",
    "index_sha256",
    "t2g_bytes",
    "t2g_sha256",
    "transcriptome_bytes",
    "transcriptome_sha256",
    "reference_hash_id",
    "kb_python_version",
    "kallisto_version",
    "bustools_version",
]

RAW_PATH_FIELDS = [
    "sample1_path",
    "sample2_path",
    "config_path",
    "summary_path",
    "output_bus_path",
    "matrix_ec_path",
    "transcripts_path",
    "whitelist_path",
    "adata_multimap_path",
    "barcodes_path",
    "genes_path",
    "gene_names_path",
    "adata_path",
    "run_info_path",
    "analysis_path",
    "index_path",
    "t2g_path",
    "transcriptome_path",
]

RUN_ROW_HASH_ARTIFACTS = (
    "config",
    "summary",
    "output_bus",
    "matrix_ec",
    "transcripts",
    "whitelist",
    "adata_multimap",
    "barcodes",
    "genes",
    "gene_names",
    "adata",
    "run_info",
    "analysis",
    "index",
    "t2g",
    "transcriptome",
)


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    if not slug:
        raise InventoryContractError(f"cannot sanitize identifier {value!r}")
    return slug


def _run_id(relative_run_dir: Path) -> str:
    return "__".join(_slug(part) for part in relative_run_dir.parts)


def _read_config(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise InventoryContractError(f"cannot read config {path}: {exc}") from exc
    if not isinstance(config, dict):
        raise InventoryContractError(f"config is not a mapping: {path}")
    return config


def _required_config_path(config: Mapping[str, Any], key: str, config_path: Path) -> Path:
    value = config.get(key)
    if not isinstance(value, str) or not value or value == "None":
        raise InventoryContractError(f"{config_path}: missing required config value {key}")
    return Path(value).expanduser().resolve()


def _sample_stem(path_value: str) -> str:
    name = Path(path_value).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"(?:_R?[12](?:_001)?)$", "", name, flags=re.IGNORECASE)


def _logical_id(sample1: str, sample2: str, config_path: Path) -> str:
    first = _sample_stem(sample1)
    second = _sample_stem(sample2)
    if first != second:
        raise InventoryContractError(
            f"{config_path}: paired FASTQ names do not resolve to one logical input "
            f"({first!r} != {second!r})"
        )
    return _slug(first)


def _classification(logical_id: str) -> tuple[str, str]:
    if logical_id in {"SRR12682296", "SRR12682297", "SRR12682298"}:
        return "ebv_control", "EBV"
    if logical_id in {"SRR6825024", "SRR6825025"}:
        return "hiv_control", "HIV-1"
    return "skin", ""


def _find_whitelist(kb_dir: Path, chemistry: str) -> Path:
    expected = kb_dir / {
        "10xv2": "10x_version2_whitelist.txt",
        "10xv3": "10x_version3_whitelist.txt",
    }.get(chemistry, "")
    if chemistry not in {"10xv2", "10xv3"}:
        raise InventoryContractError(f"unsupported chemistry {chemistry!r} in {kb_dir}")
    if not expected.is_file():
        raise InventoryContractError(
            f"{kb_dir}: missing archived {chemistry} whitelist {expected.name}"
        )
    other = list(kb_dir.glob("10x_version*_whitelist.txt"))
    if other != [expected]:
        found = ", ".join(path.name for path in sorted(other))
        raise InventoryContractError(f"{kb_dir}: expected exactly {expected.name}; found [{found}]")
    return expected


def _require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise InventoryContractError(f"missing required {label}: {path}")
    return path.resolve()


def _sha256(path: Path, cache: dict[Path, tuple[int, str]]) -> tuple[int, str]:
    resolved = path.resolve()
    if resolved not in cache:
        digest = hashlib.sha256()
        with resolved.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(chunk)
        cache[resolved] = (resolved.stat().st_size, digest.hexdigest())
    return cache[resolved]


def _reference_hash_id(hashes: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in hashes:
        digest.update(value.encode("ascii"))
        digest.update(b"\0")
    return f"sha256:{digest.hexdigest()}"


def _write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {field_name: row.get(field_name, "") for field_name in fields} for row in rows
        )


_GENE_TOTAL_RE = re.compile(
    r"^(?P<identifier>[^;]+);(?P<value>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)$"
)
_VIRUS_TOTAL_RE = re.compile(
    r"^(?P<identifier>.+) has a viral load of: "
    r"(?P<value>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?) UMIs\.$"
)
_OVERALL_TOTAL_RE = re.compile(
    r"^Total amount of viral load found: "
    r"(?P<value>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)$"
)


def _finite_float(value: str, *, source: Path, line_number: int) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise LegacyReconstructionError(
            f"{source}:{line_number}: non-finite legacy summary value {value!r}"
        )
    return parsed


def parse_legacy_summary(path: Path) -> LegacySummary:
    """Parse v2 gene, virus, and overall totals without dropping omissions."""

    gene_totals: dict[str, float] = {}
    virus_totals: dict[str, float] = {}
    total_viral_load: Optional[float] = None
    no_viral_genes = False
    in_gene_table = False
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        line = raw_line.strip()
        if line == "Gene ID; Gene Count":
            in_gene_table = True
            continue
        if in_gene_table:
            if not line:
                in_gene_table = False
                continue
            match = _GENE_TOTAL_RE.match(line)
            if not match:
                raise LegacyReconstructionError(
                    f"{path}:{line_number}: malformed legacy gene total {line!r}"
                )
            identifier = match.group("identifier")
            if identifier in gene_totals:
                raise LegacyReconstructionError(
                    f"{path}:{line_number}: duplicate legacy gene {identifier!r}"
                )
            gene_totals[identifier] = _finite_float(
                match.group("value"), source=path, line_number=line_number
            )
            continue
        virus_match = _VIRUS_TOTAL_RE.match(line)
        if virus_match:
            identifier = virus_match.group("identifier")
            if identifier in virus_totals:
                raise LegacyReconstructionError(
                    f"{path}:{line_number}: duplicate legacy virus {identifier!r}"
                )
            virus_totals[identifier] = _finite_float(
                virus_match.group("value"), source=path, line_number=line_number
            )
            continue
        overall_match = _OVERALL_TOTAL_RE.match(line)
        if overall_match:
            if total_viral_load is not None:
                raise LegacyReconstructionError(f"{path}:{line_number}: duplicate total viral load")
            total_viral_load = _finite_float(
                overall_match.group("value"), source=path, line_number=line_number
            )
            continue
        if line.startswith("No viral gene IDs found in this sample"):
            no_viral_genes = True
    return LegacySummary(
        gene_totals=gene_totals,
        virus_totals=virus_totals,
        total_viral_load=total_viral_load,
        no_viral_genes=no_viral_genes,
    )


def _comparison_rows(
    record_type: str,
    reconstructed: Mapping[str, float],
    reported: Mapping[str, float],
    *,
    tolerance: float,
) -> list[dict[str, Any]]:
    rows = []
    for identifier in sorted(set(reconstructed) | set(reported)):
        reconstructed_value = reconstructed.get(identifier)
        summary_value = reported.get(identifier)
        if reconstructed_value is None:
            status = "missing_in_reconstruction"
            delta = None
        elif summary_value is None:
            status = "missing_in_summary"
            delta = None
        else:
            delta = reconstructed_value - summary_value
            status = (
                "match"
                if math.isclose(
                    reconstructed_value,
                    summary_value,
                    rel_tol=0.0,
                    abs_tol=tolerance,
                )
                else "mismatch"
            )
        rows.append(
            {
                "record_type": record_type,
                "identifier": identifier,
                "reconstructed_value": reconstructed_value,
                "summary_value": summary_value,
                "delta": delta,
                "status": status,
            }
        )
    return rows


def reconstruct_legacy(
    h5ad_path: Path,
    summary_path: Path,
    viral_accessions_path: Path,
    *,
    virus_name_map: Optional[Mapping[str, str]] = None,
    tolerance: float = 1e-6,
) -> LegacyReconstruction:
    """Recreate v2.2 calls as ``counts_original + counts_corrected`` and audit them."""

    if tolerance < 0 or not math.isfinite(tolerance):
        raise LegacyReconstructionError(f"invalid absolute tolerance {tolerance!r}")
    import anndata as ad
    import numpy as np

    adata = ad.read_h5ad(h5ad_path)
    required_layers = ("counts_original", "counts_corrected")
    missing_layers = [name for name in required_layers if name not in adata.layers]
    if missing_layers:
        raise LegacyReconstructionError(f"{h5ad_path}: missing required layers {missing_layers}")
    combined = adata.layers["counts_original"] + adata.layers["counts_corrected"]
    if combined.shape != adata.shape:
        raise LegacyReconstructionError(
            f"{h5ad_path}: combined layer shape {combined.shape} != AnnData shape {adata.shape}"
        )
    totals = np.asarray(combined.sum(axis=0)).reshape(-1)
    if totals.shape != (adata.n_vars,) or not np.all(np.isfinite(totals)):
        raise LegacyReconstructionError(f"{h5ad_path}: invalid reconstructed gene totals")
    accessions = {
        line.strip()
        for line in viral_accessions_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    var_index = {str(gene): index for index, gene in enumerate(adata.var_names)}
    called_gene_totals = {
        gene: float(totals[var_index[gene]])
        for gene in sorted(accessions)
        if gene in var_index and float(totals[var_index[gene]]) > 1.0
    }

    if virus_name_map is None:
        from viralscan.constants import VIRUS_NAME_MAP

        virus_name_map = VIRUS_NAME_MAP
    reconstructed_virus_totals: dict[str, float] = {}
    for gene, count in called_gene_totals.items():
        virus = gene
        for prefix, candidate in virus_name_map.items():
            if prefix in gene:
                virus = candidate
                break
        reconstructed_virus_totals[virus] = reconstructed_virus_totals.get(virus, 0.0) + count
    reconstructed_total = float(sum(called_gene_totals.values()))
    summary = parse_legacy_summary(summary_path)
    comparisons = _comparison_rows(
        "gene",
        called_gene_totals,
        summary.gene_totals,
        tolerance=tolerance,
    )
    comparisons.extend(
        _comparison_rows(
            "virus",
            reconstructed_virus_totals,
            summary.virus_totals,
            tolerance=tolerance,
        )
    )
    if summary.total_viral_load is None:
        total_status = (
            "not_reported_no_calls"
            if summary.no_viral_genes and reconstructed_total == 0.0
            else "missing_in_summary"
        )
        total_delta = None
    else:
        total_delta = reconstructed_total - summary.total_viral_load
        total_status = (
            "match"
            if math.isclose(
                reconstructed_total,
                summary.total_viral_load,
                rel_tol=0.0,
                abs_tol=tolerance,
            )
            else "mismatch"
        )
    comparisons.append(
        {
            "record_type": "total",
            "identifier": "total_viral_load",
            "reconstructed_value": reconstructed_total,
            "summary_value": summary.total_viral_load,
            "delta": total_delta,
            "status": total_status,
        }
    )
    return LegacyReconstruction(
        called_gene_totals=called_gene_totals,
        reconstructed_virus_totals=reconstructed_virus_totals,
        reconstructed_total_viral_load=reconstructed_total,
        summary=summary,
        comparisons=tuple(comparisons),
    )


def _atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.with_suffix(path.suffix + ".tmp")
    staging.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staging.replace(path)


def _select_manifest_row(manifest_path: Path, run_id: str) -> dict[str, str]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        rows = [
            row for row in csv.DictReader(handle, delimiter="\t") if row.get("run_id") == run_id
        ]
    if len(rows) != 1:
        raise RunRowError(
            f"run ID {run_id!r} must select exactly one raw-manifest row; found {len(rows)}"
        )
    return rows[0]


def _manifest_hash_check(row: Mapping[str, str]) -> dict[str, Any]:
    records = []
    cache: dict[Path, tuple[int, str]] = {}
    for artifact in RUN_ROW_HASH_ARTIFACTS:
        path_value = row.get(f"{artifact}_path", "")
        expected_hash = row.get(f"{artifact}_sha256", "")
        expected_size_text = row.get(f"{artifact}_bytes", "")
        if not path_value or not expected_hash or not expected_size_text:
            records.append(
                {
                    "artifact": artifact,
                    "status": "manifest_field_missing",
                    "expected_sha256": expected_hash,
                    "observed_sha256": None,
                    "expected_bytes": expected_size_text or None,
                    "observed_bytes": None,
                }
            )
            continue
        path = Path(path_value)
        if not path.is_absolute():
            records.append(
                {
                    "artifact": artifact,
                    "status": "path_not_absolute",
                    "expected_sha256": expected_hash,
                    "observed_sha256": None,
                    "expected_bytes": expected_size_text,
                    "observed_bytes": None,
                }
            )
            continue
        if not path.is_file():
            records.append(
                {
                    "artifact": artifact,
                    "status": "missing",
                    "expected_sha256": expected_hash,
                    "observed_sha256": None,
                    "expected_bytes": expected_size_text,
                    "observed_bytes": None,
                }
            )
            continue
        try:
            observed_size, observed_hash = _sha256(path, cache)
        except OSError as error:
            records.append(
                {
                    "artifact": artifact,
                    "status": "read_error",
                    "expected_sha256": expected_hash,
                    "observed_sha256": None,
                    "expected_bytes": expected_size_text,
                    "observed_bytes": None,
                    "error": str(error),
                }
            )
            continue
        try:
            expected_size = int(expected_size_text)
        except ValueError:
            expected_size = None
        if expected_size is None:
            status = "invalid_expected_bytes"
        elif observed_size != expected_size:
            status = "size_mismatch"
        elif observed_hash != expected_hash:
            status = "sha256_mismatch"
        else:
            status = "match"
        records.append(
            {
                "artifact": artifact,
                "status": status,
                "expected_sha256": expected_hash,
                "observed_sha256": observed_hash,
                "expected_bytes": expected_size_text,
                "observed_bytes": observed_size,
            }
        )
    return {
        "schema_version": "1.0.0",
        "run_id": row.get("run_id"),
        "status": (
            "match" if all(record["status"] == "match" for record in records) else "mismatch"
        ),
        "artifacts": records,
    }


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _manifest_path(row: Mapping[str, str], field_name: str) -> Path:
    value = row.get(field_name, "")
    if not value:
        raise RunRowError(f"raw manifest is missing {field_name}")
    path = Path(value)
    if not path.is_absolute():
        raise RunRowError(f"{field_name} must be absolute in the raw manifest")
    return path.resolve()


def _validate_run_path_contract(
    row: Mapping[str, str],
    source_root: Path,
) -> Path:
    """Prove that hashed paths are exactly the paths legacy and v3 will consume."""

    source_root = source_root.expanduser().resolve()
    if not source_root.is_dir():
        raise RunRowError(f"archive source root is not a directory: {source_root}")
    config_path = _manifest_path(row, "config_path")
    sample_dir = config_path.parent
    if not _is_within(sample_dir, source_root):
        raise RunRowError(
            f"manifest sample directory {sample_dir} is outside source root {source_root}"
        )
    chemistry = row.get("chemistry", "")
    whitelist_name = {
        "10xv2": "10x_version2_whitelist.txt",
        "10xv3": "10x_version3_whitelist.txt",
    }.get(chemistry)
    if whitelist_name is None:
        raise RunRowError(f"unsupported manifest chemistry {chemistry!r}")
    kb_dir = sample_dir / "kb-python"
    counts_dir = kb_dir / "counts_unfiltered"
    expected_source_paths = {
        "config_path": sample_dir / "config.yaml",
        "summary_path": sample_dir / "summary.txt",
        "output_bus_path": kb_dir / "output.bus",
        "matrix_ec_path": kb_dir / "matrix.ec",
        "transcripts_path": kb_dir / "transcripts.txt",
        "whitelist_path": kb_dir / whitelist_name,
        "adata_multimap_path": counts_dir / "adata_multimap.h5ad",
        "barcodes_path": counts_dir / "cells_x_genes.barcodes.txt",
        "genes_path": counts_dir / "cells_x_genes.genes.txt",
        "gene_names_path": counts_dir / "cells_x_genes.genes.names.txt",
        "adata_path": counts_dir / "adata.h5ad",
        "run_info_path": kb_dir / "run_info.json",
        "analysis_path": sample_dir / "log" / "analysis.txt",
    }
    for field_name, expected in expected_source_paths.items():
        observed = _manifest_path(row, field_name)
        expected = expected.resolve()
        if observed != expected:
            raise RunRowError(
                f"{field_name} does not match the executed path: {observed} != {expected}"
            )
        if not _is_within(observed, source_root):
            raise RunRowError(f"{field_name} is outside source root {source_root}")

    config = _read_config(config_path)
    if config.get("technology") != chemistry:
        raise RunRowError(
            "manifest chemistry does not match config technology: "
            f"{chemistry!r} != {config.get('technology')!r}"
        )
    expected_reference_paths = {
        "index_path": _required_config_path(config, "index", config_path),
        "t2g_path": _required_config_path(config, "transcripts", config_path),
    }
    expected_reference_paths["transcriptome_path"] = (
        expected_reference_paths["index_path"].parent / "transcriptome.fa"
    )
    for field_name, expected in expected_reference_paths.items():
        observed = _manifest_path(row, field_name)
        expected = expected.resolve()
        if observed != expected:
            raise RunRowError(
                f"{field_name} does not match the config-derived path: {observed} != {expected}"
            )
    return sample_dir


def _write_legacy_outputs(
    output_dir: Path,
    row: Mapping[str, str],
    result: LegacyReconstruction,
) -> dict[str, Any]:
    audit_rows = [
        {
            "run_id": row["run_id"],
            "logical_id": row.get("logical_id", ""),
            **comparison,
        }
        for comparison in result.comparisons
    ]
    audit_fields = [
        "run_id",
        "logical_id",
        "record_type",
        "identifier",
        "reconstructed_value",
        "summary_value",
        "delta",
        "status",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "legacy_audit.tsv", audit_rows, audit_fields)
    payload = {
        "schema_version": "1.0.0",
        "run_id": row["run_id"],
        "logical_id": row.get("logical_id", ""),
        "called_gene_totals": dict(result.called_gene_totals),
        "reconstructed_virus_totals": dict(result.reconstructed_virus_totals),
        "reconstructed_total_viral_load": result.reconstructed_total_viral_load,
        "summary": {
            "gene_totals": dict(result.summary.gene_totals),
            "virus_totals": dict(result.summary.virus_totals),
            "total_viral_load": result.summary.total_viral_load,
            "no_viral_genes": result.summary.no_viral_genes,
        },
        "audit_status_counts": dict(
            Counter(comparison["status"] for comparison in result.comparisons)
        ),
    }
    _atomic_json(output_dir / "legacy_reconstruction.json", payload)
    _atomic_json(
        output_dir / "status.json",
        {
            "schema_version": "1.0.0",
            "run_id": row["run_id"],
            "status": "success",
            "audit_status_counts": payload["audit_status_counts"],
        },
    )
    return payload


def run_row(
    raw_manifest: Path,
    run_id: str,
    output_root: Path,
    *,
    source_root: Path,
    threads: int = 8,
) -> dict[str, Any]:
    """Run the frozen legacy and v3 comparisons for exactly one manifest row."""

    safe_run_id = _slug(run_id)
    if safe_run_id != run_id:
        raise RunRowError(f"run ID is not a sanitized deterministic ID: {run_id!r}")
    if threads < 1:
        raise RunRowError("threads must be at least one")
    row = _select_manifest_row(raw_manifest, run_id)
    sample_dir = _validate_run_path_contract(row, source_root)
    resolved_source_root = source_root.expanduser().resolve()
    row_output = (output_root.expanduser().resolve() / run_id).resolve()
    if _is_within(row_output, resolved_source_root):
        raise RunRowError(
            "output directory must be outside the read-only archive source root "
            f"{resolved_source_root}"
        )
    if row_output.exists():
        raise RunRowError(
            f"refusing to reuse existing row output; choose a clean named output root: {row_output}"
        )
    row_output.mkdir(parents=True)
    stage = "verify_before"
    errors: list[dict[str, str]] = []
    legacy_payload: Optional[dict[str, Any]] = None
    v3_argv: Optional[list[str]] = None
    try:
        before = _manifest_hash_check(row)
        _atomic_json(row_output / "hash_check_before.json", before)
        if before["status"] != "match":
            raise RunRowError("one or more current inputs do not match the raw manifest")

        stage = "legacy_reconstruction"
        legacy_result = reconstruct_legacy(
            Path(row["adata_multimap_path"]),
            Path(row["summary_path"]),
            Path(row["analysis_path"]),
        )
        legacy_payload = _write_legacy_outputs(row_output / "legacy", row, legacy_result)

        stage = "v3_benchmark"
        from scripts import benchmark_v3_multimap

        v3_argv = [
            "--sample-id",
            run_id,
            "--sample-dir",
            str(sample_dir),
            "--whitelist",
            row["whitelist_path"],
            "--t2g",
            row["t2g_path"],
            "--method",
            "host-conservative",
            "--output",
            str(row_output / "v3" / "result.json"),
            "--threads",
            str(threads),
        ]
        return_code = benchmark_v3_multimap.main(v3_argv)
        if return_code != 0:
            raise RunRowError(f"v3 benchmark returned non-zero status {return_code}")
    except Exception as error:
        errors.append(
            {
                "stage": stage,
                "error_type": type(error).__name__,
                "message": str(error),
            }
        )
    finally:
        after = _manifest_hash_check(row)
        _atomic_json(row_output / "hash_check_after.json", after)
        if after["status"] != "match":
            errors.append(
                {
                    "stage": "verify_after",
                    "error_type": "RunRowError",
                    "message": "one or more source inputs changed or no longer match",
                }
            )

    status = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "status": "failed" if errors else "success",
        "threads": threads,
        "legacy_audit_status_counts": (
            legacy_payload["audit_status_counts"] if legacy_payload else {}
        ),
        "v3_argv": v3_argv,
        "errors": errors,
    }
    _atomic_json(row_output / "status.json", status)
    if errors:
        _atomic_json(row_output / "failure.json", status)
    else:
        (row_output / "failure.json").unlink(missing_ok=True)
    return status


def _planned_rows(manifest_path: Path) -> list[dict[str, str]]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    run_ids = [row.get("run_id", "") for row in rows]
    if not rows or any(not run_id for run_id in run_ids):
        raise RunRowError("planned manifest must contain non-empty run_id values")
    duplicates = sorted(run_id for run_id, count in Counter(run_ids).items() if count > 1)
    if duplicates:
        raise RunRowError(f"planned manifest contains duplicate run IDs: {duplicates}")
    return sorted(rows, key=lambda row: row["run_id"])


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RunRowError(f"cannot read JSON object {path}: {error}") from error
    if not isinstance(payload, dict):
        raise RunRowError(f"JSON payload is not an object: {path}")
    return payload


def _sanitize_aggregate_text(value: Any) -> str:
    """Remove absolute filesystem paths from small aggregate failure tables."""

    return re.sub(r"(?<!\S)/(?:[^\s,;]+)", "<absolute-path>", str(value))


def summarize_results(
    raw_manifest: Path,
    run_root: Path,
    output_dir: Path,
) -> dict[str, str]:
    """Aggregate every planned row, including explicit and implicit failures."""

    planned = _planned_rows(raw_manifest)
    run_rows: list[dict[str, Any]] = []
    virus_rows: list[dict[str, Any]] = []
    legacy_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    for plan in planned:
        run_id = plan["run_id"]
        logical_id = plan.get("logical_id", "")
        row_dir = run_root / run_id
        status_path = row_dir / "status.json"
        status: dict[str, Any]
        if not status_path.is_file():
            status = {
                "status": "missing",
                "errors": [
                    {
                        "stage": "summarize",
                        "error_type": "MissingOutput",
                        "message": "planned row has no status.json",
                    }
                ],
            }
        else:
            try:
                status = _read_json_object(status_path)
            except RunRowError as error:
                status = {
                    "status": "invalid",
                    "errors": [
                        {
                            "stage": "summarize",
                            "error_type": type(error).__name__,
                            "message": str(error),
                        }
                    ],
                }
        row_status = str(status.get("status", "invalid"))
        legacy_total: Any = ""
        unique_mass: Any = ""
        ambiguous_mass: Any = ""
        selected_mass: Any = ""
        if row_status == "success":
            try:
                legacy_payload = _read_json_object(
                    row_dir / "legacy" / "legacy_reconstruction.json"
                )
                v3_payload = _read_json_object(row_dir / "v3" / "result.json")
                legacy_total = legacy_payload.get("reconstructed_total_viral_load", "")
                unique_mass = v3_payload.get("unique_molecule_mass", "")
                ambiguous_mass = v3_payload.get("allocated_ambiguous_mass", "")
                selected_mass = v3_payload.get("selected_matrix_mass", "")
            except RunRowError as error:
                row_status = "invalid"
                status = {
                    "status": row_status,
                    "errors": [
                        {
                            "stage": "summarize",
                            "error_type": type(error).__name__,
                            "message": str(error),
                        }
                    ],
                }
        run_rows.append(
            {
                "run_id": run_id,
                "logical_id": logical_id,
                "sample_class": plan.get("sample_class", ""),
                "expected_target": plan.get("expected_target", ""),
                "status": row_status,
                "legacy_total_viral_load": legacy_total,
                "v3_unique_molecule_mass": unique_mass,
                "v3_allocated_ambiguous_mass": ambiguous_mass,
                "v3_selected_matrix_mass": selected_mass,
            }
        )

        legacy_path = row_dir / "legacy" / "legacy_audit.tsv"
        appended_legacy = False
        if legacy_path.is_file():
            try:
                with legacy_path.open(newline="", encoding="utf-8") as handle:
                    for legacy in csv.DictReader(handle, delimiter="\t"):
                        legacy_rows.append(
                            {
                                "run_id": run_id,
                                "logical_id": logical_id,
                                "row_status": row_status,
                                "record_type": legacy.get("record_type", ""),
                                "identifier": legacy.get("identifier", ""),
                                "reconstructed_value": legacy.get("reconstructed_value", ""),
                                "summary_value": legacy.get("summary_value", ""),
                                "delta": legacy.get("delta", ""),
                                "audit_status": legacy.get("status", ""),
                            }
                        )
                        appended_legacy = True
            except OSError:
                appended_legacy = False
        if not appended_legacy:
            legacy_rows.append(
                {
                    "run_id": run_id,
                    "logical_id": logical_id,
                    "row_status": row_status,
                    "record_type": "",
                    "identifier": "",
                    "reconstructed_value": "",
                    "summary_value": "",
                    "delta": "",
                    "audit_status": "not_available",
                }
            )

        per_virus_path = row_dir / "v3" / "per_virus.tsv"
        appended_virus = False
        if per_virus_path.is_file():
            try:
                with per_virus_path.open(newline="", encoding="utf-8") as handle:
                    for virus in csv.DictReader(handle, delimiter="\t"):
                        virus_rows.append(
                            {
                                "run_id": run_id,
                                "logical_id": logical_id,
                                "row_status": row_status,
                                "virus_name": virus.get("virus_name", ""),
                                "n_features": virus.get("n_features", ""),
                                "v3_unique": virus.get("v3_unique", ""),
                                "v3_equal": virus.get("v3_equal", ""),
                                "v3_host_conservative": virus.get("v3_host_conservative", ""),
                            }
                        )
                        appended_virus = True
            except OSError:
                appended_virus = False
        if not appended_virus:
            virus_rows.append(
                {
                    "run_id": run_id,
                    "logical_id": logical_id,
                    "row_status": row_status,
                    "virus_name": "",
                    "n_features": "",
                    "v3_unique": "",
                    "v3_equal": "",
                    "v3_host_conservative": "",
                }
            )

        if row_status != "success":
            errors = status.get("errors")
            if not isinstance(errors, list) or not errors:
                errors = [
                    {
                        "stage": row_status,
                        "error_type": row_status,
                        "message": str(status.get("reason", "no failure reason recorded")),
                    }
                ]
            for error in errors:
                failure_rows.append(
                    {
                        "run_id": run_id,
                        "logical_id": logical_id,
                        "status": row_status,
                        "stage": error.get("stage", ""),
                        "error_type": error.get("error_type", ""),
                        "message": _sanitize_aggregate_text(error.get("message", "")),
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "run_metrics": output_dir / "run_metrics.tsv",
        "virus_metrics": output_dir / "virus_metrics.tsv",
        "legacy_reproduction": output_dir / "legacy_reproduction.tsv",
        "failures": output_dir / "failures.tsv",
    }
    _write_tsv(paths["run_metrics"], run_rows, list(run_rows[0]))
    _write_tsv(paths["virus_metrics"], virus_rows, list(virus_rows[0]))
    _write_tsv(paths["legacy_reproduction"], legacy_rows, list(legacy_rows[0]))
    failure_fields = [
        "run_id",
        "logical_id",
        "status",
        "stage",
        "error_type",
        "message",
    ]
    _write_tsv(paths["failures"], failure_rows, failure_fields)
    return {name: str(path) for name, path in paths.items()}


def _add_validation_error(
    errors: list[dict[str, str]],
    run_id: str,
    check: str,
    message: str,
) -> None:
    errors.append({"run_id": run_id, "check": check, "message": message})


def _finite_nonnegative_number(value: Any) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number >= 0


def _validate_v3_success(
    run_id: str,
    row_dir: Path,
    errors: list[dict[str, str]],
) -> None:
    import anndata as ad
    import h5py
    import numpy as np

    v3_dir = row_dir / "v3"
    required = [
        row_dir / "hash_check_before.json",
        row_dir / "hash_check_after.json",
        row_dir / "legacy" / "legacy_audit.tsv",
        row_dir / "legacy" / "legacy_reconstruction.json",
        v3_dir / "result.json",
        v3_dir / "status.json",
        v3_dir / "hashes.json",
        v3_dir / "count_audit.json",
        v3_dir / "adata_v3.h5ad",
        v3_dir / "per_cell.tsv",
        v3_dir / "per_virus.tsv",
    ]
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        _add_validation_error(
            errors, run_id, "required_outputs", f"missing outputs: {sorted(missing)}"
        )
        return
    for name in ("hash_check_before.json", "hash_check_after.json"):
        payload = _read_json_object(row_dir / name)
        if payload.get("status") != "match":
            _add_validation_error(errors, run_id, "source_hashes", f"{name} is not match")
    v3_status = _read_json_object(v3_dir / "status.json")
    if v3_status.get("status") != "success":
        _add_validation_error(errors, run_id, "v3_status", "v3 status is not success")
    result = _read_json_object(v3_dir / "result.json")
    if result.get("sample_id") != run_id:
        _add_validation_error(
            errors, run_id, "result_schema", "result sample_id does not match run ID"
        )
    audit = result.get("audit")
    if not isinstance(audit, dict):
        _add_validation_error(errors, run_id, "count_audit", "result audit is missing")
    else:
        fields = (
            "input_molecules",
            "unique_molecules",
            "ambiguous_molecules",
            "unresolved_molecules",
        )
        if not all(_finite_nonnegative_number(audit.get(field)) for field in fields):
            _add_validation_error(errors, run_id, "count_audit", "audit values are invalid")
        else:
            expected_input = sum(
                float(audit[field])
                for field in (
                    "unique_molecules",
                    "ambiguous_molecules",
                    "unresolved_molecules",
                )
            )
            if not math.isclose(
                float(audit["input_molecules"]),
                expected_input,
                rel_tol=0.0,
                abs_tol=1e-6,
            ):
                _add_validation_error(
                    errors,
                    run_id,
                    "molecule_conservation",
                    "input != unique + ambiguous + unresolved",
                )
    count_audit = _read_json_object(v3_dir / "count_audit.json")
    if isinstance(audit, dict):
        for field_name in (
            "input_molecules",
            "unique_molecules",
            "ambiguous_molecules",
            "unresolved_molecules",
        ):
            if count_audit.get(field_name) != audit.get(field_name):
                _add_validation_error(
                    errors,
                    run_id,
                    "count_audit",
                    f"count_audit.json disagrees on {field_name}",
                )
    endpoints = result.get("comparison_endpoints")
    if not isinstance(endpoints, dict) or set(endpoints) != {
        "v3-unique",
        "v3-equal",
        "v3-host-conservative",
    }:
        _add_validation_error(
            errors,
            run_id,
            "result_schema",
            "frozen comparison endpoints are missing or unexpected",
        )
    mass_fields = (
        "unique_molecule_mass",
        "allocated_ambiguous_mass",
        "selected_matrix_mass",
    )
    if not all(_finite_nonnegative_number(result.get(field)) for field in mass_fields):
        _add_validation_error(
            errors, run_id, "matrix_mass", "selected matrix mass values are invalid"
        )
    elif not math.isclose(
        float(result["selected_matrix_mass"]),
        float(result["unique_molecule_mass"]) + float(result["allocated_ambiguous_mass"]),
        rel_tol=0.0,
        abs_tol=1e-6,
    ):
        _add_validation_error(
            errors,
            run_id,
            "matrix_conservation",
            "selected mass != unique + allocated ambiguous",
        )

    h5ad_path = v3_dir / "adata_v3.h5ad"
    declared_shape = (int(result.get("n_cells", -1)), int(result.get("n_genes", -1)))
    with h5py.File(h5ad_path, "r") as handle:
        nodes = {"X": handle.get("X")}
        layers = handle.get("layers")
        for layer in (
            "counts_unique",
            "counts_ambiguous_allocated",
            "counts_multimap_equal",
            "counts_multimap_host_conservative",
        ):
            nodes[layer] = layers.get(layer) if layers is not None else None
        for name, node in nodes.items():
            if node is None:
                _add_validation_error(
                    errors, run_id, "matrix_schema", f"missing H5AD matrix {name}"
                )
                continue
            shape = (
                tuple(int(value) for value in node.attrs["shape"])
                if isinstance(node, h5py.Group)
                else tuple(int(value) for value in node.shape)
            )
            if shape != declared_shape:
                _add_validation_error(
                    errors,
                    run_id,
                    "matrix_alignment",
                    f"{name} shape {shape} != declared {declared_shape}",
                )
            values = node["data"][:] if isinstance(node, h5py.Group) else node[:]
            if not np.isfinite(values).all() or (values < 0).any():
                _add_validation_error(
                    errors,
                    run_id,
                    "matrix_values",
                    f"{name} contains non-finite or negative values",
                )
    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        expected_barcodes = [str(value) for value in backed.obs_names]
    finally:
        backed.file.close()
    with (v3_dir / "per_cell.tsv").open(newline="", encoding="utf-8") as handle:
        per_cell = list(csv.DictReader(handle, delimiter="\t"))
    observed_barcodes = [row.get("barcode", "") for row in per_cell]
    if observed_barcodes != expected_barcodes:
        _add_validation_error(
            errors,
            run_id,
            "barcode_alignment",
            "per-cell barcodes do not match H5AD order",
        )
    for row in per_cell:
        for metric_name in (
            "v3_unique_viral",
            "v3_equal_viral",
            "v3_host_conservative_viral",
        ):
            if not _finite_nonnegative_number(row.get(metric_name)):
                _add_validation_error(
                    errors,
                    run_id,
                    "per_cell_values",
                    f"invalid {metric_name}",
                )
                break
    with (v3_dir / "per_virus.tsv").open(newline="", encoding="utf-8") as handle:
        per_virus = list(csv.DictReader(handle, delimiter="\t"))
    for row in per_virus:
        for metric_name in (
            "n_features",
            "v3_unique",
            "v3_equal",
            "v3_host_conservative",
        ):
            if not _finite_nonnegative_number(row.get(metric_name)):
                _add_validation_error(
                    errors,
                    run_id,
                    "per_virus_values",
                    f"invalid {metric_name}",
                )
                break

    hashes = _read_json_object(v3_dir / "hashes.json")
    output_hashes = hashes.get("outputs")
    if not isinstance(output_hashes, dict):
        _add_validation_error(errors, run_id, "output_hashes", "hashes.json outputs are missing")
    else:
        cache: dict[Path, tuple[int, str]] = {}
        for name, record in sorted(output_hashes.items()):
            path = v3_dir / name
            if not path.is_file() or not isinstance(record, dict):
                _add_validation_error(
                    errors, run_id, "output_hashes", f"missing hashed output {name}"
                )
                continue
            size, digest = _sha256(path, cache)
            if size != record.get("size_bytes") or digest != record.get("sha256"):
                _add_validation_error(errors, run_id, "output_hashes", f"hash mismatch for {name}")
    with (row_dir / "legacy" / "legacy_audit.tsv").open(newline="", encoding="utf-8") as handle:
        legacy_reader = csv.DictReader(handle, delimiter="\t")
        required_legacy = {
            "record_type",
            "identifier",
            "reconstructed_value",
            "summary_value",
            "delta",
            "status",
        }
        if not required_legacy.issubset(legacy_reader.fieldnames or []):
            _add_validation_error(
                errors, run_id, "legacy_schema", "legacy audit columns are incomplete"
            )
        else:
            for row in legacy_reader:
                for field in ("reconstructed_value", "summary_value", "delta"):
                    value = row.get(field, "")
                    if value and not math.isfinite(float(value)):
                        _add_validation_error(
                            errors,
                            run_id,
                            "legacy_values",
                            f"non-finite legacy {field}",
                        )


def validate_results(raw_manifest: Path, run_root: Path) -> dict[str, Any]:
    """Validate all planned outputs without excluding failed or unsupported rows."""

    planned = _planned_rows(raw_manifest)
    errors: list[dict[str, str]] = []
    status_counts: Counter[str] = Counter()
    planned_ids = {row["run_id"] for row in planned}
    if run_root.is_dir():
        extra = sorted(
            path.name
            for path in run_root.iterdir()
            if path.is_dir() and path.name not in planned_ids
        )
        if extra:
            _add_validation_error(
                errors, "", "unplanned_outputs", f"unplanned run directories: {extra}"
            )
    for plan in planned:
        run_id = plan["run_id"]
        row_dir = run_root / run_id
        status_path = row_dir / "status.json"
        if not status_path.is_file():
            status_counts["missing"] += 1
            _add_validation_error(errors, run_id, "planned_row", "planned row has no status.json")
            continue
        try:
            status = _read_json_object(status_path)
        except RunRowError as error:
            status_counts["invalid"] += 1
            _add_validation_error(errors, run_id, "status_schema", str(error))
            continue
        row_status = str(status.get("status", "invalid"))
        status_counts[row_status] += 1
        if row_status == "success":
            if status.get("run_id") != run_id:
                _add_validation_error(errors, run_id, "status_schema", "wrapper run_id mismatch")
            try:
                _validate_v3_success(run_id, row_dir, errors)
            except (OSError, ValueError, RunRowError) as error:
                _add_validation_error(errors, run_id, "success_validation", str(error))
        elif row_status == "failed":
            if not (row_dir / "failure.json").is_file():
                _add_validation_error(
                    errors, run_id, "failure_schema", "failed row lacks failure.json"
                )
            if not isinstance(status.get("errors"), list) or not status["errors"]:
                _add_validation_error(errors, run_id, "failure_schema", "failed row lacks errors")
        elif row_status == "unsupported":
            if not status.get("reason"):
                _add_validation_error(
                    errors,
                    run_id,
                    "unsupported_schema",
                    "unsupported row lacks reason",
                )
        else:
            _add_validation_error(
                errors, run_id, "status_schema", f"unsupported status {row_status!r}"
            )
    return {
        "schema_version": "1.0.0",
        "status": "invalid" if errors else "valid",
        "planned_rows": len(planned),
        "status_counts": dict(sorted(status_counts.items())),
        "errors": errors,
    }


def _load_versions(kb_dir: Path) -> tuple[str, str, str]:
    kb_info_path = kb_dir / "kb_info.json"
    if not kb_info_path.is_file():
        return "", "", ""
    try:
        info = json.loads(kb_info_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InventoryContractError(f"cannot parse {kb_info_path}: {exc}") from exc
    kallisto = info.get("kallisto") or {}
    bustools = info.get("bustools") or {}
    return (
        str(info.get("version", "")),
        str(kallisto.get("version", "")),
        str(bustools.get("version", "")),
    )


def inventory(
    source_root: Path,
    raw_manifest: Path,
    cohort_manifest: Path,
    *,
    contract: Optional[InventoryContract] = None,
) -> list[dict[str, Any]]:
    """Inventory archived result trees and emit raw and sanitized TSV manifests."""

    contract = contract or InventoryContract()
    source_root = source_root.expanduser().resolve()
    if not source_root.is_dir():
        raise InventoryContractError(f"source root is not a directory: {source_root}")

    candidates: list[dict[str, Any]] = []
    for bus_path in sorted(source_root.glob("**/kb-python/output.bus")):
        kb_dir = bus_path.parent
        run_dir = kb_dir.parent
        config_path = _require_file(run_dir / "config.yaml", "config")
        config = _read_config(config_path)
        sample1 = config.get("sample1")
        sample2 = config.get("sample2")
        chemistry = config.get("technology")
        if not isinstance(sample1, str) or not isinstance(sample2, str):
            raise InventoryContractError(f"{config_path}: missing paired sample paths")
        if not isinstance(chemistry, str):
            raise InventoryContractError(f"{config_path}: missing technology")
        relative_run = run_dir.relative_to(source_root)
        candidates.append(
            {
                "run_id": _run_id(relative_run),
                "logical_id": _logical_id(sample1, sample2, config_path),
                "chemistry": chemistry,
                "relative_run": relative_run,
                "run_dir": run_dir,
                "kb_dir": kb_dir,
                "config": config,
                "sample1": sample1,
                "sample2": sample2,
            }
        )

    run_ids = [row["run_id"] for row in candidates]
    if len(set(run_ids)) != len(run_ids):
        raise InventoryContractError("sanitized run IDs are not unique")
    logical_counts = Counter(row["logical_id"] for row in candidates)
    chemistry_counts = Counter(row["chemistry"] for row in candidates)
    errors = []
    if len(candidates) != contract.technical_rows:
        errors.append(
            f"technical rows: expected {contract.technical_rows}, found {len(candidates)}"
        )
    if len(logical_counts) != contract.logical_inputs:
        errors.append(
            f"logical inputs: expected {contract.logical_inputs}, found {len(logical_counts)}"
        )
    for chemistry, expected in contract.chemistry_counts.items():
        found = chemistry_counts.get(chemistry, 0)
        if found != expected:
            errors.append(f"{chemistry} rows: expected {expected}, found {found}")
    unexpected_chemistries = set(chemistry_counts) - set(contract.chemistry_counts)
    if unexpected_chemistries:
        errors.append(f"unexpected chemistries: {sorted(unexpected_chemistries)}")
    if errors:
        raise InventoryContractError("frozen cohort mismatch: " + "; ".join(errors))

    digest_cache: dict[Path, tuple[int, str]] = {}
    sanitized_rows: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    for candidate in candidates:
        run_dir = candidate["run_dir"]
        kb_dir = candidate["kb_dir"]
        config = candidate["config"]
        index_path = _require_file(
            _required_config_path(config, "index", run_dir / "config.yaml"), "index"
        )
        t2g_path = _require_file(
            _required_config_path(config, "transcripts", run_dir / "config.yaml"), "t2g"
        )
        transcriptome_path = _require_file(index_path.parent / "transcriptome.fa", "transcriptome")
        paths = {
            "config": _require_file(run_dir / "config.yaml", "config"),
            "summary": _require_file(run_dir / "summary.txt", "summary"),
            "output_bus": _require_file(kb_dir / "output.bus", "BUS"),
            "matrix_ec": _require_file(kb_dir / "matrix.ec", "EC"),
            "transcripts": _require_file(kb_dir / "transcripts.txt", "transcript list"),
            "whitelist": _find_whitelist(kb_dir, candidate["chemistry"]),
            "adata_multimap": _require_file(
                kb_dir / "counts_unfiltered" / "adata_multimap.h5ad", "multimap H5AD"
            ),
            "barcodes": _require_file(
                kb_dir / "counts_unfiltered" / "cells_x_genes.barcodes.txt",
                "barcode list",
            ),
            "genes": _require_file(
                kb_dir / "counts_unfiltered" / "cells_x_genes.genes.txt",
                "gene list",
            ),
            "gene_names": _require_file(
                kb_dir / "counts_unfiltered" / "cells_x_genes.genes.names.txt",
                "gene-name list",
            ),
            "adata": _require_file(kb_dir / "counts_unfiltered" / "adata.h5ad", "original H5AD"),
            "run_info": _require_file(kb_dir / "run_info.json", "run metadata"),
            "analysis": _require_file(run_dir / "log" / "analysis.txt", "viral accession list"),
            "index": index_path,
            "t2g": t2g_path,
            "transcriptome": transcriptome_path,
        }
        hashes = {name: _sha256(path, digest_cache) for name, path in paths.items()}
        sample_class, expected_target = _classification(candidate["logical_id"])
        repeat_group = (
            f"repeat__{candidate['logical_id']}"
            if logical_counts[candidate["logical_id"]] > 1
            else ""
        )
        kb_version, kallisto_version, bustools_version = _load_versions(kb_dir)
        row: dict[str, Any] = {
            "run_id": candidate["run_id"],
            "logical_id": candidate["logical_id"],
            "technical_repeat_group": repeat_group,
            "chemistry": candidate["chemistry"],
            "sample_class": sample_class,
            "expected_target": expected_target,
            "sample1_name": Path(candidate["sample1"]).name,
            "sample2_name": Path(candidate["sample2"]).name,
            "kb_python_version": kb_version,
            "kallisto_version": kallisto_version,
            "bustools_version": bustools_version,
        }
        for name in (
            "config",
            "summary",
            "output_bus",
            "matrix_ec",
            "transcripts",
            "whitelist",
            "adata_multimap",
            "barcodes",
            "genes",
            "gene_names",
            "adata",
            "run_info",
            "analysis",
        ):
            row[f"{name}_path"] = paths[name].relative_to(source_root).as_posix()
        for name, (size, sha256) in hashes.items():
            row[f"{name}_bytes"] = size
            row[f"{name}_sha256"] = sha256
        row["reference_hash_id"] = _reference_hash_id(
            hashes[name][1] for name in ("index", "t2g", "transcriptome")
        )
        sanitized_rows.append(row)

        raw_row = dict(row)
        raw_row["sample1_path"] = candidate["sample1"]
        raw_row["sample2_path"] = candidate["sample2"]
        for name, path in paths.items():
            raw_row[f"{name}_path"] = str(path)
        raw_rows.append(raw_row)

    sanitized_rows.sort(key=lambda row: row["run_id"])
    raw_rows.sort(key=lambda row: row["run_id"])
    _write_tsv(cohort_manifest, sanitized_rows, SANITIZED_FIELDS)
    raw_fields = list(dict.fromkeys([*SANITIZED_FIELDS, *RAW_PATH_FIELDS]))
    _write_tsv(raw_manifest, raw_rows, raw_fields)
    return sanitized_rows


FRESH_SAMPLES = ("SRR12682296", "SRR12682297", "SRR12682298", "SRR6825024", "SRR6825025")
FRESH_V2_TOLERANCE = 1.0e-6  # protocol comparison.endpoints legacy-v2-reconstructed


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _fresh_record(packet_roots: Sequence[Path], task_id: str) -> tuple[dict[str, Any] | None, str]:
    """Return (record, source). Later packets win; revalidation supersedes original status."""

    for root in reversed(packet_roots):
        for kind in ("revalidation", "status"):
            path = root / kind / f"{task_id}.json"
            if path.is_file():
                return json.loads(path.read_text(encoding="utf-8")), f"{root.name}/{kind}"
    return None, "missing"


def _nested_sample_root(output: Path) -> Path:
    """Directory holding config.yaml; the CLI nests results one level below ``-o``."""

    if (output / "config.yaml").is_file():
        return output
    children = sorted(c for c in output.iterdir() if (c / "config.yaml").is_file())
    if len(children) != 1:
        raise FileNotFoundError(f"no unique config.yaml root beneath {output}")
    return children[0]


def _delta_rows(
    sample: str,
    stack: str,
    record_type: str,
    fresh: Mapping[str, float],
    archived: Mapping[str, float],
    tolerance: float,
    *,
    absent_as_zero: bool = False,
) -> list[dict[str, Any]]:
    """Per-identifier comparison. ``absent_as_zero`` is for tables that list only
    non-zero calls on one side (fresh viral_summary vs the archived all-virus table)."""

    rows: list[dict[str, Any]] = []
    for identifier in sorted(set(fresh) | set(archived)):
        f_val, a_val = fresh.get(identifier), archived.get(identifier)
        if absent_as_zero:
            f_val, a_val = f_val or 0.0, a_val or 0.0
        if f_val is None or a_val is None:
            verdict = "fresh_only" if a_val is None else "archive_only"
            delta = None
        else:
            delta = f_val - a_val
            verdict = "match" if abs(delta) <= tolerance else "mismatch"
        rows.append(
            {
                "sample_id": sample,
                "stack": stack,
                "record_type": record_type,
                "identifier": identifier,
                "fresh_value": f_val,
                "archived_value": a_val,
                "delta": delta,
                "tolerance": tolerance,
                "verdict": verdict,
                # LVC-11 requires an analyst to classify each mismatch (input, dependency,
                # nondeterminism, reference, package, unresolved); never guessed here.
                "classification": "unclassified" if verdict != "match" else "",
            }
        )
    return rows


def _fresh_rows(
    stack: str,
    record: Mapping[str, Any],
    archive_id: str,
    sample: str,
    archive: Mapping[str, list[dict[str, str]]],
    v3_tolerance: float,
) -> list[dict[str, Any]]:
    """Comparison records for one successful fresh row; raises on unreadable outputs."""

    legacy, run_metrics, virus_metrics = (
        archive["legacy"],
        archive["run_metrics"],
        archive["virus_metrics"],
    )
    if stack == "v3":
        root = _nested_sample_root(Path(record["output"]))
        summary = _read_tsv(root / "results" / "viral_summary.tsv")
        fresh = {r["virus_name"]: float(r["viral_molecules_total_est"]) for r in summary}
        archived = {
            r["virus_name"]: float(r["v3_host_conservative"])
            for r in virus_metrics
            if r["run_id"] == archive_id
        }
        return _delta_rows(
            sample, stack, "virus", fresh, archived, v3_tolerance, absent_as_zero=True
        )
    parsed = parse_legacy_summary(Path(record["v2_output_root"]) / "summary.txt")

    def archived_of(kind: str) -> dict[str, float]:
        return {
            r["identifier"]: float(r["summary_value"])
            for r in legacy
            if r["run_id"] == archive_id and r["record_type"] == kind
        }

    arch_total = {
        r["run_id"]: float(r["legacy_total_viral_load"])
        for r in run_metrics
        if r["run_id"] == archive_id
    }
    fresh_total = {archive_id: parsed.total_viral_load} if parsed.total_viral_load is not None else {}
    tol = FRESH_V2_TOLERANCE
    return [
        *_delta_rows(sample, stack, "gene", parsed.gene_totals, archived_of("gene"), tol),
        *_delta_rows(sample, stack, "virus", parsed.virus_totals, archived_of("virus"), tol),
        *_delta_rows(sample, stack, "total_viral_load", fresh_total, arch_total, tol),
    ]


def compare_fresh_vs_archive(
    packet_roots: Sequence[Path],
    archive_dir: Path,
    *,
    v3_tolerance: float = 0.0,
) -> dict[str, Any]:
    """Compare fresh v2/v3 control outputs with the archived tracked tables.

    ``packet_roots`` are searched later-wins, so attempt-3 rows override attempt-2
    rows and revalidation records override the original status. Archived side:
    ``legacy_reproduction.tsv`` / ``run_metrics.tsv`` (v2 ``legacy-v2-reconstructed``)
    and ``virus_metrics.tsv`` (v3 host-conservative product endpoint). v3-unique and
    v3-equal are not compared (the fresh run exports them only as AnnData layers).
    Failed, missing, or unreadable rows stay in the output with their reason.
    """

    archive = {
        "legacy": _read_tsv(archive_dir / "legacy_reproduction.tsv"),
        "run_metrics": _read_tsv(archive_dir / "run_metrics.tsv"),
        "virus_metrics": _read_tsv(archive_dir / "virus_metrics.tsv"),
    }
    rows: list[dict[str, Any]] = []
    run_status: list[dict[str, Any]] = []
    for sample in FRESH_SAMPLES:
        for stack in ("v2", "v3"):
            task_id = f"{stack}__{sample}"
            record, source = _fresh_record(packet_roots, task_id)
            entry: dict[str, Any] = {"task_id": task_id, "record_source": source}
            run_status.append(entry)
            if record is None or record.get("status") != "success":
                entry.update(
                    status="not_compared",
                    reason="no record" if record is None else f"row status {record.get('status')}",
                    exit_code=None if record is None else record.get("exit_code"),
                )
                continue
            try:
                rows += _fresh_rows(stack, record, f"{sample}__{sample}", sample, archive, v3_tolerance)
                entry["status"] = "compared"
            except (OSError, ValueError, KeyError) as error:  # one bad row must not sink the report
                entry.update(status="not_compared", reason=f"{type(error).__name__}: {error}")
    counts = Counter(row["verdict"] for row in rows)
    return {
        "schema_version": "1.0.0",
        "endpoints_compared": ["legacy-v2-reconstructed", "v3-host-conservative"],
        "endpoints_not_compared": ["v3-unique", "v3-equal"],
        "rows_compared": sum(1 for r in run_status if r["status"] == "compared"),
        "rows_not_compared": sum(1 for r in run_status if r["status"] != "compared"),
        "run_status": run_status,
        "verdict_counts": dict(counts),
        "records": rows,
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    inventory_parser = subparsers.add_parser("inventory", help="inventory legacy results")
    inventory_parser.add_argument("--source-root", type=Path, required=True)
    inventory_parser.add_argument("--raw-manifest", type=Path, required=True)
    inventory_parser.add_argument("--cohort-manifest", type=Path, required=True)
    run_parser = subparsers.add_parser(
        "run-row", help="run legacy and v3 comparisons for one raw-manifest row"
    )
    run_parser.add_argument("--raw-manifest", type=Path, required=True)
    run_parser.add_argument("--run-id", required=True)
    run_parser.add_argument("--source-root", type=Path, required=True)
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--threads", type=int, default=8)
    summarize_parser = subparsers.add_parser(
        "summarize", help="aggregate every planned comparison row"
    )
    summarize_parser.add_argument("--raw-manifest", type=Path, required=True)
    summarize_parser.add_argument("--run-root", type=Path, required=True)
    summarize_parser.add_argument("--output-dir", type=Path, required=True)
    validate_parser = subparsers.add_parser(
        "validate", help="validate all planned comparison outputs"
    )
    validate_parser.add_argument("--raw-manifest", type=Path, required=True)
    validate_parser.add_argument("--run-root", type=Path, required=True)
    validate_parser.add_argument("--report", type=Path, required=True)
    fresh_parser = subparsers.add_parser(
        "fresh-vs-archive", help="compare fresh v2/v3 control outputs with archived tables"
    )
    fresh_parser.add_argument(
        "--packet-root",
        type=Path,
        action="append",
        required=True,
        dest="packet_roots",
        help="repeatable; later packets override earlier ones per task row",
    )
    fresh_parser.add_argument("--archive-dir", type=Path, required=True)
    fresh_parser.add_argument("--report", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.command == "inventory":
        rows = inventory(args.source_root, args.raw_manifest, args.cohort_manifest)
        print(f"Inventoried {len(rows)} technical rows.")
        return 0
    if args.command == "run-row":
        status = run_row(
            args.raw_manifest,
            args.run_id,
            args.output_root,
            source_root=args.source_root,
            threads=args.threads,
        )
        print(json.dumps(status, indent=2, sort_keys=True))
        return 0 if status["status"] == "success" else 1
    if args.command == "summarize":
        outputs = summarize_results(args.raw_manifest, args.run_root, args.output_dir)
        print(json.dumps(outputs, indent=2, sort_keys=True))
        return 0
    if args.command == "validate":
        report = validate_results(args.raw_manifest, args.run_root)
        _atomic_json(args.report, report)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if report["status"] == "valid" else 1
    if args.command == "fresh-vs-archive":
        report = compare_fresh_vs_archive(args.packet_roots, args.archive_dir)
        _atomic_json(args.report, report)
        print(json.dumps({k: v for k, v in report.items() if k != "records"}, indent=2))
        return 0
    raise AssertionError(f"unhandled command {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
