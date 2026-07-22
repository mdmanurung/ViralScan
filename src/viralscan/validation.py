"""Environment and completed-run validation for ViralScan v3."""

from __future__ import annotations

import importlib.util
import json
import shutil
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from viralscan.run_safety import RUN_MANIFEST, sha256_file

PYTHON_REQUIREMENTS = ("anndata", "jsonschema", "numpy", "pandas", "scipy", "snakemake")
FULL_TOOLS = (
    "kb",
    "kallisto",
    "bustools",
    "STAR",
    "minimap2",
    "samtools",
    "blastn",
    "makeblastdb",
    "cd-hit-est",
    "Rscript",
    "snakemake",
)
V3_SCHEMA_PACKAGE = "viralscan.schemas.v3"
REQUIRED_V3_SCHEMAS = (
    "count_audit.schema.json",
    "evidence_manifest.schema.json",
    "h5ad_contract.json",
    "reference_manifest.schema.json",
    "run_manifest.schema.json",
    "validation_protocol.schema.json",
)


@dataclass(frozen=True)
class ValidationIssue:
    level: str
    code: str
    message: str
    path: str = ""

    def as_dict(self) -> dict[str, str]:
        return {
            "level": self.level,
            "code": self.code,
            "message": self.message,
            "path": self.path,
        }


def packaged_schema_resource(name: str) -> Any:
    """Return a required packaged v3 schema resource, failing closed if absent."""
    if name not in REQUIRED_V3_SCHEMAS:
        raise ValueError(f"Unknown ViralScan v3 schema: {name}")
    resource = resources.files(V3_SCHEMA_PACKAGE).joinpath(name)
    if not resource.is_file():
        raise FileNotFoundError(f"Required packaged ViralScan v3 schema is missing: {name}")
    return resource


def _read_packaged_schema(name: str) -> dict[str, Any]:
    resource = packaged_schema_resource(name)
    try:
        document = json.loads(resource.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Packaged ViralScan v3 schema is invalid JSON: {name}: {exc}") from exc
    if not isinstance(document, dict):
        raise ValueError(f"Packaged ViralScan v3 schema must be a JSON object: {name}")
    return document


def validate_json_schema(document: Any, schema_path: Any) -> list[ValidationIssue]:
    """Validate a JSON document against one of the shipped v3 schemas."""
    try:
        import jsonschema

        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        validator = jsonschema.Draft202012Validator(schema)
        return [
            ValidationIssue(
                "error",
                "schema_validation",
                error.message,
                str(schema_path),
            )
            for error in sorted(validator.iter_errors(document), key=lambda item: list(item.path))
        ]
    except (OSError, json.JSONDecodeError, ImportError) as exc:
        return [ValidationIssue("error", "schema_unavailable", str(exc), str(schema_path))]


def doctor_report(profile: str = "full") -> dict[str, Any]:
    if profile not in {"pip", "full"}:
        raise ValueError("doctor profile must be 'pip' or 'full'.")
    python = {name: importlib.util.find_spec(name) is not None for name in PYTHON_REQUIREMENTS}
    tools = {name: shutil.which(name) for name in FULL_TOOLS} if profile == "full" else {}
    schemas: dict[str, bool] = {}
    schema_errors: dict[str, str] = {}
    for name in REQUIRED_V3_SCHEMAS:
        try:
            _read_packaged_schema(name)
        except (OSError, ImportError, TypeError, ValueError) as exc:
            schemas[name] = False
            schema_errors[name] = str(exc)
        else:
            schemas[name] = True
    ok = all(python.values()) and all(tools.values()) and all(schemas.values())
    return {
        "profile": profile,
        "ok": ok,
        "python": python,
        "tools": tools,
        "schemas": schemas,
        "schema_errors": schema_errors,
    }


def _matrix_issues(adata: Any, path: Path) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    required = ("counts_unique", "counts_ambiguous_allocated")
    for layer in required:
        if layer not in adata.layers:
            issues.append(ValidationIssue("error", "missing_layer", layer, str(path)))
    if issues:
        return issues

    unique = adata.layers["counts_unique"]
    ambiguous = adata.layers["counts_ambiguous_allocated"]
    for name, matrix in (
        ("X", adata.X),
        ("counts_unique", unique),
        ("counts_ambiguous_allocated", ambiguous),
    ):
        data = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).reshape(-1)
        if not np.isfinite(data).all():
            issues.append(ValidationIssue("error", "non_finite", name, str(path)))
        if (data < 0).any():
            issues.append(ValidationIssue("error", "negative_mass", name, str(path)))
    delta = adata.X - (unique + ambiguous)
    delta_data = delta.data if sparse.issparse(delta) else np.asarray(delta).reshape(-1)
    if delta_data.size and float(np.max(np.abs(delta_data))) > 1e-9:
        issues.append(
            ValidationIssue("error", "x_layer_mismatch", "X != unique + ambiguous", str(path))
        )

    audit = adata.uns.get("molecule_audit")
    if not isinstance(audit, dict):
        issues.append(ValidationIssue("error", "missing_audit", "molecule_audit", str(path)))
    else:
        total = int(audit.get("input_molecules", -1))
        partition = sum(
            int(audit.get(name, -1))
            for name in ("unique_molecules", "ambiguous_molecules", "unresolved_molecules")
        )
        if total != partition:
            issues.append(
                ValidationIssue("error", "audit_conservation", f"{total} != {partition}", str(path))
            )
        allocated = float(ambiguous.sum())
        if not np.isclose(allocated, float(audit.get("allocated_ambiguous_mass", -1)), atol=1e-9):
            issues.append(
                ValidationIssue("error", "audit_layer_mismatch", str(allocated), str(path))
            )
    return issues


def validate_run(run_dir: Path, verify_inputs: bool = True) -> dict[str, Any]:
    import anndata as ad

    run_dir = run_dir.resolve()
    issues: list[ValidationIssue] = []
    manifest_path = run_dir / RUN_MANIFEST
    if not manifest_path.is_file():
        issues.append(ValidationIssue("error", "missing_manifest", RUN_MANIFEST, str(run_dir)))
        manifest: dict[str, Any] = {}
    else:
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            issues.append(
                ValidationIssue("error", "invalid_manifest", str(exc), str(manifest_path))
            )
            manifest = {}

    if manifest and manifest.get("schema_version") != "3.0.0":
        issues.append(
            ValidationIssue(
                "error", "wrong_schema", str(manifest.get("schema_version")), str(manifest_path)
            )
        )

    try:
        schema_resource = packaged_schema_resource("run_manifest.schema.json")
    except (OSError, ImportError, TypeError, ValueError) as exc:
        issues.append(
            ValidationIssue(
                "error",
                "schema_unavailable",
                str(exc),
                "run_manifest.schema.json",
            )
        )
    else:
        if manifest:
            issues.extend(validate_json_schema(manifest, schema_resource))

    if verify_inputs and manifest:
        options = manifest.get("options", {})
        for key, expected in manifest.get("input_fingerprints", {}).items():
            field, index_raw = key.split(":", 1)
            values = str(options.get(field, "")).split(",")
            index = int(index_raw)
            if index >= len(values) or not Path(values[index]).is_file():
                issues.append(
                    ValidationIssue(
                        "error", "missing_input", key, values[index] if index < len(values) else ""
                    )
                )
            elif sha256_file(Path(values[index]).resolve()) != expected:
                issues.append(
                    ValidationIssue("error", "input_fingerprint_mismatch", key, values[index])
                )

    h5ads = sorted(run_dir.rglob("adata_multimap.h5ad"))
    if not h5ads:
        issues.append(
            ValidationIssue("error", "missing_h5ad", "No adata_multimap.h5ad found", str(run_dir))
        )
    for path in h5ads:
        try:
            adata = ad.read_h5ad(path)
        except Exception as exc:
            issues.append(ValidationIssue("error", "invalid_h5ad", str(exc), str(path)))
            continue
        if adata.uns.get("count_schema_version") != "3.0.0":
            issues.append(
                ValidationIssue("error", "legacy_h5ad", "Not a v3 count schema", str(path))
            )
            continue
        issues.extend(_matrix_issues(adata, path))

    return {
        "schema_version": "3.0.0",
        "run_dir": str(run_dir),
        "ok": not any(issue.level == "error" for issue in issues),
        "h5ad_files": [str(path) for path in h5ads],
        "issues": [issue.as_dict() for issue in issues],
    }
