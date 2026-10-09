"""Environment and completed-run validation for ViralScan v3."""

from __future__ import annotations

import functools
import importlib.util
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from viralscan.run_safety import RUN_COMPLETE, RUN_MANIFEST, sha256_file

PYTHON_REQUIREMENTS = ("anndata", "jsonschema", "numpy", "pandas", "scipy")
FULL_PYTHON_REQUIREMENTS = ("snakemake",)
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
#: ``kb info`` lines naming the binaries kb itself runs: ``kallisto: 0.52.0 (/path)``.
_KB_TOOL_RE = re.compile(r"^(kallisto|bustools):\s+(\S+)\s+\((.+)\)\s*$", re.MULTILINE)


@functools.cache
def kb_tools() -> dict[str, tuple[str, str]]:
    """``{name: (version, path)}`` for the kallisto and bustools that ``kb`` runs.

    REL-16: a conda ``kallisto`` on PATH and kb's bundled one both report the
    same version, yet the conda binary segfaults on a ``kb ref`` index and the
    bundled one spins forever on a conda-built index. An index must be read by
    the binary that built it, and ``kb count`` always uses its own, so every
    direct call follows the ``kb`` on PATH, not PATH itself. Empty when ``kb``
    is missing or its output is unreadable.
    """
    kb = shutil.which("kb")
    if kb is None:
        return {}
    try:
        out = subprocess.run(
            [kb, "info"], capture_output=True, text=True, timeout=120, check=False
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return {}
    return {
        name: (version, path)
        for name, version, path in _KB_TOOL_RE.findall(out)
        if Path(path).is_file()
    }


def tool_path(name: str) -> str | None:
    """The binary ViralScan runs for *name*: kb's own kallisto/bustools, else PATH."""
    if name in ("kallisto", "bustools") and name in kb_tools():
        return kb_tools()[name][1]
    return shutil.which(name)


def tool_provenance(names: tuple[str, ...] = ("kallisto", "bustools")) -> dict[str, Any]:
    """Resolved path, reported version and SHA-256 per tool, for run manifests.

    Same-version binaries differ (REL-16), so only the hash tells a parity run
    which one actually ran. Missing tools are recorded as ``None``.
    """
    out: dict[str, Any] = {}
    for name in names:
        path = tool_path(name)
        out[name] = (
            None
            if path is None
            else {
                "path": path,
                "version": kb_tools().get(name, (None, None))[0],
                "sha256": sha256_file(Path(path)),
            }
        )
    return out


V3_SCHEMA_PACKAGE = "viralscan.schemas.v3"
REQUIRED_V3_SCHEMAS = (
    "artifact_inventory.schema.json",
    "claim_registry.schema.json",
    "count_audit.schema.json",
    "evidence_manifest.schema.json",
    "h5ad_contract.json",
    "positive_control_receipt.schema.json",
    "positive_control_report.schema.json",
    "reference_manifest.schema.json",
    "run_complete.schema.json",
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
    """Validate a JSON document against one of the shipped v3 schemas.

    ``jsonschema`` is a hard dependency, so an ImportError here is a broken
    install rather than a finding about the document, and it propagates. An
    unreadable or malformed schema still reports an issue, because a run may
    legitimately be validated against a damaged installation and the caller
    needs every other issue in the same report.
    """
    import jsonschema

    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [ValidationIssue("error", "schema_unavailable", str(exc), str(schema_path))]
    if not _is_json_schema(schema):
        # A document with no schema keywords validates everything, so wiring one
        # in would look like enforcement while checking nothing.
        return [
            ValidationIssue(
                "error",
                "not_a_json_schema",
                f"{getattr(schema_path, 'name', schema_path)} declares no JSON Schema keywords",
                str(schema_path),
            )
        ]
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


def _is_json_schema(schema: Any) -> bool:
    """True when the document actually constrains anything.

    ``h5ad_contract.json`` is a prose contract, not a JSON Schema. Handing it to
    a validator accepts every input silently, so callers must be stopped from
    treating it as one.
    """
    if not isinstance(schema, dict):
        return False
    return bool({"$schema", "type", "properties", "required", "$ref"} & set(schema))


def h5ad_contract() -> dict[str, Any]:
    """Return the packaged H5AD contract that ``_matrix_issues`` enforces."""
    return _read_packaged_schema("h5ad_contract.json")


class SchemaContractError(RuntimeError):
    """Raised when ViralScan is about to write an artifact that violates its schema."""


def require_schema_valid(document: Any, schema_name: str, path: Any = "") -> None:
    """Refuse to write a public artifact that does not satisfy its shipped schema.

    This raises rather than returning issues, unlike the validate-run path. The
    difference is whose fault the failure is: at validate-run the artifact is
    input and a malformed one is a finding to report, while here ViralScan is the
    author, so a violation is a defect in this code and writing the file anyway
    would publish it under a schema it does not meet.
    """
    issues = validate_json_schema(_plain(document), packaged_schema_resource(schema_name))
    if issues:
        detail = "; ".join(issue.message for issue in issues)
        raise SchemaContractError(f"refusing to write {path or schema_name}: {detail}")


def _plain(value: Any) -> Any:
    """Coerce numpy scalars so a document validates as plain JSON.

    anndata restores ``uns`` integers as ``np.int64``, and counters summed with
    numpy arrive the same way, neither of which jsonschema accepts as ``integer``.
    Values are unchanged; only their Python types are. Applied at both boundaries,
    because a spurious refusal to write would be a production-only failure that
    no fixture built from Python ints can reproduce.
    """
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def doctor_report(profile: str = "full") -> dict[str, Any]:
    if profile not in {"pip", "full"}:
        raise ValueError("doctor profile must be 'pip' or 'full'.")
    needed = PYTHON_REQUIREMENTS + (FULL_PYTHON_REQUIREMENTS if profile == "full" else ())
    python = {name: importlib.util.find_spec(name) is not None for name in needed}
    tools = {name: tool_path(name) for name in FULL_TOOLS} if profile == "full" else {}
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


def _matrix_issues(
    adata: Any, path: Path, contract: dict[str, Any], audit_schema: Any
) -> list[ValidationIssue]:
    """Re-derive the H5AD count contract from a persisted file.

    ``multimapping.py`` already raises on these invariants when it builds the
    matrix, so this is not a duplicate of that check: it re-establishes them from
    bytes on disk, which is what catches a truncated write, a hand-edited file, or
    an artifact produced by a different version.

    The required layers and uns keys come from the packaged contract rather than
    from literals here, so the contract file cannot drift away from what is
    enforced.
    """
    issues: list[ValidationIssue] = []
    required = tuple(contract.get("required_layers") or {})
    for layer in required:
        if layer not in adata.layers:
            issues.append(ValidationIssue("error", "missing_layer", layer, str(path)))
    for key in contract.get("required_uns") or []:
        if key not in adata.uns:
            issues.append(ValidationIssue("error", "missing_uns", key, str(path)))
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
        # The audit record is the document count_audit.schema.json governs. Until
        # this was wired up the schema shipped without ever validating anything,
        # so two of its seven required fields had no reader at all.
        issues.extend(validate_json_schema(_plain(dict(audit)), audit_schema))
        total = int(audit.get("input_molecules", -1))
        partition = sum(
            int(audit.get(name, -1))
            for name in ("unique_molecules", "ambiguous_molecules", "unresolved_molecules")
        )
        if total != partition:
            issues.append(
                ValidationIssue("error", "audit_conservation", f"{total} != {partition}", str(path))
            )
        resolved = int(audit.get("resolved_molecules", -1))
        expected_resolved = sum(
            int(audit.get(name, -1)) for name in ("unique_molecules", "ambiguous_molecules")
        )
        if resolved != expected_resolved:
            issues.append(
                ValidationIssue(
                    "error",
                    "audit_resolved_mismatch",
                    f"{resolved} != {expected_resolved}",
                    str(path),
                )
            )
        allocated = float(ambiguous.sum())
        if not np.isclose(allocated, float(audit.get("allocated_ambiguous_mass", -1)), atol=1e-9):
            issues.append(
                ValidationIssue("error", "audit_layer_mismatch", str(allocated), str(path))
            )
        # The contract's third invariant: unique and ambiguous are non-overlapping
        # partitions of the molecules, so the unique layer must carry exactly the
        # audited unique molecules and no allocated mass.
        unique_mass = float(unique.sum())
        if not np.isclose(unique_mass, float(audit.get("unique_molecules", -1)), atol=1e-9):
            issues.append(
                ValidationIssue("error", "audit_unique_mismatch", str(unique_mass), str(path))
            )
    return issues


SCHEMA_BY_ARTIFACT = {
    "reference_manifest.json": "reference_manifest.schema.json",
    "evidence_manifest.json": "evidence_manifest.schema.json",
}


def _schema_or_issue(name: str) -> tuple[Any, ValidationIssue | None]:
    """Resolve a packaged schema, or the issue explaining why it is unusable.

    ``validate_run`` reports rather than raises, so a damaged installation has to
    become a finding in the same report as everything else it managed to check.
    """
    try:
        return packaged_schema_resource(name), None
    except (OSError, ImportError, TypeError, ValueError) as exc:
        return None, ValidationIssue("error", "schema_unavailable", str(exc), name)


def _sibling_manifest_issues(run_dir: Path) -> list[ValidationIssue]:
    """Validate every schema-governed manifest found anywhere under the run.

    These two shipped schemas had no reader before: nothing validated a reference
    or evidence manifest at any boundary. They are validated where they are found
    rather than at a fixed path, because a run may contain several.
    """
    issues: list[ValidationIssue] = []
    for artifact, schema_name in SCHEMA_BY_ARTIFACT.items():
        found = sorted(run_dir.rglob(artifact))
        if not found:
            continue
        schema_resource, schema_issue = _schema_or_issue(schema_name)
        if schema_issue is not None:
            issues.append(schema_issue)
            continue
        for path in found:
            try:
                document = json.loads(path.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError) as exc:
                issues.append(ValidationIssue("error", "invalid_manifest", str(exc), str(path)))
                continue
            issues.extend(
                ValidationIssue(issue.level, issue.code, issue.message, str(path))
                for issue in validate_json_schema(document, schema_resource)
            )
    return issues


def _completion_issues(run_dir: Path, manifest: dict[str, Any]) -> list[ValidationIssue]:
    """Check ``run_complete.json`` against the manifest and the artifacts on disk.

    Required when the manifest declares ``completion_marker``; runs that predate
    the field only get a warning for a missing marker (and are still verified if
    one is present).
    """
    path = run_dir / RUN_COMPLETE
    required = manifest.get("completion_marker") is True
    if not path.is_file():
        level = "error" if required else "warning"
        return [ValidationIssue(level, "missing_completion_marker", RUN_COMPLETE, str(run_dir))]
    try:
        marker = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return [ValidationIssue("error", "invalid_completion_marker", str(exc), str(path))]
    schema, schema_issue = _schema_or_issue("run_complete.schema.json")
    if schema_issue is not None:
        return [schema_issue]
    issues = [
        ValidationIssue(i.level, i.code, i.message, str(path))
        for i in validate_json_schema(marker, schema)
    ]
    if issues:
        return issues
    if manifest and marker["run_fingerprint"] != manifest.get("run_fingerprint"):
        issues.append(
            ValidationIssue(
                "error", "completion_fingerprint_mismatch", "marker is for another run", str(path)
            )
        )
    for rel, expected in marker["artifacts"].items():
        target = (run_dir / rel).resolve()
        if run_dir not in target.parents or not target.is_file():
            issues.append(ValidationIssue("error", "completion_artifact_missing", rel, str(path)))
        elif sha256_file(target) != expected:
            issues.append(ValidationIssue("error", "completion_artifact_mismatch", rel, str(path)))
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
        # No `if manifest:` guard. An empty mapping is exactly the case the schema
        # is there to reject, and skipping it reported a clean run for a manifest
        # that was missing every required field.
        issues.extend(validate_json_schema(manifest, schema_resource))

    issues.extend(_sibling_manifest_issues(run_dir))
    issues.extend(_completion_issues(run_dir, manifest))

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

    # --no-multimapping runs have no adata_multimap.h5ad by design
    # (scripts/multimap.py writes it only when config.multimapping is on); their
    # count matrix is kb's own adata.h5ad, which is not a v3-contract file, so
    # only its existence and readability are checked.
    multimapping_enabled = bool(manifest.get("options", {}).get("multimapping", True))
    h5ad_name = "adata_multimap.h5ad" if multimapping_enabled else "adata.h5ad"
    h5ads = sorted(run_dir.rglob(h5ad_name))
    if not h5ads:
        issues.append(
            ValidationIssue("error", "missing_h5ad", f"No {h5ad_name} found", str(run_dir))
        )
    if h5ads and multimapping_enabled:
        try:
            contract = h5ad_contract()
        except (OSError, ImportError, TypeError, ValueError) as exc:
            contract = {}
            issues.append(
                ValidationIssue("error", "schema_unavailable", str(exc), "h5ad_contract.json")
            )
        audit_schema, audit_schema_issue = _schema_or_issue("count_audit.schema.json")
        if audit_schema_issue is not None:
            issues.append(audit_schema_issue)
    for path in h5ads:
        try:
            adata = ad.read_h5ad(path)
        except Exception as exc:
            issues.append(ValidationIssue("error", "invalid_h5ad", str(exc), str(path)))
            continue
        if not multimapping_enabled:
            continue
        if adata.uns.get("count_schema_version") != "3.0.0":
            issues.append(
                ValidationIssue("error", "legacy_h5ad", "Not a v3 count schema", str(path))
            )
            continue
        if not contract or audit_schema is None:
            # The contract itself is unreadable, already reported above. Checking
            # the matrix against an empty contract would report a clean file.
            continue
        issues.extend(_matrix_issues(adata, path, contract, audit_schema))

    return {
        "schema_version": "3.0.0",
        "run_dir": str(run_dir),
        "ok": not any(issue.level == "error" for issue in issues),
        "h5ad_files": [str(path) for path in h5ads],
        "issues": [issue.as_dict() for issue in issues],
    }
