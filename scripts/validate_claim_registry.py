#!/usr/bin/env python3
"""Validate the ViralScan v3 public claim registry and marker coverage."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Optional

try:
    from scripts.governance_utils import (
        git_snapshot_errors,
        institutional_path_errors,
        json_schema_errors,
        load_json,
        relative_path_error,
        repository_file_error,
        sha256_file,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from governance_utils import (  # type: ignore[no-redef]
        git_snapshot_errors,
        institutional_path_errors,
        json_schema_errors,
        load_json,
        relative_path_error,
        repository_file_error,
        sha256_file,
    )

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = REPO_ROOT / "claims" / "registry.json"
DEFAULT_SCHEMA = REPO_ROOT / "schemas" / "v3" / "claim_registry.schema.json"
DEFAULT_SCOPE = REPO_ROOT / "config" / "public_ship_scope.json"
DEFAULT_INVENTORY = REPO_ROOT / "analysis" / "v3_artifact_inventory.tsv"

MARKER = re.compile(
    r"<!-- viralscan-claim:(?P<id>[a-z0-9][a-z0-9-]*) "
    r"status=(?P<status>[a-z0-9_]+) -->"
)
MARKER_START = re.compile(r"<!--\s*viralscan-claim:")


def validate_registry_document(document: Any, repo_root: Path, schema_path: Path) -> list[str]:
    errors = json_schema_errors(document, schema_path)
    if not isinstance(document, dict) or not isinstance(document.get("claims"), list):
        return sorted(set(errors))

    claim_ids: set[str] = set()
    artifact_ids: dict[str, tuple[str, str]] = {}
    for index, claim in enumerate(document["claims"]):
        if not isinstance(claim, dict):
            continue
        claim_id = claim.get("id", "")
        if claim_id in claim_ids:
            errors.append(f"claim {index}: duplicate claim id {claim_id!r}")
        claim_ids.add(claim_id)

        expected_marker = f"<!-- viralscan-claim:{claim_id} status={claim.get('status', '')} -->"
        for source in claim.get("source_locations", []):
            if not isinstance(source, dict):
                continue
            path_value = source.get("path", "")
            problem = relative_path_error(path_value)
            if problem:
                errors.append(f"claim {claim_id!r} source {path_value!r}: {problem}")
            else:
                file_problem = repository_file_error(repo_root, path_value)
                if file_problem == "file not found":
                    errors.append(f"claim {claim_id!r}: source file not found: {path_value}")
                elif file_problem:
                    errors.append(
                        f"claim {claim_id!r}: unsafe source file {path_value}: {file_problem}"
                    )
            if source.get("marker") != expected_marker:
                errors.append(f"claim {claim_id!r}: source marker must equal {expected_marker!r}")

        for artifact in claim.get("artifacts", []):
            if not isinstance(artifact, dict):
                continue
            artifact_id = artifact.get("artifact_id", "")
            path_value = artifact.get("path", "")
            expected = artifact.get("sha256", "")
            previous = artifact_ids.get(artifact_id)
            if previous is not None and previous != (path_value, expected):
                errors.append(
                    f"artifact identity {artifact_id!r} maps to conflicting path or sha256"
                )
            artifact_ids[artifact_id] = (path_value, expected)
            problem = relative_path_error(path_value)
            if problem:
                errors.append(f"claim {claim_id!r} artifact {path_value!r}: {problem}")
                continue
            artifact_path = repo_root / path_value
            file_problem = repository_file_error(repo_root, path_value)
            if file_problem == "file not found":
                errors.append(f"claim {claim_id!r}: missing artifact {path_value}")
            elif file_problem:
                errors.append(f"claim {claim_id!r}: unsafe artifact {path_value}: {file_problem}")
            else:
                observed = sha256_file(artifact_path)
                if observed != expected:
                    errors.append(
                        f"claim {claim_id!r}: stale artifact sha256 for {path_value}: "
                        f"expected {expected}, observed {observed}"
                    )

        for field in ("input_hashes", "reference_hashes"):
            values = claim.get(field, [])
            if isinstance(values, list):
                ids = [item.get("id") for item in values if isinstance(item, dict)]
                for duplicate in sorted({item for item in ids if ids.count(item) > 1}):
                    errors.append(f"claim {claim_id!r}: duplicate {field} id {duplicate!r}")

        serialized = json.dumps(claim, sort_keys=True)
        errors.extend(institutional_path_errors(serialized, f"claim {claim_id!r}"))
        if claim.get("evidence_version") == "pre_v3" and claim.get("status") == "validated_v3":
            errors.append(f"claim {claim_id!r}: pre_v3 evidence cannot be validated_v3")
    return sorted(set(errors))


def coverage_errors(
    claims_by_id: dict[str, dict[str, Any]],
    claim_bearing: list[Path],
    repo_root: Path,
) -> list[str]:
    errors: list[str] = []
    markers_by_path: dict[str, list[tuple[str, str, str]]] = {}
    for relative in claim_bearing:
        path_value = relative.as_posix()
        problem = relative_path_error(path_value)
        if problem:
            errors.append(f"claim-bearing path {path_value!r}: {problem}")
            continue
        path = repo_root / relative
        if not path.is_file():
            errors.append(f"claim-bearing file not found: {path_value}")
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError) as exc:
            errors.append(f"cannot read claim-bearing file {path_value}: {exc}")
            continue
        canonical = list(MARKER.finditer(text))
        if len(MARKER_START.findall(text)) != len(canonical):
            errors.append(f"malformed claim marker in {path_value}")
        markers_by_path[path_value] = [
            (match.group("id"), match.group("status"), match.group(0)) for match in canonical
        ]
        for claim_id, status, _marker in markers_by_path[path_value]:
            claim = claims_by_id.get(claim_id)
            if claim is None:
                errors.append(f"unregistered claim marker {claim_id!r} in {path_value}")
            elif claim.get("status") != status:
                errors.append(
                    f"claim marker status mismatch for {claim_id!r} in {path_value}: "
                    f"marker {status!r}, registry {claim.get('status')!r}"
                )

    for claim_id, claim in claims_by_id.items():
        if not claim.get("public"):
            continue
        for source in claim.get("source_locations", []):
            if not isinstance(source, dict):
                continue
            path_value = source.get("path", "")
            marker = source.get("marker", "")
            matches = [
                found
                for _identifier, _status, found in markers_by_path.get(path_value, [])
                if found == marker
            ]
            if len(matches) != 1:
                errors.append(
                    f"claim {claim_id!r}: source marker not found exactly once in {path_value}"
                )
    return sorted(set(errors))


def _load_scope_files(path: Path) -> tuple[list[Path], list[Path], list[str]]:
    if not path.is_file():
        return [], [], [f"ship-scope config not found: {path}"]
    try:
        scope = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return [], [], [f"cannot load ship-scope config {path}: {exc}"]
    if not isinstance(scope, dict):
        return [], [], ["ship-scope config must be a JSON object"]
    errors: list[str] = []
    lists: dict[str, list[Path]] = {}
    for key in ("public_docs", "claim_bearing"):
        values = scope.get(key)
        if not isinstance(values, list) or not all(isinstance(value, str) for value in values):
            errors.append(f"ship-scope {key} must be a list of paths")
            lists[key] = []
        else:
            lists[key] = [Path(value) for value in values]
    return lists["public_docs"], lists["claim_bearing"], errors


def _inventory_errors(
    claims_by_id: dict[str, dict[str, Any]], inventory_path: Path, git_root: Path
) -> list[str]:
    if not inventory_path.is_file():
        return [f"artifact inventory not found: {inventory_path}"]
    try:
        with inventory_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        return [f"cannot read artifact inventory {inventory_path}: {exc}"]
    errors: list[str] = []
    by_id: dict[str, dict[str, str]] = {}
    by_path: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(rows, start=2):
        artifact_id = row.get("artifact_id", "")
        path_value = row.get("path", "")
        if artifact_id in by_id:
            errors.append(f"inventory row {row_number}: duplicate artifact_id {artifact_id!r}")
        else:
            by_id[artifact_id] = row
        if path_value in by_path:
            errors.append(f"inventory row {row_number}: duplicate artifact path {path_value!r}")
        else:
            by_path[path_value] = row

    for claim_id, claim in claims_by_id.items():
        linked_rows: dict[str, dict[str, str]] = {}
        artifact_ids: list[str] = []
        for artifact in claim.get("artifacts", []):
            if not isinstance(artifact, dict):
                continue
            artifact_id = artifact.get("artifact_id", "")
            artifact_ids.append(str(artifact_id))
            row = by_id.get(artifact_id)
            if row is None:
                errors.append(
                    f"claim {claim_id!r}: artifact {artifact_id!r} is absent from inventory"
                )
                continue
            linked_rows[artifact_id] = row
            if row.get("path") != artifact.get("path"):
                errors.append(f"claim {claim_id!r}: inventory path mismatch for {artifact_id!r}")
            if row.get("artifact_sha256") != artifact.get("sha256"):
                errors.append(f"claim {claim_id!r}: inventory sha256 mismatch for {artifact_id!r}")
            try:
                row_claims = json.loads(row.get("claim_ids", ""))
            except json.JSONDecodeError:
                row_claims = []
            if claim_id not in row_claims:
                errors.append(
                    f"claim {claim_id!r}: inventory artifact {artifact_id!r} lacks claim linkage"
                )
        for duplicate in sorted({value for value in artifact_ids if artifact_ids.count(value) > 1}):
            errors.append(f"claim {claim_id!r}: duplicate artifact id {duplicate!r}")

        for field in ("input_hashes", "reference_hashes"):
            for dependency in claim.get(field, []):
                if not isinstance(dependency, dict):
                    continue
                dependency_id = dependency.get("id", "")
                row = by_id.get(str(dependency_id))
                if row is None:
                    errors.append(f"claim {claim_id!r}: unknown {field} id {dependency_id!r}")
                    continue
                linked_rows[str(dependency_id)] = row
                expected = row.get("artifact_sha256", "")
                observed = dependency.get("sha256", "")
                if observed != expected:
                    errors.append(
                        f"claim {claim_id!r}: {field} digest mismatch for {dependency_id!r}: "
                        f"expected {expected}, observed {observed}"
                    )
                try:
                    row_claims = json.loads(row.get("claim_ids", ""))
                except json.JSONDecodeError:
                    row_claims = []
                if claim_id not in row_claims:
                    errors.append(
                        f"claim {claim_id!r}: inventory dependency {dependency_id!r} "
                        "lacks claim linkage"
                    )

        for source in claim.get("source_locations", []):
            if not isinstance(source, dict):
                continue
            path_value = source.get("path", "")
            row = by_path.get(str(path_value))
            if row is None:
                errors.append(
                    f"claim {claim_id!r}: source document {path_value!r} is absent from inventory"
                )
                continue
            try:
                row_claims = json.loads(row.get("claim_ids", ""))
            except json.JSONDecodeError:
                row_claims = []
            if claim_id not in row_claims:
                errors.append(
                    f"claim {claim_id!r}: inventory source document {path_value!r} "
                    "lacks claim linkage"
                )

        if (
            claim.get("status") == "validated_v3"
            and claim.get("claim_type") in {"workflow_execution", "scientific_result"}
            and not any(
                row.get("kind") in {"result", "validation_receipt"}
                for artifact_id, row in linked_rows.items()
                if artifact_id in artifact_ids
            )
        ):
            errors.append(
                f"claim {claim_id!r}: validated {claim.get('claim_type')} requires an "
                "inventoried result or validation_receipt artifact"
            )

        git_sha = claim.get("git_sha", "")
        if linked_rows and git_sha == "not_applicable":
            errors.append(
                f"claim {claim_id!r}: linked evidence requires an exact git commit snapshot"
            )
        elif isinstance(git_sha, str) and len(git_sha) == 40:
            snapshot_paths: set[tuple[str, str]] = set()
            for row in linked_rows.values():
                path_value = row.get("path", "")
                digest = row.get("artifact_sha256", "")
                key = (path_value, digest)
                if key in snapshot_paths:
                    continue
                snapshot_paths.add(key)
                errors.extend(
                    git_snapshot_errors(
                        git_root,
                        git_sha,
                        path_value,
                        digest,
                        context=f"claim {claim_id!r}",
                    )
                )
    return errors


def _private_claim_leakage_errors(
    claims_by_id: dict[str, dict[str, Any]], public_docs: list[Path], repo_root: Path
) -> list[str]:
    claim = claims_by_id.get("v3-ebv-molecule-baseline")
    if not claim or claim.get("public") or claim.get("status") == "validated_v3":
        return []

    errors: list[str] = []
    for relative in public_docs:
        path = repo_root / relative
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if "103,145,071" in text:
            errors.append(
                f"private incomplete EBV claim value 103,145,071 appears in public document "
                f"{relative.as_posix()}"
            )
        if re.search(r"retained\s+EBV\s+baseline\s+verif(?:y|ies)", text, re.IGNORECASE):
            errors.append(
                "private incomplete validation wording 'retained EBV baseline verifies' "
                f"appears in public document {relative.as_posix()}"
            )
    return errors


def validate_registry_file(
    registry_path: Optional[Path] = None,
    *,
    repo_root: Path = REPO_ROOT,
    schema_path: Optional[Path] = None,
    scope_path: Optional[Path] = None,
    inventory_path: Optional[Path] = None,
    git_root: Optional[Path] = None,
    coverage: bool = False,
) -> list[str]:
    registry_path = registry_path or repo_root / "claims/registry.json"
    schema_path = schema_path or repo_root / "schemas/v3/claim_registry.schema.json"
    scope_path = scope_path or repo_root / "config/public_ship_scope.json"
    inventory_path = inventory_path or repo_root / "analysis/v3_artifact_inventory.tsv"
    if git_root is None:
        configured_git_root = os.environ.get("VIRALSCAN_GIT_ROOT")
        git_root = Path(configured_git_root) if configured_git_root else repo_root
    if not registry_path.is_file():
        return [f"claim registry not found: {registry_path}"]
    try:
        document = load_json(registry_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return [f"cannot load claim registry {registry_path}: {exc}"]
    errors = validate_registry_document(document, repo_root, schema_path)
    claims = document.get("claims", []) if isinstance(document, dict) else []
    claims_by_id = {
        claim["id"]: claim
        for claim in claims
        if isinstance(claim, dict) and isinstance(claim.get("id"), str)
    }
    errors.extend(_inventory_errors(claims_by_id, inventory_path, git_root))
    if coverage:
        public_docs, claim_files, scope_errors = _load_scope_files(scope_path)
        errors.extend(scope_errors)
        errors.extend(coverage_errors(claims_by_id, claim_files, repo_root))
        errors.extend(_private_claim_leakage_errors(claims_by_id, public_docs, repo_root))
    return sorted(set(errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("registry", nargs="?", type=Path)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--git-root", type=Path)
    parser.add_argument("--schema", type=Path)
    parser.add_argument("--scope", type=Path)
    parser.add_argument("--inventory", type=Path)
    parser.add_argument("--coverage", action="store_true")
    args = parser.parse_args()
    errors = validate_registry_file(
        args.registry,
        repo_root=args.repo_root,
        git_root=args.git_root,
        schema_path=args.schema,
        scope_path=args.scope,
        inventory_path=args.inventory,
        coverage=args.coverage,
    )
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("claim registry validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
