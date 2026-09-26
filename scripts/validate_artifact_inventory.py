#!/usr/bin/env python3
"""Validate the claim-bearing ViralScan v3 artifact inventory."""

from __future__ import annotations

import argparse
import csv
import json
import os
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
DEFAULT_INVENTORY = REPO_ROOT / "analysis" / "v3_artifact_inventory.tsv"
DEFAULT_SCHEMA = REPO_ROOT / "schemas" / "v3" / "artifact_inventory.schema.json"
DEFAULT_REGISTRY = REPO_ROOT / "claims" / "registry.json"

INVENTORY_COLUMNS = [
    "schema_version",
    "artifact_id",
    "kind",
    "path",
    "scope",
    "software_version",
    "counting_version",
    "git_sha",
    "input_sha256s",
    "reference_sha256s",
    "artifact_sha256",
    "schema",
    "layer",
    "denominator",
    "generation_command",
    "rebuild_eligibility",
    "claim_ids",
]
JSON_COLUMNS = ("input_sha256s", "reference_sha256s", "claim_ids")


def _known_claims(path: Optional[Path]) -> tuple[set[str], dict[str, str], list[str]]:
    if path is None:
        return set(), {}, []
    if not path.is_file():
        return set(), {}, [f"claim registry not found: {path}"]
    try:
        document = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return set(), {}, [f"cannot load claim registry {path}: {exc}"]
    claims = document.get("claims", []) if isinstance(document, dict) else []
    ids: set[str] = set()
    statuses: dict[str, str] = {}
    for claim in claims:
        if isinstance(claim, dict) and isinstance(claim.get("id"), str):
            claim_id = claim["id"]
            ids.add(claim_id)
            statuses[claim_id] = str(claim.get("status", ""))
    return ids, statuses, []


def _parsed_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, Any]], list[str]]:
    parsed: list[dict[str, Any]] = []
    errors: list[str] = []
    for index, row in enumerate(rows, start=2):
        item: dict[str, Any] = dict(row)
        for column in JSON_COLUMNS:
            try:
                item[column] = json.loads(row.get(column, ""))
            except json.JSONDecodeError as exc:
                errors.append(f"row {index} {column}: invalid JSON: {exc.msg}")
                item[column] = None
        parsed.append(item)
    return parsed, errors


def validate_inventory_rows(
    rows: list[dict[str, str]],
    repo_root: Path,
    schema_path: Path,
    known_claim_ids: Optional[set[str]] = None,
    claim_statuses: Optional[dict[str, str]] = None,
    *,
    git_root: Optional[Path] = None,
) -> list[str]:
    parsed, errors = _parsed_rows(rows)
    errors.extend(json_schema_errors({"schema_version": "1.1.0", "artifacts": parsed}, schema_path))

    artifact_ids: set[str] = set()
    artifact_paths: set[str] = set()
    artifacts_by_id: dict[str, dict[str, str]] = {}
    known_claim_ids = known_claim_ids or set()
    claim_statuses = claim_statuses or {}
    git_root = git_root or repo_root
    for index, (raw, item) in enumerate(zip(rows, parsed), start=2):
        artifact_id = raw.get("artifact_id", "")
        path_value = raw.get("path", "")
        if artifact_id in artifact_ids:
            errors.append(f"row {index}: duplicate artifact_id {artifact_id!r}")
        else:
            artifacts_by_id[artifact_id] = raw
        artifact_ids.add(artifact_id)
        if path_value in artifact_paths:
            errors.append(f"row {index}: duplicate artifact path {path_value!r}")
        artifact_paths.add(path_value)

        for column in INVENTORY_COLUMNS:
            errors.extend(institutional_path_errors(raw.get(column, ""), f"row {index} {column}"))

        path_problem = relative_path_error(path_value)
        if path_problem:
            errors.append(f"row {index} path {path_value!r}: {path_problem}")
            continue
        artifact_path = repo_root / path_value
        file_problem = repository_file_error(repo_root, path_value)
        if file_problem == "file not found":
            errors.append(f"row {index}: missing artifact {path_value}")
        elif file_problem:
            errors.append(f"row {index}: unsafe artifact {path_value}: {file_problem}")
        else:
            observed = sha256_file(artifact_path)
            expected = raw.get("artifact_sha256", "")
            if observed != expected:
                errors.append(
                    f"row {index}: stale artifact sha256 for {path_value}: "
                    f"expected {expected}, observed {observed}"
                )

        git_sha = raw.get("git_sha", "")
        expected = raw.get("artifact_sha256", "")
        if len(git_sha) == 40 and len(expected) == 64:
            errors.extend(
                git_snapshot_errors(
                    git_root,
                    git_sha,
                    path_value,
                    expected,
                    context=f"row {index} artifact {artifact_id!r}",
                )
            )
        generation_command = raw.get("generation_command", "")
        if generation_command.startswith("git show "):
            expected_command = f"git show {git_sha}:{path_value}"
            if generation_command != expected_command:
                errors.append(f"row {index}: git retrieval command must equal {expected_command!r}")

        schema_value = raw.get("schema", "")
        if schema_value != "not_applicable":
            problem = relative_path_error(schema_value)
            if problem:
                errors.append(f"row {index} schema {schema_value!r}: {problem}")
            else:
                schema_problem = repository_file_error(repo_root, schema_value)
                if schema_problem == "file not found":
                    errors.append(f"row {index}: schema artifact not found: {schema_value}")
                elif schema_problem:
                    errors.append(
                        f"row {index}: unsafe schema artifact {schema_value}: {schema_problem}"
                    )

        for column in ("input_sha256s", "reference_sha256s"):
            values = item.get(column)
            if isinstance(values, list):
                ids = [value.get("id") for value in values if isinstance(value, dict)]
                for duplicate in sorted({value for value in ids if ids.count(value) > 1}):
                    errors.append(f"row {index}: duplicate {column} id {duplicate!r}")

        claim_ids = item.get("claim_ids")
        if isinstance(claim_ids, list) and known_claim_ids:
            for claim_id in claim_ids:
                if claim_id not in known_claim_ids:
                    errors.append(f"row {index}: unknown claim id {claim_id!r}")
                if claim_statuses.get(str(claim_id)) == "validated_v3" and raw.get(
                    "counting_version", ""
                ).lower().startswith("pre-v3"):
                    errors.append(
                        f"row {index}: validated_v3 claim {claim_id!r} uses pre-v3 counting"
                    )

    for index, item in enumerate(parsed, start=2):
        for column in ("input_sha256s", "reference_sha256s"):
            values = item.get(column)
            if not isinstance(values, list):
                continue
            for dependency in values:
                if not isinstance(dependency, dict):
                    continue
                dependency_id = dependency.get("id", "")
                dependency_row = artifacts_by_id.get(str(dependency_id))
                if dependency_row is None:
                    errors.append(f"row {index}: unknown {column} id {dependency_id!r}")
                    continue
                expected = dependency_row.get("artifact_sha256", "")
                observed = dependency.get("sha256", "")
                if observed != expected:
                    errors.append(
                        f"row {index}: {column} digest mismatch for {dependency_id!r}: "
                        f"expected {expected}, observed {observed}"
                    )
    return sorted(set(errors))


def _resolve_default(path: Optional[Path], repo_root: Path, relative: str) -> Path:
    return path if path is not None else repo_root / relative


def validate_inventory_file(
    inventory_path: Optional[Path] = None,
    *,
    repo_root: Path = REPO_ROOT,
    schema_path: Optional[Path] = None,
    registry_path: Optional[Path] = None,
    git_root: Optional[Path] = None,
) -> list[str]:
    inventory_path = _resolve_default(
        inventory_path, repo_root, "analysis/v3_artifact_inventory.tsv"
    )
    schema_path = _resolve_default(
        schema_path, repo_root, "schemas/v3/artifact_inventory.schema.json"
    )
    registry_path = _resolve_default(registry_path, repo_root, "claims/registry.json")
    if git_root is None:
        configured_git_root = os.environ.get("VIRALSCAN_GIT_ROOT")
        git_root = Path(configured_git_root) if configured_git_root else repo_root
    if not inventory_path.is_file():
        return [f"artifact inventory not found: {inventory_path}"]
    try:
        with inventory_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != INVENTORY_COLUMNS:
                return [
                    "artifact inventory header mismatch: expected " + "\t".join(INVENTORY_COLUMNS)
                ]
            rows = list(reader)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        return [f"cannot read artifact inventory {inventory_path}: {exc}"]
    if not rows:
        return ["artifact inventory contains no rows"]

    claim_ids, statuses, claim_errors = _known_claims(registry_path)
    return sorted(
        set(
            claim_errors
            + validate_inventory_rows(
                rows,
                repo_root,
                schema_path,
                known_claim_ids=claim_ids,
                claim_statuses=statuses,
                git_root=git_root,
            )
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", nargs="?", type=Path)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--git-root", type=Path)
    parser.add_argument("--schema", type=Path)
    parser.add_argument("--registry", type=Path)
    args = parser.parse_args()
    errors = validate_inventory_file(
        args.inventory,
        repo_root=args.repo_root,
        git_root=args.git_root,
        schema_path=args.schema,
        registry_path=args.registry,
    )
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("artifact inventory validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
