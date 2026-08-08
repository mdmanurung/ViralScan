"""Shared fail-closed helpers for repository governance validators."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
from typing import Any

from jsonschema import Draft202012Validator
from jsonschema.exceptions import SchemaError

INSTITUTIONAL_PREFIXES = tuple(
    f"/{directory}/" for directory in ("exports", "home", "lustre", "gpfs", "scratch")
)


class DuplicateKeyError(ValueError):
    """Raised when JSON repeats a mapping key."""


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateKeyError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def load_json(path: Path) -> Any:
    """Load JSON while rejecting duplicate keys."""
    with path.open(encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=_unique_object)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_path_error(value: Any) -> str:
    """Return an error for unsafe/non-canonical repository-relative paths."""
    if not isinstance(value, str) or not value:
        return "path must be a non-empty string"
    if "\\" in value:
        return "path must use POSIX separators"
    path = PurePosixPath(value)
    if path.is_absolute():
        return "path must be repository-relative"
    if any(part in {"", ".", ".."} for part in path.parts):
        return "path must be normalized and cannot contain '.' or '..'"
    return ""


def json_schema_errors(document: Any, schema_path: Path) -> list[str]:
    """Validate a document against a Draft 2020-12 schema."""
    if not schema_path.is_file():
        return [f"schema not found: {schema_path}"]
    try:
        schema = load_json(schema_path)
        Draft202012Validator.check_schema(schema)
    except (OSError, ValueError, json.JSONDecodeError, SchemaError) as exc:
        return [f"invalid schema {schema_path}: {exc}"]

    errors: list[str] = []
    validator = Draft202012Validator(schema)
    for error in sorted(validator.iter_errors(document), key=lambda item: list(item.path)):
        location = ".".join(str(part) for part in error.absolute_path) or "<root>"
        errors.append(f"schema {location}: {error.message}")
    return errors


def institutional_path_errors(value: str, context: str) -> list[str]:
    return [
        f"institutional absolute path {prefix!r} in {context}"
        for prefix in INSTITUTIONAL_PREFIXES
        if prefix in value
    ]
