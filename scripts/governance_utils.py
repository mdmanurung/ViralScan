"""Shared fail-closed helpers for repository governance validators."""

from __future__ import annotations

import hashlib
import json
import subprocess
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
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        return "path cannot contain control characters"
    if "\\" in value:
        return "path must use POSIX separators"
    path = PurePosixPath(value)
    if path.is_absolute():
        return "path must be repository-relative"
    if any(part in {"", ".", ".."} for part in path.parts):
        return "path must be normalized and cannot contain '.' or '..'"
    return ""


def repository_file_error(repo_root: Path, value: str) -> str:
    """Return an error unless *value* names a regular, in-tree, non-symlink file."""
    problem = relative_path_error(value)
    if problem:
        return problem

    root = repo_root.resolve()
    candidate = root
    for part in PurePosixPath(value).parts:
        candidate = candidate / part
        if candidate.is_symlink():
            return "path must not traverse or name a symlink"
    if not candidate.exists():
        return "file not found"
    if not candidate.is_file():
        return "path is not a regular file"
    try:
        candidate.resolve(strict=True).relative_to(root)
    except (OSError, ValueError):
        return "path escapes the repository root"
    return ""


def _git(git_root: Path, arguments: list[str]) -> subprocess.CompletedProcess[bytes]:
    """Run Git with an argument vector and captured bytes."""
    return subprocess.run(
        ["git", "-C", str(git_root), *arguments],
        check=False,
        capture_output=True,
    )


def git_snapshot_errors(
    git_root: Path,
    git_sha: str,
    path_value: str,
    expected_sha256: str,
    *,
    context: str,
) -> list[str]:
    """Validate that an exact commit contains a regular blob with the expected bytes."""
    path_problem = relative_path_error(path_value)
    if path_problem:
        return [f"{context}: unsafe git path {path_value!r}: {path_problem}"]

    try:
        object_type = _git(git_root, ["cat-file", "-t", git_sha])
    except OSError as exc:
        return [f"{context}: cannot execute git: {exc}"]
    if object_type.returncode != 0:
        return [f"{context}: git commit {git_sha!r} not found in {git_root}"]
    observed_type = object_type.stdout.decode("ascii", errors="replace").strip()
    if observed_type != "commit":
        return [
            f"{context}: git object {git_sha!r} is not a commit "
            f"(object type {observed_type or 'unknown'})"
        ]

    tree_entry = _git(git_root, ["ls-tree", "--full-tree", "-z", git_sha, "--", path_value])
    if tree_entry.returncode != 0:
        return [f"{context}: cannot inspect {path_value!r} in git commit {git_sha}"]
    entries = [entry for entry in tree_entry.stdout.split(b"\0") if entry]
    if not entries:
        return [f"{context}: path {path_value!r} is absent from git commit {git_sha}"]
    if len(entries) != 1:
        return [f"{context}: path {path_value!r} is ambiguous in git commit {git_sha}"]
    try:
        metadata, encoded_path = entries[0].split(b"\t", 1)
        mode, entry_type, object_id = metadata.decode("ascii").split()
        observed_path = encoded_path.decode("utf-8")
    except (UnicodeDecodeError, ValueError):
        return [f"{context}: malformed git tree entry for {path_value!r} at {git_sha}"]
    if observed_path != path_value:
        return [f"{context}: git tree resolved {path_value!r} as unexpected path {observed_path!r}"]
    if mode == "120000":
        return [f"{context}: path {path_value!r} is a symlink in git commit {git_sha}"]
    if entry_type != "blob" or mode not in {"100644", "100755"}:
        return [
            f"{context}: path {path_value!r} is not a regular git blob "
            f"in commit {git_sha} (mode {mode}, type {entry_type})"
        ]

    blob = _git(git_root, ["cat-file", "blob", object_id])
    if blob.returncode != 0:
        return [f"{context}: cannot read git blob {object_id} for {path_value!r}"]
    observed_sha256 = hashlib.sha256(blob.stdout).hexdigest()
    if observed_sha256 != expected_sha256:
        return [
            f"{context}: git blob sha256 mismatch for {path_value}: "
            f"expected {expected_sha256}, observed {observed_sha256} at {git_sha}"
        ]
    return []


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
