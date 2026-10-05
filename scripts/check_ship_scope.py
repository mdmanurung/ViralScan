#!/usr/bin/env python3
"""Validate ViralScan's positive wheel, sdist, and Docker-context allowlists."""

from __future__ import annotations

import argparse
import json
import re
import stat
import sys
import tarfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any, Optional

try:
    from scripts.governance_utils import (
        institutional_path_errors,
        load_json,
        relative_path_error,
        repository_file_error,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from governance_utils import (  # type: ignore[no-redef]
        institutional_path_errors,
        load_json,
        relative_path_error,
        repository_file_error,
    )

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCOPE = REPO_ROOT / "config" / "public_ship_scope.json"
ROOT_KEYS = {
    "schema_version",
    "wheel",
    "sdist",
    "docker_context",
    "public_docs",
    "claim_bearing",
}
LIST_KEYS = ("public_docs", "claim_bearing")
ARCHIVE_KEYS = ("wheel", "sdist", "docker_context")
SDIST_METADATA = {"PKG-INFO", "setup.cfg"}
SDIST_METADATA_PREFIX = "src/ViralScan.egg-info/"


def _member_error(value: Any, *, allow_dist_info: bool = False) -> str:
    normalized = value
    if isinstance(value, str) and "{dist_info}" in value:
        if not allow_dist_info or not value.startswith("{dist_info}/"):
            return "{dist_info} is permitted only as the wheel metadata root"
        if value.count("{dist_info}") != 1:
            return "wheel member contains more than one {dist_info} placeholder"
        normalized = value.replace("{dist_info}", "distribution.dist-info", 1)
    problem = relative_path_error(normalized)
    if problem:
        return problem
    if any(character in str(normalized) for character in "*?["):
        return "member must be an exact relative path, not a glob"
    return ""


def load_ship_scope(
    path: Path = DEFAULT_SCOPE, *, repo_root: Path = REPO_ROOT
) -> tuple[dict[str, Any], list[str]]:
    if not path.is_file():
        return {}, [f"ship-scope config not found: {path}"]
    try:
        scope = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, [f"cannot load ship-scope config {path}: {exc}"]
    if not isinstance(scope, dict):
        return {}, ["ship-scope config must be a JSON object"]

    errors: list[str] = []
    if set(scope) != ROOT_KEYS:
        errors.append(
            f"ship-scope keys must be exactly {sorted(ROOT_KEYS)}; observed {sorted(scope)}"
        )
    if scope.get("schema_version") != "1.1.0":
        errors.append("ship-scope schema_version must equal '1.1.0'")

    for key in ARCHIVE_KEYS:
        section = scope.get(key)
        if not isinstance(section, dict) or set(section) != {"members"}:
            errors.append(f"ship-scope {key} must contain only a members list")
            continue
        values = section.get("members")
        if not isinstance(values, list) or not all(isinstance(value, str) for value in values):
            errors.append(f"ship-scope {key}.members must be a list of strings")
            continue
        if len(values) != len(set(values)):
            errors.append(f"ship-scope {key}.members contains a duplicate")
        for value in values:
            problem = _member_error(value, allow_dist_info=key == "wheel")
            if problem:
                errors.append(f"ship-scope {key} member {value!r}: {problem}")

    for key in LIST_KEYS:
        values = scope.get(key)
        if not isinstance(values, list) or not all(isinstance(value, str) for value in values):
            errors.append(f"ship-scope {key} must be a list of strings")
            continue
        if len(values) != len(set(values)):
            errors.append(f"ship-scope {key} contains a duplicate")
        for value in values:
            problem = _member_error(value)
            if problem:
                errors.append(f"ship-scope {key} member {value!r}: {problem}")
                continue
            file_problem = repository_file_error(repo_root, value)
            if file_problem == "file not found":
                errors.append(f"ship-scope {key} file not found: {value}")
            elif file_problem:
                errors.append(f"ship-scope {key} file {value!r}: {file_problem}")
    return scope, sorted(set(errors))


def _safe_archive_name(name: str) -> tuple[str, str]:
    normalized = name.rstrip("/")
    problem = relative_path_error(normalized)
    return normalized, problem


def _text_content_errors(data: bytes, context: str) -> list[str]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        return [f"{context} is not UTF-8: {exc}"]
    return institutional_path_errors(text, context)


def _wheel_members(path: Path) -> tuple[set[str], list[str]]:
    members: set[str] = set()
    errors: list[str] = []
    entries: list[tuple[str, bytes]] = []
    raw_names: set[str] = set()
    try:
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                if info.is_dir():
                    continue
                mode = info.external_attr >> 16
                file_type = stat.S_IFMT(mode)
                if stat.S_ISLNK(mode):
                    errors.append(f"wheel contains symlink: {info.filename}")
                    continue
                if file_type not in {0, stat.S_IFREG}:
                    errors.append(f"wheel contains unsupported member type: {info.filename}")
                    continue
                name, problem = _safe_archive_name(info.filename)
                if problem:
                    errors.append(f"unsafe wheel member {info.filename!r}: {problem}")
                    continue
                if name in raw_names:
                    errors.append(f"duplicate wheel member: {name}")
                raw_names.add(name)
                try:
                    data = archive.read(info)
                except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
                    errors.append(f"cannot read wheel member {name}: {exc}")
                    continue
                entries.append((name, data))
                errors.extend(_text_content_errors(data, f"wheel member {name}"))
    except (OSError, zipfile.BadZipFile) as exc:
        return set(), [f"cannot read wheel {path}: {exc}"]

    dist_info_dirs = {
        PurePosixPath(name).parts[0]
        for name, _data in entries
        if PurePosixPath(name).parts and PurePosixPath(name).parts[0].endswith(".dist-info")
    }
    if len(dist_info_dirs) != 1:
        errors.append(
            f"wheel must contain exactly one .dist-info directory, found {sorted(dist_info_dirs)}"
        )
    for name, _data in entries:
        parts = list(PurePosixPath(name).parts)
        if parts and parts[0] in dist_info_dirs:
            parts[0] = "{dist_info}"
        normalized = PurePosixPath(*parts).as_posix()
        if normalized in members:
            errors.append(f"duplicate wheel member after normalization: {normalized}")
        members.add(normalized)
    return members, errors


def _sdist_members(path: Path) -> tuple[set[str], list[str]]:
    members: set[str] = set()
    errors: list[str] = []
    entries: list[str] = []
    try:
        with tarfile.open(path, "r:*") as archive:
            for info in archive.getmembers():
                if info.isdir():
                    continue
                if info.issym() or info.islnk():
                    errors.append(f"sdist contains link: {info.name}")
                    continue
                if not info.isfile():
                    errors.append(f"sdist contains unsupported member type: {info.name}")
                    continue
                name, problem = _safe_archive_name(info.name)
                if problem:
                    errors.append(f"unsafe sdist member {info.name!r}: {problem}")
                    continue
                extracted = archive.extractfile(info)
                if extracted is None:
                    errors.append(f"cannot read sdist member {name}")
                    continue
                try:
                    data = extracted.read()
                except OSError as exc:
                    errors.append(f"cannot read sdist member {name}: {exc}")
                    continue
                entries.append(name)
                errors.extend(_text_content_errors(data, f"sdist member {name}"))
    except (OSError, tarfile.TarError) as exc:
        return set(), [f"cannot read sdist {path}: {exc}"]

    roots = {PurePosixPath(name).parts[0] for name in entries if PurePosixPath(name).parts}
    if len(roots) != 1:
        errors.append(f"sdist must contain exactly one root directory, found {sorted(roots)}")
    for name in entries:
        parts = PurePosixPath(name).parts
        normalized = PurePosixPath(*parts[1:]).as_posix() if len(parts) > 1 else ""
        if not normalized:
            errors.append(f"sdist file is not below its root directory: {name}")
            continue
        if normalized in members:
            errors.append(f"duplicate sdist member after normalization: {normalized}")
        members.add(normalized)
    return members, errors


def validate_archive_members(path: Path, scope: dict[str, Any], target: str) -> list[str]:
    if target == "wheel":
        observed, errors = _wheel_members(path)
    elif target == "sdist":
        observed, errors = _sdist_members(path)
    else:
        return [f"unknown archive target: {target}"]
    expected = set(scope.get(target, {}).get("members", []))
    for member in sorted(observed - expected):
        errors.append(f"unexpected {target} member: {member}")
    for member in sorted(expected - observed):
        errors.append(f"missing {target} member: {member}")
    return sorted(set(errors))


def _effective_docker_context_members(root: Path) -> tuple[set[str], list[str]]:
    ignore_path = root / ".dockerignore"
    ignore_problem = repository_file_error(root, ".dockerignore")
    if ignore_problem:
        return set(), [f"invalid .dockerignore: {ignore_problem}"]
    try:
        lines = [
            line.strip()
            for line in ignore_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except (OSError, UnicodeDecodeError) as exc:
        return set(), [f"cannot read .dockerignore as UTF-8: {exc}"]
    if not lines or lines[0] != "**":
        return set(), [".dockerignore must start with the catch-all exclusion '**'"]

    members: set[str] = set()
    errors: list[str] = []
    for line_number, line in enumerate(lines[1:], start=2):
        if not line.startswith("!"):
            errors.append(f".dockerignore line {line_number} must be an exact allowlist negation")
            continue
        is_directory = line.endswith("/")
        value = line[1:].rstrip("/")
        problem = _member_error(value)
        if not value or problem:
            errors.append(f".dockerignore line {line_number} {line!r}: {problem or 'empty path'}")
            continue
        if is_directory:
            continue
        file_problem = repository_file_error(root, value)
        if file_problem:
            errors.append(f"docker context member {value!r}: {file_problem}")
            continue
        if value in members:
            errors.append(f"duplicate effective docker context member: {value}")
        members.add(value)
    return members, errors


def effective_docker_context_members(root: Path) -> set[str]:
    """Return exact regular files selected by the required positive Docker allowlist."""
    return _effective_docker_context_members(root)[0]


def ship_scope_source_paths(scope: dict[str, Any]) -> tuple[set[str], list[str]]:
    """Map every project-controlled shipped member to its repository source path."""
    sources: set[str] = set()
    errors: list[str] = []

    for member in scope.get("wheel", {}).get("members", []):
        if member.startswith("{dist_info}/"):
            continue
        if not member.startswith("viralscan/"):
            errors.append(f"wheel member has no project source mapping: {member}")
            continue
        sources.add(f"src/{member}")

    for member in scope.get("sdist", {}).get("members", []):
        if member in SDIST_METADATA or member.startswith(SDIST_METADATA_PREFIX):
            continue
        sources.add(member)

    for member in scope.get("docker_context", {}).get("members", []):
        sources.add(member)
    return sources, errors


def _repository_text_errors(repo_root: Path, relative: str) -> list[str]:
    problem = repository_file_error(repo_root, relative)
    if problem:
        return [f"ship-scope source {relative!r}: {problem}"]
    try:
        data = (repo_root / relative).read_bytes()
    except OSError as exc:
        return [f"cannot read ship-scope source {relative}: {exc}"]
    return _text_content_errors(data, f"ship-scope source {relative}")


def validate_repository_scope(
    repo_root: Path = REPO_ROOT, scope_path: Optional[Path] = None
) -> list[str]:
    scope_path = scope_path or repo_root / "config/public_ship_scope.json"
    scope, errors = load_ship_scope(scope_path, repo_root=repo_root)
    if not scope:
        return errors

    source_paths, mapping_errors = ship_scope_source_paths(scope)
    errors.extend(mapping_errors)
    for relative in sorted(source_paths):
        errors.extend(_repository_text_errors(repo_root, relative))

    sdist_members = set(scope.get("sdist", {}).get("members", []))
    public_docs = set(scope.get("public_docs", []))
    if not public_docs.issubset(sdist_members):
        errors.append(
            f"public_docs missing from sdist allowlist: {sorted(public_docs - sdist_members)}"
        )
    claim_bearing = set(scope.get("claim_bearing", []))
    if not claim_bearing.issubset(public_docs):
        errors.append(
            f"claim_bearing files must be public_docs: {sorted(claim_bearing - public_docs)}"
        )

    observed_context, context_errors = _effective_docker_context_members(repo_root)
    errors.extend(context_errors)
    for relative in sorted(observed_context - source_paths):
        errors.extend(_repository_text_errors(repo_root, relative))
    expected_context = set(scope.get("docker_context", {}).get("members", []))
    for member in sorted(observed_context - expected_context):
        errors.append(f"unexpected docker_context member: {member}")
    for member in sorted(expected_context - observed_context):
        errors.append(f"missing docker_context member: {member}")

    dockerfile = repo_root / "Dockerfile"
    try:
        docker_text = dockerfile.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        errors.append(f"cannot read Dockerfile: {exc}")
    else:
        shell_wholesale = re.compile(r"(?mi)^\s*(?:COPY|ADD)(?:\s+--[^\s]+)*\s+(?:\.|\./)\s+")
        json_wholesale = re.compile(r"(?mi)^\s*(?:COPY|ADD)\s+\[\s*[\"'](?:\.|\./)[\"']\s*,")
        if shell_wholesale.search(docker_text) or json_wholesale.search(docker_text):
            errors.append("Dockerfile must not COPY or ADD the repository wholesale")
    return sorted(set(errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--scope", type=Path)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--docker-context", type=Path, required=True)
    args = parser.parse_args()
    scope_path = args.scope or args.repo_root / "config/public_ship_scope.json"
    scope, errors = load_ship_scope(scope_path, repo_root=args.docker_context)
    if scope:
        errors.extend(validate_archive_members(args.wheel, scope, "wheel"))
        errors.extend(validate_archive_members(args.sdist, scope, "sdist"))
        errors.extend(validate_repository_scope(args.docker_context, scope_path))
    errors = sorted(set(errors))
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("ship-scope validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
