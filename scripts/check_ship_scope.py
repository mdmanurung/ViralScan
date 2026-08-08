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
    from scripts.governance_utils import load_json, relative_path_error
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from governance_utils import load_json, relative_path_error  # type: ignore[no-redef]

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCOPE = REPO_ROOT / "config" / "public_ship_scope.json"
ROOT_KEYS = {
    "schema_version",
    "wheel",
    "sdist",
    "docker_context",
    "public_docs",
    "claim_bearing",
    "text_files",
}
LIST_KEYS = ("public_docs", "claim_bearing", "text_files")
ARCHIVE_KEYS = ("wheel", "sdist", "docker_context")


def _member_error(value: Any, *, allow_dist_info: bool = False) -> str:
    if allow_dist_info and isinstance(value, str):
        value = value.replace("{dist_info}", "distribution.dist-info")
    problem = relative_path_error(value)
    if problem:
        return problem
    if any(character in str(value) for character in "*?["):
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
    if scope.get("schema_version") != "1.0.0":
        errors.append("ship-scope schema_version must equal '1.0.0'")

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
            elif not (repo_root / value).is_file():
                errors.append(f"ship-scope {key} file not found: {value}")
    return scope, sorted(set(errors))


def _safe_archive_name(name: str) -> tuple[str, str]:
    normalized = name.rstrip("/")
    problem = relative_path_error(normalized)
    return normalized, problem


def _wheel_members(path: Path) -> tuple[set[str], list[str]]:
    members: set[str] = set()
    errors: list[str] = []
    raw_names: list[str] = []
    try:
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                if info.is_dir():
                    continue
                mode = info.external_attr >> 16
                if stat.S_ISLNK(mode):
                    errors.append(f"wheel contains symlink: {info.filename}")
                    continue
                name, problem = _safe_archive_name(info.filename)
                if problem:
                    errors.append(f"unsafe wheel member {info.filename!r}: {problem}")
                    continue
                raw_names.append(name)
    except (OSError, zipfile.BadZipFile) as exc:
        return set(), [f"cannot read wheel {path}: {exc}"]

    dist_info_dirs = {
        PurePosixPath(name).parts[0]
        for name in raw_names
        if PurePosixPath(name).parts and PurePosixPath(name).parts[0].endswith(".dist-info")
    }
    if len(dist_info_dirs) != 1:
        errors.append(
            f"wheel must contain exactly one .dist-info directory, found {sorted(dist_info_dirs)}"
        )
    for name in raw_names:
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
    names: list[str] = []
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
                names.append(name)
    except (OSError, tarfile.TarError) as exc:
        return set(), [f"cannot read sdist {path}: {exc}"]
    roots = {PurePosixPath(name).parts[0] for name in names if PurePosixPath(name).parts}
    if len(roots) != 1:
        errors.append(f"sdist must contain exactly one root directory, found {sorted(roots)}")
    for name in names:
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


def effective_docker_context_members(root: Path) -> set[str]:
    ignore_path = root / ".dockerignore"
    if not ignore_path.is_file():
        return set()
    try:
        lines = [
            line.strip()
            for line in ignore_path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except (OSError, UnicodeDecodeError):
        return set()
    if not lines or lines[0] != "**":
        return set()
    members: set[str] = set()
    for line in lines[1:]:
        if not line.startswith("!"):
            return set()
        value = line[1:].rstrip("/")
        if not value or _member_error(value):
            return set()
        path = root / value
        if path.is_file():
            members.add(value)
    return members


def validate_repository_scope(
    repo_root: Path = REPO_ROOT, scope_path: Optional[Path] = None
) -> list[str]:
    scope_path = scope_path or repo_root / "config/public_ship_scope.json"
    scope, errors = load_ship_scope(scope_path, repo_root=repo_root)
    if not scope:
        return errors
    sdist_members = set(scope.get("sdist", {}).get("members", []))
    generated_sdist = {"PKG-INFO", "setup.cfg"}
    for member in sdist_members - generated_sdist:
        if member.startswith("src/ViralScan.egg-info/"):
            continue
        if not (repo_root / member).is_file():
            errors.append(f"sdist allowlist source file not found: {member}")
    for member in scope.get("wheel", {}).get("members", []):
        if member.startswith("{dist_info}/"):
            continue
        source = repo_root / "src" / member
        if not source.is_file():
            errors.append(f"wheel allowlist source file not found: src/{member}")
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
    text_files = set(scope.get("text_files", []))
    required_text = public_docs | claim_bearing
    if not required_text.issubset(text_files):
        errors.append(
            f"public text files missing from text_files: {sorted(required_text - text_files)}"
        )

    observed_context = effective_docker_context_members(repo_root)
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
        if re.search(r"(?mi)^\s*(?:COPY|ADD)\s+\.\s", docker_text):
            errors.append("Dockerfile must not COPY or ADD the repository wholesale")
    return sorted(set(errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", type=Path, default=DEFAULT_SCOPE)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--docker-context", type=Path, required=True)
    args = parser.parse_args()
    scope, errors = load_ship_scope(args.scope, repo_root=REPO_ROOT)
    if scope:
        errors.extend(validate_archive_members(args.wheel, scope, "wheel"))
        errors.extend(validate_archive_members(args.sdist, scope, "sdist"))
        errors.extend(validate_repository_scope(args.docker_context, args.scope))
    errors = sorted(set(errors))
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("ship-scope validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
