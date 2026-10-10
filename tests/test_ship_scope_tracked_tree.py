"""PKG-01: every tracked package file must be in the wheel, sdist and Docker allowlists.

The allowlists in ``config/public_ship_scope.json`` and ``.dockerignore`` are positive
(deny-by-default), so a new ``src/viralscan`` module that is not added to them silently
disappears from the built artifact or the container. Compare them with the git tree.
"""

import fnmatch
import shutil
import subprocess
from pathlib import Path

import pytest

from scripts.check_ship_scope import effective_docker_context_members, load_ship_scope

ROOT = Path(__file__).resolve().parents[1]


def _tracked_package_files() -> set[str]:
    git = shutil.which("git")
    if git is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout; cannot enumerate tracked package files")
    out = subprocess.run(
        [git, "ls-files", "src/viralscan"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout
    return {line for line in out.splitlines() if line}


def _scope() -> dict:
    scope, errors = load_ship_scope(ROOT / "config" / "public_ship_scope.json", repo_root=ROOT)
    assert errors == []
    return scope


def test_tracked_package_files_are_in_every_allowlist() -> None:
    tracked = _tracked_package_files()
    assert "src/viralscan/menu.py" in tracked  # guards against an empty/misrooted listing
    scope = _scope()
    wheel = {f"src/{m}" for m in scope["wheel"]["members"] if m.startswith("viralscan/")}
    sdist = set(scope["sdist"]["members"])
    docker = set(scope["docker_context"]["members"])
    effective = effective_docker_context_members(ROOT)
    for label, shipped in (
        ("wheel", wheel),
        ("sdist", sdist),
        ("docker_context", docker),
        (".dockerignore", effective),
    ):
        pkg = {m for m in shipped if m.startswith("src/viralscan/")}
        assert sorted(tracked - pkg) == [], f"tracked but not shipped in {label}"
        assert sorted(pkg - tracked) == [], f"shipped in {label} but not tracked"


def test_non_python_package_files_are_covered_by_package_data() -> None:
    """The wheel builder only picks up non-.py files named in ``package-data``."""
    from setuptools.config.pyprojecttoml import read_configuration

    cfg = read_configuration(str(ROOT / "pyproject.toml"), expand=False)
    package_data = cfg["tool"]["setuptools"]["package-data"]
    for rel in sorted(_tracked_package_files()):
        path = Path(rel).relative_to("src")
        if path.suffix == ".py":
            continue
        # package-data globs are relative to the package directory; any ancestor package counts.
        covered = any(
            fnmatch.fnmatch(path.relative_to(*path.parts[:depth]).as_posix(), pattern)
            for depth in range(1, len(path.parts))
            for pattern in package_data.get(".".join(path.parts[:depth]), [])
        )
        assert covered, f"{rel} is not matched by [tool.setuptools.package-data]"
