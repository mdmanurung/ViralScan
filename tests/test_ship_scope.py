import json
import tarfile
import zipfile
from pathlib import Path

from scripts.check_ship_scope import (
    effective_docker_context_members,
    load_ship_scope,
    validate_archive_members,
    validate_repository_scope,
)

ROOT = Path(__file__).resolve().parents[1]


def _scope() -> dict:
    return {
        "schema_version": "1.0.0",
        "wheel": {
            "members": [
                "viralscan/__init__.py",
                "{dist_info}/METADATA",
            ],
        },
        "sdist": {
            "members": ["PKG-INFO", "README.md", "src/viralscan/__init__.py"],
        },
        "docker_context": {
            "members": [".dockerignore", "Dockerfile", "environment.yml"],
        },
        "public_docs": ["README.md"],
        "claim_bearing": ["README.md"],
        "text_files": ["README.md", "Dockerfile", "environment.yml"],
    }


def test_scope_rejects_duplicates_wildcards_and_nonmembers(tmp_path: Path) -> None:
    scope = _scope()
    scope["wheel"]["members"] = ["viralscan/*.py", "viralscan/*.py"]
    path = tmp_path / "scope.json"
    path.write_text(json.dumps(scope), encoding="utf-8")
    errors = load_ship_scope(path, repo_root=tmp_path)[1]
    assert any("exact relative path" in error for error in errors)
    assert any("duplicate" in error for error in errors)


def test_wheel_and_sdist_reject_unexpected_members(tmp_path: Path) -> None:
    wheel = tmp_path / "ViralScan.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("viralscan/__init__.py", "")
        archive.writestr("viralscan-3.0.0.dist-info/METADATA", "")
        archive.writestr("private.txt", "")
    assert any(
        "unexpected wheel member" in error
        for error in validate_archive_members(wheel, _scope(), "wheel")
    )

    sdist = tmp_path / "ViralScan.tar.gz"
    source = tmp_path / "private.txt"
    source.write_text("private", encoding="utf-8")
    with tarfile.open(sdist, "w:gz") as archive:
        archive.add(source, arcname="ViralScan-3.0.0/private.txt")
    assert any(
        "unexpected sdist member" in error
        for error in validate_archive_members(sdist, _scope(), "sdist")
    )


def test_effective_docker_context_is_positive_allowlist(tmp_path: Path) -> None:
    for name in ("Dockerfile", "environment.yml", "private.txt"):
        (tmp_path / name).write_text(name, encoding="utf-8")
    (tmp_path / ".dockerignore").write_text(
        "**\n!.dockerignore\n!Dockerfile\n!environment.yml\n", encoding="utf-8"
    )
    assert effective_docker_context_members(tmp_path) == {
        ".dockerignore",
        "Dockerfile",
        "environment.yml",
    }


def test_repository_ship_scope_is_internally_consistent() -> None:
    assert validate_repository_scope(ROOT) == []
