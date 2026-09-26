import csv
import hashlib
import json
import subprocess
from pathlib import Path

from scripts.validate_artifact_inventory import (
    INVENTORY_COLUMNS,
    validate_inventory_file,
    validate_inventory_rows,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas" / "v3" / "artifact_inventory.schema.json"


def _row(
    path: Path,
    *,
    artifact_id: str = "fixture",
    relative_path: str = "artifact.txt",
) -> dict[str, str]:
    return {
        "schema_version": "1.1.0",
        "artifact_id": artifact_id,
        "kind": "input",
        "path": relative_path,
        "scope": "public",
        "software_version": "ViralScan 3.0.0.dev0",
        "counting_version": "v3 molecule contract",
        "git_sha": "a" * 40,
        "input_sha256s": "[]",
        "reference_sha256s": "[]",
        "artifact_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "schema": "not_applicable",
        "layer": "not_applicable",
        "denominator": "one fixture",
        "generation_command": "printf evidence",
        "rebuild_eligibility": "rebuildable",
        "claim_ids": "[]",
    }


def _commit(repo: Path) -> str:
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "ViralScan tests"], check=True)
    subprocess.run(
        ["git", "-C", str(repo), "config", "user.email", "tests@viralscan.invalid"],
        check=True,
    )
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "fixture"], check=True)
    return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()


def test_rejects_incomplete_inventory_rows(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    row = _row(artifact)
    row["generation_command"] = ""
    errors = validate_inventory_rows([row], tmp_path, SCHEMA)
    assert any("generation_command" in error for error in errors)


def test_rejects_duplicate_artifact_identity_and_path(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    row = _row(artifact)
    errors = validate_inventory_rows([row, dict(row)], tmp_path, SCHEMA)
    assert any("duplicate artifact_id" in error for error in errors)
    assert any("duplicate artifact path" in error for error in errors)


def test_rejects_missing_stale_and_institutional_artifacts(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    row = _row(artifact)
    artifact.write_text("changed\n", encoding="utf-8")
    errors = validate_inventory_rows([row], tmp_path, SCHEMA)
    assert any("stale artifact sha256" in error for error in errors)

    row["path"] = "missing.txt"
    errors = validate_inventory_rows([row], tmp_path, SCHEMA)
    assert any("missing artifact" in error for error in errors)

    row["generation_command"] = "tool --input /exports/private/file"
    errors = validate_inventory_rows([row], tmp_path, SCHEMA)
    assert any("institutional absolute path" in error for error in errors)


def test_rejects_incomplete_tsv_header(tmp_path: Path) -> None:
    inventory = tmp_path / "inventory.tsv"
    with inventory.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_COLUMNS[:-1], delimiter="\t")
        writer.writeheader()
    errors = validate_inventory_file(inventory, repo_root=tmp_path, schema_path=SCHEMA)
    assert any("header" in error for error in errors)


def test_rejects_missing_and_non_commit_git_objects(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    commit = _commit(tmp_path)
    row = _row(artifact)

    row["git_sha"] = "a" * 40
    errors = validate_inventory_rows([row], tmp_path, SCHEMA)
    assert any("git commit" in error and "not found" in error for error in errors)

    blob = subprocess.check_output(
        ["git", "-C", str(tmp_path), "rev-parse", f"{commit}:artifact.txt"], text=True
    ).strip()
    row["git_sha"] = blob
    errors = validate_inventory_rows([row], tmp_path, SCHEMA)
    assert any("not a commit" in error for error in errors)


def test_rejects_absent_symlink_and_blob_mismatch_git_paths(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("original\n", encoding="utf-8")
    target = tmp_path / "target.txt"
    target.write_text("target\n", encoding="utf-8")
    link = tmp_path / "link.txt"
    link.symlink_to("target.txt")
    commit = _commit(tmp_path)

    absent = tmp_path / "absent.txt"
    absent.write_text("added after commit\n", encoding="utf-8")
    absent_row = _row(absent, artifact_id="absent", relative_path="absent.txt")
    absent_row["git_sha"] = commit
    errors = validate_inventory_rows([absent_row], tmp_path, SCHEMA)
    assert any("absent from git commit" in error for error in errors)

    link_row = _row(link, artifact_id="link", relative_path="link.txt")
    link_row["git_sha"] = commit
    errors = validate_inventory_rows([link_row], tmp_path, SCHEMA)
    assert any("symlink" in error for error in errors)

    artifact.write_text("changed after commit\n", encoding="utf-8")
    stale_snapshot = _row(artifact)
    stale_snapshot["git_sha"] = commit
    errors = validate_inventory_rows([stale_snapshot], tmp_path, SCHEMA)
    assert any("git blob sha256 mismatch" in error for error in errors)


def test_rejects_unknown_and_digest_mismatched_dependencies(tmp_path: Path) -> None:
    first = tmp_path / "first.txt"
    second = tmp_path / "second.txt"
    first.write_text("first\n", encoding="utf-8")
    second.write_text("second\n", encoding="utf-8")
    commit = _commit(tmp_path)
    first_row = _row(first, artifact_id="first", relative_path="first.txt")
    second_row = _row(second, artifact_id="second", relative_path="second.txt")
    first_row["git_sha"] = commit
    second_row["git_sha"] = commit

    first_row["input_sha256s"] = json.dumps(
        [{"id": "unknown", "sha256": "0" * 64}], separators=(",", ":")
    )
    errors = validate_inventory_rows([first_row, second_row], tmp_path, SCHEMA)
    assert any("unknown input_sha256s id 'unknown'" in error for error in errors)

    first_row["input_sha256s"] = json.dumps(
        [{"id": "second", "sha256": "0" * 64}], separators=(",", ":")
    )
    errors = validate_inventory_rows([first_row, second_row], tmp_path, SCHEMA)
    assert any("input_sha256s digest mismatch for 'second'" in error for error in errors)


def test_validation_receipt_is_an_allowed_inventoried_artifact_kind(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("immutable successful execution receipt\n", encoding="utf-8")
    commit = _commit(tmp_path)
    row = _row(artifact)
    row["kind"] = "validation_receipt"
    row["git_sha"] = commit

    assert validate_inventory_rows([row], tmp_path, SCHEMA) == []


def test_repository_artifact_inventory_is_valid() -> None:
    assert validate_inventory_file() == []
