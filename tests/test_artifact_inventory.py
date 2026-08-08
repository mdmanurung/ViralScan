import csv
import hashlib
from pathlib import Path

from scripts.validate_artifact_inventory import (
    INVENTORY_COLUMNS,
    validate_inventory_file,
    validate_inventory_rows,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas" / "v3" / "artifact_inventory.schema.json"


def _row(path: Path) -> dict[str, str]:
    return {
        "schema_version": "1.0.0",
        "artifact_id": "fixture",
        "kind": "input",
        "path": "artifact.txt",
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


def test_repository_artifact_inventory_is_valid() -> None:
    assert validate_inventory_file() == []
