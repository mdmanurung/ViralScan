import hashlib
import json
from pathlib import Path

from scripts.validate_claim_registry import (
    coverage_errors,
    validate_registry_document,
    validate_registry_file,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schemas" / "v3" / "claim_registry.schema.json"


def _claim(artifact: Path) -> dict:
    digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
    return {
        "id": "test-claim",
        "statement": "The fixture supports the test contract.",
        "public": True,
        "status": "validated_v3",
        "evidence_version": "v3",
        "source_locations": [
            {
                "path": "claim.md",
                "marker": "<!-- viralscan-claim:test-claim status=validated_v3 -->",
            }
        ],
        "artifacts": [
            {
                "artifact_id": "fixture",
                "path": "artifact.txt",
                "sha256": digest,
            }
        ],
        "git_sha": "a" * 40,
        "input_hashes": [{"id": "fixture", "sha256": digest}],
        "reference_hashes": [],
        "schema": "not_applicable",
        "layer": "not_applicable",
        "denominator": "one fixture",
        "generation_command": "pytest tests/test_claim_registry.py",
        "validation_scope": "synthetic test fixture",
    }


def _registry(claim: dict) -> dict:
    return {
        "schema_version": "1.0.0",
        "release_line": "3.0.0.dev0",
        "policy": "Only validated_v3 claims are eligible as validated v3 claims.",
        "claims": [claim],
    }


def test_rejects_incomplete_claim_metadata(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    claim = _claim(artifact)
    del claim["generation_command"]
    errors = validate_registry_document(_registry(claim), tmp_path, SCHEMA)
    assert any("generation_command" in error for error in errors)


def test_rejects_missing_and_stale_claim_artifacts(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    claim = _claim(artifact)
    artifact.write_text("changed\n", encoding="utf-8")
    errors = validate_registry_document(_registry(claim), tmp_path, SCHEMA)
    assert any("stale artifact sha256" in error for error in errors)

    claim["artifacts"][0]["path"] = "missing.txt"
    errors = validate_registry_document(_registry(claim), tmp_path, SCHEMA)
    assert any("missing artifact" in error for error in errors)


def test_rejects_duplicate_claim_ids_and_validated_legacy_evidence(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    claim = _claim(artifact)
    duplicate = json.loads(json.dumps(claim))
    errors = validate_registry_document(
        {**_registry(claim), "claims": [claim, duplicate]}, tmp_path, SCHEMA
    )
    assert any("duplicate claim id" in error for error in errors)

    claim["evidence_version"] = "pre_v3"
    errors = validate_registry_document(_registry(claim), tmp_path, SCHEMA)
    assert any("pre_v3" in error for error in errors)


def test_coverage_rejects_unregistered_and_missing_markers(tmp_path: Path) -> None:
    (tmp_path / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    public = tmp_path / "claim.md"
    public.write_text(
        "<!-- viralscan-claim:not-registered status=validated_v3 -->\n",
        encoding="utf-8",
    )
    errors = coverage_errors(
        {"test-claim": _claim(tmp_path / "artifact.txt")},
        [Path("claim.md")],
        tmp_path,
    )
    assert any("unregistered claim marker" in error for error in errors)
    assert any("source marker not found" in error for error in errors)


def test_repository_claim_registry_is_valid_and_covered() -> None:
    assert validate_registry_file(coverage=True) == []
