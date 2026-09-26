import csv
import hashlib
import json
from pathlib import Path

from scripts.validate_artifact_inventory import INVENTORY_COLUMNS
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
        "claim_type": "software_invariant",
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
        "schema_version": "1.1.0",
        "release_line": "3.0.0.dev0",
        "policy": "Only validated_v3 claims are eligible as validated v3 claims.",
        "claims": [claim],
    }


def _inventory_row(artifact: Path, *, kind: str = "command") -> dict[str, str]:
    digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
    return {
        "schema_version": "1.1.0",
        "artifact_id": "fixture",
        "kind": kind,
        "path": "artifact.txt",
        "scope": "public",
        "software_version": "ViralScan 3.0.0.dev0",
        "counting_version": "v3 molecule contract",
        "git_sha": "a" * 40,
        "input_sha256s": "[]",
        "reference_sha256s": "[]",
        "artifact_sha256": digest,
        "schema": "not_applicable",
        "layer": "not_applicable",
        "denominator": "one fixture",
        "generation_command": "pytest tests/test_claim_registry.py",
        "rebuild_eligibility": "rebuildable",
        "claim_ids": '["test-claim"]',
    }


def _write_inventory(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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


def test_claim_type_is_required(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    claim = _claim(artifact)
    del claim["claim_type"]
    errors = validate_registry_document(_registry(claim), tmp_path, SCHEMA)
    assert any("claim_type" in error and "required" in error for error in errors)


def test_validated_workflow_requires_result_or_validation_receipt(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("test source only\n", encoding="utf-8")
    (tmp_path / "claim.md").write_text(
        "<!-- viralscan-claim:test-claim status=validated_v3 -->\n", encoding="utf-8"
    )
    claim = _claim(artifact)
    claim["claim_type"] = "workflow_execution"
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(_registry(claim)), encoding="utf-8")
    inventory_path = tmp_path / "inventory.tsv"
    _write_inventory(inventory_path, [_inventory_row(artifact, kind="command")])

    errors = validate_registry_file(
        registry_path,
        repo_root=tmp_path,
        schema_path=SCHEMA,
        inventory_path=inventory_path,
    )
    assert any("result or validation_receipt" in error for error in errors)

    _write_inventory(inventory_path, [_inventory_row(artifact, kind="validation_receipt")])
    errors = validate_registry_file(
        registry_path,
        repo_root=tmp_path,
        schema_path=SCHEMA,
        inventory_path=inventory_path,
    )
    assert not any("result or validation_receipt" in error for error in errors)


def test_claim_source_requires_reciprocal_inventory_linkage(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    (tmp_path / "claim.md").write_text(
        "<!-- viralscan-claim:test-claim status=validated_v3 -->\n", encoding="utf-8"
    )
    claim = _claim(artifact)
    claim["claim_type"] = "software_invariant"
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(_registry(claim)), encoding="utf-8")
    inventory_path = tmp_path / "inventory.tsv"
    _write_inventory(inventory_path, [_inventory_row(artifact)])

    errors = validate_registry_file(
        registry_path,
        repo_root=tmp_path,
        schema_path=SCHEMA,
        inventory_path=inventory_path,
    )
    assert any("source document 'claim.md' is absent from inventory" in error for error in errors)


def test_claim_dependencies_and_git_snapshot_are_validated(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("evidence\n", encoding="utf-8")
    (tmp_path / "claim.md").write_text(
        "<!-- viralscan-claim:test-claim status=validated_v3 -->\n", encoding="utf-8"
    )
    claim = _claim(artifact)
    claim["claim_type"] = "software_invariant"
    claim["input_hashes"] = [{"id": "unknown", "sha256": "0" * 64}]
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(_registry(claim)), encoding="utf-8")
    inventory_path = tmp_path / "inventory.tsv"
    _write_inventory(inventory_path, [_inventory_row(artifact)])

    errors = validate_registry_file(
        registry_path,
        repo_root=tmp_path,
        schema_path=SCHEMA,
        inventory_path=inventory_path,
    )
    assert any("unknown input_hashes id 'unknown'" in error for error in errors)
    assert any("git commit" in error and "not found" in error for error in errors)

    claim["input_hashes"] = [{"id": "fixture", "sha256": "0" * 64}]
    registry_path.write_text(json.dumps(_registry(claim)), encoding="utf-8")
    errors = validate_registry_file(
        registry_path,
        repo_root=tmp_path,
        schema_path=SCHEMA,
        inventory_path=inventory_path,
    )
    assert any("input_hashes digest mismatch for 'fixture'" in error for error in errors)


def test_private_incomplete_ebv_validation_cannot_leak_into_public_docs(tmp_path: Path) -> None:
    public_doc = tmp_path / "public.md"
    public_doc.write_text(
        "The retained EBV baseline verifies conservation on 103,145,071 BUS records.\n",
        encoding="utf-8",
    )
    claim = {
        **_claim(public_doc),
        "id": "v3-ebv-molecule-baseline",
        "claim_type": "scientific_result",
        "public": False,
        "status": "provenance_incomplete",
        "source_locations": [],
        "artifacts": [],
        "git_sha": "not_applicable",
        "input_hashes": [],
        "reason": "No complete immutable provenance is registered.",
    }
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(_registry(claim)), encoding="utf-8")
    inventory_path = tmp_path / "inventory.tsv"
    _write_inventory(inventory_path, [])
    scope = {
        "schema_version": "1.1.0",
        "wheel": {"members": []},
        "sdist": {"members": ["public.md"]},
        "docker_context": {"members": []},
        "public_docs": ["public.md"],
        "claim_bearing": [],
    }
    scope_path = tmp_path / "scope.json"
    scope_path.write_text(json.dumps(scope), encoding="utf-8")

    errors = validate_registry_file(
        registry_path,
        repo_root=tmp_path,
        schema_path=SCHEMA,
        scope_path=scope_path,
        inventory_path=inventory_path,
        coverage=True,
    )
    assert any("103,145,071" in error for error in errors)
    assert any("retained EBV baseline verifies" in error for error in errors)


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
