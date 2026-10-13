"""SCI-05 round 12 protocol-side regressions (DEV-027): M1, M5, M6, B1/M7, m1, m2."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from scripts.validate_v3_protocol import (
    DEFAULT_LEDGER,
    frozen_inputs_sha256,
    load_yaml,
    pre_freeze_inventory_sha256,
    schema_sha256,
    validate_protocol,
    validate_protocol_file,
    verify_inventory_evidence,
    workflow_matrix_sha256,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = ROOT / "analysis" / "v3_validation" / "protocol.yaml"
SCHEMA_PATH = ROOT / "schemas" / "v3" / "validation_protocol.schema.json"
PACKAGED_SCHEMA = ROOT / "src" / "viralscan" / "schemas" / "v3" / "validation_protocol.schema.json"
V03 = "influenza_gse336913"


@pytest.fixture
def protocol() -> dict:
    return load_yaml(PROTOCOL_PATH)


@pytest.fixture
def schema() -> dict:
    return load_yaml(SCHEMA_PATH)


@pytest.fixture
def ledger() -> dict:
    return load_yaml(DEFAULT_LEDGER)


def _dataset(document: dict, dataset_id: str) -> dict:
    return next(item for item in document["datasets"] if item["id"] == dataset_id)


def _redigest(document: dict) -> None:
    document["frozen_inputs"]["contract_sha256"] = frozen_inputs_sha256(document)


def _redigest_matrix(document: dict) -> None:
    matrix = document["workflow_matrix"]
    matrix["contract_sha256"] = workflow_matrix_sha256(matrix)


def test_the_committed_protocol_validates_with_its_ledger_and_evidence() -> None:
    assert validate_protocol_file() == []


# --- M1: the reviewer's tamper must fail --------------------------------------------------


def _flip_role(document: dict) -> None:
    _dataset(document, "hhv6b_srr20710641")["analysis_role"] = "training-and-holdout"


def _reword_planting_rule(document: dict) -> None:
    rule = _dataset(document, "synthetic_factorial")["planting_rule"]
    assert "pbmc_10x_healthy_v3 is excluded as a" in rule
    _dataset(document, "synthetic_factorial")["planting_rule"] = rule.replace(
        "pbmc_10x_healthy_v3 is excluded as a", "pbmc_10x_healthy_v3 is a permitted"
    )


def _relabel(document: dict) -> None:
    _dataset(document, "ebv_srr12682296")["label"] = "EBV scored positive"


def _bump_version(document: dict) -> None:
    document["protocol_version"] = "0.3.1"


def _open_the_gate(document: dict) -> None:
    document["freeze"]["training_execution_allowed"] = True


def _rederive_seeds(document: dict) -> None:
    document["seeds"]["derivation"] += " Also hash the moon phase."


TAMPERS = [
    _flip_role,
    _reword_planting_rule,
    _relabel,
    _bump_version,
    _open_the_gate,
    _rederive_seeds,
]


@pytest.mark.parametrize("tamper", TAMPERS, ids=lambda fn: fn.__name__)
def test_a_tamper_outside_the_old_digest_scope_is_detected(
    protocol: dict, schema: dict, ledger: dict, tamper
) -> None:
    tampered = deepcopy(protocol)
    tamper(tampered)

    errors = validate_protocol(tampered, schema, ledger=ledger)

    assert any("frozen_inputs contract_sha256 does not match" in error for error in errors)


@pytest.mark.parametrize("tamper", TAMPERS, ids=lambda fn: fn.__name__)
def test_a_redigested_tamper_is_an_undocumented_amendment(
    protocol: dict, schema: dict, ledger: dict, tamper
) -> None:
    tampered = deepcopy(protocol)
    tamper(tampered)
    _redigest(tampered)

    errors = validate_protocol(tampered, schema, ledger=ledger)

    assert any("undocumented amendment" in error for error in errors)


def test_an_asset_pending_reason_rewrite_is_not_an_amendment(
    protocol: dict, schema: dict, ledger: dict
) -> None:
    edited = deepcopy(protocol)
    _dataset(edited, "kshv_gse190558")["assets"][0]["pending_reason"] = "Reworded prose only."

    assert validate_protocol(edited, schema, ledger=ledger) == []


# --- M5: workflow-row coverage -----------------------------------------------------------


def test_every_required_dataset_has_a_workflow_row_in_the_committed_protocol(
    protocol: dict,
) -> None:
    covered = {
        dataset_id
        for workflow in protocol["workflow_matrix"]["workflows"]
        for dataset_id in workflow["dataset_ids"]
    }
    required = {item["id"] for item in protocol["datasets"] if item["required"]}
    assert required <= covered
    assert (
        "synthetic_mixed_host_virus" in protocol["workflow_matrix"]["workflows"][0]["dataset_ids"]
    )


def test_dropping_a_required_dataset_row_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    for workflow in tampered["workflow_matrix"]["workflows"]:
        workflow["dataset_ids"] = [
            item for item in workflow["dataset_ids"] if item != "synthetic_mixed_host_virus"
        ]
    tampered["workflow_matrix"]["expected_row_count"] -= 1
    _redigest_matrix(tampered)

    errors = validate_protocol(tampered, schema)

    assert any(
        "required dataset 'synthetic_mixed_host_virus' is referenced by no workflow row" in e
        for e in errors
    )


def test_a_context_row_without_rows_must_be_explicitly_marked(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    del _dataset(tampered, V03)["no_workflow_rows_reason"]

    errors = validate_protocol(tampered, schema)

    assert any(f"dataset '{V03}' has no workflow row and no explicit" in e for e in errors)


def test_a_stale_marker_on_a_dataset_with_rows_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    _dataset(tampered, "ebv_srr12682296")["no_workflow_rows_reason"] = "stale"

    errors = validate_protocol(tampered, schema)

    assert any("declares no_workflow_rows_reason but is referenced by" in e for e in errors)


def test_a_required_dataset_cannot_hide_behind_the_marker(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    _dataset(tampered, V03)["required"] = True

    errors = validate_protocol(tampered, schema)

    assert any(f"required dataset '{V03}'" in e for e in errors)


# --- B1 / M7 / Q80: V03 is a pending context row, not a scored T2 row -------------------


def test_v03_is_a_pending_context_row(protocol: dict) -> None:
    v03 = _dataset(protocol, V03)
    assert v03["analysis_role"] == "context-only"
    assert v03["truth_status"] == "contextual"
    assert v03["class"] == "public-positive"
    assert v03["required"] is False
    assert all(asset["digest_status"] == "pending" for asset in v03["assets"])
    for workflow in protocol["workflow_matrix"]["workflows"]:
        assert V03 not in workflow["dataset_ids"]
    reason = v03["no_workflow_rows_reason"]
    for needed in ("PRE-01", "5-prime read model", "read structure only", "counts-only"):
        assert needed in reason


def test_no_t2_endpoint_or_scored_role_was_introduced(protocol: dict, schema: dict) -> None:
    ids = {item["id"] for item in protocol["endpoints"]}
    # E11 is the descriptive real-data endpoint (DEV-029); no scored T2 endpoint exists.
    assert {item for item in ids if item.startswith("E11")} <= {"E11_real_data_descriptive"}
    denominators = {item["id"] for item in protocol["harmonization"]["denominators"]}
    assert not any(item.startswith("D26") for item in denominators)
    truth_values = schema["$defs"]["dataset"]["properties"]["truth_status"]["enum"]
    assert "design-known" not in truth_values


def test_dev_027_covers_every_digest_scope_it_moves(ledger: dict) -> None:
    record = next(item for item in ledger["deviations"] if item["deviation_id"] == "DEV-027")
    assert record["digest_scope"] == "harmonization.contract_sha256"
    assert set(record["additional_digest_changes"]) == {
        "workflow_matrix.contract_sha256",
        "partitions.contract_sha256",
        "frozen_inputs.contract_sha256",
    }
    assert record["git_sha"] == "88f02a1"
    assert record["outcome_triggered"] is False
    assert {c["deviation_id"] for c in record["field_corrections"]} == {
        f"DEV-0{n}" for n in range(20, 27)
    }
    assert {c["field"] for c in record["field_corrections"]} == {"outcome_triggered_justification"}
    assert record["user_decision_needed"]


# --- M6: outcome definition and pre-freeze inventory --------------------------------------


def _with_entry(protocol: dict, **entry) -> dict:
    tampered = deepcopy(protocol)
    section = tampered["pre_freeze_runs"]
    section["inventory"].append(
        {
            "id": "PFR-99",
            "kind": "outcome-generating",
            "period": "2026-10",
            "protocol_datasets": [],
            "other_datasets": [],
            "summary": "synthetic regression entry",
            "evidence": [],
            **entry,
        }
    )
    section["inventory_sha256"] = pre_freeze_inventory_sha256(section)
    return tampered


def test_an_outcome_run_on_a_scored_dataset_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = _with_entry(protocol, protocol_datasets=["pbmc_10x_healthy_v3"])

    errors = validate_protocol(tampered, schema)

    assert any("whose analysis_role 'external-evaluation' is scored" in e for e in errors)


def test_an_outcome_run_on_a_pending_dataset_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = _with_entry(protocol, protocol_datasets=[V03])

    errors = validate_protocol(tampered, schema)

    assert any("declared pending and must stay unseen" in e for e in errors)


def test_a_download_of_a_pending_dataset_is_not_an_outcome_run(
    protocol: dict, schema: dict
) -> None:
    tampered = _with_entry(protocol, kind="data-acquisition", protocol_datasets=[V03])

    errors = validate_protocol(tampered, schema)

    assert not any("pre-freeze" in e for e in errors)


def test_the_committed_inventory_has_no_outcome_run_on_a_scored_dataset(protocol: dict) -> None:
    datasets = {item["id"]: item for item in protocol["datasets"]}
    for entry in protocol["pre_freeze_runs"]["inventory"]:
        if entry["kind"] != "outcome-generating":
            continue
        for dataset_id in entry["protocol_datasets"]:
            assert datasets[dataset_id]["analysis_role"] == "context-only", entry["id"]
            assert datasets[dataset_id]["analysis_role"] != "training-and-holdout"


def test_an_edited_inventory_fails_its_own_hash(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["pre_freeze_runs"]["inventory"][0]["summary"] += " (edited)"

    errors = validate_protocol(tampered, schema)

    assert any("inventory_sha256 does not match the inventory" in e for e in errors)


def test_repo_evidence_that_no_longer_hashes_is_reported(protocol: dict) -> None:
    assert verify_inventory_evidence(protocol) == []
    tampered = deepcopy(protocol)
    evidence = next(
        item
        for entry in tampered["pre_freeze_runs"]["inventory"]
        for item in entry["evidence"]
        if item["location"] == "repo"
    )
    evidence["sha256"] = "0" * 64

    errors = verify_inventory_evidence(tampered)

    assert any("no longer hashes to its recorded sha256" in e for e in errors)


def test_the_inventory_is_covered_by_the_frozen_inputs_digest(protocol: dict) -> None:
    tampered = deepcopy(protocol)
    del tampered["pre_freeze_runs"]["inventory"][-1]

    assert frozen_inputs_sha256(tampered) != protocol["frozen_inputs"]["contract_sha256"]


# --- code-side M4 request: B15 structured, and inside the digest ----------------------------


def test_b15_empty_droplet_values_are_structured_and_match_the_prose(protocol: dict) -> None:
    structure = protocol["partitions"]["generated_sample_structure"]
    assert structure["empty_barcodes_per_sample"] == 3000
    assert structure["empty_median_umi"] == 50
    prose = protocol["harmonization"]["cell_universe"]["synthetic_rule"]
    assert "3,000 empty barcodes" in prose
    assert "median 50" in prose


def test_the_generator_structure_is_covered_by_the_frozen_inputs_digest(
    protocol: dict, schema: dict, ledger: dict
) -> None:
    tampered = deepcopy(protocol)
    tampered["partitions"]["generated_sample_structure"]["empty_barcodes_per_sample"] = 10
    assert frozen_inputs_sha256(tampered) != protocol["frozen_inputs"]["contract_sha256"]
    errors = validate_protocol(tampered, schema, ledger=ledger)
    assert any("partitions contract_sha256 does not match" in e for e in errors)


# --- m2: the schema is pinned -------------------------------------------------------------


def test_both_schema_copies_are_identical_and_match_the_pin(protocol: dict) -> None:
    assert PACKAGED_SCHEMA.read_bytes() == SCHEMA_PATH.read_bytes()
    assert schema_sha256(load_yaml(PACKAGED_SCHEMA)) == protocol["frozen_inputs"]["schema_sha256"]


def test_a_widened_schema_enum_is_detected(protocol: dict, schema: dict, ledger: dict) -> None:
    widened = deepcopy(schema)
    widened["$defs"]["dataset"]["properties"]["analysis_role"]["enum"].append("scored-t2")

    errors = validate_protocol(protocol, widened, ledger=ledger)

    assert any("frozen_inputs.schema_sha256" in e for e in errors)


# --- m1: the append-only anchor is not vacuous on a push to main --------------------------


def test_ci_compares_the_ledger_with_the_previous_head_on_push() -> None:
    workflow = yaml.safe_load((ROOT / ".github" / "workflows" / "ci.yml").read_text())
    steps = workflow["jobs"]["lint"]["steps"]
    step = next(item for item in steps if item.get("name") == "deviation ledger is append-only")
    assert "github.event.before" in step["env"]["LEDGER_BASE"]
    assert "origin/main" in step["env"]["LEDGER_BASE"]
    assert '"$LEDGER_BASE"' in step["run"]
    assert "--revision origin/main" not in step["run"]
