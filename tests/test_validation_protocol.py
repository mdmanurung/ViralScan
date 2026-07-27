from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from scripts.validate_v3_protocol import (
    DEFAULT_LEDGER,
    calibration_sha256,
    failure_reporting_sha256,
    frozen_inputs_sha256,
    harmonization_dependencies_sha256,
    harmonization_sha256,
    load_yaml,
    partitions_sha256,
    require_execution_allowed,
    validate_amendment_ledger,
    validate_protocol,
    validate_protocol_file,
    workflow_matrix_sha256,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = ROOT / "analysis" / "v3_validation" / "protocol.yaml"
SCHEMA_PATH = ROOT / "schemas" / "v3" / "validation_protocol.schema.json"


@pytest.fixture
def protocol() -> dict:
    return load_yaml(PROTOCOL_PATH)


@pytest.fixture
def schema() -> dict:
    return load_yaml(SCHEMA_PATH)


def test_committed_protocol_is_a_valid_non_executable_draft(protocol: dict, schema: dict) -> None:
    assert validate_protocol(protocol, schema) == []
    assert protocol["status"] == "draft"
    assert protocol["execution_readiness"]["training_allowed"] is False
    assert protocol["execution_readiness"]["holdout_allowed"] is False
    assert protocol["execution_readiness"]["training_blockers"]
    assert protocol["execution_readiness"]["holdout_blockers"]


def test_draft_cannot_pass_the_frozen_execution_gate(protocol: dict, schema: dict) -> None:
    errors = validate_protocol(protocol, schema, require_frozen=True)
    assert "protocol is not frozen" in errors
    assert "training execution blockers remain" in errors
    assert any("unfrozen asset" in error for error in errors)


def test_holdout_has_a_stricter_gate_than_training(protocol: dict, schema: dict) -> None:
    errors = validate_protocol(protocol, schema, phase="holdout")
    assert "holdout execution is not allowed" in errors
    assert "holdout execution blockers remain" in errors
    assert "holdout requires frozen training_results_sha256" in errors
    assert "holdout requires frozen thresholds_sha256" in errors


def test_schema_rejects_unknown_top_level_keys(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    broken["outcome_selected_barcodes"] = True
    errors = validate_protocol(broken, schema)
    assert any("Additional properties are not allowed" in error for error in errors)


def test_schema_rejects_deleting_planned_freeze_sections(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    del broken["planned_freeze_sections"]
    errors = validate_protocol(broken, schema, phase="training")
    assert any("'planned_freeze_sections' is a required property" in error for error in errors)


def test_schema_rejects_deleting_harmonization_contract(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    del broken["harmonization"]
    errors = validate_protocol(broken, schema)
    assert any("'harmonization' is a required property" in error for error in errors)


def test_semantic_validation_rejects_duplicate_ids(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    broken["datasets"].append(deepcopy(broken["datasets"][0]))
    errors = validate_protocol(broken, schema)
    assert "datasets contains duplicate id 'synthetic_factorial'" in errors


def test_semantic_validation_rejects_wrong_chemistry_geometry(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    broken["supported_scope"]["chemistries"][0]["umi_length"] = 12
    errors = validate_protocol(broken, schema)
    assert "chemistry '10xv2' geometry (16, 12) does not match (16, 10)" in errors


def test_verified_assets_require_a_sha256(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    asset = broken["datasets"][0]["assets"][0]
    asset["digest_status"] = "verified"
    asset["sha256"] = None
    errors = validate_protocol(broken, schema)
    assert any("is not of type 'string'" in error for error in errors)


def test_duplicate_yaml_keys_fail_closed(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.yaml"
    path.write_text("schema_version: '3.0.0'\nschema_version: '3.0.0'\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate YAML key"):
        load_yaml(path)


def test_public_control_geometry_and_recorded_hashes(protocol: dict) -> None:
    datasets = {item["id"]: item for item in protocol["datasets"]}
    expected = {
        "hhv6b_srr20710641": (
            "10xv3",
            "a2ade8f967394938c49f33dba837ba40d0bfb79702249d5235e3421b2899b866",
        ),
        "ebv_srr12682296": (
            "10xv2",
            "f3a4929c6d3dcd2d8d7b86918705c9ba687415d25179fc53d265476cf2e4d412",
        ),
        "hsv1_srr8315713": (
            "drop-seq",
            "7b4aefc295c4971c08725ed739046c5fa5601b2ff50a505134572dc934291fe7",
        ),
    }
    for dataset_id, (chemistry, r1_sha256) in expected.items():
        dataset = datasets[dataset_id]
        assert dataset["chemistries"] == [chemistry]
        assert dataset["assets"][0]["sha256"] == r1_sha256
        assert dataset["assets"][0]["digest_status"] == "recorded-reverify"


def test_protocol_contains_no_institutional_absolute_paths() -> None:
    text = PROTOCOL_PATH.read_text(encoding="utf-8")
    assert "/exports/" not in text
    assert "/home/" not in text


def test_estimation_objectives_are_not_mislabeled_as_falsifiable(protocol: dict) -> None:
    hypotheses = {item["id"]: item for item in protocol["hypotheses"]}
    assert hypotheses["H3_positive_recovery"]["class"] == "estimation"
    assert hypotheses["H7_workflow_tradeoff"]["class"] == "estimation"


def test_sci02_harmonization_decisions_are_frozen(protocol: dict) -> None:
    harmonization = protocol["harmonization"]
    assert harmonization["status"] == "frozen"
    assert harmonization["outcome_blinding"]["outcome_independent"] is True
    assert harmonization["contract_sha256"] == harmonization_sha256(harmonization)
    assert harmonization["dependent_fields_sha256"] == (harmonization_dependencies_sha256(protocol))
    assert "viral-count-matrix" in harmonization["outcome_blinding"]["prohibited_inputs"]
    for section_name in (
        "cell_universe",
        "feature_universe",
        "count_layers_denominators",
    ):
        assert protocol["planned_freeze_sections"][section_name]["status"] == "frozen"
    assert "harmonization" not in {
        item["id"] for item in protocol["execution_readiness"]["training_blockers"]
    }


def test_schema_rejects_outcome_selected_primary_cell_source(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    broken["harmonization"]["cell_universe"]["source_precedence"][2] = (
        "viral-positive-barcode-intersection"
    )
    errors = validate_protocol(broken, schema)
    assert any(
        "schema harmonization.cell_universe.source_precedence.2" in error
        and "'shared-host-only-emptydrops' was expected" in error
        for error in errors
    )


def test_schema_rejects_observed_feature_universe(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    broken["harmonization"]["feature_universe"]["source"] = (
        "intersection-of-nonzero-observed-features"
    )
    errors = validate_protocol(broken, schema)
    assert any(
        "schema harmonization.feature_universe.source" in error
        and "frozen-reference-and-tool-capability-manifests" in error
        for error in errors
    )


def test_cross_tool_primary_requires_integer_gene_unique_molecules(
    protocol: dict, schema: dict
) -> None:
    primary = protocol["harmonization"]["count_layers"]["cross_tool_primary"]
    assert primary["viralscan_source"] == "counts_unique"
    assert primary["unit"] == "integer-unique-molecule"
    assert primary["adapter_unit"] == "molecule"
    assert primary["adapter_ambiguity_rule"] == "gene-unique-only"

    broken = deepcopy(protocol)
    broken["harmonization"]["count_layers"]["cross_tool_primary"]["unit"] = (
        "fractional-umi-estimate"
    )
    errors = validate_protocol(broken, schema)
    assert any(
        "schema harmonization.count_layers.cross_tool_primary.unit" in error
        and "'integer-unique-molecule' was expected" in error
        for error in errors
    )


def test_semantic_validation_rejects_cross_section_layer_drift(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    broken["supported_scope"]["cross_tool_primary_layer"] = "counts_unique"
    broken["harmonization"]["count_layers"]["cross_tool_primary"]["viralscan_source"] = "X"
    errors = validate_protocol(broken, schema)
    assert any(
        "schema harmonization.count_layers.cross_tool_primary.viralscan_source" in error
        for error in errors
    )


def test_semantic_validation_rejects_duplicate_denominator_ids(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    denominators = broken["harmonization"]["denominators"]
    denominators[-1]["id"] = denominators[0]["id"]
    errors = validate_protocol(broken, schema)
    assert "harmonization denominators contain duplicate ids" in errors
    assert any("denominator registry mismatch" in error for error in errors)


def test_semantic_validation_rejects_positive_rate_denominator_drift(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    denominator = next(
        item
        for item in broken["harmonization"]["denominators"]
        if item["id"] == "D12_positive_cell_rate"
    )
    denominator["denominator"] = "method-positive-cells-only"
    broken["harmonization"]["contract_sha256"] = harmonization_sha256(broken["harmonization"])
    amended_schema = deepcopy(schema)
    amended_schema["$defs"]["harmonization"]["properties"]["contract_sha256"]["const"] = broken[
        "harmonization"
    ]["contract_sha256"]
    errors = validate_protocol(broken, amended_schema)
    assert any("denominator 'D12_positive_cell_rate' does not match" in error for error in errors)


def test_semantic_validation_rejects_sci02_freeze_marker_drift(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    broken["planned_freeze_sections"]["cell_universe"]["status"] = "pending"
    errors = validate_protocol(broken, schema)
    assert "frozen harmonization requires planned section 'cell_universe' to be frozen" in errors


def test_semantic_validation_rejects_stale_harmonization_blocker(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    broken["execution_readiness"]["training_blockers"].append(
        {
            "id": "harmonization",
            "description": "stale",
            "resolution_task": "SCI-02",
        }
    )
    errors = validate_protocol(broken, schema)
    assert "frozen harmonization cannot remain a training blocker" in errors


def test_cell_anchor_contract_has_no_silent_fallback(protocol: dict) -> None:
    cells = protocol["harmonization"]["cell_universe"]
    emptydrops = cells["host_only_emptydrops"]
    assert emptydrops == {
        "caller": "DropletUtils::emptyDrops",
        "matrix": "raw-host-only-gene-molecule-matrix",
        "fdr": 0.01,
        "lower": 100,
        "niters": 10000,
        "seed_source": "seeds.cell_calling",
        "fallback": "forbidden",
        "version_requirement": "exact-version-and-container-digest-in-SCI-04-workflow-row",
    }
    assert cells["common_across_workflow_rows"] is True
    assert cells["missing_anchor_barcode_policy"].endswith("otherwise-row-failure")


def test_barcode_contract_is_library_scoped_and_collision_fatal(protocol: dict) -> None:
    normalization = protocol["harmonization"]["barcode_normalization"]
    assert normalization["scope"] == "within-library-only"
    assert normalization["cross_library_matching"] == "forbidden"
    assert normalization["canonical_alphabet_regex"] == "^[ACGT]+$"
    assert normalization["chemistry_length_source"] == "supported-scope-chemistry"
    assert normalization["collision_policy"].startswith("fail-dataset-preflight")


def test_missing_workflow_outputs_remain_null_failures(protocol: dict) -> None:
    policy = protocol["harmonization"]["row_failure_policy"]
    assert policy["retain_all_prespecified_rows"] is True
    assert policy["include_in_failure_denominator"] is True
    assert policy["impute_accuracy"] is False
    assert policy["noncomplete_metric_value"] is None
    assert set(policy["row_execution_statuses"]) == {
        "not_started",
        "complete",
        "failed",
        "timeout",
        "oom",
        "invalid_input",
        "unsupported",
    }
    assert set(policy["endpoint_comparability_statuses"]) == {
        "comparable",
        "incomparable",
        "not_applicable",
    }


def test_harmonization_audits_include_identity_hash_columns(protocol: dict) -> None:
    artifacts = {
        item["path"]: set(item["required_columns"])
        for item in protocol["harmonization"]["audit_artifacts"]
    }
    assert {
        "source_sha256",
        "canonical_key",
        "universe_role",
    } <= artifacts["analysis/v3_validation/generated/cell_universe_manifest.tsv"]
    assert {
        "canonical_gene_id",
        "canonical_virus_id",
        "sequence_sha256",
        "mapping_cardinality",
    } <= artifacts["analysis/v3_validation/generated/feature_universe.tsv"]
    assert {
        "source_artifact_sha256",
        "feature_manifest_sha256",
        "umi_collapse_semantics",
        "ambiguity_rule",
    } <= artifacts["analysis/v3_validation/generated/count_layer_map.tsv"]
    assert {
        "cell_universe_sha256",
        "numerator_feature_universe_sha256",
        "denominator_feature_universe_sha256",
        "truth_manifest_sha256",
        "layer_map_sha256",
    } <= artifacts["analysis/v3_validation/generated/denominator_audit.tsv"]


def test_cell_recovery_has_distinct_precision_recall_and_auprc_populations(
    protocol: dict,
) -> None:
    denominators = {item["id"]: item for item in protocol["harmonization"]["denominators"]}
    assert denominators["D4_infected_cell_precision"]["denominator"].startswith(
        "all called-positive cells"
    )
    assert denominators["D5_infected_cell_recall"]["denominator"].startswith(
        "all exact truth-infected cells"
    )
    assert denominators["D6_infected_cell_auprc"]["denominator"].startswith(
        "every truth-labelled barcode"
    )
    assert (
        denominators["D6_infected_cell_auprc"]["count_layer"]
        == "continuous-cell-score-id-frozen-in-SCI-03"
    )
    assert protocol["harmonization"]["derived_metric_rules"]["cell_f1"].startswith(
        "two-times-precision"
    )


def test_every_endpoint_has_a_frozen_denominator_mapping(protocol: dict) -> None:
    endpoint_ids = {item["id"] for item in protocol["endpoints"]}
    mapping = protocol["harmonization"]["endpoint_denominator_map"]
    denominator_ids = {item["id"] for item in protocol["harmonization"]["denominators"]}
    assert set(mapping) == endpoint_ids
    assert all(set(ids) <= denominator_ids for ids in mapping.values())
    assert {
        "D7_exact_negative_cell_false_calls",
        "D8_exact_negative_sample_false_calls",
        "D9_exact_negative_molecule_false_calls",
    } <= set(mapping["E2_negative_false_calls"])
    assert {
        "D15_two_step_truth_read_loss",
        "D16_two_step_truth_molecule_loss",
    } == set(mapping["E7_two_step_loss"])


def test_presumed_negatives_are_not_used_as_false_positive_truth(protocol: dict) -> None:
    denominators = {item["id"]: item for item in protocol["harmonization"]["denominators"]}
    exact = denominators["D7_exact_negative_cell_false_calls"]
    contextual = denominators["D10_presumed_negative_unexpected_calls"]
    assert exact["class"] == "primary"
    assert "exact synthetic" in exact["population"]
    assert contextual["class"] == "secondary"
    assert "presumed-negative" in contextual["population"]
    assert "false" not in contextual["numerator"]


def test_molecule_truth_projection_scores_corrected_keys_once(protocol: dict) -> None:
    projection = protocol["harmonization"]["molecule_truth_projection"]
    assert projection["corrected_key"].endswith("corrected-UMI")
    assert projection["correction_contract"] == (
        "outcome-independent-canonical-truth-normalizer-frozen-by-VAL-09"
    )
    assert projection["gene_projection"] == "ECs-to-distinct-compatible-genes"
    assert projection["collision_rule"] == (
        "one-corrected-key-contributes-at-most-one-truth-molecule"
    )
    assert set(projection["truth_classes"]) == {"unique", "ambiguous", "unresolved"}
    assert "false negative" in projection["adapter_rule"]


def test_sibling_confusion_keeps_molecule_and_cell_units_separate(
    protocol: dict,
) -> None:
    denominators = {item["id"]: item for item in protocol["harmonization"]["denominators"]}
    molecule = denominators["D13_sibling_molecule_confusion"]
    cell = denominators["D24_sibling_cell_confusion"]
    assert molecule["count_layer"] == "counts_unique-plus-molecule-ambiguity-audit"
    assert "canonical planted target molecule" in molecule["denominator"]
    assert cell["count_layer"].startswith("viralscan-product-evidence-X")
    assert "truth-labelled anchor cell" in cell["denominator"]
    assert set(protocol["harmonization"]["endpoint_denominator_map"]["E5_sibling_confusion"]) == {
        molecule["id"],
        cell["id"],
    }


def test_filter_boundary_loss_uses_exact_lineage_not_downstream_counts(
    protocol: dict,
) -> None:
    denominators = {item["id"]: item for item in protocol["harmonization"]["denominators"]}
    molecule_loss = denominators["D16_two_step_truth_molecule_loss"]
    assert molecule_loss["count_layer"] == "exact-read-lineage-molecule-survival"
    assert "no planted target fragment survives" in molecule_loss["numerator"]


def test_unattempted_rows_and_not_applicable_endpoints_have_explicit_denominators(
    protocol: dict,
) -> None:
    denominators = {item["id"]: item for item in protocol["harmonization"]["denominators"]}
    assert "including unattempted" in denominators["D1_count_invariant_runs"]["denominator"]
    assert (
        "excluding predeclared not-applicable"
        in denominators["D23_endpoint_incomparability_rate"]["denominator"]
    )


def test_frozen_contract_digest_rejects_normative_prose_drift(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    broken["harmonization"]["cell_universe"]["public_data_rule"] = "select viral-positive barcodes"
    errors = validate_protocol(broken, schema)
    assert "harmonization contract_sha256 does not match the canonical SCI-02 contract" in errors


def test_dependency_digest_rejects_endpoint_prose_that_selects_positive_cells(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    endpoint = next(item for item in broken["endpoints"] if item["id"] == "E4_cell_recovery")
    endpoint["reporting"] = "Compute only over ViralScan-positive barcodes."
    errors = validate_protocol(broken, schema)
    assert (
        "harmonization dependent_fields_sha256 does not match "
        "SCI-02-dependent protocol fields" in errors
    )


def test_dependency_digest_rejects_outcome_based_dataset_exclusion(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    dataset = next(item for item in broken["datasets"] if item["id"] == "synthetic_host_only")
    dataset["exclusion_rule"] = "Exclude replicates containing probable viral calls."
    errors = validate_protocol(broken, schema)
    assert any("dependent_fields_sha256" in error for error in errors)


def test_dependency_digest_rejects_hypothesis_outcome_selection(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    hypothesis = next(item for item in broken["hypotheses"] if item["id"] == "H3_positive_recovery")
    hypothesis["statement"] = "Select cells after viral outcomes are observed."
    errors = validate_protocol(broken, schema)
    assert any("dependent_fields_sha256" in error for error in errors)


def test_frozen_schema_rejects_rehashed_harmonization_amendment(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    broken["harmonization"]["feature_universe"]["primary_rule"] = (
        "keep only nonzero observed viral features"
    )
    broken["harmonization"]["contract_sha256"] = harmonization_sha256(broken["harmonization"])
    # Read the pinned digest from the schema rather than hard-coding it, so an
    # amendment recorded in the ledger updates one place, not two.
    pinned = schema["$defs"]["harmonization"]["properties"]["contract_sha256"]["const"]
    errors = validate_protocol(broken, schema)
    assert any(
        "schema harmonization.contract_sha256" in error and pinned in error for error in errors
    )


def test_semantic_validation_rejects_rehashed_audit_column_loss(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    artifact = next(
        item
        for item in broken["harmonization"]["audit_artifacts"]
        if item["path"].endswith("feature_universe.tsv")
    )
    artifact["required_columns"].remove("canonical_gene_id")
    artifact["required_columns"].append("filler")
    broken["harmonization"]["contract_sha256"] = harmonization_sha256(broken["harmonization"])
    amended_schema = deepcopy(schema)
    amended_schema["$defs"]["harmonization"]["properties"]["contract_sha256"]["const"] = broken[
        "harmonization"
    ]["contract_sha256"]
    errors = validate_protocol(broken, amended_schema)
    assert "harmonization audit columns do not match the frozen SCI-02 registry" in errors


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("hypotheses",), None),
        (("endpoints",), None),
        (("harmonization", "denominators"), None),
        (("harmonization", "audit_artifacts"), None),
        (("execution_readiness", "training_blockers"), None),
        (("planned_freeze_sections",), None),
    ],
)
def test_malformed_documents_return_schema_errors_without_crashing(
    protocol: dict,
    schema: dict,
    path: tuple[str, ...],
    value: object,
) -> None:
    broken = deepcopy(protocol)
    target = broken
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    errors = validate_protocol(broken, schema, phase="training")
    assert errors
    assert all(error.startswith("schema ") for error in errors)


def test_semantic_validation_rejects_duplicate_nested_asset_ids(
    protocol: dict, schema: dict
) -> None:
    broken = deepcopy(protocol)
    assets = broken["datasets"][0]["assets"]
    assets.append(deepcopy(assets[0]))
    errors = validate_protocol(broken, schema)
    assert (
        "datasets.synthetic_factorial.assets contains duplicate id 'synthetic_truth_manifest'"
        in errors
    )


def test_semantic_validation_rejects_duplicate_blocker_ids(protocol: dict, schema: dict) -> None:
    broken = deepcopy(protocol)
    blockers = broken["execution_readiness"]["training_blockers"]
    duplicated_id = blockers[0]["id"]
    blockers.append(deepcopy(blockers[0]))
    errors = validate_protocol(broken, schema)
    assert (
        f"execution_readiness.training_blockers contains duplicate id {duplicated_id!r}" in errors
    )


def test_unhashable_yaml_mapping_key_fails_closed(tmp_path: Path) -> None:
    bad_protocol = tmp_path / "unhashable.yaml"
    bad_protocol.write_text("? [a, b]\n: value\n", encoding="utf-8")
    errors = validate_protocol_file(bad_protocol, SCHEMA_PATH)
    assert len(errors) == 1
    assert "unhashable YAML mapping key" in errors[0]


def test_invalid_custom_schema_fails_closed(tmp_path: Path) -> None:
    bad_schema = tmp_path / "bad-schema.yaml"
    bad_schema.write_text("type: definitely-not-a-json-schema-type\n", encoding="utf-8")
    errors = validate_protocol_file(PROTOCOL_PATH, bad_schema)
    assert len(errors) == 1
    assert errors[0].startswith("invalid validation schema:")


def test_sci03_sections_are_drafted_but_not_frozen(protocol: dict) -> None:
    """SCI-05 round 1 rejected these freezes; they must not claim frozen status."""
    partitions = protocol["partitions"]
    calibration = protocol["calibration"]

    assert partitions["status"] == "pending"
    assert calibration["status"] == "pending"
    assert partitions["frozen_at"] is None
    assert calibration["frozen_at"] is None
    assert partitions["contract_sha256"] == partitions_sha256(partitions)
    assert calibration["contract_sha256"] == calibration_sha256(calibration)

    planned = protocol["planned_freeze_sections"]
    assert planned["partitions"]["status"] == "pending"
    assert planned["calibration_metrics_lod"]["status"] == "pending"

    blocker_ids = {b["id"] for b in protocol["execution_readiness"]["training_blockers"]}
    assert "partitions_metrics" in blocker_ids
    assert "data_hashes" in blocker_ids


def test_a_frozen_partition_would_require_frozen_stratification_factors(
    protocol: dict,
) -> None:
    """SCI-05 F1: the split strata cannot be computed while factor levels are empty."""
    factor_status = {f["id"]: f for f in protocol["factors"]}
    unfrozen = [
        name
        for name in protocol["partitions"]["stratification_factors"]
        if factor_status[name].get("status") != "frozen"
    ]

    assert unfrozen, "if every factor is frozen, partitions may be frozen again"
    assert protocol["partitions"]["status"] == "pending"


def test_partitions_allocate_whole_biological_samples(protocol: dict) -> None:
    partitions = protocol["partitions"]
    assert "biological sample" in partitions["unit"].lower()
    assert 0 < partitions["holdout_fraction"] < 1
    assert partitions["holdout_evaluations_allowed"] == 1
    prohibitions = " ".join(partitions["leakage_prohibitions"]).lower()
    for leak in ("template", "locus", "molecule", "cell barcode"):
        assert leak in prohibitions, leak


def test_uncertainty_resamples_samples_not_cells(protocol: dict) -> None:
    uncertainty = protocol["calibration"]["uncertainty"]
    assert "sample" in uncertainty["resampling_unit"].lower()
    prohibited = " ".join(uncertainty["prohibited_units"]).lower()
    assert "cells treated as independent" in prohibited
    assert "molecules treated as independent" in prohibited
    assert uncertainty["replicates"] >= 1000


def test_limit_of_detection_forbids_extrapolation(protocol: dict) -> None:
    lod = protocol["calibration"]["limit_of_detection"]
    assert lod["extrapolation"].lower().startswith("prohibited")
    assert "0.95" in lod["reported_quantity"]


def test_threshold_search_is_training_only(protocol: dict) -> None:
    search = protocol["calibration"]["threshold_search"]
    assert "training" in search["data"].lower()
    assert "holdout" in search["data"].lower()
    prohibited = " ".join(search["prohibited"]).lower()
    assert "holdout" in prohibited


def test_tampering_with_the_partitions_contract_is_detected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["partitions"]["holdout_fraction"] = 0.5

    errors = validate_protocol(tampered, schema)

    assert any("partitions contract_sha256" in error for error in errors)


def test_tampering_with_the_calibration_contract_is_detected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["calibration"]["uncertainty"]["replicates"] = 10

    errors = validate_protocol(tampered, schema)

    assert any("calibration contract_sha256" in error for error in errors)


def test_frozen_partitions_require_frozen_seeds(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["partitions"]["status"] = "frozen"
    tampered["partitions"]["frozen_at"] = "2026-07-27"
    tampered["seeds"]["split"] = {
        "status": "pending",
        "value": None,
        "resolution_task": "SCI-03",
        "freeze_required": True,
        "pending_reason": "unfrozen for this test",
    }
    tampered["partitions"]["contract_sha256"] = partitions_sha256(tampered["partitions"])
    tampered["planned_freeze_sections"]["partitions"]["status"] = "frozen"

    errors = validate_protocol(tampered, schema)

    assert "frozen partitions require seeds.split to be frozen" in errors


def test_planned_section_cannot_be_frozen_before_its_contract(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["planned_freeze_sections"]["partitions"]["status"] = "frozen"

    errors = validate_protocol(tampered, schema)

    assert "planned section 'partitions' cannot be frozen before partitions" in errors


def test_calibration_metrics_must_reference_declared_endpoints(
    protocol: dict, schema: dict
) -> None:
    tampered = deepcopy(protocol)
    tampered["calibration"]["metrics"][0]["endpoint_id"] = "E99_does_not_exist"
    tampered["calibration"]["contract_sha256"] = calibration_sha256(tampered["calibration"])

    errors = validate_protocol(tampered, schema)

    assert any("references unknown endpoint" in error for error in errors)


def test_stratification_factors_must_be_declared_factors(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["partitions"]["stratification_factors"] = ["not_a_declared_factor"]
    tampered["partitions"]["contract_sha256"] = partitions_sha256(tampered["partitions"])

    errors = validate_protocol(tampered, schema)

    assert any("is not a declared factor" in error for error in errors)


def test_schema_rejects_deleting_the_partitions_contract(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    del tampered["partitions"]["contract_sha256"]

    errors = validate_protocol(tampered, schema)

    assert any("contract_sha256" in error for error in errors)


def test_packaged_schema_matches_the_canonical_schema() -> None:
    packaged = ROOT / "src" / "viralscan" / "schemas" / "v3" / "validation_protocol.schema.json"
    assert packaged.read_text(encoding="utf-8") == SCHEMA_PATH.read_text(encoding="utf-8")


def test_sci04_sections_are_drafted_but_not_frozen(protocol: dict) -> None:
    """SCI-05 round 1 rejected these freezes; they must not claim frozen status."""
    matrix = protocol["workflow_matrix"]
    reporting = protocol["failure_and_deviation_reporting"]

    assert matrix["status"] == "pending"
    assert reporting["status"] == "pending"
    assert matrix["contract_sha256"] == workflow_matrix_sha256(matrix)
    assert reporting["contract_sha256"] == failure_reporting_sha256(reporting)

    planned = protocol["planned_freeze_sections"]
    assert planned["workflow_matrix"]["status"] == "pending"
    assert planned["failure_and_deviation_reporting"]["status"] == "pending"

    blocker_ids = {b["id"] for b in protocol["execution_readiness"]["training_blockers"]}
    assert "workflow_rows" in blocker_ids
    assert "tool_environments" in blocker_ids
    assert "independent_review_round_2" in blocker_ids


def test_every_comparator_named_by_the_plan_has_a_workflow(protocol: dict) -> None:
    tools = {w["tool"] for w in protocol["workflow_matrix"]["workflows"]}
    assert tools == {
        "viralscan",
        "starsolo",
        "traditional-host-subtraction",
        "venus",
        "viral-track",
        "virtus",
        "viralscan-2.2.0",
    }


def test_kallisto_two_step_is_explicitly_excluded_with_a_revisit_condition(
    protocol: dict,
) -> None:
    excluded = {w["id"]: w for w in protocol["workflow_matrix"]["excluded_workflows"]}
    assert "X_kallisto_two_step" in excluded
    assert "exact fragment lineage" in excluded["X_kallisto_two_step"]["revisit_condition"]


def test_dedicated_comparators_cover_the_three_required_positives(protocol: dict) -> None:
    required = {"hhv6b_srr20710641", "ebv_srr12682296", "hsv1_srr8315713"}
    for workflow in protocol["workflow_matrix"]["workflows"]:
        if workflow["tool"] in {"venus", "viral-track", "virtus"}:
            assert required.issubset(set(workflow["dataset_ids"])), workflow["id"]
            # F5: the only exact-truth dataset must reach them too.
            assert "synthetic_factorial" in workflow["dataset_ids"], workflow["id"]


def test_row_count_drift_is_detected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["workflow_matrix"]["workflows"][0]["dataset_ids"].pop()
    tampered["workflow_matrix"]["contract_sha256"] = workflow_matrix_sha256(
        tampered["workflow_matrix"]
    )

    errors = validate_protocol(tampered, schema)

    assert any("does not match the" in error and "enumerated rows" in error for error in errors)


def test_workflow_referencing_an_unknown_dataset_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    tampered["workflow_matrix"]["workflows"][0]["dataset_ids"][0] = "not_a_dataset"
    tampered["workflow_matrix"]["contract_sha256"] = workflow_matrix_sha256(
        tampered["workflow_matrix"]
    )

    errors = validate_protocol(tampered, schema)

    assert any("references unknown dataset" in error for error in errors)


def test_workflow_referencing_an_unknown_reference_is_rejected(
    protocol: dict, schema: dict
) -> None:
    tampered = deepcopy(protocol)
    tampered["workflow_matrix"]["workflows"][0]["reference_id"] = "not_a_reference"
    tampered["workflow_matrix"]["contract_sha256"] = workflow_matrix_sha256(
        tampered["workflow_matrix"]
    )

    errors = validate_protocol(tampered, schema)

    assert any("references unknown reference" in error for error in errors)


def test_unpinned_environments_block_execution_but_not_the_freeze(
    protocol: dict, schema: dict
) -> None:
    """The matrix may freeze which rows exist before REL-03 supplies digests."""
    assert validate_protocol(protocol, schema) == []

    errors = validate_protocol(protocol, schema, phase="training")

    assert "workflow environment pinning is not frozen" in errors
    assert any("no pinned tool version or container digest" in error for error in errors)


def test_tampering_with_the_workflow_matrix_contract_is_detected(
    protocol: dict, schema: dict
) -> None:
    tampered = deepcopy(protocol)
    tampered["workflow_matrix"]["primary_comparison_rules"].append("anything goes")

    errors = validate_protocol(tampered, schema)

    assert any("workflow_matrix contract_sha256" in error for error in errors)


def test_failure_reporting_forbids_silent_row_loss(protocol: dict) -> None:
    reporting = protocol["failure_and_deviation_reporting"]
    assert reporting["row_failure_policy_ref"] == "harmonization.row_failure_policy"
    prohibited = " ".join(reporting["prohibited"]).lower()
    assert "deleting or omitting a planned row" in prohibited
    assert "imputing an accuracy value" in prohibited
    assert "outcome_triggered" in reporting["deviation_record"]["required_fields"]


def test_every_frozen_section_amendment_has_a_deviation_record() -> None:
    """The amendment_rule requires a ledger entry for any post-freeze change."""
    ledger = load_yaml(ROOT / "analysis" / "v3_validation" / "deviations.yaml")
    protocol = load_yaml(PROTOCOL_PATH)
    required = set(
        protocol["failure_and_deviation_reporting"]["deviation_record"]["required_fields"]
    )

    assert ledger["deviations"], "ledger must not be empty once a frozen section has changed"
    for record in ledger["deviations"]:
        missing = required - set(record)
        assert not missing, f"{record.get('deviation_id')} omits {sorted(missing)}"


def test_the_latest_ledger_record_matches_the_current_matrix_digest(protocol: dict) -> None:
    """The ledger is append-only, so the newest record for a section must be current."""
    ledger = load_yaml(ROOT / "analysis" / "v3_validation" / "deviations.yaml")
    matrix_records = [r for r in ledger["deviations"] if "workflow_matrix" in r["protocol_section"]]

    assert matrix_records
    latest = matrix_records[-1]
    assert latest["protocol_sha256_after"] == protocol["workflow_matrix"]["contract_sha256"]
    assert latest["protocol_sha256_after"] != latest["protocol_sha256_before"]
    assert all(r["outcome_triggered"] is False for r in ledger["deviations"])


def test_the_virtus_addition_is_recorded_and_not_outcome_triggered() -> None:
    ledger = load_yaml(ROOT / "analysis" / "v3_validation" / "deviations.yaml")
    record = next(r for r in ledger["deviations"] if r["deviation_id"] == "DEV-001")

    assert record["protocol_section"] == "workflow_matrix"
    assert record["outcome_triggered"] is False
    assert any("virtus" in row for row in record["affected_workflow_row_ids"])


def test_every_recorded_digest_change_is_reproducible(protocol: dict) -> None:
    """A ledger claiming an 'after' digest nobody can recompute is worthless."""
    ledger = load_yaml(ROOT / "analysis" / "v3_validation" / "deviations.yaml")
    live = {
        "harmonization.contract_sha256": harmonization_sha256(protocol["harmonization"]),
        "frozen_inputs.contract_sha256": frozen_inputs_sha256(protocol),
        "partitions.contract_sha256": partitions_sha256(protocol["partitions"]),
        "calibration.contract_sha256": calibration_sha256(protocol["calibration"]),
        "workflow_matrix.contract_sha256": workflow_matrix_sha256(protocol["workflow_matrix"]),
        "failure_and_deviation_reporting.contract_sha256": failure_reporting_sha256(
            protocol["failure_and_deviation_reporting"]
        ),
    }
    latest = ledger["deviations"][-1]
    for scope, change in (latest.get("additional_digest_changes") or {}).items():
        assert change["after"] == live[scope], scope
    assert latest["protocol_sha256_after"] == live[latest["digest_scope"]]


def test_comparators_span_more_than_one_architecture(protocol: dict) -> None:
    """A single shared architecture cannot distinguish tool-specific from general results."""
    modes = {
        w["mode"] for w in protocol["workflow_matrix"]["workflows"] if w["role"] == "comparator"
    }
    assert len(modes) >= 3
    assert "comparator_architecture_note" in protocol["workflow_matrix"]


@pytest.fixture
def ledger() -> dict:
    return load_yaml(DEFAULT_LEDGER)


def test_the_committed_protocol_and_ledger_agree(protocol: dict, ledger: dict) -> None:
    assert validate_amendment_ledger(protocol, ledger) == []


def test_a_frozen_section_without_a_baseline_record_is_rejected(protocol: dict) -> None:
    """The genesis case: harmonization froze before the ledger existed."""
    errors = validate_amendment_ledger(protocol, {"deviations": []})

    assert any("no deviation-ledger baseline record" in error for error in errors)


def test_a_silently_rehashed_frozen_section_is_detected(protocol: dict, ledger: dict) -> None:
    """F7: editing a frozen section and recomputing its own digest used to pass."""
    tampered = deepcopy(protocol)
    tampered["harmonization"]["decision_timing"] = "after-outcomes-are-inspected"
    tampered["harmonization"]["contract_sha256"] = harmonization_sha256(tampered["harmonization"])

    errors = validate_amendment_ledger(tampered, ledger)

    assert any("undocumented amendment" in error for error in errors)


def test_a_broken_ledger_chain_is_detected(protocol: dict, ledger: dict) -> None:
    tampered_ledger = deepcopy(ledger)
    matrix_records = [
        r for r in tampered_ledger["deviations"] if "workflow_matrix" in r["protocol_section"]
    ]
    assert len(matrix_records) >= 2
    matrix_records[-1]["protocol_sha256_before"] = "f" * 64
    frozen = deepcopy(protocol)
    frozen["workflow_matrix"]["status"] = "frozen"
    frozen["workflow_matrix"]["frozen_at"] = "2026-07-27"
    frozen["workflow_matrix"]["contract_sha256"] = workflow_matrix_sha256(frozen["workflow_matrix"])

    errors = validate_amendment_ledger(frozen, tampered_ledger)

    assert any("is not continuous" in error for error in errors)


def test_a_pending_section_is_not_bound_by_the_amendment_rule(protocol: dict) -> None:
    """Only sections currently claiming frozen must be accounted for."""
    assert protocol["partitions"]["status"] == "pending"

    errors = validate_amendment_ledger(protocol, {"deviations": []})

    assert not any("'partitions'" in error for error in errors)


def test_validate_protocol_file_fails_closed_without_a_ledger(tmp_path: Path) -> None:
    errors = validate_protocol_file(
        PROTOCOL_PATH, SCHEMA_PATH, ledger_path=tmp_path / "absent.yaml"
    )

    assert len(errors) == 1
    assert "deviation ledger not found" in errors[0]


def test_validate_protocol_file_checks_the_real_ledger() -> None:
    """End-to-end: the committed protocol against the committed ledger."""
    assert validate_protocol_file(PROTOCOL_PATH, SCHEMA_PATH) == []


def test_execution_gate_blocks_training_today() -> None:
    """Blockers are open by design, so the gate must refuse to let rows execute."""
    with pytest.raises(SystemExit) as excinfo:
        require_execution_allowed("training")

    assert "execution gate failed" in str(excinfo.value)


def test_execution_gate_permits_the_draft_phase() -> None:
    """The draft gate is what CI checks; it must pass on the committed protocol."""
    require_execution_allowed("draft")


def test_the_execution_gate_obligation_is_declared_in_the_protocol(protocol: dict) -> None:
    """F8: the guard exists, so the protocol must say where it has to be called."""
    gate = protocol["workflow_matrix"]["environment_pinning"]["execution_gate"]

    assert gate["guard"] == "scripts.validate_v3_protocol.require_execution_allowed"
    assert set(gate["required_call_sites"]) == {
        "manifest_build_time",
        "per_row_pre_execution",
    }
    assert gate["status"] == "unimplemented"

    blocker_ids = {b["id"] for b in protocol["execution_readiness"]["training_blockers"]}
    assert "execution_gate_unimplemented" in blocker_ids


def test_ci_runs_the_draft_gate_and_not_the_training_gate() -> None:
    """A --phase training job in CI would be permanently red by design."""
    ci = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    # Only executed lines count; the comment above the step names the phases
    # precisely so a future editor does not re-add them.
    commands = [
        line.strip()
        for line in ci.splitlines()
        if "validate_v3_protocol.py" in line and not line.strip().startswith("#")
    ]

    assert commands, "CI must run the protocol validator"
    assert any(command.endswith("validate_v3_protocol.py") for command in commands)
    assert not any("--phase" in command for command in commands)


def test_frozen_partitions_requires_factors_with_non_empty_levels(
    protocol: dict, schema: dict
) -> None:
    """F1: round 1 froze partitions while four of five factors had empty levels."""
    tampered = deepcopy(protocol)
    tampered["partitions"]["status"] = "frozen"
    tampered["partitions"]["frozen_at"] = "2026-07-27"
    tampered["partitions"]["contract_sha256"] = partitions_sha256(tampered["partitions"])
    tampered["planned_freeze_sections"]["partitions"]["status"] = "frozen"

    errors = validate_protocol(tampered, schema)

    assert any("to be frozen with non-empty levels" in error for error in errors), (
        "a declared-but-empty factor must not satisfy a partitions freeze"
    )


def test_frozen_calibration_requires_enough_abundance_levels(protocol: dict, schema: dict) -> None:
    """F2: a probit LOD is not identifiable with too few abundance levels."""
    tampered = deepcopy(protocol)
    tampered["calibration"]["status"] = "frozen"
    tampered["calibration"]["frozen_at"] = "2026-07-27"
    tampered["calibration"]["contract_sha256"] = calibration_sha256(tampered["calibration"])
    tampered["planned_freeze_sections"]["calibration_metrics_lod"]["status"] = "frozen"
    for factor in tampered["factors"]:
        if factor["id"] == "viral_abundance":
            factor["status"] = "frozen"
            factor["levels"] = [1, 10, 100]  # three: below the threshold

    errors = validate_protocol(tampered, schema)

    assert any("at least 5 levels" in error for error in errors)


def test_enough_abundance_levels_satisfies_the_calibration_precondition(
    protocol: dict, schema: dict
) -> None:
    """The rule must be satisfiable, not merely blocking."""
    tampered = deepcopy(protocol)
    tampered["calibration"]["status"] = "frozen"
    tampered["calibration"]["frozen_at"] = "2026-07-27"
    tampered["calibration"]["contract_sha256"] = calibration_sha256(tampered["calibration"])
    tampered["planned_freeze_sections"]["calibration_metrics_lod"]["status"] = "frozen"
    for factor in tampered["factors"]:
        if factor["id"] == "viral_abundance":
            factor["status"] = "frozen"
            factor["levels"] = [1, 10, 100, 1000, 10000]

    errors = validate_protocol(tampered, schema)

    assert not any("viral_abundance" in error for error in errors)


def test_the_dormant_preconditions_do_not_fire_today(protocol: dict, schema: dict) -> None:
    """Both sections are pending, so neither new rule may affect current validation."""
    assert protocol["partitions"]["status"] == "pending"
    assert protocol["calibration"]["status"] == "pending"

    assert validate_protocol(protocol, schema) == []


def test_partitions_apportionment_is_not_ceiling_biased(protocol: dict) -> None:
    """F11: per-stratum ceiling rounding gives 50% holdout at n=2 against a 0.3 target."""
    partitions = protocol["partitions"]

    assert "largest-remainder" in partitions["algorithm"]
    assert "rounding up" not in partitions["algorithm"]
    assert "apportionment_rationale" in partitions


def test_partitions_sort_key_is_pinned(protocol: dict) -> None:
    """F12: sample_2 and sample_10 order differently under ASCII and natural sort."""
    sort_key = protocol["partitions"]["sort_key"].lower()

    assert "byte-wise ascii" in sort_key
    assert "case sensitive" in sort_key
    assert "no locale" in sort_key


def test_partitions_split_scope_is_explicit(protocol: dict) -> None:
    """F22: single-sample datasets cannot be split and must be named exempt."""
    partitions = protocol["partitions"]

    assert partitions["applies_to_dataset_roles"] == ["training-and-holdout"]
    assert "exempt_dataset_roles_rationale" in partitions
    assert "unsplittable_stratum_reporting" in partitions


def test_reference_construction_is_a_declared_leakage_path(protocol: dict) -> None:
    """F17: deciding which genes are quantifiable can carry holdout information."""
    prohibitions = " ".join(protocol["partitions"]["leakage_prohibitions"]).lower()

    assert "d-list" in prohibitions
    assert "feature universe" in prohibitions
    assert "training-eligible" in prohibitions


def test_an_empty_eligible_grid_has_a_declared_outcome(protocol: dict) -> None:
    """F3: the procedure had an undefined branch on a plausible input."""
    policy = protocol["calibration"]["threshold_search"]["no_eligible_grid_point_policy"].lower()

    assert "calibration has failed" in policy
    assert "never be relaxed" in policy
    assert "never be widened" in policy


def test_the_tie_breaker_names_its_standard_error(protocol: dict) -> None:
    """F13: three plausible estimators give three different frozen thresholds."""
    text = protocol["calibration"]["tie_breaker_standard_error"].lower()

    assert "bootstrap" in text
    assert "not a binomial" in text
    assert "not a cross-validation" in text


def test_the_bootstrap_refuses_degenerate_small_samples(protocol: dict) -> None:
    """F4: BCa on n=1 returns a zero-width interval that reads as certainty."""
    uncertainty = protocol["calibration"]["uncertainty"]

    assert uncertainty["minimum_samples_for_interval"] >= 3
    policy = uncertainty["below_minimum_policy"].lower()
    assert "not-estimable" in policy
    assert "zero-width" in policy
    assert "stratification_scope" in uncertainty


def test_failed_rows_cannot_flatter_the_negative_false_call_rate(protocol: dict) -> None:
    """F9: crashed rows were being counted as clean negative conditions."""
    metrics = {m["id"]: m for m in protocol["calibration"]["metrics"]}
    rule = metrics["M2_negative_false_call_rate"]["denominator_rule"].lower()

    assert "complete, comparable" in rule
    assert "technical-noncompletion" in rule


def test_two_step_loss_denominates_against_truth(protocol: dict) -> None:
    """F10: a measured denominator shrinks with the numerator and hides loss."""
    metrics = {m["id"]: m for m in protocol["calibration"]["metrics"]}
    rule = metrics["M7_two_step_loss"]["denominator_rule"].lower()

    assert "truth-eligible" in rule


def test_cross_tool_parity_is_symmetric_and_defined_at_zero(protocol: dict) -> None:
    """F20: denominating on the comparator alone is undefined at zero counts."""
    metrics = {m["id"]: m for m in protocol["calibration"]["metrics"]}
    parity = metrics["M9_cross_tool_parity"]

    assert "symmetric" in parity["definition"].lower()
    assert "not-applicable" in parity["denominator_rule"].lower()


def test_dedicated_comparators_run_both_reference_arms(protocol: dict) -> None:
    """F14: one arm cannot separate 'tool is worse' from 'reference differed'."""
    arms: dict[str, set[str]] = {}
    for workflow in protocol["workflow_matrix"]["workflows"]:
        if workflow["tool"] in {"venus", "viral-track", "virtus"}:
            arms.setdefault(workflow["tool"], set()).add(workflow["reference_resolution"])

    assert set(arms) == {"venus", "viral-track", "virtus"}
    for tool, resolutions in arms.items():
        assert resolutions == {"native-published", "matched-accession-index"}, tool


def test_native_arms_use_a_native_reference_not_the_curated_index(protocol: dict) -> None:
    """A native arm pinned to the curated index would not be a native arm."""
    references = {r["id"]: r for r in protocol["references"]}
    for workflow in protocol["workflow_matrix"]["workflows"]:
        if workflow.get("reference_resolution") == "native-published":
            assert workflow["reference_id"] != "curated_human_virus", workflow["id"]
            assert references[workflow["reference_id"]]["profile"] == "native-published"
        elif workflow.get("reference_resolution") == "matched-accession-index":
            assert workflow["reference_id"] == "curated_human_virus", workflow["id"]


def test_the_matrix_declares_its_fairness_disclosures(protocol: dict) -> None:
    """F6 and F15: both asymmetries favour ViralScan and must be disclosed."""
    matrix = protocol["workflow_matrix"]

    breadth = matrix["breadth_disclosure_rule"].lower()
    assert "no comparator" in breadth
    assert "do not yet exist" in breadth

    tuning = matrix["comparator_tuning_asymmetry"].lower()
    assert "published defaults" in tuning
    assert "counts_unique" in tuning
    assert any("untuned" in risk.lower() for risk in protocol["principal_risks"])


def test_the_tool_environments_blocker_names_every_comparator(protocol: dict, schema: dict) -> None:
    """F24: the text already went stale once when VIRTUS was added."""
    assert validate_protocol(protocol, schema) == []

    tampered = deepcopy(protocol)
    for blocker in tampered["execution_readiness"]["training_blockers"]:
        if blocker["id"] == "tool_environments":
            blocker["description"] = "starsolo and venus versions are pending."

    errors = validate_protocol(tampered, schema)

    assert any("omits comparator" in error for error in errors)


def test_digest_scope_is_a_required_deviation_field(protocol: dict) -> None:
    required = protocol["failure_and_deviation_reporting"]["deviation_record"]["required_fields"]

    assert "digest_scope" in required


def test_rewriting_a_frozen_seed_is_detected(protocol: dict, schema: dict) -> None:
    """R2-F1, reproduced by the round-2 reviewer: this used to produce no errors."""
    tampered = deepcopy(protocol)
    assert tampered["seeds"]["split"]["status"] == "frozen"
    tampered["seeds"]["split"]["value"] = 999999999

    errors = validate_protocol(tampered, schema)

    assert any("frozen_inputs contract_sha256" in error for error in errors)


def test_widening_a_frozen_factor_is_detected(protocol: dict, schema: dict) -> None:
    """R2-F1: frozen factor levels determine the strata and must be tamper-evident."""
    tampered = deepcopy(protocol)
    widened = False
    for factor in tampered["factors"]:
        if factor.get("status") == "frozen":
            factor["levels"] = list(factor.get("levels") or []) + ["made-up-level"]
            widened = True
    assert widened, "the fixture must contain at least one frozen factor"

    errors = validate_protocol(tampered, schema)

    assert any("frozen_inputs contract_sha256" in error for error in errors)


def test_frozen_inputs_covers_every_frozen_seed_and_factor(protocol: dict) -> None:
    frozen_seeds = {
        name
        for name, seed in protocol["seeds"].items()
        if isinstance(seed, dict) and seed.get("status") == "frozen"
    }
    assert frozen_seeds, "seeds are frozen today, so the digest must cover them"
    assert protocol["frozen_inputs"]["status"] == "frozen"
    assert protocol["frozen_inputs"]["contract_sha256"] == frozen_inputs_sha256(protocol)
    assert set(protocol["frozen_inputs"]["covers"]) == {"seeds", "factors"}


def test_pending_seeds_may_still_change_freely(protocol: dict, schema: dict) -> None:
    """The rule binds frozen entries only; VAL-01 must still be able to land."""
    tampered = deepcopy(protocol)
    assert tampered["seeds"]["generation"]["status"] == "pending"
    tampered["seeds"]["generation"]["pending_reason"] = "reworded during VAL-01 design"

    errors = validate_protocol(tampered, schema)

    assert not any("frozen_inputs" in error for error in errors)


def test_frozen_inputs_states_what_it_does_not_cover(protocol: dict) -> None:
    """An integrity rail that oversells its scope is worse than a narrow one."""
    limitation = protocol["frozen_inputs"]["limitation"].lower()

    assert "pending" in limitation
    assert "verified" in limitation


def test_the_primary_accuracy_claim_names_its_arm(protocol: dict) -> None:
    """R2-F2: rule 2 demanded identical sequences, which the native arm violates."""
    matrix = protocol["workflow_matrix"]
    rule_two = matrix["primary_comparison_rules"][1].lower()

    assert "matched-accession-index" in rule_two
    assert "native-published arms deliberately violate it" in rule_two

    scope = matrix["arm_claim_scope"].lower()
    assert "matched-accession-index arm carries every primary" in scope
    assert "reported separately" in scope


def test_every_arm_is_paired_with_its_counterpart(protocol: dict) -> None:
    """R2-F3: an unpaired arm cannot separate a tool difference from a reference one."""
    workflows = {w["id"]: w for w in protocol["workflow_matrix"]["workflows"]}
    natives = {
        i for i, w in workflows.items() if w.get("reference_resolution") == "native-published"
    }
    matched = {
        i
        for i, w in workflows.items()
        if w.get("reference_resolution") == "matched-accession-index"
    }

    assert natives and len(natives) == len(matched)
    for matched_id in matched:
        counterpart = workflows[matched_id]["native_counterpart"]
        assert counterpart in natives
        assert workflows[counterpart]["tool"] == workflows[matched_id]["tool"]


def test_an_unpaired_native_arm_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    for workflow in tampered["workflow_matrix"]["workflows"]:
        if workflow.get("reference_resolution") == "matched-accession-index":
            workflow["native_counterpart"] = "W7a_venus_native"
    tampered["workflow_matrix"]["contract_sha256"] = workflow_matrix_sha256(
        tampered["workflow_matrix"]
    )

    errors = validate_protocol(tampered, schema)

    assert any("has no matched counterpart" in error for error in errors)
    assert any("is claimed by both" in error for error in errors)


def test_a_cross_tool_arm_pairing_is_rejected(protocol: dict, schema: dict) -> None:
    tampered = deepcopy(protocol)
    for workflow in tampered["workflow_matrix"]["workflows"]:
        if workflow["id"] == "W7b_venus_matched":
            workflow["native_counterpart"] = "W9a_virtus_native"
    tampered["workflow_matrix"]["contract_sha256"] = workflow_matrix_sha256(
        tampered["workflow_matrix"]
    )

    errors = validate_protocol(tampered, schema)

    assert any("of a different tool" in error for error in errors)


def test_the_accession_linkage_admits_it_is_not_yet_checkable(protocol: dict) -> None:
    """Declaring a linkage that cannot be verified yet must say so, not imply it holds."""
    rule = protocol["workflow_matrix"]["accession_linkage_rule"].lower()

    assert "null" in rule
    assert "not yet checkable" in rule
    assert "rel-03" in rule


def test_the_predecessor_is_a_scored_comparator(protocol: dict) -> None:
    """A benchmark justifying a rewrite must include the thing being rewritten."""
    legacy = [w for w in protocol["workflow_matrix"]["workflows"] if w["tool"] == "viralscan-2.2.0"]

    assert len(legacy) == 2, "both reference arms, as for every dedicated comparator"
    for arm in legacy:
        assert "synthetic_factorial" in arm["dataset_ids"], arm["id"]
        assert "synthetic_host_only" in arm["dataset_ids"], arm["id"]
        assert "synthetic_host_homology" in arm["dataset_ids"], arm["id"]


def test_an_improvement_claim_may_not_rest_on_counts(protocol: dict) -> None:
    """The two versions report different units, so magnitude carries no verdict."""
    rule = protocol["workflow_matrix"]["predecessor_comparison_rule"].lower()

    assert "different units" in rule
    assert "definitional" in rule
    assert "precision, recall" in rule
    assert "general superiority" in rule
