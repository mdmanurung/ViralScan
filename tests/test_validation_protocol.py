from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from scripts.validate_v3_protocol import (
    calibration_sha256,
    failure_reporting_sha256,
    harmonization_dependencies_sha256,
    harmonization_sha256,
    load_yaml,
    partitions_sha256,
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
    errors = validate_protocol(broken, schema)
    assert any(
        "schema harmonization.contract_sha256" in error
        and "699854b71169222e74d26c2119f1c9d7742d5962ac4ed8c02dfa9f288e3339dc" in error
        for error in errors
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
