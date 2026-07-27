#!/usr/bin/env python3
"""Validate the repo's ViralScan v3 scientific protocol."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator, FormatChecker
from jsonschema.exceptions import SchemaError

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = REPO_ROOT / "analysis" / "v3_validation" / "protocol.yaml"
DEFAULT_SCHEMA = REPO_ROOT / "schemas" / "v3" / "validation_protocol.schema.json"
DEFAULT_LEDGER = REPO_ROOT / "analysis" / "v3_validation" / "deviations.yaml"

# Minimum distinct viral-abundance levels per chemistry. SCI-05 round 1 (F2)
# gave four to five as the identifiability floor; round 5 (R5-F4) raised it to
# seven so the grid can bracket the detection knee with levels to spare, since a
# floor-value grid fails if a single level turns out uninformative.
MINIMUM_ABUNDANCE_LEVELS = 7


class UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that rejects duplicate mapping keys."""


def _construct_unique_mapping(
    loader: UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False
) -> dict[Any, Any]:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            if key in mapping:
                raise ValueError(f"duplicate YAML key: {key!r}")
            mapping[key] = loader.construct_object(value_node, deep=deep)
        except TypeError as exc:
            raise ValueError(f"unhashable YAML mapping key: {key!r}") from exc
    return mapping


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML mapping while rejecting duplicate keys."""
    with path.open(encoding="utf-8") as handle:
        document = yaml.load(handle, Loader=UniqueKeyLoader)
    if not isinstance(document, dict):
        raise ValueError("protocol must be a YAML mapping")
    return document


def _duplicate_ids(items: Any, section: str) -> list[str]:
    if not isinstance(items, list):
        return []
    ids = [item.get("id") for item in items if isinstance(item, dict)]
    return [
        f"{section} contains duplicate id {item_id!r}"
        for item_id in sorted(set(ids))
        if ids.count(item_id) > 1
    ]


def _canonical_sha256(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def harmonization_sha256(harmonization: dict[str, Any]) -> str:
    """Return the canonical SCI-02 contract digest, excluding its digest field."""
    payload = {key: value for key, value in harmonization.items() if key != "contract_sha256"}
    return _canonical_sha256(payload)


def partitions_sha256(partitions: dict[str, Any]) -> str:
    """Return the canonical SCI-03 partition-contract digest, excluding its digest field."""
    payload = {key: value for key, value in partitions.items() if key != "contract_sha256"}
    return _canonical_sha256(payload)


def calibration_sha256(calibration: dict[str, Any]) -> str:
    """Return the canonical SCI-03 calibration-contract digest, excluding its digest field."""
    payload = {key: value for key, value in calibration.items() if key != "contract_sha256"}
    return _canonical_sha256(payload)


def _validate_sci03(document: dict[str, Any]) -> list[str]:
    """Check SCI-03 decisions that JSON Schema cannot express across sections."""
    errors: list[str] = []
    partitions = document.get("partitions")
    calibration = document.get("calibration")
    if not isinstance(partitions, dict) or not isinstance(calibration, dict):
        return ["partitions and calibration must be objects"]

    planned = document.get("planned_freeze_sections", {})
    pairs = (
        ("partitions", partitions, "partitions"),
        ("calibration_metrics_lod", calibration, "calibration"),
    )
    for planned_name, section, section_name in pairs:
        section_frozen = section.get("status") == "frozen"
        planned_section = planned.get(planned_name, {}) if isinstance(planned, dict) else {}
        planned_status = (
            planned_section.get("status") if isinstance(planned_section, dict) else None
        )
        if section_frozen and planned_status != "frozen":
            errors.append(
                f"frozen {section_name} requires planned section {planned_name!r} to be frozen"
            )
        if not section_frozen and planned_status == "frozen":
            errors.append(
                f"planned section {planned_name!r} cannot be frozen before {section_name}"
            )
        if section_frozen and not section.get("frozen_at"):
            errors.append(f"frozen {section_name} requires frozen_at")

    # R2-F1: the ledger check only runs when a ledger is supplied, but a rewritten
    # frozen seed must fail every path, including in-memory validation.
    frozen_inputs = document.get("frozen_inputs")
    if isinstance(frozen_inputs, dict) and frozen_inputs.get("status") == "frozen":
        if frozen_inputs.get("contract_sha256") != frozen_inputs_sha256(document):
            errors.append(
                "frozen_inputs contract_sha256 does not match the frozen seed and "
                "factor values it covers"
            )

    if partitions.get("contract_sha256") != partitions_sha256(partitions):
        errors.append("partitions contract_sha256 does not match the canonical SCI-03 contract")
    if calibration.get("contract_sha256") != calibration_sha256(calibration):
        errors.append("calibration contract_sha256 does not match the canonical SCI-03 contract")

    # A frozen split is unreproducible without its seed, and an unfrozen seed
    # after freeze would let the partition be reselected post hoc.
    seeds = document.get("seeds", {})
    if partitions.get("status") == "frozen" and isinstance(seeds, dict):
        for seed_name in ("root", "split", "cell_calling", "evidence_sampling", "bootstrap"):
            seed = seeds.get(seed_name, {})
            if not isinstance(seed, dict) or seed.get("status") != "frozen":
                errors.append(f"frozen partitions require seeds.{seed_name} to be frozen")

    # Every declared metric must resolve to a real endpoint.
    endpoint_ids = {
        item.get("id") for item in document.get("endpoints", []) if isinstance(item, dict)
    }
    for metric in calibration.get("metrics", []):
        if not isinstance(metric, dict):
            continue
        endpoint_id = metric.get("endpoint_id")
        if endpoint_id not in endpoint_ids:
            errors.append(
                f"calibration metric {metric.get('id')!r} references unknown endpoint "
                f"{endpoint_id!r}"
            )
    lod_endpoint = calibration.get("limit_of_detection", {})
    if isinstance(lod_endpoint, dict):
        endpoint_id = lod_endpoint.get("endpoint_id")
        if endpoint_id not in endpoint_ids:
            errors.append(f"limit_of_detection references unknown endpoint {endpoint_id!r}")

    # Stratification factors must be declared factors.
    factors_by_id = {
        item.get("id"): item for item in document.get("factors", []) if isinstance(item, dict)
    }
    factor_ids = set(factors_by_id)
    for factor in partitions.get("stratification_factors", []):
        if factor not in factor_ids:
            errors.append(f"partitions stratification factor {factor!r} is not a declared factor")

    # A declared factor is not a usable one. SCI-05 round 1 (F1) found partitions
    # frozen while four of five stratification factors had empty level lists, so
    # the stratum cross-product could not be computed at all. Being declared was
    # the only thing previously checked.
    if partitions.get("status") == "frozen":
        for factor_id in partitions.get("stratification_factors", []):
            factor = factors_by_id.get(factor_id)
            if not isinstance(factor, dict):
                continue
            if factor.get("status") != "frozen" or not factor.get("levels"):
                errors.append(
                    f"frozen partitions requires stratification factor {factor_id!r} "
                    "to be frozen with non-empty levels"
                )

    # F2: a probit limit-of-detection curve is not identifiable without enough
    # abundance levels spanning sub-detection to saturating detection. Five is the
    # conservative end of the reviewer's range, leaving room for one level to land
    # uninformatively without collapsing the fit.
    if calibration.get("status") == "frozen":
        abundance = factors_by_id.get("viral_abundance")
        levels = abundance.get("levels") if isinstance(abundance, dict) else None
        if not isinstance(abundance, dict) or abundance.get("status") != "frozen":
            errors.append(
                "frozen calibration requires factors.viral_abundance to be frozen "
                "for probit limit-of-detection identifiability"
            )
        elif len(levels or []) < MINIMUM_ABUNDANCE_LEVELS:
            errors.append(
                f"frozen calibration requires factors.viral_abundance to declare at least "
                f"{MINIMUM_ABUNDANCE_LEVELS} levels for probit limit-of-detection "
                f"identifiability, found {len(levels or [])}"
            )
    for factor in calibration.get("uncertainty", {}).get("stratified_by", []):
        if factor not in factor_ids:
            errors.append(f"uncertainty stratification factor {factor!r} is not a declared factor")

    return errors


def workflow_matrix_sha256(workflow_matrix: dict[str, Any]) -> str:
    """Return the canonical SCI-04 workflow-matrix digest, excluding its digest field."""
    payload = {key: value for key, value in workflow_matrix.items() if key != "contract_sha256"}
    return _canonical_sha256(payload)


def failure_reporting_sha256(failure_reporting: dict[str, Any]) -> str:
    """Return the canonical SCI-04 failure-reporting digest, excluding its digest field."""
    payload = {key: value for key, value in failure_reporting.items() if key != "contract_sha256"}
    return _canonical_sha256(payload)


def _validate_sci04(document: dict[str, Any], phase: str) -> list[str]:
    """Check SCI-04 decisions that JSON Schema cannot express across sections."""
    errors: list[str] = []
    matrix = document.get("workflow_matrix")
    reporting = document.get("failure_and_deviation_reporting")
    if not isinstance(matrix, dict) or not isinstance(reporting, dict):
        return ["workflow_matrix and failure_and_deviation_reporting must be objects"]

    planned = document.get("planned_freeze_sections", {})
    pairs = (
        ("workflow_matrix", matrix, "workflow_matrix"),
        ("failure_and_deviation_reporting", reporting, "failure_and_deviation_reporting"),
    )
    for planned_name, section, section_name in pairs:
        section_frozen = section.get("status") == "frozen"
        planned_section = planned.get(planned_name, {}) if isinstance(planned, dict) else {}
        planned_status = (
            planned_section.get("status") if isinstance(planned_section, dict) else None
        )
        if section_frozen and planned_status != "frozen":
            errors.append(
                f"frozen {section_name} requires planned section {planned_name!r} to be frozen"
            )
        if not section_frozen and planned_status == "frozen":
            errors.append(
                f"planned section {planned_name!r} cannot be frozen before {section_name}"
            )
        if section_frozen and not section.get("frozen_at"):
            errors.append(f"frozen {section_name} requires frozen_at")

    if matrix.get("contract_sha256") != workflow_matrix_sha256(matrix):
        errors.append(
            "workflow_matrix contract_sha256 does not match the canonical SCI-04 contract"
        )
    if reporting.get("contract_sha256") != failure_reporting_sha256(reporting):
        errors.append(
            "failure_and_deviation_reporting contract_sha256 does not match "
            "the canonical SCI-04 contract"
        )

    workflows = matrix.get("workflows", [])
    errors.extend(_duplicate_ids(workflows, "workflow_matrix.workflows"))

    dataset_ids = {
        item.get("id") for item in document.get("datasets", []) if isinstance(item, dict)
    }
    reference_ids = {
        item.get("id") for item in document.get("references", []) if isinstance(item, dict)
    }
    declared_rows = 0
    for workflow in workflows:
        if not isinstance(workflow, dict):
            continue
        workflow_id = workflow.get("id")
        if workflow.get("reference_id") not in reference_ids:
            errors.append(
                f"workflow {workflow_id!r} references unknown reference "
                f"{workflow.get('reference_id')!r}"
            )
        rows = workflow.get("dataset_ids", [])
        declared_rows += len(rows)
        for dataset_id in rows:
            if dataset_id not in dataset_ids:
                errors.append(f"workflow {workflow_id!r} references unknown dataset {dataset_id!r}")
        if len(set(rows)) != len(rows):
            errors.append(f"workflow {workflow_id!r} repeats a dataset id")

    # A drifting row count is how a silently dropped comparison hides.
    if matrix.get("expected_row_count") != declared_rows:
        errors.append(
            f"workflow_matrix expected_row_count {matrix.get('expected_row_count')!r} "
            f"does not match the {declared_rows} enumerated rows"
        )

    # R2-F2 and R2-F3: a two-arm comparator is only meaningful if the arms are
    # actually paired. An unpaired native arm has nothing to be compared against,
    # and a matched arm whose counterpart does not resolve cannot attribute a
    # difference to the reference rather than to the tool.
    by_id = {w.get("id"): w for w in workflows if isinstance(w, dict)}
    native_arms = {
        w.get("id"): w
        for w in workflows
        if isinstance(w, dict) and w.get("reference_resolution") == "native-published"
    }
    claimed_counterparts: dict[str, str] = {}
    for workflow in workflows:
        if not isinstance(workflow, dict):
            continue
        if workflow.get("reference_resolution") != "matched-accession-index":
            continue
        workflow_id = workflow.get("id")
        counterpart = workflow.get("native_counterpart")
        if not counterpart:
            errors.append(f"matched arm {workflow_id!r} declares no native_counterpart")
            continue
        if counterpart not in native_arms:
            errors.append(
                f"matched arm {workflow_id!r} names {counterpart!r}, which is not a "
                "native-published arm"
            )
            continue
        if by_id[counterpart].get("tool") != workflow.get("tool"):
            errors.append(
                f"matched arm {workflow_id!r} is paired with {counterpart!r} of a different tool"
            )
        if counterpart in claimed_counterparts:
            errors.append(
                f"native arm {counterpart!r} is claimed by both "
                f"{claimed_counterparts[counterpart]!r} and {workflow_id!r}"
            )
        claimed_counterparts[counterpart] = workflow_id
    for native_id in native_arms:
        if native_id not in claimed_counterparts:
            errors.append(f"native arm {native_id!r} has no matched counterpart")

    # The blocker naming the unpinned tools must name all of them. It has already
    # gone stale once: VIRTUS was added without updating the text. A prose
    # instruction to keep it current is not a rail, so tie it to the matrix.
    comparator_tools = {
        workflow.get("tool")
        for workflow in workflows
        if isinstance(workflow, dict) and workflow.get("role") == "comparator"
    }
    readiness = document.get("execution_readiness", {})
    blockers = readiness.get("training_blockers", []) if isinstance(readiness, dict) else []
    for blocker in blockers:
        if not isinstance(blocker, dict) or blocker.get("id") != "tool_environments":
            continue
        description = (blocker.get("description") or "").lower()
        missing = sorted(
            tool for tool in comparator_tools if tool and tool.lower() not in description
        )
        if missing:
            errors.append(
                "tool_environments blocker description omits comparator "
                f"{'tools' if len(missing) > 1 else 'tool'} {', '.join(repr(t) for t in missing)}"
            )

    # The matrix may freeze which rows exist before REL-03 supplies versions and
    # digests, but no row may execute against an unpinned environment.
    if phase in {"training", "holdout"}:
        pinning = matrix.get("environment_pinning", {})
        if isinstance(pinning, dict) and pinning.get("status") != "frozen":
            errors.append("workflow environment pinning is not frozen")
        for workflow in workflows:
            if not isinstance(workflow, dict):
                continue
            if not workflow.get("tool_version") or not workflow.get("container_digest"):
                errors.append(
                    f"workflow {workflow.get('id')!r} has no pinned tool version or container digest"
                )

    return errors


def harmonization_dependencies_sha256(document: dict[str, Any]) -> str:
    """Hash protocol prose whose meaning must not contradict SCI-02."""
    datasets = []
    for dataset in document.get("datasets", []):
        if not isinstance(dataset, dict):
            continue
        datasets.append(
            {
                key: dataset.get(key)
                for key in (
                    "id",
                    "class",
                    "truth_status",
                    "expected_targets",
                    "required",
                    "exclusion_rule",
                )
            }
        )
    payload = {
        "supported_scope": document.get("supported_scope"),
        "hypotheses": document.get("hypotheses"),
        "endpoints": document.get("endpoints"),
        "datasets": datasets,
        "exclusions": document.get("exclusions"),
    }
    return _canonical_sha256(payload)


def _validate_harmonization(document: dict[str, Any]) -> list[str]:
    """Check SCI-02 decisions that JSON Schema cannot express across sections."""
    errors: list[str] = []
    harmonization = document.get("harmonization")
    if not isinstance(harmonization, dict):
        return errors

    denominator_ids = [
        item.get("id") for item in harmonization.get("denominators", []) if isinstance(item, dict)
    ]
    expected_denominator_ids = {
        "D1_count_invariant_runs",
        "D2_unique_molecule_precision",
        "D3_unique_molecule_recall",
        "D4_infected_cell_precision",
        "D5_infected_cell_recall",
        "D6_infected_cell_auprc",
        "D7_exact_negative_cell_false_calls",
        "D8_exact_negative_sample_false_calls",
        "D9_exact_negative_molecule_false_calls",
        "D10_presumed_negative_unexpected_calls",
        "D11_empty_droplet_false_calls",
        "D12_positive_cell_rate",
        "D13_sibling_molecule_confusion",
        "D14_host_virus_allocation",
        "D15_two_step_truth_read_loss",
        "D16_two_step_truth_molecule_loss",
        "D17_lod_detection_probability",
        "D18_tier_specificity",
        "D19_tier_ppv",
        "D20_cross_tool_absolute_unique",
        "D21_secondary_burden_per_10000",
        "D22_workflow_technical_failure_rate",
        "D23_endpoint_incomparability_rate",
        "D24_sibling_cell_confusion",
    }
    if len(denominator_ids) != len(set(denominator_ids)):
        errors.append("harmonization denominators contain duplicate ids")
    if set(denominator_ids) != expected_denominator_ids:
        missing = sorted(expected_denominator_ids - set(denominator_ids))
        extra = sorted(set(denominator_ids) - expected_denominator_ids)
        errors.append(
            f"harmonization denominator registry mismatch; missing={missing}, extra={extra}"
        )

    artifact_paths = [
        item.get("path")
        for item in harmonization.get("audit_artifacts", [])
        if isinstance(item, dict)
    ]
    expected_artifact_paths = {
        "analysis/v3_validation/generated/cell_universe_manifest.tsv",
        "analysis/v3_validation/generated/barcode_normalization_audit.tsv",
        "analysis/v3_validation/generated/feature_universe.tsv",
        "analysis/v3_validation/generated/count_layer_map.tsv",
        "analysis/v3_validation/generated/denominator_audit.tsv",
        "analysis/v3_validation/generated/workflow_failures.tsv",
        "analysis/v3_validation/generated/endpoint_comparability.tsv",
    }
    if len(artifact_paths) != len(set(artifact_paths)):
        errors.append("harmonization audit artifacts contain duplicate paths")
    if set(artifact_paths) != expected_artifact_paths:
        missing = sorted(expected_artifact_paths - set(artifact_paths))
        extra = sorted(set(artifact_paths) - expected_artifact_paths)
        errors.append(
            f"harmonization audit artifact registry mismatch; missing={missing}, extra={extra}"
        )

    observed_digest = harmonization.get("contract_sha256")
    expected_digest = harmonization_sha256(harmonization)
    if observed_digest != expected_digest:
        errors.append("harmonization contract_sha256 does not match the canonical SCI-02 contract")
    observed_dependencies_digest = harmonization.get("dependent_fields_sha256")
    expected_dependencies_digest = harmonization_dependencies_sha256(document)
    if observed_dependencies_digest != expected_dependencies_digest:
        errors.append(
            "harmonization dependent_fields_sha256 does not match SCI-02-dependent protocol fields"
        )

    planned = document.get("planned_freeze_sections", {})
    sci02_sections = ("cell_universe", "feature_universe", "count_layers_denominators")
    is_frozen = harmonization.get("status") == "frozen"
    if isinstance(planned, dict):
        for section_name in sci02_sections:
            section = planned.get(section_name, {})
            section_status = section.get("status") if isinstance(section, dict) else None
            if is_frozen and section_status != "frozen":
                errors.append(
                    f"frozen harmonization requires planned section {section_name!r} to be frozen"
                )
            if not is_frozen and section_status == "frozen":
                errors.append(
                    f"planned section {section_name!r} cannot be frozen before harmonization"
                )

    readiness = document.get("execution_readiness", {})
    blockers = readiness.get("training_blockers", []) if isinstance(readiness, dict) else []
    blocker_ids = {item.get("id") for item in blockers if isinstance(item, dict) and item.get("id")}
    if is_frozen and "harmonization" in blocker_ids:
        errors.append("frozen harmonization cannot remain a training blocker")
    if not is_frozen and "harmonization" not in blocker_ids:
        errors.append("unfrozen harmonization must remain a training blocker")

    cell_universe = harmonization.get("cell_universe", {})
    if not isinstance(cell_universe, dict):
        cell_universe = {}
    feature_universe = harmonization.get("feature_universe", {})
    if not isinstance(feature_universe, dict):
        feature_universe = {}
    count_layers = harmonization.get("count_layers", {})
    cross_tool = (
        count_layers.get("cross_tool_primary", {}) if isinstance(count_layers, dict) else {}
    )
    scope = document.get("supported_scope", {})
    if not isinstance(scope, dict):
        scope = {}
    if isinstance(cross_tool, dict):
        if cross_tool.get("viralscan_source") != scope.get("cross_tool_primary_layer"):
            errors.append(
                "cross-tool primary ViralScan layer must match supported_scope.cross_tool_primary_layer"
            )
        if cross_tool.get("cell_universe") != cell_universe.get("id"):
            errors.append("cross-tool primary layer references the wrong cell universe")
        if cross_tool.get("feature_universe") != feature_universe.get("id"):
            errors.append("cross-tool primary layer references the wrong feature universe")

    denominators = {
        item.get("id"): item
        for item in harmonization.get("denominators", [])
        if isinstance(item, dict) and item.get("id")
    }
    expected_denominator_fields = {
        "D1_count_invariant_runs": (
            "every prespecified count-producing workflow row whose tool declares a "
            "molecule-conservation invariant, including unattempted and noncomplete such rows",
            "complete run matrix and audit",
            "selected-method-X-and-count-audit",
        ),
        "D2_unique_molecule_precision": (
            "all reported unique target molecules including assignments outside truth as false positives",
            "frozen shared cells and common target genes",
            "counts_unique",
        ),
        "D3_unique_molecule_recall": (
            "all eligible planted target molecules on the same frozen cells and genes including unresolved or lost truth as false negatives",
            "frozen shared cells and common target genes",
            "counts_unique",
        ),
        "D4_infected_cell_precision": (
            "all called-positive cells at that tier within the frozen shared anchor",
            "shared-called-cell-anchor",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D5_infected_cell_recall": (
            "all exact truth-infected cells within the frozen shared anchor",
            "shared-called-cell-anchor",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D6_infected_cell_auprc": (
            "every truth-labelled barcode in the frozen shared called-cell anchor",
            "shared-called-cell-anchor",
            "continuous-cell-score-id-frozen-in-SCI-03",
        ),
        "D7_exact_negative_cell_false_calls": (
            "every barcode in each exact-negative frozen shared called-cell anchor",
            "exact synthetic host-only and planted host-homology negatives only",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D8_exact_negative_sample_false_calls": (
            "every exact-negative biological sample in the frozen partition",
            "exact synthetic host-only and planted host-homology negative samples",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D9_exact_negative_molecule_false_calls": (
            "all selected-method host-plus-virus molecules on the same frozen cells and common total-expression genes",
            "exact synthetic host-only and planted host-homology negatives only",
            "viralscan-product-evidence-X-host-conservative",
        ),
        "D10_presumed_negative_unexpected_calls": (
            "every barcode in each presumed-negative frozen shared called-cell anchor",
            "observational or reagent presumed-negative controls",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D11_empty_droplet_false_calls": (
            "every prespecified frozen empty-droplet barcode",
            "exact synthetic empty-droplet evaluation universe",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D12_positive_cell_rate": (
            "every barcode in the frozen shared called-cell anchor",
            "shared-called-cell-anchor",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D13_sibling_molecule_confusion": (
            "every eligible canonical planted target molecule for the prespecified virus pair",
            "exact synthetic sibling-virus strata on frozen cells and common target genes",
            "counts_unique-plus-molecule-ambiguity-audit",
        ),
        "D14_host_virus_allocation": (
            "every planted mixed host-virus corrected CB-UMI truth key",
            "exact synthetic mixed host-virus molecule truth",
            "selected-method-X-and-count-audit",
        ),
        "D15_two_step_truth_read_loss": (
            "every eligible planted target fragment entering that boundary",
            "exact paired read truth for each matched combined and two-step run",
            "exact-read-lineage",
        ),
        "D16_two_step_truth_molecule_loss": (
            "every eligible planted corrected CB-UMI target molecule entering that boundary",
            "exact molecule truth for each matched combined and two-step run",
            "exact-read-lineage-molecule-survival",
        ),
        "D17_lod_detection_probability": (
            "every planned valid biological replicate-virus row at each abundance level within a chemistry",
            "frozen synthetic LOD grid per chemistry, with homology reported as a covariate "
            "rather than as a separate fit, and technical failures retained separately",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D18_tier_specificity": (
            "all exact truth-negative anchor cells equal to true negatives plus false positives",
            "exact synthetic truth-labelled shared cell anchors",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D19_tier_ppv": (
            "all tier-positive exact truth-labelled anchor cells equal to true positives plus false positives",
            "exact synthetic truth-labelled shared cell anchors",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
        "D20_cross_tool_absolute_unique": (
            "not-applicable-absolute-count",
            "frozen shared cells and common target genes",
            "counts_unique-or-documented-gene-unique-integer-equivalent",
        ),
        "D21_secondary_burden_per_10000": (
            "all molecules from that same layer on the frozen common host-plus-virus genes",
            "frozen shared cells and common total-expression genes",
            "same-layer-for-numerator-and-denominator",
        ),
        "D22_workflow_technical_failure_rate": (
            "every prespecified workflow row",
            "complete workflow matrix",
            "not-applicable",
        ),
        "D23_endpoint_incomparability_rate": (
            "every prespecified applicable row-endpoint pair excluding predeclared not-applicable pairs",
            "complete workflow-by-endpoint matrix",
            "not-applicable",
        ),
        "D24_sibling_cell_confusion": (
            "every truth-labelled anchor cell for the prespecified virus pair",
            "exact synthetic sibling-virus strata on the frozen shared cell anchor",
            "viralscan-product-evidence-X-host-conservative-plus-read-qc",
        ),
    }
    for denominator_id, expected_fields in expected_denominator_fields.items():
        item = denominators.get(denominator_id, {})
        observed_fields = (
            item.get("denominator"),
            item.get("population"),
            item.get("count_layer"),
        )
        if observed_fields != expected_fields:
            errors.append(
                f"denominator {denominator_id!r} does not match its frozen "
                "denominator, population, and count-layer contract"
            )

    expected_endpoint_map = {
        "E1_count_invariants": {"D1_count_invariant_runs"},
        "E2_negative_false_calls": {
            "D7_exact_negative_cell_false_calls",
            "D8_exact_negative_sample_false_calls",
            "D9_exact_negative_molecule_false_calls",
            "D11_empty_droplet_false_calls",
        },
        "E3_molecule_recovery": {
            "D2_unique_molecule_precision",
            "D3_unique_molecule_recall",
        },
        "E4_cell_recovery": {
            "D4_infected_cell_precision",
            "D5_infected_cell_recall",
            "D6_infected_cell_auprc",
        },
        "E5_sibling_confusion": {
            "D13_sibling_molecule_confusion",
            "D24_sibling_cell_confusion",
        },
        "E6_host_virus_allocation": {"D14_host_virus_allocation"},
        "E7_two_step_loss": {
            "D15_two_step_truth_read_loss",
            "D16_two_step_truth_molecule_loss",
        },
        "E8_limit_of_detection": {"D17_lod_detection_probability"},
        "E9_tier_calibration": {"D18_tier_specificity", "D19_tier_ppv"},
        "E10_cross_tool_parity": {
            "D20_cross_tool_absolute_unique",
            "D21_secondary_burden_per_10000",
            "D22_workflow_technical_failure_rate",
            "D23_endpoint_incomparability_rate",
        },
    }
    endpoint_map = harmonization.get("endpoint_denominator_map", {})
    if isinstance(endpoint_map, dict):
        observed_endpoint_map = {
            endpoint_id: set(ids) if isinstance(ids, list) else set()
            for endpoint_id, ids in endpoint_map.items()
        }
        if observed_endpoint_map != expected_endpoint_map:
            errors.append("endpoint_denominator_map does not match the frozen SCI-02 registry")

    expected_audit_columns = {
        "analysis/v3_validation/generated/cell_universe_manifest.tsv": {
            "dataset_id",
            "biological_sample_id",
            "library_id",
            "lane_id",
            "chemistry",
            "original_barcode",
            "cb_sequence",
            "gem_group",
            "canonical_key",
            "source",
            "source_sha256",
            "universe_role",
            "membership",
            "anchor_method",
            "anchor_config_sha256",
            "host_matrix_sha256",
            "whitelist_sha256",
            "environment_sha256",
        },
        "analysis/v3_validation/generated/barcode_normalization_audit.tsv": {
            "dataset_id",
            "library_id",
            "original_barcode",
            "normalized_sequence",
            "gem_group",
            "canonical_key",
            "collision",
        },
        "analysis/v3_validation/generated/feature_universe.tsv": {
            "workflow_row_id",
            "source_feature_id",
            "canonical_gene_id",
            "canonical_virus_id",
            "feature_type",
            "accession_version",
            "sequence_sha256",
            "mapping_cardinality",
            "quantifiable",
            "included_primary",
            "reason",
        },
        "analysis/v3_validation/generated/count_layer_map.tsv": {
            "workflow_row_id",
            "endpoint_id",
            "source_artifact",
            "source_artifact_sha256",
            "source_layer",
            "unit",
            "integer_status",
            "umi_collapse_semantics",
            "ambiguity_rule",
            "ambiguity_model",
            "barcode_domain",
            "cell_universe_sha256",
            "feature_manifest_sha256",
            "feature_universe_sha256",
        },
        "analysis/v3_validation/generated/denominator_audit.tsv": {
            "workflow_row_id",
            "endpoint_id",
            "numerator_name",
            "numerator_value",
            "numerator_feature_universe_sha256",
            "denominator_name",
            "denominator_value",
            "denominator_feature_universe_sha256",
            "cell_universe_sha256",
            "truth_manifest_sha256",
            "layer_map_sha256",
        },
        "analysis/v3_validation/generated/workflow_failures.tsv": {
            "workflow_row_id",
            "attempt_id",
            "status",
            "stage",
            "exit_code",
            "reason",
            "log_path",
        },
        # R6-F5: the not-applicable semantics and the endpoint-incomparability
        # rate had no source column until this table was declared.
        "analysis/v3_validation/generated/endpoint_comparability.tsv": {
            "workflow_row_id",
            "workflow_id",
            "tool",
            "endpoint_id",
            "count_layer",
            "comparability_status",
            "reason",
            "adapter_id",
            "adapter_declared_before_outcomes",
        },
    }
    observed_audit_columns = {
        item.get("path"): set(item.get("required_columns", []))
        for item in harmonization.get("audit_artifacts", [])
        if isinstance(item, dict)
    }
    if observed_audit_columns != expected_audit_columns:
        errors.append("harmonization audit columns do not match the frozen SCI-02 registry")

    return errors


def frozen_inputs_sha256(document: dict[str, Any]) -> str:
    """Digest every already-frozen seed value and factor level set.

    Seeds and factors are not sections with their own status and digest, but they
    are load-bearing for reproducibility: a frozen split seed that can be silently
    rewritten makes the partition unreproducible, and the protocol's own rule that
    seeds are never selected after viewing outcomes would have no evidence behind
    it. SCI-05 round 2 (R2-F1) reproduced exactly that. This digest is computed
    from the document rather than from a section, so it covers both.
    """
    seeds = document.get("seeds", {})
    frozen_seeds = {
        name: value.get("value")
        for name, value in (seeds.items() if isinstance(seeds, dict) else [])
        if isinstance(value, dict) and value.get("status") == "frozen"
    }
    frozen_factors = {
        item.get("id"): item.get("levels")
        for item in document.get("factors", [])
        if isinstance(item, dict) and item.get("status") == "frozen"
    }

    # R3-F2: dataset and reference identity sat outside every digest, so an SRA
    # accession could be swapped after the partition was drawn with no error.
    # harmonization_dependencies_sha256 hashes only a partial dataset projection.
    def _identity(container: str) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for item in document.get(container, []):
            if not isinstance(item, dict):
                continue
            out[item.get("id")] = {
                "accessions": item.get("accessions"),
                "chemistries": item.get("chemistries"),
                "source": item.get("source"),
                "profile": item.get("profile"),
                "assets": [
                    {
                        "id": asset.get("id"),
                        "locator": asset.get("locator"),
                        "digest_status": asset.get("digest_status"),
                        "sha256": asset.get("sha256"),
                    }
                    for asset in item.get("assets", [])
                    if isinstance(asset, dict)
                ],
            }
        return out

    # The section's own declarations must be inside its digest, or the list of
    # records it declares could be trimmed without changing anything. Its
    # contract_sha256 is excluded to avoid circularity.
    section = document.get("frozen_inputs")
    self_declaration = (
        {key: value for key, value in section.items() if key != "contract_sha256"}
        if isinstance(section, dict)
        else None
    )

    return _canonical_sha256(
        {
            "self": self_declaration,
            "seeds": frozen_seeds,
            "factors": frozen_factors,
            "datasets": _identity("datasets"),
            "references": _identity("references"),
            # R3-F3: blockers could be deleted and training_allowed flipped with
            # no error at the draft phase CI actually runs.
            "execution_readiness": document.get("execution_readiness"),
        }
    )


# Digests computed from the whole document rather than from one section.
_DOCUMENT_DIGEST_FNS: dict[str, Any] = {
    "frozen_inputs": frozen_inputs_sha256,
}

_SECTION_DIGEST_FNS: dict[str, Any] = {
    "harmonization": harmonization_sha256,
    "partitions": partitions_sha256,
    "calibration": calibration_sha256,
    "workflow_matrix": workflow_matrix_sha256,
    "failure_and_deviation_reporting": failure_reporting_sha256,
}


def _ledger_chain_for_section(ledger: dict[str, Any], section: str) -> list[tuple[str, str]]:
    """Return this section's (before, after) digest pairs in ledger order."""
    scope = f"{section}.contract_sha256"
    chain: list[tuple[str, str]] = []
    for record in ledger.get("deviations", []):
        if not isinstance(record, dict):
            continue
        if record.get("digest_scope") == scope:
            chain.append(
                (record.get("protocol_sha256_before"), record.get("protocol_sha256_after"))
            )
        extra = record.get("additional_digest_changes") or {}
        if isinstance(extra, dict) and scope in extra:
            change = extra[scope] or {}
            chain.append((change.get("before"), change.get("after")))
    return chain


def ledger_record_sha256(record: dict[str, Any], previous: str) -> str:
    """Hash one ledger record together with its predecessor's hash.

    Chaining makes an *inconsistent* edit detectable: editing record N changes its
    own hash and invalidates the chain from N+1 onward, so history cannot be
    rewritten by touching a single record.

    It does not make the ledger append-only in fact. An author who re-chains the
    whole file after editing produces a self-consistent ledger that this check
    accepts; SCI-05 round 4 (R4-F1) reproduced exactly that. The external anchor
    is `scripts/check_ledger_append_only.py`, which compares against committed
    history rather than against the file itself.
    """
    payload = {key: value for key, value in record.items() if key != "record_sha256"}
    return _canonical_sha256({"previous": previous, "record": payload})


def validate_ledger_integrity(ledger: dict[str, Any]) -> list[str]:
    """Verify the ledger's own record chain.

    SCI-05 round 3 (R3-F1) reproduced an edit-in-place that left no trace: the
    amendment rail rested on a ledger whose history could be silently rewritten.
    """
    errors: list[str] = []
    previous = ""
    for index, record in enumerate(ledger.get("deviations", [])):
        if not isinstance(record, dict):
            errors.append(f"deviation ledger record {index} is not a mapping")
            return errors
        record_id = record.get("deviation_id", index)
        claimed = record.get("record_sha256")
        if not claimed:
            errors.append(f"deviation record {record_id!r} has no record_sha256")
            return errors
        expected = ledger_record_sha256(record, previous)
        if claimed != expected:
            errors.append(
                f"deviation record {record_id!r} breaks the ledger chain; it was "
                "edited in place rather than superseded by a new record"
            )
            return errors
        previous = claimed
    return errors


def validate_amendment_ledger(document: dict[str, Any], ledger: dict[str, Any]) -> list[str]:
    """Check that every frozen section's digest is accounted for in the ledger.

    This detects an *undocumented* amendment, not an *illegitimate* one. An author
    who edits a frozen section and appends a matching record still passes; that is
    authorial honesty, not tamper-evidence. Anchoring on git would have the same
    ceiling, since git also trusts the working tree at commit time.
    """
    errors: list[str] = []
    # The chain must hold before any digest claim in it can be trusted.
    errors.extend(validate_ledger_integrity(ledger))
    if errors:
        return errors

    # R4-F3: deleting a record and re-chaining left no trace, because nothing
    # declared which records ought to exist. The expected ids live inside
    # frozen_inputs, which is itself digested, so a deletion must now also edit
    # and re-digest the protocol rather than only the ledger.
    frozen_inputs = document.get("frozen_inputs")
    if isinstance(frozen_inputs, dict) and frozen_inputs.get("status") == "frozen":
        declared = list(frozen_inputs.get("recorded_deviations") or [])
        present = [
            record.get("deviation_id")
            for record in ledger.get("deviations", [])
            if isinstance(record, dict)
        ]
        missing = [item for item in declared if item not in present]
        if missing:
            errors.append("deviation ledger is missing declared records: " + ", ".join(missing))
        if present[: len(declared)] != declared:
            errors.append("deviation ledger reorders or replaces declared records")
    covered = {name: (fn, False) for name, fn in _SECTION_DIGEST_FNS.items()}
    covered.update({name: (fn, True) for name, fn in _DOCUMENT_DIGEST_FNS.items()})
    for name, (digest_fn, from_document) in covered.items():
        section = document.get(name)
        if not isinstance(section, dict) or section.get("status") != "frozen":
            # The amendment rule only binds sections that currently claim frozen.
            continue
        current = section.get("contract_sha256")
        chain = _ledger_chain_for_section(ledger, name)
        if not chain:
            errors.append(
                f"frozen section {name!r} has no deviation-ledger baseline record; "
                "add a genesis record before the amendment rule can be enforced"
            )
            continue
        for index in range(1, len(chain)):
            if chain[index][0] != chain[index - 1][1]:
                errors.append(f"deviation ledger for {name!r} is not continuous at record {index}")
        if chain[-1][1] != current:
            errors.append(
                f"frozen section {name!r} digest does not match the latest ledger "
                "record; undocumented amendment"
            )
        recomputed = digest_fn(document) if from_document else digest_fn(section)
        if recomputed != current:
            errors.append(f"frozen section {name!r} does not hash to its claimed digest")
    return errors


def validate_protocol(
    document: dict[str, Any],
    schema: dict[str, Any],
    *,
    require_frozen: bool = False,
    phase: str = "draft",
    ledger: dict[str, Any] | None = None,
) -> list[str]:
    """Return structural and cross-reference validation errors.

    ``ledger`` is optional so that in-memory callers may validate a synthetic
    document without one. ``validate_protocol_file`` always supplies it.
    """
    if phase not in {"draft", "training", "holdout"}:
        raise ValueError("phase must be draft, training, or holdout")
    if require_frozen and phase == "draft":
        phase = "training"
    try:
        Draft202012Validator.check_schema(schema)
    except SchemaError as exc:
        return [f"invalid validation schema: {exc.message}"]
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors = [
        f"schema {'.'.join(str(part) for part in error.absolute_path) or '<root>'}: {error.message}"
        for error in sorted(
            validator.iter_errors(document), key=lambda item: list(item.absolute_path)
        )
    ]

    # Structural errors can make later cross-reference walks unsafe. Fail closed
    # with the schema diagnostics before applying semantic checks.
    if errors:
        return sorted(set(errors))

    for section in (
        "hypotheses",
        "endpoints",
        "datasets",
        "factors",
        "references",
        "exclusions",
    ):
        errors.extend(_duplicate_ids(document.get(section), section))
    for container_name in ("datasets", "references"):
        for container in document.get(container_name, []):
            if isinstance(container, dict):
                errors.extend(
                    _duplicate_ids(
                        container.get("assets"),
                        f"{container_name}.{container.get('id')}.assets",
                    )
                )
    readiness_for_ids = document.get("execution_readiness", {})
    if isinstance(readiness_for_ids, dict):
        errors.extend(
            _duplicate_ids(
                readiness_for_ids.get("training_blockers"),
                "execution_readiness.training_blockers",
            )
        )
        errors.extend(
            _duplicate_ids(
                readiness_for_ids.get("holdout_blockers"),
                "execution_readiness.holdout_blockers",
            )
        )

    scope = document.get("supported_scope", {})
    if isinstance(scope, dict):
        errors.extend(_duplicate_ids(scope.get("chemistries"), "supported_scope.chemistries"))
        expected_geometry = {
            "10xv2": (16, 10),
            "10xv3": (16, 12),
            "drop-seq": (12, 8),
        }
        for chemistry in scope.get("chemistries", []):
            if not isinstance(chemistry, dict) or chemistry.get("id") not in expected_geometry:
                continue
            observed = (chemistry.get("cell_barcode_length"), chemistry.get("umi_length"))
            expected = expected_geometry[chemistry["id"]]
            if observed != expected:
                errors.append(
                    f"chemistry {chemistry['id']!r} geometry {observed} does not match {expected}"
                )

    endpoint_ids = {
        item.get("id") for item in document.get("endpoints", []) if isinstance(item, dict)
    }
    for hypothesis in document.get("hypotheses", []):
        if not isinstance(hypothesis, dict):
            continue
        for endpoint_id in hypothesis.get("endpoint_ids", []):
            if endpoint_id not in endpoint_ids:
                errors.append(
                    f"hypothesis {hypothesis.get('id')!r} references unknown endpoint {endpoint_id!r}"
                )

    seed_values = []
    seeds = document.get("seeds", {})
    if isinstance(seeds, dict):
        for name, decision in seeds.items():
            if name == "derivation" or not isinstance(decision, dict):
                continue
            value = decision.get("value")
            if isinstance(value, int):
                seed_values.append((name, value))
    duplicated_seeds = sorted(
        {value for _, value in seed_values if [x[1] for x in seed_values].count(value) > 1}
    )
    if duplicated_seeds:
        errors.append(f"frozen stage seeds must be distinct; duplicates: {duplicated_seeds}")

    errors.extend(_validate_harmonization(document))
    errors.extend(_validate_sci03(document))
    errors.extend(_validate_sci04(document, phase))
    if ledger is not None:
        errors.extend(validate_amendment_ledger(document, ledger))

    if document.get("status") == "frozen" or phase != "draft":
        if document.get("status") != "frozen":
            errors.append("protocol is not frozen")
        readiness = document.get("execution_readiness", {})
        if not readiness.get("training_allowed"):
            errors.append("training execution is not allowed")
        if readiness.get("training_blockers"):
            errors.append("training execution blockers remain")
        freeze = document.get("freeze", {})
        if not freeze.get("training_execution_allowed"):
            errors.append("freeze record does not allow training execution")

        for container_name in ("datasets", "references"):
            for container in document.get(container_name, []):
                if not isinstance(container, dict):
                    continue
                for asset in container.get("assets", []):
                    if (
                        isinstance(asset, dict)
                        and asset.get("freeze_required")
                        and asset.get("digest_status") != "verified"
                    ):
                        errors.append(
                            f"{container_name} {container.get('id')!r} has unfrozen asset {asset.get('id')!r}"
                        )
        for factor in document.get("factors", []):
            if isinstance(factor, dict) and factor.get("status") != "frozen":
                errors.append(f"factor {factor.get('id')!r} is not frozen")
        if isinstance(seeds, dict):
            for name, decision in seeds.items():
                if (
                    name != "derivation"
                    and isinstance(decision, dict)
                    and decision.get("status") != "frozen"
                ):
                    errors.append(f"seed {name!r} is not frozen")
        for name, section in document.get("planned_freeze_sections", {}).items():
            if (
                isinstance(section, dict)
                and section.get("freeze_required")
                and section.get("status") != "frozen"
            ):
                errors.append(f"planned section {name!r} is not frozen")

        if phase == "holdout":
            if not readiness.get("holdout_allowed"):
                errors.append("holdout execution is not allowed")
            if readiness.get("holdout_blockers"):
                errors.append("holdout execution blockers remain")
            if not freeze.get("holdout_execution_allowed"):
                errors.append("freeze record does not allow holdout execution")
            for field in ("training_results_sha256", "thresholds_sha256"):
                value = freeze.get(field)
                if not isinstance(value, str) or len(value) != 64:
                    errors.append(f"holdout requires frozen {field}")

    return sorted(set(errors))


def validate_protocol_file(
    protocol_path: Path = DEFAULT_PROTOCOL,
    schema_path: Path = DEFAULT_SCHEMA,
    *,
    require_frozen: bool = False,
    phase: str = "draft",
    ledger_path: Path = DEFAULT_LEDGER,
) -> list[str]:
    """Load and validate a protocol file; missing files fail closed.

    The ledger is required here rather than optional: the committed protocol is
    always checked against the committed amendment ledger.
    """
    if not protocol_path.is_file():
        return [f"protocol not found: {protocol_path}"]
    if not schema_path.is_file():
        return [f"schema not found: {schema_path}"]
    if not ledger_path.is_file():
        return [f"deviation ledger not found: {ledger_path}"]
    try:
        document = load_yaml(protocol_path)
        schema = load_yaml(schema_path)
        ledger = load_yaml(ledger_path)
        return validate_protocol(
            document,
            schema,
            require_frozen=require_frozen,
            phase=phase,
            ledger=ledger,
        )
    except (OSError, TypeError, ValueError, yaml.YAMLError) as exc:
        return [str(exc)]


def require_execution_allowed(
    phase: str,
    *,
    protocol_path: Path = DEFAULT_PROTOCOL,
    schema_path: Path = DEFAULT_SCHEMA,
    ledger_path: Path = DEFAULT_LEDGER,
) -> None:
    """Fail closed unless the protocol permits an outcome-generating run at ``phase``.

    This is the single reusable choke point for any future SCI-04 workflow-row
    executor. It must be called both when a run manifest is built and immediately
    before each row's subprocess executes; see
    ``workflow_matrix.environment_pinning.execution_gate`` in the protocol.

    It is deliberately opt-in: nothing in the ordinary ``viralscan`` CLI calls it,
    because the research protocol is not shipped in the installed package. There
    is no bypass flag — an escape hatch on an integrity rail is that rail's own
    failure mode.
    """
    errors = validate_protocol_file(
        protocol_path,
        schema_path,
        phase=phase,
        ledger_path=ledger_path,
    )
    if errors:
        joined = "\n".join(f"  - {error}" for error in errors)
        raise SystemExit(f"protocol execution gate failed for phase {phase!r}:\n{joined}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protocol", nargs="?", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument(
        "--require-frozen",
        action="store_true",
        help="deprecated alias for --phase training",
    )
    parser.add_argument(
        "--phase",
        choices=("draft", "training", "holdout"),
        default="draft",
        help="validation gate to enforce (default: draft structure only)",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    errors = validate_protocol_file(
        args.protocol,
        args.schema,
        require_frozen=args.require_frozen,
        phase=args.phase,
    )
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"valid: {args.protocol}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
