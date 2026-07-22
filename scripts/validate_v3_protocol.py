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
            "every prespecified count-producing workflow row including unattempted and noncomplete rows",
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
            "every planned valid biological replicate-virus row at each abundance chemistry and homology stratum",
            "frozen synthetic LOD grid with technical failures retained separately",
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
    }
    observed_audit_columns = {
        item.get("path"): set(item.get("required_columns", []))
        for item in harmonization.get("audit_artifacts", [])
        if isinstance(item, dict)
    }
    if observed_audit_columns != expected_audit_columns:
        errors.append("harmonization audit columns do not match the frozen SCI-02 registry")

    return errors


def validate_protocol(
    document: dict[str, Any],
    schema: dict[str, Any],
    *,
    require_frozen: bool = False,
    phase: str = "draft",
) -> list[str]:
    """Return structural and cross-reference validation errors."""
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
) -> list[str]:
    """Load and validate a protocol file; missing files fail closed."""
    if not protocol_path.is_file():
        return [f"protocol not found: {protocol_path}"]
    if not schema_path.is_file():
        return [f"schema not found: {schema_path}"]
    try:
        document = load_yaml(protocol_path)
        schema = load_yaml(schema_path)
        return validate_protocol(
            document,
            schema,
            require_frozen=require_frozen,
            phase=phase,
        )
    except (OSError, TypeError, ValueError, yaml.YAMLError) as exc:
        return [str(exc)]


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
