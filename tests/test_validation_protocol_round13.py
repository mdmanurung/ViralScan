"""SCI-05 round 13 protocol-side regressions (DEV-029): M6, M7, BM1, BM3, BM4, minors."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.validate_v3_protocol import DEFAULT_LEDGER, load_yaml

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = ROOT / "analysis" / "v3_validation" / "protocol.yaml"


@pytest.fixture(scope="module")
def protocol() -> dict:
    return load_yaml(PROTOCOL_PATH)


def test_g3_freezes_probable_only_and_names_what_is_out(protocol: dict) -> None:
    scope = protocol["calibration"]["freeze_scope"]
    for needle in ("probable tier only", "cumulative", "AUPRC", "support fraction", "M8"):
        assert needle in scope
    blockers = {b["id"] for b in protocol["execution_readiness"]["training_blockers"]}
    assert "tier_cutoff_measures" in blockers


def test_m8_nominal_confidence_metric_is_gone(protocol: dict) -> None:
    ids = {m["id"] for m in protocol["calibration"]["metrics"]}
    assert not any(i.startswith("M8") for i in ids)


def test_real_data_denominators_feed_a_descriptive_endpoint(protocol: dict) -> None:
    mapping = protocol["harmonization"]["endpoint_denominator_map"]
    assert set(mapping["E11_real_data_descriptive"]) == {
        "D10_presumed_negative_unexpected_calls",
        "D12_positive_cell_rate",
    }
    mapped = {d for ids in mapping.values() for d in ids}
    for dataset_denominator in ("D10_presumed_negative_unexpected_calls", "D12_positive_cell_rate"):
        assert dataset_denominator in mapped
    hypotheses = {h["id"]: h for h in protocol["hypotheses"]}
    assert hypotheses["H8_real_data_descriptive"]["endpoint_ids"] == ["E11_real_data_descriptive"]
    endpoints = {e["id"]: e for e in protocol["endpoints"]}
    assert "No accuracy" in endpoints["E11_real_data_descriptive"]["acceptance_rule"]
    assert "recovery metrics" in endpoints["E10_cross_tool_parity"]["reporting"]
    assert "E11" in endpoints["E10_cross_tool_parity"]["reporting"]


def test_anellovirus_gate_no_longer_implies_a_sensitivity_claim(protocol: dict) -> None:
    rule = next(d for d in protocol["datasets"] if d["id"] == "anellovirus_orthogonal_positive")[
        "exclusion_rule"
    ]
    assert "No real-sensitivity claim is made whatever" in rule


def test_unknown_origin_is_never_a_false_positive(protocol: dict) -> None:
    rules = protocol["harmonization"]["derived_metric_rules"]
    assert "never as a false positive" in rules["unknown_origin"]
    assert rules["f1"].startswith("two-times-true-positives")
    denominators = {d["id"]: d for d in protocol["harmonization"]["denominators"]}
    for key in ("D2_unique_molecule_precision", "D4_infected_cell_precision"):
        assert "unknown origin" in denominators[key]["denominator"]
    assert "unknown-origin" in " ".join(protocol["failure_and_deviation_reporting"]["prohibited"])


def test_degenerate_bootstrap_has_a_reason_coded_policy(protocol: dict) -> None:
    policy = protocol["calibration"]["uncertainty"]["degenerate_interval_policy"]
    for needle in ("constant-bootstrap-distribution", "undefined-statistic", "zero-width"):
        assert needle in policy


def test_f1_objective_freezes_unit_strata_and_pooling(protocol: dict) -> None:
    objective = protocol["calibration"]["threshold_search"]["objective"]
    for needle in ("(cell, target", "micro pooling", "synthetic_factorial"):
        assert needle in objective


def test_synthetic_factorial_states_the_per_library_crossing(protocol: dict) -> None:
    rule = next(d for d in protocol["datasets"] if d["id"] == "synthetic_factorial")[
        "planting_rule"
    ]
    assert "full abundance by host-virus-homology grid" in rule
    assert "low_complexity" in rule


def test_dev_029_is_ledgered_and_declared(protocol: dict) -> None:
    ledger = load_yaml(DEFAULT_LEDGER)
    record = ledger["deviations"][-1]
    assert record["deviation_id"] == "DEV-029"
    assert "DEV-029" in protocol["frozen_inputs"]["recorded_deviations"]
    assert record["outcome_triggered"] is False
    assert record["digest_scope"] == "harmonization.contract_sha256"


def test_g3_attestation_wording_and_known_limitation(protocol: dict) -> None:
    rule = protocol["pre_freeze_runs"]["outcome_definition"]["scored_dataset_rule"]
    assert "no recorded pre-freeze outcome on scored or protected pending-row datasets" in rule
    risks = " ".join(protocol["principal_risks"])
    assert "protocol v2" in risks and "29 percent" in risks and "42 percent" in risks


def test_tier_cutoff_measures_are_defined_and_still_block_training(protocol: dict) -> None:
    """Q135: the three measures are written; the implementation is pending, so training stays blocked."""
    scope = protocol["calibration"]["freeze_scope"]
    for needle in (
        "m(c,v) >= M and F(s,v) >= f and S(s,v) >= s*",
        "at least 1 molecule",
        "virus_best, host_best or tie",
        "zero denominator makes S not-estimable",
        "broadcast to every cell",
    ):
        assert needle in " ".join(scope.split())
    blocker = next(
        b
        for b in protocol["execution_readiness"]["training_blockers"]
        if b["id"] == "tier_cutoff_measures"
    )
    assert "implementation is pending" in " ".join(blocker["description"].split())
