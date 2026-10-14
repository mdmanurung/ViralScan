"""SCI-05 round 15 (DEV-032): Q139 a bracket fallback, N2 to N6 and the Biomni LOD overflow, one test each."""

import copy
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import check_g3_diff as g3  # noqa: E402
import v3_lod as lod  # noqa: E402
import v3_score as score  # noqa: E402
import validate_v3_protocol as gate  # noqa: E402

LEVELS = [1, 3, 10, 30, 100, 300, 1000]


@pytest.fixture(scope="module")
def protocol() -> dict:
    return gate.load_yaml(gate.DEFAULT_PROTOCOL)


def _text(protocol, *path) -> str:
    node = protocol
    for key in path:
        node = node[key]
    return " ".join(node.split())


# ── Q139 a: model-free bracket ────────────────────────────────────────────────


def _rows(k_by_level, n=4):
    return [(lv, int(i < k)) for lv, k in k_by_level.items() for i in range(n)]


def test_bracket_runs_from_highest_level_below_95_to_lowest_level_at_or_above():
    b = lod.lod_bracket(_rows({1: 0, 3: 1, 10: 3, 30: 4, 100: 4}, n=4))
    assert (b["status"], b["lower"], b["upper"], b["non_monotone"]) == ("bracketed", 10, 30, False)
    # 19 of 20 is exactly 95 percent and reaches it; 18 of 20 does not
    twenty = [(10, int(i < 19)) for i in range(20)] + [(1, 0)] * 20
    assert lod.lod_bracket(twenty)["upper"] == 10
    assert lod.lod_bracket([(10, int(i < 18)) for i in range(20)] + [(1, 0)] * 20)["upper"] is None


def test_bracket_reports_only_a_lower_bound_when_no_level_reaches_95():
    b = lod.lod_bracket(_rows({1: 0, 10: 1, 100: 3}))
    assert (b["status"], b["lower"], b["upper"]) == ("above_range", 100, None)


def test_bracket_reports_only_an_upper_bound_when_every_level_reaches_95_and_flags_a_dip():
    b = lod.lod_bracket(_rows({1: 4, 10: 4}))
    assert (b["status"], b["lower"], b["upper"]) == ("below_range", None, 1)
    dip = lod.lod_bracket(_rows({1: 4, 10: 0, 100: 4}))  # 95 percent at 1, below at 10
    assert (dip["lower"], dip["upper"], dip["non_monotone"]) == (10, 100, True)
    with pytest.raises(ValueError):
        lod.lod_bracket([])


def _probe_units(lod95, slope, seed, n_units=10, per_level=2):
    """The round-15 probe: 10 libraries x 7 levels x 2 samples per level (140 rows), probit truth."""
    rng = np.random.default_rng(1000 + seed)
    l50 = math.log10(lod95) - 1.645 / slope
    return {
        f"u{u}": [
            (float(x), int(rng.random() < stats.norm.cdf(slope * (math.log10(x) - l50))))
            for x in LEVELS
            for _ in range(per_level)
        ]
        for u in range(n_units)
    }


@pytest.mark.parametrize("lod95,slope,seed", [(100, 5, 0), (30, 5, 1), (30, 2, 2)])
def test_a_not_estimable_interval_carries_its_reason_and_a_bracket_that_matches_the_data(
    lod95, slope, seed
):
    units = _probe_units(lod95, slope, seed)
    out = lod.lod_interval(units, seed=seed, strata=dict.fromkeys(units, "all"), replicates=300)
    if out["interval_status"] == "ok":  # slope 2 can be estimable; the fallback is then absent
        assert "bracket" not in out
        return
    assert out["reason_code"] and "bracket" in out
    b, rows = out["bracket"], [r for u in units.values() for r in u]
    for entry in b["levels"]:
        at = [d for a, d in rows if a == entry["level"]]
        assert (entry["detected"], entry["n"]) == (sum(at), len(at)) == (sum(at), 20)
    if b["lower"] is not None and b["upper"] is not None:
        assert b["lower"] < b["upper"]


def test_the_review_probe_slope_5_is_not_estimable_with_a_bracket_on_every_seed():
    for seed in range(3):
        units = _probe_units(100, 5, seed)
        out = lod.lod_interval(units, seed=seed, strata=dict.fromkeys(units, "a"), replicates=200)
        assert out["interval_status"] == "not_estimable" and out["bracket"]["levels"]


def test_an_ok_interval_has_no_bracket_and_a_failed_point_fit_still_gets_one():
    units = {f"s{lv}_{i}": [(lv, int(i % 3 != 0 or lv > 10))] for lv in LEVELS for i in range(12)}
    sep = lod.lod_interval(
        {k: v for k, v in units.items()}, seed=1, strata=dict.fromkeys(units, "a"), replicates=50
    )
    assert sep["interval_status"] in {"ok", "not_estimable"}
    separated = {f"u{i}": [(1, 0), (10, 0), (100, 1), (1000, 1)] for i in range(4)}
    out = lod.lod_interval(separated, seed=1, strata=dict.fromkeys(separated, "a"), replicates=20)
    assert out["reason_code"] == "point-estimate-failed" and out["point"]["status"] == "separation"
    assert out["bracket"]["status"] == "bracketed"
    assert (out["bracket"]["lower"], out["bracket"]["upper"]) == (10, 100)


def test_the_policy_text_states_the_bracket_and_the_new_reason_codes(protocol):
    fit = _text(protocol, "calibration", "limit_of_detection", "failed_fit_policy")
    for needle in (
        "model-free bracket",
        "highest tested level",
        "only the lower bound",
        "point-estimate-failed",
        "failed-fit",
        "nonfinite_lod",
        "no fit",
    ):
        assert needle in fit, needle
    interval = _text(protocol, "calibration", "uncertainty", "degenerate_interval_policy")
    assert "point-estimate-failed" in interval and "failed-fit" in interval
    e8 = next(e for e in protocol["endpoints"] if e["id"] == "E8_limit_of_detection")
    assert "model-free bracket" in e8["reporting"] and "not estimable" in e8["reporting"]


# ── Biomni residual Major: log-space range check ─────────────────────────────


def _tiny_slope(monkeypatch, only_after=0):
    real, calls = lod.fit_probit, {"n": 0}

    def fake(abundance, detected):
        calls["n"] += 1
        if calls["n"] <= only_after:
            return real(abundance, detected)
        return dict(status="ok", intercept=-5.0, slope=1e-5, n=len(abundance))

    monkeypatch.setattr(lod, "fit_probit", fake)


def test_a_tiny_positive_slope_is_not_reached_instead_of_overflowing(monkeypatch):
    _tiny_slope(monkeypatch)
    out = lod.lod95([1, 10, 100], [0, 1, 1])  # 10 ** 5e5 used to raise OverflowError
    assert out["status"] == "not_reached" and out["lod"] is None
    monkeypatch.setattr(
        lod, "fit_probit", lambda a, d: dict(status="ok", intercept=5.0, slope=1e-5, n=3)
    )
    assert lod.lod95([1, 10, 100], [0, 1, 1])["status"] == "below_range"
    monkeypatch.setattr(
        lod, "fit_probit", lambda a, d: dict(status="ok", intercept=math.nan, slope=1.0, n=3)
    )
    assert lod.lod95([1, 10, 100], [0, 1, 1])["status"] == "nonfinite_lod"


def test_the_overflow_is_a_reason_coded_result_on_the_point_and_the_bootstrap_paths(monkeypatch):
    units = {f"u{i}": [(1, 0), (10, i % 2), (100, int(i < 3)), (1000, 1)] for i in range(4)}
    strata = dict.fromkeys(units, "a")
    _tiny_slope(monkeypatch, only_after=1)  # point fit real, every replicate tiny slope
    boot = lod.lod_interval(units, seed=1, strata=strata, replicates=5)
    assert boot["interval_status"] == "not_estimable" and boot["reason_code"] == "lod-not-reached"
    assert boot["failed_fit_statuses"] == {"not_reached": 5}
    _tiny_slope(monkeypatch)  # point fit itself overflows
    point = lod.lod_interval(units, seed=1, strata=strata, replicates=5)
    assert point["reason_code"] == "point-estimate-failed"
    assert point["point"]["status"] == "not_reached" and "bracket" in point


# ── N4: strata required ───────────────────────────────────────────────────────


def test_lod_interval_requires_strata_for_every_unit():
    units = {f"u{i}": [(1, 0), (10, 1), (100, 1)] for i in range(3)}
    with pytest.raises(TypeError):
        lod.lod_interval(units, seed=1)  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="strata has no label"):
        lod.lod_interval(units, seed=1, strata={"u0": "a"}, replicates=5)


# ── N2: rowless exact library refused; blocker named ─────────────────────────


def _host(dataset="synthetic_host_only"):
    return dict(
        dataset_id=dataset, origin="unknown", planted_virus_id="", cell_barcode="h", umi="h"
    )


def test_an_exact_library_without_truth_rows_is_refused_with_a_reason_code(protocol):
    bare = dict(sample_id="s1", dataset="synthetic_host_only", tier="T1", called={}, truth_rows=[])
    with pytest.raises(ValueError, match=score.EXACT_WITHOUT_ROWS):
        score.score_by_truth_status([bare], protocol=protocol)
    missing_key = {k: v for k, v in bare.items() if k != "truth_rows"}
    with pytest.raises(ValueError, match="exact-dataset-without-truth-rows"):
        score.score_by_truth_status([missing_key], protocol=protocol)
    # rows bind it; a presumed-negative library may still be rowless (descriptive only)
    backed = dict(bare, truth_rows=[_host()])
    pbmc = dict(sample_id="p", dataset="pbmc_10x_healthy_v3", tier="T4", called={}, truth_rows=[])
    out = score.score_by_truth_status([backed, pbmc], protocol=protocol)
    assert (
        out["exact"]["T1"]["unconditional"]["n"] == 1
        and "pbmc_10x_healthy_v3" in out["descriptive"]
    )


def test_run_provenance_binding_is_a_named_training_blocker_with_all_three_bindings(protocol):
    blockers = {b["id"]: b for b in protocol["execution_readiness"]["training_blockers"]}
    text = " ".join(blockers["run_provenance_binding"]["description"].split())
    for needle in (
        "bind the library's dataset to the run manifest",
        "exact-dataset-without-truth-rows",
        "cell-calling method per library",
        "hash the GTF or panel behind D1",
    ):
        assert needle in text, needle
    assert blockers["run_provenance_binding"]["resolution_task"] == "SCI-04"
    assert protocol["execution_readiness"]["training_allowed"] is False


# ── N5: scorer record hygiene ─────────────────────────────────────────────────


def test_precision_recall_f1_reports_the_unfiltered_recall_denominator():
    out = score.precision_recall_f1({"c1", "c2", "c3", "c4"}, {"c1", "c2"}, twin_excluded={"c1"})
    assert (out["tp"], out["fn"]) == (1, 2)  # twin-filtered
    assert (out["n_recalled"], out["n_truth"], out["recall"]) == (2, 4, 0.5)


def test_an_exact_dataset_requires_tier_t1(protocol):
    lib = dict(
        sample_id="s1", dataset="synthetic_host_only", tier="T3", called={}, truth_rows=[_host()]
    )
    with pytest.raises(ValueError, match="tier must be T1"):
        score.score_by_truth_status([lib], protocol=protocol)


# ── N6: G3 signing diff check ─────────────────────────────────────────────────


def _g3_edit(protocol):
    """The edits a G3 signing commit makes, applied to a copy of the reviewed draft."""
    doc = copy.deepcopy(protocol)
    for f in doc["factors"]:
        if f.get("status") == "pending":
            f["status"] = "frozen"
            f.pop("pending_reason", None)
    doc["seeds"]["generation"].update(status="frozen", value=20260727005)
    for s in ("partitions", "calibration", "defaults_selection", "workflow_matrix"):
        doc[s].update(status="frozen", frozen_at="2026-10-15T00:00:00Z")
    doc["calibration"]["contract_sha256"] = "0" * 64
    for s in ("partitions", "calibration_metrics_lod", "workflow_matrix", "independent_review"):
        doc["planned_freeze_sections"][s]["status"] = "frozen"
    doc["freeze"].update(
        freeze_manifest="analysis/v3_validation/freeze_sidecar.yaml", frozen_at="x"
    )
    doc["pre_freeze_runs"]["inventory_sha256"] = "1" * 64
    doc["frozen_inputs"]["recorded_deviations"].append("DEV-033")
    doc["frozen_inputs"]["contract_sha256"] = "2" * 64
    doc["execution_readiness"]["training_blockers"] = [
        b
        for b in doc["execution_readiness"]["training_blockers"]
        if b["id"] not in {"partitions_metrics", "independent_review_round_2", "independent_review"}
    ]
    return doc


def test_the_g3_diff_check_accepts_exactly_the_freeze_edits(protocol):
    assert g3.protocol_violations(protocol, protocol) == []
    assert g3.protocol_violations(protocol, _g3_edit(protocol)) == []


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.update(status="frozen"),
        lambda d: d["freeze"].update(training_execution_allowed=True),
        lambda d: d["execution_readiness"].update(training_allowed=True),
        lambda d: d["endpoints"][0].update({"class": "changed"}),
        lambda d: d["calibration"]["limit_of_detection"].update(failed_fit_policy="changed"),
        lambda d: d["frozen_inputs"].update(schema_sha256="3" * 64),
        lambda d: d["execution_readiness"].update(
            training_blockers=[
                b
                for b in d["execution_readiness"]["training_blockers"]
                if b["id"] != "run_provenance_binding"
            ]
        ),
        lambda d: d["frozen_inputs"]["recorded_deviations"].pop(0),
    ],
)
def test_the_g3_diff_check_rejects_anything_else(protocol, mutate):
    doc = _g3_edit(protocol)
    mutate(doc)
    assert g3.protocol_violations(protocol, doc)


def test_the_g3_diff_check_rejects_files_outside_the_signing_set():
    assert g3.file_violations(["analysis/v3_validation/protocol.yaml", "PLAN.md"]) == []
    assert g3.file_violations(["scripts/v3_lod.py", "schemas/v3/validation_protocol.schema.json"])
