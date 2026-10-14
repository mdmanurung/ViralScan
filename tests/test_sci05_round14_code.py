"""SCI-05 round 14 (DEV-031): one regression per finding, written to fail on the reviewed tree."""

from __future__ import annotations

import copy
import inspect
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_lod as lod  # noqa: E402
import v3_score as score  # noqa: E402
import validate_v3_protocol as gate  # noqa: E402


@pytest.fixture(scope="module")
def protocol() -> dict:
    return gate.load_yaml(gate.DEFAULT_PROTOCOL)


# ── D6 is not applicable, never failed ────────────────────────────────────────


def test_d6_is_not_applicable_at_g3_and_its_row_is_not_failed(protocol):
    d6 = next(
        d for d in protocol["harmonization"]["denominators"] if d["id"] == "D6_infected_cell_auprc"
    )
    assert d6["count_layer"] == "deferred-no-continuous-score-id-at-G3"
    assert "not_applicable" in d6["missing_output"] and "row-not-failed" in d6["missing_output"]
    assert "row-failed" not in d6["missing_output"].replace("row-not-failed", "")
    assert "D6" in " ".join(protocol["calibration"]["freeze_scope"].split())


# ── recall runs over all planted truth (Q122) ────────────────────────────────


def test_recall_ignores_the_twin_exclusion_but_precision_and_f1_take_it():
    truth = {"c1", "c2", "c3", "c4"}
    out = score.precision_recall_f1(truth, {"c1", "c2"}, twin_excluded={"c1"})
    assert out["recall"] == 2 / 4  # not 1/3: c1 stays in the truth
    assert (out["tp"], out["fp"], out["fn"]) == (1, 0, 2)  # the F1 pool drops c1 on both sides
    assert out["f1"] == pytest.approx(2 / 4) and out["precision"] == 1.0


def test_recall_ignores_the_unknown_origin_rule():
    truth = {("c1", "u1", "A"), ("c2", "u2", "A")}
    out = score.precision_recall_f1(
        truth,
        truth | {("c3", "u3", "B")},
        real_background=True,
        twin_excluded={"c2"},
        cell_of=lambda m: m[0],
    )
    assert out["recall"] == 1.0 and out["n_unknown_origin"] == 1


# ── the production scorer binds to the verified protocol ─────────────────────


def test_score_production_has_no_protocol_override_and_needs_a_run_config():
    params = inspect.signature(score.score_production).parameters
    assert "protocol" not in params and params["run_config"].kind is inspect.Parameter.KEYWORD_ONLY


def test_score_production_refuses_a_run_that_did_not_use_emptydrops(monkeypatch):
    monkeypatch.setattr(score, "require_protocol_run_config", lambda *a, **k: None)
    for caller in ("external", "knee", None):
        with pytest.raises(SystemExit, match="emptydrops"):
            score.score_production("training", [], run_config=SimpleNamespace(cell_calling=caller))


def _row(sample="s1", dataset="synthetic_host_only", **over):
    base = dict(
        dataset_id=dataset,
        biological_sample_id=sample,
        background_library_id=sample,
        twin_sample_id="",
        origin="cell",
        planted_virus_id="HSV1",
    )
    return base | over


def test_a_librarys_dataset_is_read_from_its_truth_rows_not_from_a_flag(protocol):
    lib = dict(
        sample_id="s1", dataset="pbmc_10x_healthy_v3", tier="T1", called={}, truth_rows=[_row()]
    )
    with pytest.raises(ValueError, match="contradicts its truth rows"):
        score.score_by_truth_status([lib], protocol=protocol)
    del lib["dataset"]  # rows alone decide
    assert score.library_dataset(lib) == "synthetic_host_only"


def test_the_negative_only_scorer_refuses_planted_positives_and_says_so(protocol):
    lib = dict(sample_id="s1", tier="T1", called={}, truth_rows=[_row()])
    with pytest.raises(ValueError, match="negative-only scorer"):
        score.score_by_truth_status([lib], protocol=protocol)
    host = _row(sample="s2", origin="unknown", planted_virus_id="")
    empty = dict(
        sample_id="s2", dataset="synthetic_host_only", tier="T1", called={}, truth_rows=[host]
    )
    assert (
        score.score_by_truth_status([empty], protocol=protocol)["domain"]
        == "negative_libraries_only"
    )


def test_real_background_is_read_from_the_manifest_not_passed():
    synthetic = [_row()]
    real = [_row(twin_sample_id="s1_twin")]
    assert not score.real_background_from_manifest(synthetic)
    assert score.real_background_from_manifest(real)
    assert score.real_background_from_manifest([_row(background_library_id="libA")])
    row = {**_row(), "is_planted_molecule": "true", "cell_barcode": "c", "umi": "u"}
    called = {("c", "u", "HSV1"), ("c2", "u2", "HSV1")}
    assert score.score_planted_calls([row], called)["fp"] == 1  # exact truth: an extra call is a FP
    row["twin_sample_id"] = "s1_twin"
    out = score.score_planted_calls([row], called)
    assert (out["fp"], out["n_unknown_origin"]) == (0, 1)  # real background: unknown origin


def test_the_ungated_report_does_not_call_a_presumed_label_a_false_call():
    lib = dict(sample_id="s1", tier="T4", called={"EBV": {("c", "u")}}, truth_rows=[])
    text = score.format_negative_report(score.score_negative_libraries([lib]))
    assert "unexpected call" in text and "false call" not in text.replace("not a false call", "")


# ── --production guards ───────────────────────────────────────────────────────


def test_production_refuses_synthetic_factorial_without_a_real_background(protocol):
    cfg = dict(
        n_cells=2000,
        reads_per_cell=25000,
        n_empty=3000,
        split_seed=None,
        abundance=gt.factor_levels(protocol)["viral_abundance"],
        chemistries=gt.factor_levels(protocol)["chemistry"],
        datasets=sorted(gt.DATASETS),
        covariates={},
        replicates=4,
    )
    assert any(
        "synthetic_factorial needs a real" in e for e in gt.production_design_errors(cfg, protocol)
    )


def test_real_background_anchor_is_a_named_blocker_apart_from_the_simulator(protocol):
    blockers = {b["id"]: b for b in protocol["execution_readiness"]["training_blockers"]}
    assert "real_background_anchor" in blockers and "simulator" in blockers
    text = " ".join(blockers["real_background_anchor"]["description"].split())
    assert "EmptyDrops" in text and "B15" in text and "donor-library list" in text


# ── digest enforcement does not depend on frozen_inputs.status ───────────────


def _schema():
    return json.loads(gate.DEFAULT_SCHEMA.read_text())


def test_frozen_inputs_cannot_be_downgraded_to_pending(protocol):
    doc = copy.deepcopy(protocol)
    doc["frozen_inputs"]["status"] = "pending"
    errors = gate.validate_protocol(doc, _schema(), ledger=gate.load_yaml(gate.DEFAULT_LEDGER))
    assert any("frozen_inputs" in e for e in errors)


def test_the_digest_is_recomputed_and_the_ledger_bound_even_when_pending(protocol):
    doc = copy.deepcopy(protocol)
    doc["frozen_inputs"]["status"] = "pending"
    doc["seeds"]["split"]["value"] = doc["seeds"]["split"]["value"] + 1
    ledger = gate.load_yaml(gate.DEFAULT_LEDGER)
    assert any("does not match the frozen seed" in e for e in gate._validate_sci03(doc))
    doc2 = copy.deepcopy(protocol)
    doc2["frozen_inputs"]["status"] = "pending"
    assert any("cannot be downgraded" in e for e in gate.validate_amendment_ledger(doc2, ledger))
    doc3 = copy.deepcopy(protocol)
    doc3["frozen_inputs"]["status"] = "pending"
    doc3["frozen_inputs"]["recorded_deviations"] = doc3["frozen_inputs"]["recorded_deviations"][:-1]
    assert any(
        "has records not declared" in e for e in gate.validate_amendment_ledger(doc3, ledger)
    )


def test_both_schema_copies_agree_on_the_downgrade_and_frozen_at_rules():
    one = (ROOT / "schemas/v3/validation_protocol.schema.json").read_bytes()
    two = (ROOT / "src/viralscan/schemas/v3/validation_protocol.schema.json").read_bytes()
    assert one == two
    schema = json.loads(one)
    assert schema["$defs"]["frozenInputs"]["properties"]["status"] == {"const": "frozen"}


# ── f1_interval: frozen stratification, units indivisible ────────────────────


def test_f1_interval_requires_strata_and_resamples_within_them():
    units = {f"x{i}": (1, 0, 0) for i in range(3)} | {f"y{i}": (0, 1, 1) for i in range(3)}
    with pytest.raises(TypeError):
        score.f1_interval(units, seed=1, replicates=20)  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="stratum"):
        score.f1_interval(units, seed=1, replicates=20, strata={"x0": "a"})
    strata = {k: k[0] for k in units}
    out = score.f1_interval(units, seed=1, replicates=200, strata=strata)
    # within-stratum resampling keeps three x units and three y units in every replicate, so the
    # pooled counts never move; an unstratified bootstrap would vary
    assert out["n_strata"] == 2 and out["reason_code"] == "constant-bootstrap-distribution"
    free = score.f1_interval(units, seed=1, replicates=200, strata=dict.fromkeys(units, "all"))
    assert free["interval_status"] == "ok"


def test_protocol_states_the_unit_level_stratification(protocol):
    rule = " ".join(protocol["calibration"]["uncertainty"]["stratification_unit_rule"].split())
    assert "indivisible" in rule and "chemistry alone" in rule and "full abundance grid" in rule


# ── BCa / LOD failure policy ──────────────────────────────────────────────────


def test_a_replicate_that_does_not_reach_the_lod_is_reason_coded_not_infinite(monkeypatch):
    rows = [
        (level, int(i < k))
        for level, k in [(10, 2), (100, 5), (1000, 9), (10000, 10)]
        for i in range(10)
    ]
    real = lod._estimate
    calls = {"n": 0}

    def estimate(units, keys):
        calls["n"] += 1
        out = real(units, keys)
        return out if calls["n"] == 1 else dict(out, status="not_reached", lod=None)

    monkeypatch.setattr(lod, "_estimate", estimate)
    result = lod.lod_interval(
        {f"u{i}": rows for i in range(3)},
        seed=1,
        strata=dict.fromkeys(["u0", "u1", "u2"], "all"),
        replicates=5,
    )
    assert result["reason_code"] == "lod-not-reached" and result["failed_fit_statuses"] == {
        "not_reached": 5
    }
    assert lod._comparable({"status": "not_reached", "lod": None}) is None


def test_bca_has_no_silent_zero_acceleration_and_checks_its_bounds():
    constant_jack = lod.bca_interval(0.5, [0.1, 0.4, 0.6, 0.9, 0.3], [0.5, 0.5, 0.5])
    assert constant_jack["reason_code"] == "acceleration-undefined"
    assert (
        lod.bca_interval(0.5, [0.1, math.inf, 0.6], [0.4, 0.5, 0.6])["reason_code"]
        == "nonfinite-statistic"
    )
    assert (
        lod.bca_interval(0.5, [0.1, None, 0.6], [0.4, 0.5, 0.6])["reason_code"]
        == "undefined-statistic"
    )
    ok = lod.bca_interval(0.5, [0.1, 0.4, 0.6, 0.9, 0.3, 0.7], [0.4, 0.5, 0.7, 0.45])
    assert ok["interval_status"] == "ok" and ok["interval"][0] < ok["interval"][1]
    assert all(math.isfinite(b) for b in ok["interval"])


def test_every_bca_reason_code_is_in_the_protocol_vocabulary(protocol):
    policy = " ".join(protocol["calibration"]["uncertainty"]["degenerate_interval_policy"].split())
    for code in (
        "constant-bootstrap-distribution",
        "undefined-statistic",
        "nonfinite-statistic",
        "below-minimum-units",
        "bias-correction-undefined",
        "acceleration-undefined",
        "invalid-interval",
        "lod-not-reached",
    ):
        assert code in policy, code


# ── probit quasi-separation (the agent's probe) ──────────────────────────────


@pytest.mark.parametrize("detected_at_10", [2, 3, 1])
def test_a_tie_of_detected_and_missed_at_the_knee_is_quasi_separation(detected_at_10):
    levels = [1, 3, 10, 30, 100]
    xs = [x for x in levels for _ in range(4)]
    k = {1: 0, 3: 0, 10: detected_at_10, 30: 4, 100: 4}
    ys = [int(i < k[x]) for x in levels for i in range(4)]
    fit = lod.fit_probit(xs, ys)
    assert fit["status"] == "quasi_separation" and fit["slope"] is None
    assert lod.lod95(xs, ys)["lod"] is None


def test_overlapping_data_still_fits_and_failed_fit_policy_is_frozen(protocol):
    fit = lod.fit_probit([1, 1, 10, 10, 100, 100], [0, 1, 0, 1, 1, 1])
    assert fit["status"] == "ok"
    policy = " ".join(protocol["calibration"]["limit_of_detection"]["failed_fit_policy"].split())
    assert "quasi_separation" in policy and "no share of failed replicates" in policy


# ── Q135 S population ─────────────────────────────────────────────────────────


def test_s_population_is_restricted_to_molecules_aligned_to_an_index_gene_of_v(protocol):
    scope = " ".join(protocol["calibration"]["freeze_scope"].split())
    assert "at least one alignment to an index gene of v" in scope
    assert "host-only molecules never fill the denominator" in scope
