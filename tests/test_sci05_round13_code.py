"""SCI-05 round 13, code side: one regression per Major (M1-M5), written to fail on the reviewed tree."""

from __future__ import annotations

import copy
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_partition as part  # noqa: E402
import v3_score as score  # noqa: E402
import validate_v3_protocol as gate  # noqa: E402

from tests._truth_panel_real import write_background  # noqa: E402
from tests.test_generate_truth_panel import salt_file  # noqa: E402
from tests.test_generate_truth_panel_m04 import D1_ARGS  # noqa: E402


@pytest.fixture
def protocol() -> dict:
    return gate.load_yaml(gate.DEFAULT_PROTOCOL)


# ── M1: the run gate is not open at "draft" ───────────────────────────────────


@pytest.mark.parametrize("phase", ["draft", "", "validation"])
def test_m1_every_gate_refuses_a_non_execution_phase(phase):
    with pytest.raises(SystemExit, match="refused phase"):
        gate.require_execution_allowed(phase)
    with pytest.raises(SystemExit, match="refused phase"):
        score.score_production(phase, [], run_config=object())
    with pytest.raises(SystemExit, match="refused phase"):
        gate.require_protocol_run_config(object(), phase)


# ── M2: gating fields are inside frozen_inputs_sha256 ─────────────────────────


def _mutations(doc: dict) -> dict:
    ref = doc["references"][0]
    asset = next(a for a in ref["assets"] if "freeze_required" in a)
    return {
        "reference host_release": lambda d: d["references"][0].__setitem__("host_release", "x"),
        "reference viral_panel": lambda d: d["references"][0].__setitem__("viral_panel", "x"),
        "reference purpose": lambda d: d["references"][0].__setitem__("purpose", "x"),
        "asset freeze_required": lambda d: next(
            a for a in d["references"][0]["assets"] if "freeze_required" in a
        ).__setitem__("freeze_required", not asset["freeze_required"]),
        "planned section status": lambda d: next(
            iter(d["planned_freeze_sections"].values())
        ).__setitem__("status", "tampered"),
        "planned section freeze_required": lambda d: d["planned_freeze_sections"][
            "independent_review"
        ].__setitem__("freeze_required", False),
        "seed freeze_required": lambda d: next(
            v for v in d["seeds"].values() if isinstance(v, dict) and "freeze_required" in v
        ).__setitem__("freeze_required", False),
    }


def test_m2_each_execution_gating_field_moves_the_frozen_inputs_digest(protocol):
    base = gate.frozen_inputs_sha256(protocol)
    for name, mutate in _mutations(protocol).items():
        doc = copy.deepcopy(protocol)
        mutate(doc)
        assert gate.frozen_inputs_sha256(doc) != base, name


def test_m2_asset_pending_reason_stays_outside_the_digest(protocol):
    doc = copy.deepcopy(protocol)
    for ref in doc["references"]:
        for asset in ref["assets"]:
            asset["pending_reason"] = "reworded"
    assert gate.frozen_inputs_sha256(doc) == gate.frozen_inputs_sha256(protocol)


# ── M3: apportion implements the frozen algorithm ─────────────────────────────


def test_m3_apportion_matches_the_frozen_text_worked_examples():
    # quota = size * f, floors, then largest remainders up to ceil(N * f) seats
    assert part.apportion({"a": 4, "b": 4, "c": 3}, 0.3) == {"a": 2, "b": 1, "c": 1}
    seats = part.apportion(
        {f"s{i:02d}": 4 for i in range(42)}, 0.3
    )  # the protocol's production shape
    assert sum(seats.values()) == 51  # ceil(50.4), not round(50.4) = 50
    assert set(seats.values()) == {1, 2}
    # exact decimal arithmetic: 10 * 0.3 is 3 seats, not ceil(3.0000000000000004) = 4
    assert part.apportion({"a": 5, "b": 5}, 0.3) == {"a": 2, "b": 1}


def test_m3_a_single_library_chemistry_is_training():
    assert set(part.split_libraries({"10xv3": ["only"]}, 1, 0.3).values()) == {part.TRAINING}


# ── M4 and M5: --production holds the frozen design ───────────────────────────


def full_cfg(protocol: dict, **over) -> dict:
    levels = gt.factor_levels(protocol)
    cfg = dict(
        n_cells=2000,
        reads_per_cell=25000,
        n_empty=3000,
        split_seed=None,
        abundance=levels["viral_abundance"],
        chemistries=levels["chemistry"],
        datasets=sorted(gt.DATASETS),
        covariates={f: "cycle" for f, _, _ in gt.COVARIATES},
        replicates=4,
    )
    return cfg | over


def test_m5_a_complete_design_passes_and_each_gap_is_refused(protocol):
    # r14 M5: a complete design on a synthetic background is refused (synthetic_factorial needs
    # a real one), and a real background is refused too (blocker real_background_anchor), so no
    # production design passes today; every gap adds its own error on top of the background one.
    base = gt.production_design_errors(full_cfg(protocol), protocol)
    assert len(base) == 1 and "synthetic_factorial needs a real" in base[0]
    real = gt.production_design_errors(full_cfg(protocol, background_r1="x"), protocol)
    assert len(real) == 1 and "real background" in real[0]
    gaps = {
        "abundance": ["10"],
        "chemistries": ["10xv3"],
        "datasets": ["synthetic_factorial"],
        "replicates": 2,
        "covariates": {f: "0.01" for f, _, _ in gt.COVARIATES},
    }
    for key, value in gaps.items():
        errors = gt.production_design_errors(full_cfg(protocol, **{key: value}), protocol)
        assert len(errors) > len(base), key


def test_m5_the_cli_default_abundance_is_not_a_production_design(protocol):
    cfg = gt.parse_args(["--generation-seed", "1", "--out", "x", "--production"])
    assert any("viral_abundance" in e for e in gt.production_design_errors(cfg, protocol))


def test_m4_production_refuses_a_real_background(protocol):
    cfg = full_cfg(protocol, background_r1="lib_R1.fastq.gz")
    assert any("real background" in e for e in gt.production_design_errors(cfg, protocol))


def test_m5_production_generate_runs_the_design_check(tmp_path, monkeypatch, protocol):
    monkeypatch.setattr(gt, "require_execution_allowed", lambda *a, **k: None)
    frozen = copy.deepcopy(protocol)
    for factor in frozen["factors"]:
        factor["status"] = "frozen"
    monkeypatch.setattr(gt, "load_protocol", lambda: frozen)
    cfg = gt.parse_args(
        [
            "--generation-seed",
            "1",
            "--out",
            str(tmp_path / "p"),
            "--production",
            "--replicates",
            "4",
            "--n-cells",
            "2000",
            "--n-empty",
            "3000",
            "--reads-per-cell",
            "25000",
        ]
    )  # sizes are frozen-correct; only the design is short (abundance 10 only)
    with pytest.raises(ValueError, match="viral_abundance"):
        gt.generate(cfg)


def test_m5_donor_library_split_reaches_the_real_background_generator(tmp_path):
    background = write_background(tmp_path / "in")
    base = [*D1_ARGS, "--background-r1", str(background[0]), "--background-r2", str(background[1])]
    base += ["--split-seed", "7"]
    first = gt.generate(gt.parse_args([*base, "--out", str(tmp_path / "plain")]))
    assert {d["partition"] for d in first["samples"].values()} == {part.TRAINING}  # no donor list
    library = first["inputs"]["background"]["library_id"]
    donors = next(
        d
        for d in ([library, f"other_{k}"] for k in range(200))  # 2 libraries: 1 holdout seat
        if part.split_libraries({"10xv3": d}, 7, 0.3)[library] == part.HOLDOUT
    )
    out, hold = tmp_path / "out", tmp_path / "hold"
    manifest = gt.generate(
        gt.parse_args(
            [
                *base,
                "--out",
                str(out),
                "--donor-libraries",
                *donors,
                "--holdout-dir",
                str(hold),
                "--holdout-salt-file",
                str(salt_file(tmp_path / "salt")),
            ]
        )
    )
    assert manifest["samples"] == {} and len(manifest["holdout_samples"]) == len(first["samples"])
    with pytest.raises(ValueError, match="must list this library"):
        gt.generate(
            gt.parse_args([*base, "--out", str(tmp_path / "bad"), "--donor-libraries", "x"])
        )


# ── Biomni round 13: BM1, BM3, BM4, BM5, BM6 ──────────────────────────────────

import math  # noqa: E402

import v3_lod as lod  # noqa: E402


def test_bm1_unknown_origin_calls_are_not_false_positives_on_a_real_background():
    truth = {("c1", "u1", "A"), ("c2", "u2", "A")}
    called = truth | {("c3", "u3", "A"), ("c4", "u4", "B"), ("c5", "u5", "A")}
    synthetic = score.precision_recall_f1(truth, called)
    assert (synthetic["fp"], synthetic["n_unknown_origin"]) == (3, 0)
    real = score.precision_recall_f1(
        truth, called, real_background=True, known_false={("c4", "u4", "B")}
    )
    assert (real["tp"], real["fp"], real["n_unknown_origin"]) == (2, 1, 2)
    assert real["precision"] == pytest.approx(2 / 3) and real["recall"] == 1.0
    assert real["f1"] == pytest.approx(4 / 5)  # 2TP / (2TP + FP + FN)
    none_known = score.precision_recall_f1(truth, called, real_background=True)
    assert none_known["fp"] == 0 and none_known["precision"] == 1.0


def test_bm1_twin_called_cells_leave_the_e3_e4_denominators():
    truth = {("c1", "u1", "A"), ("c2", "u2", "A")}
    called = truth | {("c2", "u9", "A"), ("c3", "u3", "A")}
    out = score.precision_recall_f1(truth, called, twin_excluded={"c2"}, cell_of=lambda m: m[0])
    assert (out["tp"], out["fp"], out["fn"], out["n_twin_excluded"]) == (1, 1, 0, 1)
    cells = score.precision_recall_f1({"c1", "c2"}, {"c1", "c2", "c3"}, twin_excluded={"c2"})
    assert (cells["tp"], cells["fp"], cells["recall"]) == (1, 1, 1.0)


def test_bm3_f1_is_the_count_form_and_precision_stays_undefined_separately():
    nothing_called = score.precision_recall_f1({"a", "b"}, set())
    assert nothing_called["precision"] is None and nothing_called["f1"] == 0.0
    assert score.precision_recall_f1(set(), set())["f1"] is None
    assert score.f1_from_counts(3, 1, 2) == pytest.approx(6 / 9)


def test_bm3_the_bootstrap_uses_the_same_statistic():
    units = {"a": (8, 2, 0), "b": (5, 0, 5), "c": (6, 1, 1), "d": (9, 0, 1)}
    out = score.f1_interval(units, seed=3, replicates=300, strata=dict.fromkeys(units, "x"))
    assert out["point"] == pytest.approx(score.f1_from_counts(28, 3, 7))
    assert (
        out["interval_status"] == "ok" and out["interval"][0] <= out["point"] <= out["interval"][1]
    )
    assert (
        score.f1_interval({"a": (1, 0, 0)}, seed=3, strata={"a": "x"})["reason_code"]
        == "below-minimum-units"
    )


def test_bm4_constant_undefined_and_nonfinite_results_have_no_interval():
    # all-zero false calls: every replicate and jackknife value is 0
    zero = lod.bca_interval(0.0, [0.0] * 200, [0.0] * 5)
    assert zero["interval"] is None and zero["reason_code"] == "constant-bootstrap-distribution"
    assert (
        lod.bca_interval(math.nan, [0.1, 0.2], [0.1, 0.2])["reason_code"] == "nonfinite-statistic"
    )
    assert (
        lod.bca_interval(0.1, [0.1, math.nan], [0.1, 0.2])["reason_code"] == "nonfinite-statistic"
    )
    assert (
        lod.bca_interval(0.1, [0.1, 0.2], [math.inf, 0.2])["reason_code"] == "nonfinite-statistic"
    )
    # a degenerate F1 (nothing ever called or planted) is non-finite, not a zero-width interval
    abc = dict.fromkeys("abc", "x")
    degenerate = score.f1_interval({k: (0, 0, 0) for k in "abc"}, seed=1, replicates=50, strata=abc)
    assert degenerate["point"] is None and degenerate["reason_code"] == "undefined-statistic"
    constant = score.f1_interval({k: (4, 1, 1) for k in "abc"}, seed=1, replicates=50, strata=abc)
    assert constant["reason_code"] == "constant-bootstrap-distribution"
    assert constant["n_units"] == 3
    # the LOD interval has the same guard: identical units give an identical LOD every replicate
    rows = _lod_rows()
    result = lod.lod_interval(
        {f"u{i}": rows for i in range(3)},
        seed=1,
        strata=dict.fromkeys(["u0", "u1", "u2"], "all"),
        replicates=40,
    )
    assert result["interval"] is None and result["reason_code"] == "constant-bootstrap-distribution"


def _lod_rows():
    """One unit whose pooled probit fit is ok with LOD95 about 2,250 inside [10, 10000]."""
    return [
        (level, int(i < k))
        for level, k in [(10, 2), (100, 5), (1000, 9), (10000, 10)]
        for i in range(10)
    ]


def _library(sample_id, dataset, called):
    host = dict(
        dataset_id=dataset,
        biological_sample_id=sample_id,
        origin="unknown",
        planted_virus_id="",
        cell_barcode="h",
        umi="h",
    )  # a manifest-backed library has rows (every read is recorded, host included)
    return dict(sample_id=sample_id, dataset=dataset, tier="T1", called=called, truth_rows=[host])


def test_bm5_truth_status_comes_from_the_frozen_registry(protocol):
    status = score.truth_status_by_dataset(protocol)
    assert status["pbmc_10x_healthy_v3"] == "presumed-negative"
    assert status["synthetic_host_only"] == "exact"
    out = score.score_by_truth_status(
        [
            _library("s1", "synthetic_host_only", {}),
            _library("s2", "synthetic_host_only", {"HSV1": {("c", "u")}}),
            _library("p1", "pbmc_10x_healthy_v3", {"EBV": {("c", "u")}}),
            _library("p2", "pbmc_10x_healthy_v3", {}),
        ],
        protocol=protocol,
    )
    assert out["exact"]["T1"]["unconditional"]["false_calls"] == 1  # only the exact libraries
    report = out["descriptive"]["pbmc_10x_healthy_v3"]
    assert (
        report["report"] == "descriptive_agreement"
        and report["truth_status"] == "presumed-negative"
    )
    assert (report["n_libraries"], report["n_with_unexpected_calls"]) == (2, 1)
    assert not {"rate", "false_calls", "specificity"} & set(report)
    with pytest.raises(ValueError, match="not in the frozen registry"):
        score.score_by_truth_status([_library("x", "made_up", {})], protocol=protocol)


def test_bm5_the_gated_entry_routes_through_the_registry(monkeypatch, protocol):
    monkeypatch.setattr(score, "require_protocol_run_config", lambda *a, **k: None)
    out = score.score_production(
        "training",
        [_library("p", "pbmc_10x_healthy_v3", {})],
        run_config=SimpleNamespace(cell_calling="emptydrops"),
    )
    assert out["exact"] == {} and "pbmc_10x_healthy_v3" in out["descriptive"]


def test_bm6_a_failed_optimizer_is_reported_not_used(monkeypatch):
    from scipy import optimize

    real = optimize.minimize

    def failing(*args, **kwargs):
        fit = real(*args, **kwargs)
        fit.success = False
        fit.jac = fit.jac + 1.0  # not stationary
        return fit

    monkeypatch.setattr(lod.optimize, "minimize", failing)
    fit = lod.fit_probit([10, 10, 100, 100, 1000, 1000], [0, 1, 0, 1, 1, 1])
    assert fit["status"] == "not_converged"
    assert lod.lod95([10, 10, 100, 100, 1000, 1000], [0, 1, 0, 1, 1, 1])["lod"] is None


def test_bm6_failed_replicate_fits_keep_their_reason_in_the_interval(monkeypatch):
    real = lod._estimate
    calls = {"n": 0}

    def estimate(units, keys):
        calls["n"] += 1
        out = real(units, keys)
        return out if calls["n"] == 1 else dict(out, status="not_converged", lod=None)

    monkeypatch.setattr(lod, "_estimate", estimate)
    rows = _lod_rows()
    result = lod.lod_interval(
        {f"u{i}": rows for i in range(3)},
        seed=1,
        strata=dict.fromkeys(["u0", "u1", "u2"], "all"),
        replicates=5,
    )
    assert result["interval_status"] == "not_estimable" and result["reason_code"] == "failed-fit"
    assert result["failed_fit_statuses"] == {"not_converged": 5}


# ── minor: the freeze sidecar is opened and its pins re-hashed ───────────────


def test_minor_freeze_sidecar_is_opened_and_rehashed(tmp_path, protocol):
    import hashlib

    import yaml

    proto = tmp_path / "protocol.yaml"
    proto.write_text("a: 1\n")
    results, thresholds = tmp_path / "results.json", tmp_path / "thresholds.json"
    results.write_text("r")
    thresholds.write_text("t")
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    sidecar = tmp_path / "side.yaml"
    body = {
        "protocol_file_sha256": sha(proto),
        "artifacts": {
            "training_results": {"path": "results.json", "sha256": sha(results)},
            "thresholds": {"path": "thresholds.json", "sha256": sha(thresholds)},
        },
    }
    sidecar.write_text(yaml.safe_dump(body))
    doc = {
        "freeze": {
            "freeze_manifest": "side.yaml",
            "training_results_sha256": sha(results),
            "thresholds_sha256": sha(thresholds),
        }
    }

    def check(phase="holdout"):
        return gate.freeze_sidecar_errors(doc, proto, phase, repo_root=tmp_path)

    assert check() == [] and gate.freeze_sidecar_errors({"freeze": {}}, proto, "holdout") == []
    results.write_text("tampered")
    assert any("training_results" in e for e in check())
    assert check("training") == []  # artefacts exist only once training is done
    results.write_text("r")
    proto.write_text("a: 2\n")
    assert any("protocol_file_sha256" in e for e in check("training"))
    sidecar.unlink()
    assert any("does not exist" in e for e in check("training"))


def test_minor_m2_d25_is_checked_against_its_frozen_contract(protocol):
    doc = copy.deepcopy(protocol)
    for item in doc["harmonization"]["denominators"]:
        if item["id"] == "D25_absent_sibling_false_calls":
            item["denominator"] = "a narrower population"
    errors = gate.validate_protocol(
        doc, gate.load_yaml(gate.DEFAULT_SCHEMA), ledger=gate.load_yaml(gate.DEFAULT_LEDGER)
    )
    assert any("D25_absent_sibling_false_calls" in e and "frozen denominator" in e for e in errors)


def test_minor_m5_an_undeclared_appended_ledger_record_is_refused(protocol):
    ledger = gate.load_yaml(gate.DEFAULT_LEDGER)
    last = ledger["deviations"][-1]
    extra = {**last, "deviation_id": "DEV-999"}
    extra["record_sha256"] = gate.ledger_record_sha256(extra, last["record_sha256"])
    ledger["deviations"].append(extra)
    errors = gate.validate_amendment_ledger(protocol, ledger)
    assert any("not declared" in e and "DEV-999" in e for e in errors)
