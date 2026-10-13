"""SCI-05 round 12: does the code honour what the protocol froze? (M2, M3, M4, m5).

Reproductions of the reviewer's findings, written so each fails on the pre-fix tree:
the shipped EmptyDrops seed (100) is not ``seeds.cell_calling``; frozen keys had no consumer;
``require_execution_allowed`` had no call site; ``--production`` ignored the frozen panel sizes;
``knee`` stayed selectable.
"""

from __future__ import annotations

import copy
import importlib
import inspect
import re
import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import v3_score as score  # noqa: E402
import validate_v3_protocol as gate  # noqa: E402

from viralscan.defaults import DEFAULTS  # noqa: E402
from viralscan.runconfig import RunConfig  # noqa: E402

BLOCK = ("harmonization", "cell_universe", "host_only_emptydrops")


@pytest.fixture
def protocol() -> dict:
    return gate.load_yaml(gate.DEFAULT_PROTOCOL)


def frozen_config(protocol: dict, **over) -> RunConfig:
    block = protocol["harmonization"]["cell_universe"]["host_only_emptydrops"]
    fields = {
        "cell_calling": "emptydrops",
        "emptydrops_seed": protocol["seeds"]["cell_calling"]["value"],
        "emptydrops_fdr": block["fdr"],
        "emptydrops_lower": block["lower"],
        "emptydrops_niters": block["niters"],
    }
    return replace(RunConfig(), **{**fields, **over})


# ── M2: the EmptyDrops seed ──────────────────────────────────────────────────


def test_m2_shipped_default_seed_is_not_the_frozen_seed_and_is_refused(protocol):
    """The reviewer's reproduction: ``RunConfig().emptydrops_seed`` is 100, the anchor is 20260727002."""
    assert RunConfig().emptydrops_seed == 100  # the non-protocol CLI default is untouched
    cfg = replace(RunConfig(), cell_calling="emptydrops")
    errors = gate.run_config_errors(cfg, protocol)
    assert any("emptydrops_seed=100" in e and "20260727002" in e for e in errors)


def test_m2_frozen_values_pass(protocol):
    assert gate.run_config_errors(frozen_config(protocol), protocol) == []


@pytest.mark.parametrize(
    ("field", "value"),
    [("emptydrops_fdr", 0.05), ("emptydrops_lower", 200), ("emptydrops_niters", 5000)],
)
def test_m2_every_emptydrops_parameter_is_held(protocol, field, value):
    errors = gate.run_config_errors(frozen_config(protocol, **{field: value}), protocol)
    assert len(errors) == 1 and field in errors[0]


def test_m2_shipped_defaults_for_the_other_emptydrops_values_equal_the_protocol(protocol):
    block = protocol["harmonization"]["cell_universe"]["host_only_emptydrops"]
    assert DEFAULTS["emptydrops_fdr"] == block["fdr"]
    assert DEFAULTS["emptydrops_lower"] == block["lower"]
    assert DEFAULTS["emptydrops_niters"] == block["niters"]


def test_m2_seed_source_is_followed_not_hard_coded(protocol):
    moved = copy.deepcopy(protocol)
    moved["harmonization"]["cell_universe"]["host_only_emptydrops"]["seed_source"] = "seeds.split"
    errors = gate.run_config_errors(frozen_config(protocol), moved)
    assert any("20260727001" in e for e in errors)


def test_m2_a_pending_seed_cannot_be_used(protocol):
    with pytest.raises(ValueError, match="not a frozen integer seed"):
        gate.frozen_seed(protocol, "generation")


def test_m2_guard_runs_the_execution_gate_first(protocol):
    # The shipped protocol is a draft, so the gate itself refuses before any config is looked at.
    with pytest.raises(SystemExit, match="protocol execution gate failed"):
        gate.require_protocol_run_config(frozen_config(protocol), "training")


def test_m2_guard_refuses_a_default_config_once_the_gate_opens(monkeypatch):
    monkeypatch.setattr(gate, "require_execution_allowed", lambda *a, **k: None)
    with pytest.raises(SystemExit, match="emptydrops_seed=100"):
        gate.require_protocol_run_config(
            replace(RunConfig(), cell_calling="emptydrops"), "training"
        )


def test_m2_guard_accepts_the_frozen_config_once_the_gate_opens(monkeypatch, protocol):
    monkeypatch.setattr(gate, "require_execution_allowed", lambda *a, **k: None)
    gate.require_protocol_run_config(frozen_config(protocol), "training")


# ── m5: knee (and the silent fallbacks) are not protocol callers ─────────────


@pytest.mark.parametrize("caller", ["knee", "none", "auto"])
def test_m5_only_emptydrops_or_external_is_allowed(protocol, caller):
    errors = gate.run_config_errors(frozen_config(protocol, cell_calling=caller), protocol)
    assert any(f"cell_calling={caller!r}" in e for e in errors)


def test_m5_external_is_allowed_and_needs_no_seed(protocol):
    cfg = replace(RunConfig(), cell_calling="external")
    assert gate.run_config_errors(cfg, protocol) == []


def test_m5_the_cli_still_offers_knee_for_non_protocol_use():
    # PLAN SW-23 / DEF-00 A7: the option stays, warned, outside protocol runs.
    from viralscan.defaults import CELL_CALLING_METHODS

    assert "knee" in CELL_CALLING_METHODS


# ── M3: every frozen key has a consumer or an explicit not_yet_consumed ──────


def test_m3_registry_covers_every_frozen_numeric_value(protocol):
    assert gate.frozen_key_coverage_errors(protocol) == []


def test_m3_a_new_frozen_value_without_an_entry_fails(protocol):
    tampered = copy.deepcopy(protocol)
    tampered["seeds"]["new_stage"] = {"status": "frozen", "value": 7}
    tampered["partitions"]["generated_sample_structure"]["new_count"] = 3
    tampered["calibration"]["new_numeric"] = 1
    errors = gate.frozen_key_coverage_errors(tampered)
    assert any("seeds.new_stage.value" in e for e in errors)
    assert any("calibration.new_numeric" in e for e in errors)
    # a child of a registered prefix is covered by that prefix
    assert not any("generated_sample_structure.new_count" in e for e in errors)


def test_m3_a_stale_registry_key_fails(protocol, monkeypatch):
    monkeypatch.setitem(gate.FROZEN_KEY_CONSUMERS, "seeds.gone", "consumed: v3_lod.REPLICATES")
    assert any(
        "seeds.gone is not in the protocol" in e for e in gate.frozen_key_coverage_errors(protocol)
    )


def test_m3_a_malformed_entry_fails(protocol, monkeypatch):
    monkeypatch.setitem(gate.FROZEN_KEY_CONSUMERS, "seeds.root", "someday")
    assert any("malformed" in e for e in gate.frozen_key_coverage_errors(protocol))


def _resolve(target: str):
    module, _, attr = target.rpartition(".")
    return getattr(importlib.import_module(module), attr)


def test_m3_every_named_consumer_exists_and_reads_the_key(protocol):
    for key, entry in gate.FROZEN_KEY_CONSUMERS.items():
        kind, _, target = entry.partition(": ")
        if kind == "not_yet_consumed":
            continue
        obj = _resolve(target.split(" ")[0])
        if kind == "consumed":
            leaf = key.rsplit(".", 1)[-1]
            assert leaf in inspect.getsource(obj), f"{target} does not mention {leaf}"
        else:  # pinned: a constant that must equal the protocol value
            node = protocol
            for part in key.split("."):
                node = node[part]
            assert obj == node, f"{target}={obj!r} differs from protocol {key}={node!r}"


def test_m3_lod_constants_equal_the_protocol(protocol):
    import v3_lod

    uncertainty = protocol["calibration"]["uncertainty"]
    assert uncertainty["replicates"] == v3_lod.REPLICATES
    assert uncertainty["minimum_samples_for_interval"] == v3_lod.MINIMUM_UNITS


# ── M4: the gate has call sites, and --production holds the frozen panel sizes ─


SCRIPTS = ROOT / "scripts"
#: scripts that touch the v3 protocol, the scorer or the generator
PRODUCER = re.compile(
    r"import v3_score|generate_truth_panel|import v3_lod|--production|v3_validation|score\.PROTOCOL"
)
#: ungated on purpose: file -> why
EXEMPT = {
    "check_ledger_append_only.py": "read-only ledger check",
    "gen_workflow_comparison.py": "renders protocol text into docs; runs nothing",
    "val08_golden.py": "VAL-08 exact golden run; its generator command never passes --production",
    "validate_v3_protocol.py": "defines the guard",
}
GATE_CALLS = ("require_execution_allowed(", "require_protocol_run_config(")


def test_m4_every_protocol_producer_calls_the_gate_or_is_exempt():
    seen = {}
    for path in sorted(SCRIPTS.glob("*.py")):
        text = path.read_text()
        if PRODUCER.search(text):
            seen[path.name] = text
    gated = {n for n, t in seen.items() if any(c in t for c in GATE_CALLS)}
    assert {"generate_truth_panel.py", "v3_score.py"} <= gated
    ungated = set(seen) - gated
    assert ungated <= set(EXEMPT) | {"validate_v3_protocol.py"}, (
        f"producers with no gate call and no exemption: {sorted(ungated - set(EXEMPT))}"
    )
    assert set(EXEMPT) <= set(seen), "stale exemption"
    for name in set(EXEMPT) - {"validate_v3_protocol.py"}:
        assert "--production" not in seen[name], f"{name} must not run a production panel ungated"


def test_m4_production_runs_the_full_gate_before_anything_else(tmp_path):
    # The draft protocol fails the gate for status/freeze reasons, which the old two-line guard never saw.
    with pytest.raises(ValueError, match=r"production run refused: protocol execution gate failed"):
        gt.generate(
            gt.parse_args(["--generation-seed", "1", "--out", str(tmp_path / "p"), "--production"])
        )
    assert not (tmp_path / "p").exists()


def _production_cfg(tmp_path, *extra):
    return gt.parse_args(
        ["--generation-seed", "1", "--out", str(tmp_path / "p"), "--production", *extra]
    )


def test_m4_production_refuses_the_cli_default_panel_sizes(tmp_path, monkeypatch, protocol):
    monkeypatch.setattr(gt, "require_execution_allowed", lambda *a, **k: None)
    frozen = copy.deepcopy(protocol)
    for factor in frozen["factors"]:
        factor["status"] = "frozen"
    monkeypatch.setattr(gt, "load_protocol", lambda: frozen)
    with pytest.raises(ValueError) as exc:
        gt.generate(_production_cfg(tmp_path, "--replicates", "4"))
    message = str(exc.value)
    assert "n_cells=50 differs from the frozen cells_per_sample=2000" in message
    assert "reads_per_cell=200 differs from the frozen reads_per_cell=25000" in message


def test_m4_structure_sizes_are_read_from_the_protocol_not_hard_coded(protocol):
    cfg = {"n_cells": 2000, "reads_per_cell": 25000, "n_empty": 3000, "split_seed": None}
    cfg["background_r1"] = "lib_R1.fastq.gz"
    assert gt.production_structure_errors(cfg, protocol) == []
    changed = copy.deepcopy(protocol)
    changed["partitions"]["generated_sample_structure"]["cells_per_sample"] = 1234
    assert any("cells_per_sample=1234" in e for e in gt.production_structure_errors(cfg, changed))


def test_m4_synthetic_background_fails_closed_while_empty_barcodes_are_prose_only(protocol):
    # DEV-027 added the structured fields to the shipped protocol, so strip them to
    # reproduce the prose-only state this guard must refuse.
    prose_only = copy.deepcopy(protocol)
    structure = prose_only["partitions"]["generated_sample_structure"]
    structure.pop("empty_barcodes_per_sample", None)
    structure.pop("empty_median_umi", None)
    cfg = {"n_cells": 2000, "reads_per_cell": 25000, "n_empty": 3000, "split_seed": None}
    errors = gt.production_structure_errors(cfg, prose_only)
    assert errors and "empty_barcodes_per_sample" in errors[0]


def test_m4_structured_empty_barcodes_are_enforced_when_the_protocol_has_them(protocol):
    structured = copy.deepcopy(protocol)
    structured["partitions"]["generated_sample_structure"]["empty_barcodes_per_sample"] = 3000
    cfg = {"n_cells": 2000, "reads_per_cell": 25000, "n_empty": 20, "split_seed": None}
    assert any("n_empty=20" in e for e in gt.production_structure_errors(cfg, structured))
    cfg["n_empty"] = 3000
    assert gt.production_structure_errors(cfg, structured) == []


def test_m4_production_holds_the_frozen_split_seed(protocol):
    cfg = {"n_cells": 2000, "reads_per_cell": 25000, "n_empty": 3000, "split_seed": 5}
    cfg["background_r1"] = "x"
    assert any("split_seed=5" in e for e in gt.production_structure_errors(cfg, protocol))


def test_m4_scorer_entry_is_gated():
    with pytest.raises(SystemExit, match="protocol execution gate failed"):
        score.score_production("training", [])


def test_m4_scorer_entry_scores_when_the_gate_opens(monkeypatch):
    monkeypatch.setattr(score, "require_execution_allowed", lambda *a, **k: None)
    assert score.score_production("training", []) == {}
