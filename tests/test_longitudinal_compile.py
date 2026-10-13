"""scripts/longitudinal_compile.py (LONG-01 compiler, LONG-02 contamination demotion): synthetic fixtures."""

import csv
import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "longitudinal_compile", Path(__file__).parent.parent / "scripts" / "longitudinal_compile.py"
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)

CFG = {**mod.DEFAULTS, "burden_fold_change": 2.0}


def lib(lid, subject, tp, umis=1_000_000, run="RUN1", lane="RUN1.1", **kw):
    return {
        "library_id": lid,
        "subject_id": subject,
        "timepoint": tp,
        "host_umis": str(umis),
        "n_cells": "1000",
        "run_id": run,
        "lane": lane,
        "timepoint_order": str(tp).lstrip("T"),
        **{k: str(v) for k, v in kw.items()},
    }


def cand(lid, mol, level="L2", virus="taxid:10376", flags="", n_loci=3, thr="abc", **kw):
    return {
        "library_id": lid,
        "virus_key": virus,
        "level": level,
        "flags": flags,
        "unique_molecules": str(mol),
        "n_loci": str(n_loci),
        "thresholds_sha256": thr,
        **{k: str(v) for k, v in kw.items()},
    }


def run(cands, libs, loci=None, cfg=None):
    return mod.compile_longitudinal(cands, libs, loci, cfg or CFG)


def events(cands, libs, loci=None, cfg=None):
    return {
        (e["subject_id"], e["from_timepoint"], e["to_timepoint"]): e
        for e in run(cands, libs, loci, cfg)[1]
    }


def test_new_detection_requires_testable_earlier_visit():
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2")]
    ev = events([cand("a2", 50)], libs)[("S1", "T1", "T2")]
    assert ev["event"] == "new_detection" and ev["to_level"] == "L2"
    # earlier library without host depth cannot show absence
    libs[0]["host_umis"] = ""
    ev = events([cand("a2", 50)], libs)[("S1", "T1", "T2")]
    assert ev["event"] == "not_testable" and "no_host_depth" in ev["event_reason"]


def test_cleared_needs_testable_later_visit():
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2", umis=100)]
    c = [cand("a1", 50)]
    assert events(c, libs)[("S1", "T1", "T2")]["event"] == "cleared"
    cfg = {**CFG, "min_host_umis": 10_000}
    ev = events(c, libs, cfg=cfg)[("S1", "T1", "T2")]
    assert ev["event"] == "not_testable" and "below_min_host_umis" in ev["event_reason"]


def test_burden_is_depth_normalised_and_threshold_comes_from_config():
    # raw counts equal, host depth halves -> burden doubles
    libs = [lib("a1", "S1", "T1", umis=2_000_000), lib("a2", "S1", "T2", umis=1_000_000)]
    c = [cand("a1", 100), cand("a2", 100)]
    ev = events(c, libs)[("S1", "T1", "T2")]
    assert ev["event"] == "burden_increase" and ev["fold_change"] == pytest.approx(2.0)
    assert (
        events(c, libs, cfg={**CFG, "burden_fold_change": 3.0})[("S1", "T1", "T2")]["event"]
        == "persistent"
    )
    # no baked-in number: unset threshold -> not_testable, not a guess
    ev = events(c, libs, cfg={**CFG, "burden_fold_change": None})[("S1", "T1", "T2")]
    assert ev["event"] == "not_testable" and ev["event_reason"] == "burden_fold_change_unset"


def test_burden_decrease_and_persistent():
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2"), lib("a3", "S1", "T3")]
    ev = events([cand("a1", 100), cand("a2", 40), cand("a3", 50)], libs)
    assert ev[("S1", "T1", "T2")]["event"] == "burden_decrease"
    assert ev[("S1", "T2", "T3")]["event"] == "persistent"


def test_visits_are_ordered_by_timepoint_not_row_order_and_libraries_pool():
    libs = [lib("b", "S1", "T10"), lib("a", "S1", "T2"), lib("a_frac", "S1", "T2")]
    c = [cand("b", 90), cand("a", 20), cand("a_frac", 20)]
    traj, evs, _ = run(c, libs)
    assert [t["timepoint"] for t in traj] == ["T2", "T10"]  # numeric, not lexicographic
    assert traj[0]["unique_molecules"] == 40 and traj[0]["host_depth"] == 2_000_000
    assert evs[0]["from_timepoint"] == "T2" and evs[0]["to_timepoint"] == "T10"


def test_unorderable_subject_is_skipped_and_pooled_subjects_excluded():
    libs = [
        lib("a", "S1", "pre", timepoint_order=""),
        lib("b", "S1", "post", timepoint_order=""),
        lib("p", "S2;S3", "T1"),
        lib("u", "UNKNOWN", "T1"),
    ]
    traj, evs, man = run([cand("b", 5)], libs)
    assert traj == [] and evs == []
    assert man["skipped_subjects"] == [{"subject_id": "S1", "reason": "timepoint_order_unknown"}]
    assert man["excluded_pooled_or_unknown_subject"] == ["p", "u"]


def test_threshold_digest_mismatch_and_depth_basis_mismatch_are_not_testable():
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2")]
    ev = events([cand("a1", 50, thr="x"), cand("a2", 50, thr="y")], libs)[("S1", "T1", "T2")]
    assert ev["event"] == "not_testable" and "thresholds_sha256" in ev["event_reason"]
    libs[1]["host_umis"], libs[1]["host_reads"] = "", "1000000"
    ev = events([cand("a1", 50), cand("a2", 50)], libs)[("S1", "T1", "T2")]
    assert ev["event"] == "not_testable" and "depth_basis" in ev["event_reason"]


def test_min_detected_level_and_flag_gating():
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2")]
    # an L0 row is a candidate, not a detection, at the default min level
    assert events([cand("a2", 5, level="L0")], libs) == {}
    cfg = {**CFG, "min_detected_level": "L0"}
    assert (
        events([cand("a2", 5, level="L0")], libs, cfg=cfg)[("S1", "T1", "T2")]["event"]
        == "new_detection"
    )
    # a library carrying a configured not-testable flag cannot support an absence call
    ev = events([cand("a1", 5), cand("a2", 5, flags="not_testable")], libs)[("S1", "T1", "T2")]
    assert ev["event"] == "not_testable" and ev["event_reason"].startswith(
        "later_visit_not_testable"
    )


@pytest.mark.parametrize(
    "level,flags,eligible,block",
    [
        ("L2", "", True, ""),
        ("L3", "", True, ""),
        ("L1", "", False, "level_below_L2"),
        ("L2", "germline_integration_compatible", False, "flag:germline_integration_compatible"),
        ("L3", "screening_only;other", False, "flag:screening_only"),
    ],
)
def test_l4_eligibility(level, flags, eligible, block):
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2")]
    loci = [{"library_id": "a2", "virus_key": "taxid:10376", "molecule_key": "m1"}]
    ev = events([cand("a2", 50, level=level, flags=flags)], libs, loci)[("S1", "T1", "T2")]
    assert ev["event"] == "new_detection"
    assert ev["l4_eligible"] is eligible and ev["l4_block_reason"] == block


def test_l4_is_never_granted_to_persistent_or_unchecked_contamination():
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2")]
    ev = events([cand("a1", 50), cand("a2", 50)], libs)[("S1", "T1", "T2")]
    assert ev["event"] == "persistent" and ev["l4_eligible"] is False
    # no key table: the contamination rule could not run, so the emergence cannot be L4-eligible
    ev = events([cand("a2", 50)], libs)[("S1", "T1", "T2")]
    assert (
        ev["event"] == "new_detection" and ev["contamination_check"] == "not_testable:no_key_table"
    )
    assert ev["l4_eligible"] is False and "contamination_check" in ev["l4_block_reason"]


def _two_subjects(other_run="RUN1", other_lane="RUN1.1"):
    libs = [
        lib("a1", "S1", "T1"),
        lib("a2", "S1", "T2"),
        lib("b1", "S2", "T1", run=other_run, lane=other_lane),
    ]
    c = [cand("a2", 50), cand("b1", 80)]
    return libs, c


def _loci(a2_keys, b1_keys, kind="molecule_key"):
    return [
        {"library_id": lid, "virus_key": "taxid:10376", kind: k}
        for lid, ks in (("a2", a2_keys), ("b1", b1_keys))
        for k in ks
    ]


def test_contamination_demotes_new_detection_on_shared_molecules_same_run():
    libs, c = _two_subjects()
    ev = events(c, libs, _loci(["m1", "m2"], ["m2", "m3"]))[("S1", "T1", "T2")]
    assert ev["event"] == "demoted_contamination" and ev["original_event"] == "new_detection"
    assert (
        ev["contamination_check"] == "demoted"
        and "a2~b1(S2):molecule_keyx1" in ev["contamination_evidence"]
    )
    assert ev["l4_eligible"] is False


def test_contamination_demotes_burden_increase():
    libs, _ = _two_subjects()
    c = [cand("a1", 10), cand("a2", 100), cand("b1", 80)]
    loci = _loci(["m1"], ["m1"]) + [
        {"library_id": "a1", "virus_key": "taxid:10376", "molecule_key": "z"}
    ]
    ev = events(c, libs, loci)[("S1", "T1", "T2")]
    assert ev["event"] == "demoted_contamination" and ev["original_event"] == "burden_increase"


def test_contamination_same_lane_different_run_or_no_overlap_is_clear():
    libs, c = _two_subjects(other_run="RUN2", other_lane="RUN2.1")
    assert events(c, libs, _loci(["m1"], ["m1"]))[("S1", "T1", "T2")]["event"] == "new_detection"
    libs, c = _two_subjects()
    ev = events(c, libs, _loci(["m1"], ["m9"]))[("S1", "T1", "T2")]
    assert (
        ev["event"] == "new_detection"
        and ev["contamination_check"] == "clear"
        and ev["l4_eligible"] is True
    )


def test_contamination_lane_match_across_unknown_run_and_same_subject_ignored():
    libs, c = _two_subjects(other_run="UNKNOWN", other_lane="RUN1.1")
    assert (
        events(c, libs, _loci(["m1"], ["m1"]))[("S1", "T1", "T2")]["event"]
        == "demoted_contamination"
    )
    # same subject, same run: not another subject, so not contamination
    libs = [lib("a1", "S1", "T1"), lib("a2", "S1", "T2"), lib("a3", "S1", "T3")]
    loci = [
        {"library_id": lid, "virus_key": "taxid:10376", "molecule_key": "m1"}
        for lid in ("a2", "a3")
    ]
    ev = events([cand("a2", 50), cand("a3", 50)], libs, loci)[("S1", "T1", "T2")]
    assert ev["event"] == "new_detection"


def test_contamination_match_type_is_configurable_and_other_virus_does_not_count():
    libs, c = _two_subjects()
    loci = _loci(["g1"], ["g1"], kind="locus_id")
    assert events(c, libs, loci)[("S1", "T1", "T2")]["event"] == "demoted_contamination"
    cfg = {**CFG, "contamination_match": ["molecule_key"]}
    assert events(c, libs, loci, cfg)[("S1", "T1", "T2")]["event"] == "new_detection"
    other = [{**r, "virus_key": "taxid:1"} if r["library_id"] == "b1" else r for r in loci]
    assert events(c, libs, other)[("S1", "T1", "T2")]["contamination_check"] == "clear"


def test_contamination_not_testable_without_keys_for_the_library():
    libs, c = _two_subjects()
    loci = [{"library_id": "b1", "virus_key": "taxid:10376", "molecule_key": "m1"}]
    ev = events(c, libs, loci)[("S1", "T1", "T2")]
    assert (
        ev["event"] == "new_detection"
        and ev["contamination_check"] == "not_testable:no_keys_for_library"
    )
    assert ev["l4_eligible"] is False


def test_deterministic_and_input_order_independent():
    libs, c = _two_subjects()
    loci = _loci(["m1"], ["m1"])
    a = run(c, libs, loci)
    b = run(list(reversed(c)), list(reversed(libs)), list(reversed(loci)))
    assert a == b


def test_config_validation(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("burden_fold_chang: 2\n")
    with pytest.raises(SystemExit):
        mod.load_config(bad)
    bad.write_text("burden_fold_change: 0.5\n")
    with pytest.raises(SystemExit):
        mod.load_config(bad)
    assert mod.load_config(None)["burden_fold_change"] is None


def test_cli_writes_tables_and_manifest_and_tolerates_extra_columns(tmp_path):
    libs, c = _two_subjects()
    c = [
        {**r, "level_reason": "x", "read_virus_best_fraction": "0.9", "extra_future_column": "z"}
        for r in c
    ]
    loci = _loci(["m1"], ["m9"])

    def put(name, rows):
        p = tmp_path / name
        cols = list(dict.fromkeys(k for r in rows for k in r))
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(rows)
        return str(p)

    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("burden_fold_change: 2.0\n")
    out = tmp_path / "out"
    mod.main(
        [
            "--candidates",
            put("c.tsv", c),
            "--libraries",
            put("l.tsv", libs),
            "--loci",
            put("k.tsv", loci),
            "--config",
            str(cfg),
            "--out-dir",
            str(out),
        ]
    )
    rows = list(csv.DictReader((out / "longitudinal_events.tsv").open(), delimiter="\t"))
    assert [r["event"] for r in rows] == ["new_detection"] and rows[0]["l4_eligible"] == "True"
    assert list(rows[0]) == mod.EVENT_COLS
    traj = list(csv.DictReader((out / "longitudinal_trajectories.tsv").open(), delimiter="\t"))
    assert list(traj[0]) == mod.TRAJ_COLS and [
        t["detected"] for t in traj if t["subject_id"] == "S1"
    ] == ["False", "True"]
    man = json.loads((out / "longitudinal_manifest.json").read_text())
    assert (
        man["event_counts"] == {"new_detection": 1} and man["config"]["burden_fold_change"] == 2.0
    )
