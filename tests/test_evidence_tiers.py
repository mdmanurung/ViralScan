"""TIER-01: the probable/strong tier code against the frozen calibration-grid contract (synthetic inputs only)."""

from __future__ import annotations

import itertools

import pytest

from viralscan import evidence_tiers as et
from viralscan.evidence_tiers import Calibration, Cutoffs, GridPoint, TierEvidence


@pytest.fixture(scope="module")
def grid() -> et.Grid:
    # the frozen 45-point grid (calibration.threshold_search.grid); consistency with the protocol is tested in the reproducibility repo
    return et.Grid((1, 2, 3, 5, 10), (1, 2, 3), (0.5, 0.75, 0.9))


def _results(grid, f1=lambda c: 0.5, se=lambda c: 0.0, bad=lambda c: False, homology=lambda c: 0):
    return [GridPoint(c, f1(c), se(c), 1 if bad(c) else 0, homology(c)) for c in grid.points()]


# --- the grid is the protocol's, and only the protocol's ---------------------------------------


def test_grid_is_the_frozen_45_points(grid):
    assert grid.molecules == (1, 2, 3, 5, 10)
    assert grid.features == (1, 2, 3)
    assert grid.support == (0.5, 0.75, 0.9)
    assert len(grid.points()) == 45 == len(set(grid.points()))


def test_grid_parser_needs_the_protocol_key():
    with pytest.raises(ValueError, match="no calibration.threshold_search.grid"):
        et.grid_from_protocol({})


# --- selection: eligibility ---------------------------------------------------------------------


def test_best_f1_wins_without_ties(grid):
    target = Cutoffs(3, 2, 0.75)
    cal = et.select_cutoffs(_results(grid, f1=lambda c: 0.9 if c == target else 0.5), grid)
    assert (cal.status, cal.probable, cal.strong) == ("selected", target, None)
    assert cal.n_eligible == 45 and cal.best_f1 == 0.9


def test_a_violating_point_is_ineligible_whatever_its_f1(grid):
    top = Cutoffs(1, 1, 0.5)
    cal = et.select_cutoffs(
        _results(grid, f1=lambda c: 1.0 if c == top else 0.4, bad=lambda c: c == top), grid
    )
    assert cal.probable != top and cal.n_eligible == 44


def test_either_negative_class_makes_a_point_ineligible(grid):
    top = Cutoffs(1, 1, 0.5)
    f1 = lambda c: 1.0 if c == top else 0.4  # noqa: E731
    host = et.select_cutoffs(_results(grid, f1=f1, bad=lambda c: c == top), grid)
    homology = et.select_cutoffs(
        _results(grid, f1=f1, homology=lambda c: 3 if c == top else 0), grid
    )
    assert host.probable != top and homology.probable != top


def test_no_eligible_point_fails_and_nothing_is_relaxed(grid):
    cal = et.select_cutoffs(_results(grid, bad=lambda c: True), grid)
    assert cal.status == "failed" and cal.probable is None and cal.strong is None
    # the most conservative corner is not rescued either
    assert et.assign_tier(TierEvidence(99, 9, 1.0), cal) == "no_tier"


def test_missing_negative_count_refuses_rather_than_reading_zero(grid):
    points = _results(grid)
    points[0] = GridPoint(points[0].cutoffs, 0.5, 0.0, None, 0)
    with pytest.raises(ValueError, match="both negative-class call counts"):
        et.select_cutoffs(points, grid)


@pytest.mark.parametrize("field", ["f1", "se"])
@pytest.mark.parametrize("value", [float("nan"), -0.1, float("inf")])
def test_non_finite_or_out_of_range_inputs_refuse(grid, field, value):
    points = _results(grid)
    p = points[3]
    points[3] = GridPoint(
        p.cutoffs, value if field == "f1" else p.f1, value if field == "se" else p.se, 0, 0
    )
    with pytest.raises(ValueError):
        et.select_cutoffs(points, grid)


def test_grid_must_be_complete_distinct_and_closed(grid):
    full = _results(grid)
    with pytest.raises(ValueError, match="without a training result"):
        et.select_cutoffs(full[1:], grid)
    with pytest.raises(ValueError, match="twice"):
        et.select_cutoffs(full + full[:1], grid)
    with pytest.raises(ValueError, match="not a grid point"):
        et.select_cutoffs(full + [GridPoint(Cutoffs(4, 1, 0.5), 1.0, 0.0, 0, 0)], grid)


# --- selection: tie handling ---------------------------------------------------------------------


def test_within_one_se_the_most_conservative_point_wins(grid):
    best = Cutoffs(1, 1, 0.5)
    near = Cutoffs(5, 2, 0.75)  # 0.04 below best, inside SE 0.05
    far = Cutoffs(10, 3, 0.9)  # 0.2 below best, outside
    f1 = {best: 0.8, near: 0.76, far: 0.6}
    cal = et.select_cutoffs(
        _results(grid, f1=lambda c: f1.get(c, 0.1), se=lambda c: 0.05 if c == best else 0.0), grid
    )
    assert cal.probable == near and cal.best_f1 == 0.8 and cal.best_se == 0.05


def test_the_boundary_of_one_se_is_inclusive(grid):
    best, edge = Cutoffs(1, 1, 0.5), Cutoffs(10, 1, 0.5)
    f1 = {best: 0.75, edge: 0.5}
    cal = et.select_cutoffs(
        _results(grid, f1=lambda c: f1.get(c, 0.1), se=lambda c: 0.25 if c == best else 0.0), grid
    )
    assert cal.probable == edge


def test_conservative_order_is_molecules_then_features_then_support(grid):
    # every eligible point ties on F1 and SE: the lexicographic maximum of the cutoff tuple is the corner
    assert et.select_cutoffs(_results(grid), grid).probable == Cutoffs(10, 3, 0.9)
    # molecules dominate features and support
    ok = {Cutoffs(10, 1, 0.5), Cutoffs(5, 3, 0.9)}
    cal = et.select_cutoffs(_results(grid, bad=lambda c: c not in ok), grid)
    assert cal.probable == Cutoffs(10, 1, 0.5)
    # then features dominate support
    ok = {Cutoffs(5, 3, 0.5), Cutoffs(5, 2, 0.9)}
    cal = et.select_cutoffs(_results(grid, bad=lambda c: c not in ok), grid)
    assert cal.probable == Cutoffs(5, 3, 0.5)
    # then the larger support fraction
    ok = {Cutoffs(5, 3, 0.5), Cutoffs(5, 3, 0.75)}
    assert et.select_cutoffs(_results(grid, bad=lambda c: c not in ok), grid).probable == Cutoffs(
        5, 3, 0.75
    )


def test_the_selection_does_not_depend_on_input_order(grid):
    f1 = {c: ((i * 5) % 7) / 10 for i, c in enumerate(grid.points())}
    points = _results(grid, f1=lambda c: f1[c], se=lambda c: 0.1)
    expected = et.select_cutoffs(points, grid)
    assert et.select_cutoffs(reversed(points), grid) == expected


def test_a_tied_best_f1_with_zero_se_keeps_every_tied_point_in_the_set(grid):
    f1 = {Cutoffs(2, 1, 0.5): 0.9, Cutoffs(3, 1, 0.5): 0.9}
    cal = et.select_cutoffs(_results(grid, f1=lambda c: f1.get(c, 0.2)), grid)
    assert cal.probable == Cutoffs(3, 1, 0.5)


# --- supplied calibration result ------------------------------------------------------------------


def test_round_trip_through_a_plain_dict(grid):
    cal = et.select_cutoffs(
        _results(grid, f1=lambda c: 0.9 if c == Cutoffs(2, 2, 0.75) else 0.3), grid
    )
    assert et.calibration_from_dict(cal.to_dict(), grid) == cal


def test_a_supplied_point_off_the_grid_is_refused(grid):
    bad = {"status": "selected", "probable": {"molecules": 4, "features": 1, "support": 0.5}}
    with pytest.raises(ValueError, match="not a grid point"):
        et.calibration_from_dict(bad, grid)
    bad = {"status": "selected", "probable": {"molecules": 1, "features": 1, "support": 0.5},
           "strong": {"molecules": 4, "features": 3, "support": 0.9}}  # fmt: skip
    with pytest.raises(ValueError, match="not a grid point"):
        et.calibration_from_dict(bad, grid)


def test_strong_must_be_nested_above_probable(grid):
    probable = {"molecules": 5, "features": 2, "support": 0.75}
    for strong in (
        {"molecules": 3, "features": 3, "support": 0.9},
        {"molecules": 10, "features": 1, "support": 0.9},
        {"molecules": 10, "features": 3, "support": 0.5},
    ):
        with pytest.raises(ValueError, match="nested"):
            et.calibration_from_dict(
                {"status": "selected", "probable": probable, "strong": strong}, grid
            )
    ok = {"status": "selected", "probable": probable, "strong": probable}
    assert et.calibration_from_dict(ok, grid).strong == Cutoffs(5, 2, 0.75)


@pytest.mark.parametrize(
    "bad",
    [
        {},
        {"status": "ok"},
        {"status": "selected"},
        {"status": "selected", "probable": {"molecules": 1}},
    ],
)
def test_malformed_results_are_refused(grid, bad):
    with pytest.raises(ValueError):
        et.calibration_from_dict(bad, grid)


def test_a_failed_result_has_no_cutoffs(grid):
    cal = et.calibration_from_dict(
        {"status": "failed", "probable": {"molecules": 1, "features": 1, "support": 0.5}}, grid
    )
    assert cal == Calibration("failed")


# --- tier assignment -----------------------------------------------------------------------------


CAL = Calibration("selected", Cutoffs(3, 2, 0.75), Cutoffs(10, 3, 0.9))
PROBABLE_ONLY = Calibration("selected", Cutoffs(3, 2, 0.75))


def test_without_a_calibration_no_evidence_reaches_a_tier():
    for ev in (TierEvidence(1, 1, 1.0), TierEvidence(10**6, 99, 1.0), TierEvidence(0, 0, None)):
        assert et.assign_tier(ev, None) == "uncalibrated"
    assert {"probable", "strong"}.isdisjoint({et.assign_tier(TierEvidence(99, 9, 1.0), None)})


def test_probable_and_strong_follow_the_cutoffs_inclusively():
    assert et.assign_tier(TierEvidence(3, 2, 0.75), CAL) == "probable"  # exactly at the cutoffs
    assert et.assign_tier(TierEvidence(10, 3, 0.9), CAL) == "strong"
    assert et.assign_tier(TierEvidence(10, 3, 0.89), CAL) == "probable"
    for below in (TierEvidence(2, 2, 0.75), TierEvidence(3, 1, 0.75), TierEvidence(3, 2, 0.74)):
        assert et.assign_tier(below, CAL) == "candidate"


def test_unmeasured_support_never_reaches_a_tier_and_no_molecule_is_not_detected():
    assert et.assign_tier(TierEvidence(50, 3, None), CAL) == "candidate"
    assert et.assign_tier(TierEvidence(0, 0, None), CAL) == "not_detected"


def test_without_a_strong_point_strong_is_never_emitted():
    assert et.assign_tier(TierEvidence(10**6, 99, 1.0), PROBABLE_ONLY) == "probable"


def test_tiers_are_nested_over_the_whole_grid(grid):
    pts = grid.points()
    rank = {"not_detected": 0, "candidate": 1, "probable": 2, "strong": 3}
    evidence = [
        TierEvidence(m, f, x) for m in (0, 1, 5, 10) for f in (0, 1, 3) for x in (0.4, 0.75, 1.0)
    ]
    for p, s in itertools.product(pts, pts):
        if any(a < b for a, b in zip(s, p, strict=True)):
            continue  # not nested: calibration_from_dict refuses it
        for ev in evidence:
            tier = et.assign_tier(ev, Calibration("selected", p, s))
            # adding a strong point never removes a probable call, and a strong call is also a probable one
            assert rank[tier] >= rank[et.assign_tier(ev, Calibration("selected", p))]
            if tier == "strong":
                assert et.assign_tier(ev, Calibration("selected", s)) == "probable"


def test_a_planted_homology_pattern_cannot_reach_a_tier_under_a_failed_calibration():
    # REF-09's shape: reads that pass every cutoff still get no tier when calibration failed or is absent
    reads = TierEvidence(10, 3, 1.0)
    assert et.assign_tier(reads, Calibration("failed")) == "no_tier"
    assert et.assign_tier(reads, None) == "uncalibrated"


# --- the study-level mapping ----------------------------------------------------------------------


def test_the_mapping_table_is_the_decided_one():
    assert et.STUDY_LEVEL_TO_PACKAGE_TIER == {
        "L0": "candidate",
        "L1": "probable",
        "L2": "strong",
        "L3": "strong",
    }
    assert set(et.STUDY_LEVEL_TO_PACKAGE_TIER.values()) <= set(et.TIERS)


def test_unknown_levels_have_no_analogue():
    for level in ("L4", "L5", "", "1"):
        with pytest.raises(ValueError, match="no package-tier analogue"):
            et.package_tier_analogue(level)
