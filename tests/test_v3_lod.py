"""LOD estimator fixtures (PLAN VAL-06): hand-computed probit, no-extrapolation, BCa rules. No outcomes."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import v3_lod as lod  # noqa: E402

LEVELS = [1, 3, 10, 30, 100, 300, 1000]


def exact_rows(intercept, slope, levels, n=20000):
    """Rows whose detection rate at each level is the probit probability (rounded to 1/n)."""
    xs, ys = [], []
    for level in levels:
        p = stats.norm.cdf(intercept + slope * math.log10(level))
        k = round(p * n)
        xs += [level] * n
        ys += [1] * k + [0] * (n - k)
    return xs, ys


def test_the_mle_recovers_a_curve_that_is_exactly_probit():
    # eta = -4 + 1.5 * log10(x): LOD95 = 10 ** ((z_0.95 + 4) / 1.5) = 5795 inside [10, 10000]
    xs, ys = exact_rows(-4.0, 1.5, [10, 100, 1000, 10000])
    fit = lod.fit_probit(xs, ys)
    assert fit["status"] == "ok"
    assert fit["intercept"] == pytest.approx(-4.0, abs=0.01) and fit["slope"] == pytest.approx(
        1.5, abs=0.01
    )
    result = lod.lod95(xs, ys)
    expected = 10 ** ((stats.norm.ppf(0.95) + 4.0) / 1.5)
    assert result["status"] == "ok" and result["lod"] == pytest.approx(expected, rel=0.01)
    assert (result["lowest_tested"], result["highest_tested"]) == (10, 10000)


def test_two_level_saturated_fit_is_hand_computable():
    # 20/100 at 10 and 80/100 at 1000: slope = (z(.8) - z(.2)) / 2, intercept = z(.2) - slope
    xs = [10] * 100 + [1000] * 100
    ys = [1] * 20 + [0] * 80 + [1] * 80 + [0] * 20
    fit = lod.fit_probit(xs, ys)
    slope = (stats.norm.ppf(0.8) - stats.norm.ppf(0.2)) / 2
    assert fit["slope"] == pytest.approx(slope, abs=1e-4)
    assert fit["intercept"] == pytest.approx(stats.norm.ppf(0.2) - slope, abs=1e-4)


def test_no_extrapolation_above_or_below_the_tested_range():
    xs, ys = exact_rows(-4.0, 1.5, [10, 100, 1000])  # LOD95 = 5795 > 1000
    high = lod.lod95(xs, ys)
    assert (
        high["status"] == "not_reached" and high["lod"] is None and high["highest_tested"] == 1000
    )
    xs, ys = exact_rows(1.0, 1.5, [100, 1000, 10000])  # LOD95 far below 100
    low = lod.lod95(xs, ys)
    assert low["status"] == "below_range" and low["lod"] is None and low["lowest_tested"] == 100


@pytest.mark.parametrize(
    ("xs", "ys", "status"),
    [
        ([10, 10, 10], [0, 1, 1], "one_level"),
        ([10, 100, 1000], [0, 0, 0], "no_variation"),
        ([10, 100, 1000], [1, 1, 1], "no_variation"),
        ([10, 10, 100, 100, 1000, 1000], [0, 0, 0, 1, 1, 1], "ok"),  # overlap at 100: finite MLE
        ([10, 10, 100, 100, 1000, 1000], [0, 0, 0, 0, 1, 1], "separation"),
        ([10, 10, 100, 100, 1000, 1000], [1, 1, 1, 0, 0, 0], "non_positive_slope"),
    ],
)
def test_statuses_without_a_finite_positive_slope_are_named(xs, ys, status):
    result = lod.fit_probit(xs, ys)
    assert result["status"] == status
    if status != "ok":
        assert lod.lod95(xs, ys)["lod"] is None


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError, match="positive"):
        lod.fit_probit([0, 10], [0, 1])
    with pytest.raises(ValueError, match="differ"):
        lod.fit_probit([1, 10], [1])
    with pytest.raises(ValueError, match="binary"):
        lod.fit_probit([1, 10], [2, 0])


def simulate(seed=1, per_level=12, intercept=-2.2, slope=1.3):
    rng = np.random.default_rng(seed)
    units = {}
    for level in LEVELS:
        for i in range(per_level):
            p = stats.norm.cdf(intercept + slope * math.log10(level))
            units[f"s{level}_{i}"] = [(level, int(rng.random() < p))]
    return units


def test_interval_is_reproducible_and_brackets_the_point_estimate():
    units = simulate()
    strata = {k: k.split("_")[0] for k in units}
    one = lod.lod_interval(units, seed=7, strata=strata, replicates=300)
    two = lod.lod_interval(units, seed=7, strata=strata, replicates=300)
    assert one == two and one["interval_status"] == "ok"
    low, high = one["interval"]
    assert low <= one["point"]["lod"] <= high
    assert (one["n_units"], one["replicates"], one["seed"]) == (len(units), 300, 7)
    assert (
        lod.lod_interval(units, seed=8, strata=strata, replicates=300)["interval"]
        != one["interval"]
    )


def test_fewer_than_three_units_gives_a_point_estimate_and_no_interval():
    units = {"a": [(10, 0), (100, 1), (1000, 1)], "b": [(10, 0), (100, 0), (1000, 1)]}
    result = lod.lod_interval(units, seed=1, replicates=50)
    assert result["interval"] is None and result["interval_status"] == "not_estimable"
    assert result["point"]["status"] in {"ok", "not_reached", "below_range", "separation"}


def test_a_unit_is_resampled_whole(monkeypatch):
    """Everything nested in a unit (all its rows) enters a replicate together, never row by row."""
    units = {f"u{i}": [(10, i % 2), (100, 1), (1000, 1), (30, i % 3 == 0)] for i in range(6)}
    seen = []
    original = lod._estimate

    def spy(units_, keys):
        seen.append(list(keys))
        return original(units_, keys)

    monkeypatch.setattr(lod, "_estimate", spy)
    lod.lod_interval(units, seed=3, replicates=20)
    assert all(len(draw) in (6, 5) and set(draw) <= set(units) for draw in seen)
    assert any(len(set(draw)) < len(draw) for draw in seen)  # resampled with replacement


def test_an_inestimable_replicate_withholds_the_interval_with_a_reason():
    # separation appears in some bootstrap draws of this sparse design
    units = {f"u{i}": [(10, 0), (1000, 1)] for i in range(2)} | {
        "x": [(10, 1), (1000, 1)],
        "y": [(10, 0), (1000, 0)],
    }
    result = lod.lod_interval(units, seed=1, replicates=200)
    assert result["interval"] is None and result["reason"]
