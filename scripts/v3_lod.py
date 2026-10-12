#!/usr/bin/env python3
"""Limit-of-detection estimator of the v3 protocol (PLAN VAL-06, ``limit_of_detection`` + ``uncertainty``).

Pure functions, no outcome and no data: ``fit_probit`` is the probit regression of binary per-sample
detection on log10 planted abundance, ``lod95`` the abundance at which the fitted probability reaches
0.95, and ``lod_interval`` the bias-corrected and accelerated (BCa) bootstrap interval over the
resampling unit (``partitions.unit``: a generated sample, or a whole donor library with everything
nested in it), with the protocol's rules written out:

* one fit per virus and chemistry, the caller selects the rows;
* no extrapolation: a fitted LOD above the highest tested level is reported ``not_reached`` (with that
  level), one below the lowest tested level ``below_range``, never a number;
* fewer than ``minimum_samples_for_interval`` units: the point estimate only, interval ``not_estimable``;
* separation (every non-detected sample below every detected one) has no finite probit MLE: reported
  ``separation``, not regularised, so the choice stays the user's;
* the interval needs every bootstrap replicate and every leave-one-unit-out fit to be estimable;
  otherwise it is ``not_estimable`` and says why. No share of failed replicates is tolerated or
  excused, so no threshold is introduced here.

``replicates`` defaults to the protocol's 2000 and the interval to 95 percent; ``seed`` has no default
because ``seeds.bootstrap`` is frozen at G3.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from scipy import optimize, special, stats

LOD_LEVEL = 0.95  # protocol: LOD95
REPLICATES = 2000  # protocol: uncertainty.replicates
MINIMUM_UNITS = 3  # protocol: uncertainty.minimum_samples_for_interval
CONFIDENCE = 0.95  # protocol: uncertainty.interval


def fit_probit(abundance: Sequence[float], detected: Sequence[int | bool]) -> dict[str, Any]:
    """Probit MLE of ``P(detected) = Phi(intercept + slope * log10(abundance))``.

    ``status`` is ``ok`` or the reason there is no finite MLE: ``one_level`` (fewer than two distinct
    abundances), ``no_variation`` (all or none detected), ``separation``, ``non_positive_slope``.
    """
    x = np.asarray(abundance, dtype=float)
    y = np.asarray([1 if d else 0 for d in detected], dtype=float)
    if len(x) != len(y):
        raise ValueError("abundance and detected differ in length")
    if len(x) == 0 or not np.all(np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("abundances must be positive and finite")
    if any(d not in (0, 1) for d in detected):
        raise ValueError("detection must be binary")
    result: dict[str, Any] = dict(status="ok", intercept=None, slope=None, n=len(x))
    if len(set(x)) < 2:
        return dict(result, status="one_level")
    if y.min() == y.max():
        return dict(result, status="no_variation")
    if x[y == 0].max() < x[y == 1].min():
        return dict(result, status="separation")
    z = np.log10(x)

    def negative_log_likelihood(beta: np.ndarray) -> float:
        eta = beta[0] + beta[1] * z
        return -float(np.sum(y * special.log_ndtr(eta) + (1 - y) * special.log_ndtr(-eta)))

    start = np.array([0.0, 1.0])
    fit = optimize.minimize(negative_log_likelihood, start, method="BFGS")
    intercept, slope = (float(v) for v in fit.x)
    if not np.all(np.isfinite(fit.x)):
        return dict(result, status="separation")
    if slope <= 0:
        return dict(result, status="non_positive_slope", intercept=intercept, slope=slope)
    return dict(result, intercept=intercept, slope=slope)


def lod95(
    abundance: Sequence[float], detected: Sequence[int | bool], level: float = LOD_LEVEL
) -> dict[str, Any]:
    """LOD at ``level`` with the no-extrapolation rule; ``lod`` is a number only when ``status == ok``."""
    fit = fit_probit(abundance, detected)
    low, high = float(min(abundance)), float(max(abundance))
    out: dict[str, Any] = dict(fit, lod=None, lowest_tested=low, highest_tested=high)
    if fit["status"] != "ok":
        return out
    value = 10 ** ((float(stats.norm.ppf(level)) - fit["intercept"]) / fit["slope"])
    if value > high:
        out["status"] = "not_reached"
    elif value < low:
        out["status"] = "below_range"
    else:
        out["lod"] = value
    return out


def _estimate(
    units: Mapping[Any, Sequence[tuple[float, int]]], keys: Sequence[Any]
) -> dict[str, Any]:
    rows = [row for key in keys for row in units[key]]
    return lod95([a for a, _ in rows], [d for _, d in rows])


def _comparable(result: Mapping[str, Any]) -> float | None:
    """Ordered value of a replicate: the LOD, ``inf`` if not reached; ``None`` if not estimable."""
    if result["status"] == "ok":
        return float(result["lod"])
    if result["status"] == "not_reached":
        return math.inf
    return None


def lod_interval(
    units: Mapping[Any, Sequence[tuple[float, int]]],
    *,
    seed: int,
    strata: Mapping[Any, Any] | None = None,
    replicates: int = REPLICATES,
    confidence: float = CONFIDENCE,
    minimum_units: int = MINIMUM_UNITS,
) -> dict[str, Any]:
    """LOD95 of the pooled units and its BCa interval from a unit-level (cluster) bootstrap.

    ``units`` maps a resampling unit to its ``(abundance, detected)`` rows, so everything nested in a
    unit is resampled together. ``strata`` (unit -> label) resamples units within a stratum, keeping the
    design's chemistry/abundance structure; default one stratum.
    """
    keys = sorted(units, key=str)
    point = _estimate(units, keys)
    out: dict[str, Any] = dict(
        point=point, n_units=len(keys), replicates=replicates, seed=seed, interval=None
    )
    if len(keys) < minimum_units:
        return dict(out, interval_status="not_estimable", reason="fewer than the minimum units")
    if point["status"] != "ok":
        return dict(
            out, interval_status="not_estimable", reason=f"point estimate {point['status']}"
        )
    groups: dict[Any, list[Any]] = {}
    for key in keys:
        groups.setdefault(None if strata is None else strata[key], []).append(key)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        draw = [
            members[i]
            for members in groups.values()
            for i in rng.integers(len(members), size=len(members))
        ]
        values.append(_comparable(_estimate(units, draw)))
    if any(v is None for v in values):
        return dict(
            out, interval_status="not_estimable", reason="a bootstrap replicate was not estimable"
        )
    jack = [_comparable(_estimate(units, [k for k in keys if k != drop])) for drop in keys]
    if any(v is None or not math.isfinite(v) for v in jack):
        return dict(
            out,
            interval_status="not_estimable",
            reason="a leave-one-unit-out fit was not estimable",
        )
    boot = np.sort(np.asarray(values, dtype=float))
    theta = point["lod"]
    below = float(np.mean(boot < theta) + 0.5 * np.mean(boot == theta))
    if not 0 < below < 1:
        return dict(out, interval_status="not_estimable", reason="bias correction is undefined")
    z0 = float(stats.norm.ppf(below))
    centred = np.mean(jack) - np.asarray(jack, dtype=float)
    denominator = 6 * float(np.sum(centred**2)) ** 1.5
    acceleration = float(np.sum(centred**3)) / denominator if denominator else 0.0
    bounds = []
    for tail in ((1 - confidence) / 2, 1 - (1 - confidence) / 2):
        zt = float(stats.norm.ppf(tail))
        adjusted = float(stats.norm.cdf(z0 + (z0 + zt) / (1 - acceleration * (z0 + zt))))
        bounds.append(adjusted)
    # order statistics, widening outward, so an infinite (not reached) replicate cannot make a NaN
    lower = float(np.quantile(boot, bounds[0], method="lower"))
    upper = float(np.quantile(boot, bounds[1], method="higher"))
    return dict(
        out,
        interval=(lower, upper),
        interval_status="ok",
        bias_correction=z0,
        acceleration=acceleration,
        adjusted_levels=tuple(bounds),
    )
