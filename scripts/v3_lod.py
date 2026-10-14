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
* separation (every non-detected sample below every detected one) and quasi-separation (they share
  their boundary abundance) have no finite probit MLE: reported ``separation`` or
  ``quasi_separation``, not regularised, so the choice stays the user's;
* the interval needs every bootstrap replicate and every leave-one-unit-out fit to be estimable;
  otherwise it is ``not_estimable`` with reason ``failed-fit`` (or ``lod-not-reached`` when every
  failing replicate only failed to reach the LOD) and the failed statuses are listed. No share of
  failed replicates is tolerated or excused, so no threshold is introduced here. Reason codes
  are the protocol's hyphenated ones (``calibration.uncertainty.degenerate_interval_policy``);
* whenever the interval is not estimable (any reason code, a failed point fit included) the result
  also carries ``bracket``, the model-free fallback of :func:`lod_bracket` (SCI-05 r15 N1, Q139 a):
  the observed detection per tested level, no fit.

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
GRADIENT_TOLERANCE = (
    1e-5  # max per-sample |d NLL / d beta| accepted as stationary when BFGS flags failure
)


def _not_estimable(code: str, reason: str) -> dict[str, Any]:
    return dict(interval=None, interval_status="not_estimable", reason_code=code, reason=reason)


def _kind(value: float | None) -> str | None:
    """``undefined`` for ``None``, ``nonfinite`` for NaN or an infinity, else ``None``."""
    if value is None:
        return "undefined"
    return None if math.isfinite(value) else "nonfinite"


def bca_interval(
    point: float | None,
    boot: Sequence[float | None],
    jack: Sequence[float | None],
    *,
    confidence: float = CONFIDENCE,
) -> dict[str, Any]:
    """BCa limits of a statistic from its bootstrap and leave-one-unit-out values.

    Follows the protocol's ``calibration.uncertainty.degenerate_interval_policy`` (SCI-05 r14 M5, V1,
    V2). Returns ``interval_status`` ``ok`` with ``interval``, or ``not_estimable`` with one
    ``reason_code`` (hyphenated, as the protocol writes them):

    * ``undefined-statistic``: the point, a bootstrap or a jackknife value is ``None`` (undefined);
    * ``nonfinite-statistic``: any of them is NaN or an infinity (an infinite LOD replicate is
      not an ordered "not reached" here: the LOD caller reports that case itself);
    * ``constant-bootstrap-distribution``: every replicate is equal; a zero-width interval would read
      as perfect precision;
    * ``bias-correction-undefined``: no replicate lies on one side of the point;
    * ``acceleration-undefined``: the jackknife is constant (zero denominator, never a silent
      zero acceleration) or the adjustment ``1 - a (z0 + z)`` is not positive;
    * ``invalid-interval``: a bound is not finite, or the bounds are not ordered with positive width.
    """
    values = [point, *boot, *jack]
    kinds = {_kind(v) for v in values} - {None}
    if len(boot) == 0 or len(jack) == 0:
        kinds.add("undefined")
    if "undefined" in kinds:
        return _not_estimable(
            "undefined-statistic",
            "the statistic is undefined at the point, a replicate or a leave-one-out fit",
        )
    if kinds:
        return _not_estimable(
            "nonfinite-statistic", "a point, bootstrap or jackknife value is not finite"
        )
    assert point is not None
    boot_a, jack_a = np.sort(np.asarray(boot, dtype=float)), np.asarray(jack, dtype=float)
    if boot_a[0] == boot_a[-1]:
        return _not_estimable(
            "constant-bootstrap-distribution",
            "every bootstrap replicate is equal; a zero-width interval would be false precision",
        )
    below = float(np.mean(boot_a < point) + 0.5 * np.mean(boot_a == point))
    if not 0 < below < 1:
        return _not_estimable("bias-correction-undefined", "bias correction is undefined")
    z0 = float(stats.norm.ppf(below))
    centred = np.mean(jack_a) - jack_a
    denominator = 6 * float(np.sum(centred**2)) ** 1.5
    if not denominator > 0:
        return _not_estimable(
            "acceleration-undefined", "the jackknife is constant, so the acceleration is undefined"
        )
    acceleration = float(np.sum(centred**3)) / denominator
    bounds = []
    for tail in ((1 - confidence) / 2, 1 - (1 - confidence) / 2):
        zt = float(stats.norm.ppf(tail))
        adjust = 1 - acceleration * (z0 + zt)
        if not adjust > 0:
            return _not_estimable(
                "acceleration-undefined",
                "the acceleration adjustment 1 - a(z0 + z) is not positive",
            )
        bounds.append(float(stats.norm.cdf(z0 + (z0 + zt) / adjust)))
    interval = (
        float(np.quantile(boot_a, bounds[0], method="lower")),
        float(np.quantile(boot_a, bounds[1], method="higher")),
    )
    if not (all(math.isfinite(b) for b in (*bounds, *interval)) and interval[0] < interval[1]):
        return _not_estimable(
            "invalid-interval", "the interval bounds are not finite, ordered and of positive width"
        )
    return dict(
        interval=interval,
        interval_status="ok",
        bias_correction=z0,
        acceleration=acceleration,
        adjusted_levels=tuple(bounds),
    )


def fit_probit(abundance: Sequence[float], detected: Sequence[int | bool]) -> dict[str, Any]:
    """Probit MLE of ``P(detected) = Phi(intercept + slope * log10(abundance))``.

    ``status`` is ``ok`` or the reason there is no finite MLE: ``one_level`` (fewer than two distinct
    abundances), ``no_variation`` (all or none detected), ``separation`` (complete: every non-detected sample
    is below every detected one), ``quasi_separation`` (the two groups share their boundary
    abundance, so no finite MLE exists either and the optimizer's stopping slope is arbitrary; SCI-05
    r14 P1), ``not_converged`` (the optimizer failed and the gradient is not near zero, so
    the numbers are not an MLE), ``non_positive_slope``. Only ``ok`` carries an intercept and slope
    that may be used; SCI-05 r13 BM6: the optimizer's ``success`` flag used to be ignored.
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
    # Complete (strict) or quasi-complete (the groups touch at one abundance) separation: no finite
    # MLE, and at a tie the optimizer stops at an arbitrary slope, so the LOD would depend on the
    # tolerance (SCI-05 r14 P1). In one dimension this check is exact, not a heuristic.
    if x[y == 0].max() < x[y == 1].min():
        return dict(result, status="separation")
    if x[y == 0].max() == x[y == 1].min():
        return dict(result, status="quasi_separation")
    z = np.log10(x)

    def negative_log_likelihood(beta: np.ndarray) -> float:
        eta = beta[0] + beta[1] * z
        return -float(np.sum(y * special.log_ndtr(eta) + (1 - y) * special.log_ndtr(-eta)))

    start = np.array([0.0, 1.0])
    fit = optimize.minimize(negative_log_likelihood, start, method="BFGS")
    intercept, slope = (float(v) for v in fit.x)
    if not np.all(np.isfinite(fit.x)):
        return dict(result, status="separation")
    # BFGS often reports "precision loss" at a true optimum, so a failed flag is accepted only
    # when the gradient confirms a stationary point.
    if not fit.success and not np.max(np.abs(fit.jac)) / len(x) < GRADIENT_TOLERANCE:
        return dict(result, status="not_converged", intercept=intercept, slope=slope)
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
    # Classify in LOG space before exponentiating: a tiny positive slope makes 10 ** x overflow
    # (SCI-05 r15, Biomni residual Major). A non-finite log LOD is its own failed status, never a number.
    with np.errstate(all="ignore"):
        log_value = (float(stats.norm.ppf(level)) - fit["intercept"]) / fit["slope"]
    if math.isnan(log_value):
        out["status"] = "nonfinite_lod"
    elif log_value > math.log10(high):
        out["status"] = "not_reached"
    elif log_value < math.log10(low):
        out["status"] = "below_range"
    else:
        out["lod"] = 10**log_value
    return out


def lod_bracket(
    rows: Sequence[tuple[float, int | bool]], level: float = LOD_LEVEL
) -> dict[str, Any]:
    """Model-free bracket of the LOD from the observed detection per tested level (Q139 a).

    ``rows`` are the ``(abundance, detected)`` rows of one virus and chemistry. ``lower`` is the highest
    tested level whose empirical detection is below ``level``; ``upper`` the lowest tested level
    above it whose detection is at least ``level``. No fit, so no separation or bootstrap can fail it.
    ``status``: ``bracketed`` (both ends), ``above_range`` (no level above ``lower`` reaches
    ``level``: only ``lower`` is reported, LOD95 > ``lower``, the top tested level) or ``below_range``
    (every level reaches it: only ``upper``, LOD95 <= the lowest tested level). ``non_monotone`` flags
    a level at or above ``level`` that lies below ``lower``; the bracket then spans the dip.
    """
    counts: dict[float, list[int]] = {}
    for abundance, hit in rows:
        counts.setdefault(float(abundance), []).append(1 if hit else 0)
    if not counts:
        raise ValueError("no rows to bracket")
    levels = [
        dict(level=a, detected=sum(v), n=len(v), fraction=sum(v) / len(v))
        for a, v in sorted(counts.items())
    ]
    # tolerance: 19/20 must count as 0.95 whatever the float rounding
    reached = [r["detected"] >= level * r["n"] - 1e-9 for r in levels]
    lower = max((r["level"] for r, ok in zip(levels, reached, strict=True) if not ok), default=None)
    upper = min(
        (
            r["level"]
            for r, ok in zip(levels, reached, strict=True)
            if ok and (lower is None or r["level"] > lower)
        ),
        default=None,
    )
    status = "below_range" if lower is None else "bracketed" if upper is not None else "above_range"
    non_monotone = lower is not None and any(
        ok and r["level"] < lower for r, ok in zip(levels, reached, strict=True)
    )
    return dict(
        status=status,
        lower=lower,
        upper=upper,
        level=level,
        non_monotone=non_monotone,
        levels=levels,
    )


def _estimate(
    units: Mapping[Any, Sequence[tuple[float, int]]], keys: Sequence[Any]
) -> dict[str, Any]:
    rows = [row for key in keys for row in units[key]]
    return lod95([a for a, _ in rows], [d for _, d in rows])


def _comparable(result: Mapping[str, Any]) -> float | None:
    """The LOD of an ``ok`` fit; ``None`` for every other status, ``not_reached`` included.

    No infinity stands in for "not reached" (SCI-05 r14 M5): the protocol's failure policy makes
    a non-finite replicate a reason-coded not-estimable interval, not an ordered tail value.
    """
    return float(result["lod"]) if result["status"] == "ok" else None


def lod_interval(
    units: Mapping[Any, Sequence[tuple[float, int]]],
    *,
    seed: int,
    strata: Mapping[Any, Any],
    replicates: int = REPLICATES,
    confidence: float = CONFIDENCE,
    minimum_units: int = MINIMUM_UNITS,
) -> dict[str, Any]:
    """LOD95 of the pooled units and its BCa interval from a unit-level (cluster) bootstrap.

    ``units`` maps a resampling unit to its ``(abundance, detected)`` rows, so everything nested in a
    unit is resampled together. ``strata`` (unit -> label) is required, as for ``f1_interval``
    (SCI-05 r15 N4): units are resampled within a stratum, keeping the design's per-stratum counts
    (``partitions``/``uncertainty.stratification_unit_rule``); pass one label for a single stratum.
    When ``interval_status`` is not ``ok`` the result carries ``bracket`` (:func:`lod_bracket`).
    """
    out = _lod_interval(
        units,
        seed=seed,
        strata=strata,
        replicates=replicates,
        confidence=confidence,
        minimum_units=minimum_units,
    )
    if out["interval_status"] != "ok":
        out["bracket"] = lod_bracket([row for rows in units.values() for row in rows])
    return out


def _lod_interval(
    units: Mapping[Any, Sequence[tuple[float, int]]],
    *,
    seed: int,
    strata: Mapping[Any, Any],
    replicates: int,
    confidence: float,
    minimum_units: int,
) -> dict[str, Any]:
    keys = sorted(units, key=str)
    if missing := [k for k in keys if k not in strata]:
        raise ValueError(f"strata has no label for units {missing[:3]}")
    point = _estimate(units, keys)
    out: dict[str, Any] = dict(
        point=point, n_units=len(keys), replicates=replicates, seed=seed, interval=None
    )
    if len(keys) < minimum_units:
        return dict(
            out,
            interval_status="not_estimable",
            reason_code="below-minimum-units",
            reason="fewer than the minimum units",
        )
    if point["status"] != "ok":
        return dict(
            out,
            interval_status="not_estimable",
            reason_code="point-estimate-failed",
            reason=f"point estimate {point['status']}",
        )
    groups: dict[Any, list[Any]] = {}
    for key in keys:
        groups.setdefault(strata[key], []).append(key)
    rng = np.random.default_rng(seed)
    results = []
    for _ in range(replicates):
        draw = [
            members[i]
            for members in groups.values()
            for i in rng.integers(len(members), size=len(members))
        ]
        results.append(_estimate(units, draw))
    jacks = [_estimate(units, [k for k in keys if k != drop]) for drop in keys]
    values, jack = [_comparable(r) for r in results], [_comparable(r) for r in jacks]
    # BM6: a failed fit keeps its reason code (separation, not_reached, ...) in the result
    for label, group, found in (
        ("a bootstrap replicate", results, values),
        ("a leave-one-unit-out fit", jacks, jack),
    ):
        failed = _status_counts(r for r, v in zip(group, found, strict=True) if v is None)
        if failed:
            only_unreached = set(failed) == {"not_reached"}
            return dict(
                out,
                interval_status="not_estimable",
                reason_code="lod-not-reached" if only_unreached else "failed-fit",
                reason=f"{label} was not estimable",
                failed_fit_statuses=failed,
            )
    return dict(out, **bca_interval(point["lod"], values, jack, confidence=confidence))


def _status_counts(results: Any) -> dict[str, int]:
    counts: dict[str, int] = {}
    for result in results:
        counts[result["status"]] = counts.get(result["status"], 0) + 1
    return counts
