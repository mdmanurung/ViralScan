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
GRADIENT_TOLERANCE = (
    1e-5  # max per-sample |d NLL / d beta| accepted as stationary when BFGS flags failure
)


def bca_interval(
    point: float,
    boot: Sequence[float],
    jack: Sequence[float],
    *,
    confidence: float = CONFIDENCE,
) -> dict[str, Any]:
    """BCa limits of a statistic from its bootstrap and leave-one-unit-out values (SCI-05 r13 BM4).

    Returns ``interval_status`` ``ok`` with ``interval`` or ``not_estimable`` with a ``reason_code``:

    * ``nonfinite``: the point, a jackknife value or a bootstrap value is NaN (an infinite
      bootstrap value is a legitimate ordered "not reached" and is allowed), or the jackknife is
      infinite;
    * ``constant_bootstrap``: every replicate is equal (all-zero false calls, say). A zero-width
      interval would read as perfect precision, so it is withheld;
    * ``bias_correction_undefined``: no replicate lies on one side of the point.
    """
    boot_a, jack_a = np.asarray(boot, dtype=float), np.asarray(jack, dtype=float)
    if (
        not math.isfinite(point)
        or np.isnan(boot_a).any()
        or not np.isfinite(jack_a).all()
        or len(boot_a) == 0
    ):
        return dict(
            interval=None,
            interval_status="not_estimable",
            reason_code="nonfinite",
            reason="a point, bootstrap or jackknife value is not finite",
        )
    boot_a = np.sort(boot_a)
    if boot_a[0] == boot_a[-1]:
        return dict(
            interval=None,
            interval_status="not_estimable",
            reason_code="constant_bootstrap",
            reason="every bootstrap replicate is equal; a zero-width interval would be false precision",
        )
    below = float(np.mean(boot_a < point) + 0.5 * np.mean(boot_a == point))
    if not 0 < below < 1:
        return dict(
            interval=None,
            interval_status="not_estimable",
            reason_code="bias_correction_undefined",
            reason="bias correction is undefined",
        )
    z0 = float(stats.norm.ppf(below))
    centred = np.mean(jack_a) - jack_a
    denominator = 6 * float(np.sum(centred**2)) ** 1.5
    acceleration = float(np.sum(centred**3)) / denominator if denominator else 0.0
    bounds = []
    for tail in ((1 - confidence) / 2, 1 - (1 - confidence) / 2):
        zt = float(stats.norm.ppf(tail))
        bounds.append(float(stats.norm.cdf(z0 + (z0 + zt) / (1 - acceleration * (z0 + zt)))))
    # order statistics, widening outward, so an infinite (not reached) replicate cannot make a NaN
    return dict(
        interval=(
            float(np.quantile(boot_a, bounds[0], method="lower")),
            float(np.quantile(boot_a, bounds[1], method="higher")),
        ),
        interval_status="ok",
        bias_correction=z0,
        acceleration=acceleration,
        adjusted_levels=tuple(bounds),
    )


def fit_probit(abundance: Sequence[float], detected: Sequence[int | bool]) -> dict[str, Any]:
    """Probit MLE of ``P(detected) = Phi(intercept + slope * log10(abundance))``.

    ``status`` is ``ok`` or the reason there is no finite MLE: ``one_level`` (fewer than two distinct
    abundances), ``no_variation`` (all or none detected), ``separation`` (complete: every non-detected sample
    is below every detected one. Quasi-complete separation, where they share a boundary abundance,
    also has no MLE but BFGS stops at a finite slope there, so it is not detected; the interval's
    bootstrap then usually fails on it. Upgrade: a profile-likelihood check), ``not_converged`` (the optimizer failed and the gradient is not near zero, so
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
    if x[y == 0].max() < x[y == 1].min():  # ponytail: complete separation only, see docstring
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
        return dict(
            out,
            interval_status="not_estimable",
            reason_code="fewer_than_minimum_units",
            reason="fewer than the minimum units",
        )
    if point["status"] != "ok":
        return dict(
            out,
            interval_status="not_estimable",
            reason_code="point_estimate_failed",
            reason=f"point estimate {point['status']}",
        )
    groups: dict[Any, list[Any]] = {}
    for key in keys:
        groups.setdefault(None if strata is None else strata[key], []).append(key)
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
    # BM6: a failed fit keeps its reason code (separation, not_converged, ...) in the result
    failed_bootstrap = _status_counts(r for r, v in zip(results, values, strict=True) if v is None)
    failed_jackknife = _status_counts(
        r for r, v in zip(jacks, jack, strict=True) if v is None or not math.isfinite(v)
    )
    if failed_bootstrap:
        return dict(
            out,
            interval_status="not_estimable",
            reason_code="failed_fit",
            reason="a bootstrap replicate was not estimable",
            failed_fit_statuses=failed_bootstrap,
        )
    if failed_jackknife:
        return dict(
            out,
            interval_status="not_estimable",
            reason_code="failed_fit",
            reason="a leave-one-unit-out fit was not estimable",
            failed_fit_statuses=failed_jackknife,
        )
    return dict(out, **bca_interval(point["lod"], values, jack, confidence=confidence))


def _status_counts(results: Any) -> dict[str, int]:
    counts: dict[str, int] = {}
    for result in results:
        counts[result["status"]] = counts.get(result["status"], 0) + 1
    return counts
