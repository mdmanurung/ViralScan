"""Probable/strong evidence tiers, built to the frozen calibration-grid contract (PLAN ``TIER-01``).

Pure functions, no file I/O, no data and no calibration run. The contract is
``calibration.threshold_search`` of ``analysis/v3_validation/protocol.yaml``:

* the grid is the 45 cutoff triples ``(molecules per cell, features per sample, exact-fragment support
  fraction)``; :func:`grid_from_protocol` reads it from the already-loaded protocol, nothing here hard-codes a
  value, and no cutoff off the grid is ever adopted;
* a grid point is **eligible** only with zero probable-or-strong calls in both synthetic host-only and planted
  host-homology training negatives (:func:`select_cutoffs`);
* among eligible points within one standard error of the best training F1, the **most conservative** wins: the
  largest molecule cutoff, then the largest feature cutoff, then the largest support fraction;
* no eligible point means calibration **failed** and the tier is ``no_tier`` (UC-C6); the constraint is not
  relaxed and the grid is not widened.

Without a supplied calibration result :func:`assign_tier` returns ``uncalibrated``, never ``probable`` or
``strong``. The caller computes each grid point's F1 and its bootstrap standard error (the protocol's
``calibration.uncertainty`` resampling-unit bootstrap); this module only applies the selection rule.

Open contract items (see PLAN ``TIER-01``): the grid defines only *probable-or-stronger*, so ``strong`` is an
optional, caller-supplied grid point nested above ``probable``; the three measures are computed by
:mod:`viralscan.tier_measures` (PLAN ``TIER-02``); comparisons are inclusive (``>=``).
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, NamedTuple

from viralscan.utils import INFECTED_RELATIVE_TOLERANCE

TIER_UNCALIBRATED = "uncalibrated"  # no calibration result was supplied
TIER_NONE = "no_tier"  # calibration ran and failed (no eligible grid point)
TIER_NOT_DETECTED = "not_detected"
TIER_CANDIDATE = "candidate"  # calibrated, but the probable cutoffs are not met
TIER_PROBABLE = "probable"
TIER_STRONG = "strong"
TIERS = (
    TIER_UNCALIBRATED,
    TIER_NONE,
    TIER_NOT_DETECTED,
    TIER_CANDIDATE,
    TIER_PROBABLE,
    TIER_STRONG,
)

STATUS_SELECTED = "selected"
STATUS_FAILED = "failed"

#: The one explicit mapping between the package tiers and the study claim levels (``adjudicate_candidates.py``
#: ``level``). An analogue by position on the ladder, not an equivalence: a package tier is a calibrated
#: molecule-count rule and carries no read-level evidence; the study level adds read-level evidence (virus-best
#: fraction, complex-body fraction, sibling bleed) on top. L3 adds a gene-programme model that has no package
#: counterpart, so it maps to the highest package tier. L4 and L5 are never emitted.
STUDY_LEVEL_TO_PACKAGE_TIER: dict[str, str] = {
    "L0": TIER_CANDIDATE,
    "L1": TIER_PROBABLE,
    "L2": TIER_STRONG,
    "L3": TIER_STRONG,
}

_GRID_TEXT = re.compile(
    r"molecules per cell in \{([^{}]*)\}.*?features per sample in \{([^{}]*)\}.*?"
    r"support fraction in \{([^{}]*)\}",
    re.DOTALL,
)


class Cutoffs(NamedTuple):
    """One grid point; tuple order is the contract's order."""

    molecules: int
    features: int
    support: float


@dataclass(frozen=True)
class Grid:
    molecules: tuple[int, ...]
    features: tuple[int, ...]
    support: tuple[float, ...]

    def points(self) -> list[Cutoffs]:
        return [
            Cutoffs(m, f, s) for m in self.molecules for f in self.features for s in self.support
        ]

    def __contains__(self, point: object) -> bool:
        return point in self.points()


@dataclass(frozen=True)
class GridPoint:
    """Training result for one grid point, computed by the caller."""

    cutoffs: Cutoffs
    f1: float
    se: float  # resampling-unit bootstrap SE of ``f1`` (calibration.tie_breaker_standard_error)
    false_calls_host_only: int | None  # probable-or-strong calls in host-only negatives
    false_calls_homology: int | None  # probable-or-strong calls in planted host-homology negatives


@dataclass(frozen=True)
class Calibration:
    status: str
    probable: Cutoffs | None = None
    strong: Cutoffs | None = None
    n_eligible: int = 0
    best_f1: float | None = None
    best_se: float | None = None

    def to_dict(self) -> dict[str, Any]:
        def cut(c: Cutoffs | None) -> dict[str, Any] | None:
            return None if c is None else c._asdict()

        return {
            "status": self.status,
            "probable": cut(self.probable),
            "strong": cut(self.strong),
            "n_eligible": self.n_eligible,
            "best_f1": self.best_f1,
            "best_se": self.best_se,
        }


@dataclass(frozen=True)
class TierEvidence:
    """Measured evidence for one (sample, virus). ``None`` = not measured, which can never reach a tier."""

    max_molecules_per_cell: (
        float  # supporting molecules in the best called cell (fractions kept, Q135)
    )
    n_features: int  # distinct supporting features in the sample
    support_fraction: float | None  # exact-fragment support fraction


def _values(text: str, kind: type) -> tuple[Any, ...]:
    out = tuple(kind(part.strip()) for part in text.split(","))
    if not out or list(out) != sorted(set(out)) or any(v <= 0 for v in out):
        raise ValueError(f"grid values must be positive and strictly increasing: {text!r}")
    return out


def grid_from_protocol(protocol: Mapping[str, Any]) -> Grid:
    """Parse ``calibration.threshold_search.grid``; refuse any text that is not exactly the three sets."""
    try:
        text = protocol["calibration"]["threshold_search"]["grid"]
    except (KeyError, TypeError) as exc:
        raise ValueError("protocol has no calibration.threshold_search.grid") from exc
    match = _GRID_TEXT.search(" ".join(str(text).split()))
    if match is None or len(re.findall(r"\{[^{}]*\}", str(text))) != 3:
        raise ValueError("calibration.threshold_search.grid is not the three declared sets")
    support = _values(match.group(3), float)
    if support[-1] > 1:
        raise ValueError("support fractions must lie in (0, 1]")
    return Grid(_values(match.group(1), int), _values(match.group(2), int), support)


def _check(point: Cutoffs, grid: Grid, what: str) -> Cutoffs:
    point = Cutoffs(*point)
    if point not in grid:
        raise ValueError(
            f"{what} {tuple(point)} is not a grid point; no value outside the grid may be adopted"
        )
    return point


def select_cutoffs(points: Iterable[GridPoint], grid: Grid) -> Calibration:
    """Apply the frozen selection rule to the training results of every grid point.

    Refuses (``ValueError``) a missing, duplicated or off-grid point, a non-finite F1 or SE, and a missing
    negative-class count: an unmeasured negative class is not a zero.
    """
    seen: dict[Cutoffs, GridPoint] = {}
    for p in points:
        c = _check(p.cutoffs, grid, "cutoff")
        if c in seen:
            raise ValueError(f"grid point {tuple(c)} given twice")
        if not (math.isfinite(p.f1) and 0 <= p.f1 <= 1 and math.isfinite(p.se) and p.se >= 0):
            raise ValueError(
                f"grid point {tuple(c)}: F1 must lie in [0, 1] and the SE be finite and >= 0"
            )
        if p.false_calls_host_only is None or p.false_calls_homology is None:
            raise ValueError(f"grid point {tuple(c)}: both negative-class call counts are required")
        if p.false_calls_host_only < 0 or p.false_calls_homology < 0:
            raise ValueError(f"grid point {tuple(c)}: negative call counts")
        seen[c] = p
    missing = [tuple(c) for c in grid.points() if c not in seen]
    if missing:
        raise ValueError(f"grid points without a training result: {missing[:3]}")
    eligible = [
        p for p in seen.values() if p.false_calls_host_only == 0 and p.false_calls_homology == 0
    ]
    if not eligible:
        return Calibration(STATUS_FAILED)
    best = max(
        eligible, key=lambda p: (p.f1, p.cutoffs)
    )  # tuple only fixes which of equal-F1 points is "best"
    within = [p for p in eligible if p.f1 >= best.f1 - best.se]
    # The most conservative point is the lexicographic maximum of the cutoff tuple. Grid points are distinct, so the
    # contract's final "ascending lexicographic" tie-break can never fire.
    chosen = max(within, key=lambda p: p.cutoffs)
    return Calibration(STATUS_SELECTED, chosen.cutoffs, None, len(eligible), best.f1, best.se)


def calibration_from_dict(result: Mapping[str, Any], grid: Grid) -> Calibration:
    """Validate a supplied calibration result: its points must be on the grid, ``strong`` nested above ``probable``."""
    status = result.get("status")
    if status == STATUS_FAILED:
        return Calibration(STATUS_FAILED)
    if status != STATUS_SELECTED:
        raise ValueError(
            f"calibration status must be {STATUS_SELECTED!r} or {STATUS_FAILED!r}, got {status!r}"
        )

    def cut(key: str) -> Cutoffs | None:
        raw = result.get(key)
        if raw is None:
            return None
        try:
            return _check(Cutoffs(raw["molecules"], raw["features"], raw["support"]), grid, key)
        except (KeyError, TypeError) as exc:
            raise ValueError(f"calibration {key} needs molecules, features and support") from exc

    probable, strong = cut("probable"), cut("strong")
    if probable is None:
        raise ValueError("a selected calibration needs a probable point")
    if strong is not None and any(s < p for s, p in zip(strong, probable, strict=True)):
        raise ValueError(
            "strong must be at least as strict as probable on every cutoff (nested tiers)"
        )
    return Calibration(
        STATUS_SELECTED,
        probable,
        strong,
        int(result.get("n_eligible", 0)),
        result.get("best_f1"),
        result.get("best_se"),
    )


def _meets(evidence: TierEvidence, c: Cutoffs) -> bool:
    return (
        evidence.support_fraction is not None
        and evidence.max_molecules_per_cell >= c.molecules * (1.0 - INFECTED_RELATIVE_TOLERANCE)
        and evidence.n_features >= c.features
        and evidence.support_fraction >= c.support
    )


def assign_tier(evidence: TierEvidence, calibration: Calibration | None) -> str:
    """The tier of one (sample, virus); probable and strong need a selected calibration, nothing else gives them."""
    if calibration is None:
        return TIER_UNCALIBRATED
    if calibration.status != STATUS_SELECTED or calibration.probable is None:
        return TIER_NONE
    if evidence.max_molecules_per_cell <= 0:
        return TIER_NOT_DETECTED
    if calibration.strong is not None and _meets(evidence, calibration.strong):
        return TIER_STRONG
    return TIER_PROBABLE if _meets(evidence, calibration.probable) else TIER_CANDIDATE


def package_tier_analogue(level: str) -> str:
    """The package tier a study claim level is the analogue of (``STUDY_LEVEL_TO_PACKAGE_TIER``)."""
    try:
        return STUDY_LEVEL_TO_PACKAGE_TIER[level]
    except KeyError:
        raise ValueError(f"no package-tier analogue for study level {level!r}") from None
