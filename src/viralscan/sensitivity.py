"""Detection sensitivity: what a negative ViralScan result does and does not mean.

The question this module exists to answer is the one every viral-detection
pipeline eventually has to face: *we called nothing — is there genuinely nothing
there, or did we not look hard enough?* Those two states are indistinguishable
from a zero, so a zero has to be reported together with the sampling limit that
produced it.

Three terms decide whether a virus that **is** present gets reported:

1. **Depth.** Molecules arrive as a thinning Poisson process, so a sample with
   :math:`D` molecules carrying true viral abundance :math:`a` (fraction of the
   library) yields :math:`\\lambda = D \\cdot a` expected viral molecules, and
   :math:`P(\\text{observe } 0) = e^{-\\lambda}`. The abundance resolved with
   95 % probability needs :math:`-\\ln(0.05) = 2.996` *expected observed*
   molecules, so the **depth-only LOD95** is :math:`2.996 / D` — about one
   viral molecule per :math:`D/3` host molecules. With capture :math:`c < 1`
   only a fraction :math:`c` of true molecules is observed, so the *required
   true* molecules are :math:`2.996 / c` (5.991 at :math:`c = 0.5`); the
   record keeps the two counts in separate fields.

2. **k-mer capture.** Pseudoalignment needs an *exact* :math:`k`-mer match
   (``k = 31``). A fragment of length :math:`L` at per-base divergence
   :math:`d` survives only if it holds a run of :math:`k` substitution-free
   bases. Overlapping windows share bases, so they are not independent; the
   exact probability (:func:`fragment_capture_exact`, a dynamic programme) is
   0.71 at 5 % divergence, 0.25 at 10 %, 0.063 at 15 %, 0.013 at 20 % and
   0.0003 at 30 % for 90 bp. This is a **substitution-only heuristic** (i.i.d.
   substitutions, one target): it ignores indels, sequencing error and
   competing index members, and real recall can sit well below it (Luebbert et
   al. 2025 Fig 1c, read from the figure, +/-3 points). It never counts as a
   measured capture. **Capture is the term that decides whether a negative is
   informative**, and it is invisible in a count matrix.

3. **Allocation survival.** The default ``host-conservative`` multimap method
   credits host-virus-ambiguous molecules *zero* to the virus
   (``multimapping.build_layers``). Measured on the bundled EBV LCL run that
   costs 1.03 % of viral molecules; against host endogenous viral elements it is
   the whole ballgame.

Only term 1 is measurable from inside a run. This module therefore reports the
depth-only LOD as a *floor*, labels it as such, and refuses to present a
negative as an absence. :func:`sensitivity_record`
carry an explicit measured ``capture`` when a caller has one (a spike-in
control, or read-level alignment evidence from ``viralscan evidence``).

Validated empirically by molecule-level downsampling of the bundled EBV LCL run
(55,266,624 quantified molecules, 1,636,934 of them EBV): EBV stayed at
P(detect)=1.0000 down to 1,270 downsampled reads and first reached 0 at 127,
consistent with the Poisson floor rather than with any capture loss.

Note the depth term is *quantified molecules*, not reads. For that run the two
differ by 1.9x (103,145,071 pseudoaligned reads vs 55,266,624 molecules), and
substituting reads for molecules would understate the LOD95 by the same factor.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

#: kallisto/bustools default and the value ``kb count`` used for every run in
#: this repo's benchmarks (see kb-python/run_info.json: "k-mer length": 31).
DEFAULT_K = 31

#: ln(1 / 0.95) = 2.9957. The expected molecule count at which a Poisson
#: process is observed with 95 % probability.
LOD95_MOLECULES = -math.log(0.05)

#: Depth is the term a run can measure; capture is the term it cannot. For
#: reference, these are the LOD95 values (viral UMI per 10k host UMI) the three
#: real covid PBMC configurations and a routine 10x run produced:
#:
#:   0.0005  bundled EBV LCL, 55,266,624 molecules    -> informative
#:           (its 103,145,071 *reads* would give 0.0003 — wrong unit)
#:   0.0014  covid, no host filter, 21,613,840        -> informative
#:   0.0036  covid, D-list masking, 8,404,326         -> informative
#:   0.0056  covid, STAR host filter, 5,308,302       -> informative
#:   0.0060  10x 5k cells x 1,000 UMI                 -> informative
#:   0.0300  10x 5k cells x 200 UMI (shallow)         -> adequate
#:
#: Every one of those is below ~1 viral molecule per cell in a standard 10x
#: run. That is the point: **depth is not what limits these runs.** The covid
#: samples called SARS-CoV-2 = 0 at 21.6M quantified molecules, where the depth
#: limit was ample — the reference and the host cross-talk were the limits. So
#: the depth bands are set to say "depth is sufficient", and the record's
#: ``informative_negative`` additionally requires a *measured* capture term.
DEFAULT_INFORMATIVE_LOD_PER_10K = 0.01
ADEQUATE_LOD_PER_10K = 0.1


def poisson_zero_probability(expected: float) -> float:
    """P(observe 0 molecules) when the expected count is ``expected``."""
    if expected <= 0:
        return 1.0
    return math.exp(-expected)


def poisson_detection_probability(expected: float, threshold: int) -> float:
    """P(observe at least ``threshold`` molecules) for a Poisson(``expected``)."""
    if expected <= 0:
        return 0.0 if threshold > 0 else 1.0
    if threshold <= 0:
        return 1.0
    # Sum the lower tail explicitly: 1 - CDF(threshold-1). Each term is
    # evaluated in log space: the naive recurrence anchored at exp(-expected)
    # underflows to 0.0 for expected >~ 745, and 0 * (expected/k) stays 0, so
    # every term would collapse and the tail would wrongly return 1.0.
    cdf = 0.0
    log_expected = math.log(expected)
    for k in range(threshold):
        cdf += math.exp(-expected + k * log_expected - math.lgamma(k + 1))
        if cdf >= 1.0:
            return 0.0
    return max(0.0, min(1.0, 1.0 - cdf))


def expected_viral_molecules(depth: float, abundance: float, capture: float = 1.0) -> float:
    """Expected observed viral molecules.

    Parameters
    ----------
    depth:
        Total molecules actually quantified in the run (sum of the count
        matrix), *not* raw reads. Only quantified molecules can be detected.
    abundance:
        True viral fraction of the library, in [0, 1].
    capture:
        Fraction of true viral molecules that survive k-mer capture
        (see :func:`fragment_capture_exact`). Defaults to 1.0, i.e. the optimistic
        assumption that the query matches the reference exactly.
    """
    if depth <= 0 or abundance <= 0 or capture <= 0:
        return 0.0
    return depth * abundance * capture


def lod95(
    depth: float,
    capture: float = 1.0,
    threshold: int = 1,
) -> float:
    """Lowest true abundance detectable with 95 % probability.

    Returns a *fraction of the library* (not molecules). Multiply by 1e4 for
    the conventional "viral UMI per 10k host UMI" unit.

    ``capture`` scales the requirement: a reference that only captures a
    fraction of the query needs proportionally more true molecules.

    The linear-in-``threshold`` form (``LOD95_MOLECULES * threshold``) is a
    conservative approximation: the exact Poisson LOD is sublinear in the
    threshold (~1.3-2x lower for large thresholds), so a reported LOD is never
    over-optimistic. It is exact at the default ``threshold=1``.
    """
    if depth <= 0 or capture <= 0:
        return math.inf
    return LOD95_MOLECULES * threshold / (depth * min(1.0, capture))


def lod95_per_10k(depth: float, capture: float = 1.0, threshold: int = 1) -> float:
    """LOD95 expressed as viral UMI per 10,000 host UMI."""
    value = lod95(depth, capture=capture, threshold=threshold)
    if math.isinf(value):
        return math.inf
    return value * 1e4


def minimum_molecules_for_detection(probability: float = 0.95, threshold: int = 1) -> float:
    """Expected molecule count giving ``probability`` of clearing ``threshold``."""
    if not 0.0 < probability < 1.0:
        raise ValueError(f"probability must be in (0, 1), got {probability!r}")
    if threshold <= 0:
        return 0.0
    return -math.log(1.0 - probability) * threshold


def fragment_capture(
    divergence: float,
    read_length: int = 90,
    k: int = DEFAULT_K,
) -> float:
    """Deprecated alias of :func:`fragment_capture_exact` (SENS-CORR-02).

    This used to multiply ``(1 - (1-d)^k)`` across the ``L-k+1`` overlapping
    windows as if independent, which overstated capture (0.90 at 10 % divergence
    where the exact value is 0.25). It now delegates to the exact
    substitution-only model, so every caller moves to the corrected curve.
    A **substitution-only heuristic**: it is never a measured capture and never
    sets ``capture_measured``.
    """
    return fragment_capture_exact(divergence, read_length=read_length, k=k)


def fragment_capture_exact(
    divergence: float,
    read_length: int = 90,
    k: int = DEFAULT_K,
) -> float:
    """P(at least one run of ``k`` substitution-free bases in ``read_length`` positions).

    The **substitution-only heuristic** behind every capture-from-divergence
    number in ViralScan. Exact under independent per-base substitutions at rate
    ``divergence`` against one target; not a bound on real capture (indels,
    sequencing error and index competition push it down, clustered divergence
    pushes it up), and never a measurement.

    Computed by a dynamic programme over positions whose state is the current
    clean-run length (0..k-1). Reaching run length ``k`` is absorbing, and its probability is
    accumulated directly rather than computed as ``1 - sum(dp)``, so tiny
    probabilities do not cancel to zero. O(read_length * k).

    A fragment shorter than ``k`` holds no k-mer and returns 0.0. Endpoint
    ``divergence == 1`` is not supported.
    """
    if not (isinstance(divergence, (int, float)) and 0.0 <= divergence < 1.0):
        raise ValueError(f"divergence must be a finite number in [0, 1), got {divergence!r}")
    if k < 1 or read_length < 0:
        raise ValueError(
            f"need k >= 1 and read_length >= 0, got k={k!r}, read_length={read_length!r}"
        )
    if read_length < k:
        return 0.0
    clean = 1.0 - divergence
    run = [0.0] * k  # run[r] = P(current clean run is r, k not yet reached)
    run[0] = 1.0
    found = 0.0
    for _ in range(read_length):
        found += clean * run[-1]
        run = [divergence * sum(run), *(clean * r for r in run[:-1])]
    return min(1.0, found)


class CaptureScope(str, Enum):
    """What a measured positive-control capture is entitled to speak for."""

    EXACT_SEQUENCE = "exact_sequence"  # one pinned accession/sequence only
    VIRUS_KEY = "virus_key"  # one virus group, given approved transfer calibration
    PANEL_MECHANICS = "panel_mechanics"  # pipeline recovery check; certifies no virus


@dataclass(frozen=True)
class PositiveControl:
    """A measured control: ``capture`` (observed/expected), its ``scope`` and ``target``.

    ``target`` is the exact-sequence ID or the virus key the control was measured
    on; callers resolve it through the virus identity table, never by display name.
    """

    capture: float | None
    scope: CaptureScope = CaptureScope.PANEL_MECHANICS
    target: str | None = None


def scoped_capture(
    control: PositiveControl | None,
    virus_key: str,
    scope: CaptureScope,
) -> float | None:
    """Capture to apply to ``virus_key`` under claim ``scope``, or None if out of scope.

    None means "no capture term for this row" -- the caller reports the
    unmeasured depth-only floor and must not certify a negative, rather than
    borrowing 1.0 or another virus's value. ``PANEL_MECHANICS`` controls and
    controls whose scope or target differ from the row's return None. A
    measurement that is missing, non-finite, <= 0 or over-recovered (> 1, so not
    a clean bound on loss) is also None.
    """
    if control is None or control.scope is not scope or scope is CaptureScope.PANEL_MECHANICS:
        return None
    if control.target is None or control.target != virus_key:
        return None
    c = control.capture
    if c is None or not math.isfinite(c) or not 0.0 < c <= 1.0:
        return None
    return float(c)


def capture_cliff_text(read_length: int = 90, k: int = DEFAULT_K) -> str:
    """Prose for the divergence cliff, computed from the exact model so it cannot go stale."""
    pts = ", ".join(
        f"{fragment_capture_exact(d, read_length, k):.2g} at {d:.0%}" for d in (0.10, 0.15, 0.20)
    )
    return (
        f"substitution-only heuristic, {read_length} bp, k={k}: {pts} divergence; "
        "it ignores indels, sequencing error and index competition, so it is not a measurement"
    )


#: Interpretation bands for a run's LOD95, as viral UMI per 10k host UMI.
#: Chosen so "informative" means the run can resolve a burden small enough to
#: be biologically plausible, not merely non-zero. See
#: :data:`DEFAULT_INFORMATIVE_LOD_PER_10K` for the measured reference values.
LOD_BANDS: tuple[tuple[float, str], ...] = (
    (DEFAULT_INFORMATIVE_LOD_PER_10K, "informative"),
    (ADEQUATE_LOD_PER_10K, "adequate"),
    (1.0, "shallow"),
    (math.inf, "insufficient-depth"),
)


def classify_lod(lod_per_10k: float) -> str:
    """Bucket a LOD95 (viral UMI per 10k host UMI) into an interpretation band."""
    for ceiling, label in LOD_BANDS:
        if lod_per_10k <= ceiling:
            return label
    return "insufficient-depth"


#: Why a row's capture term is, or is not, a measurement. One vocabulary shared by
#: the TSV (``capture_status``), the prose and the HTML, so none can invent its own.
CAPTURE_MEASURED = "measured"
CAPTURE_NOT_MEASURED = "not-measured"  # no control configured (or no reason given)

#: Reasons a row is not an informative negative (``negative_blockers``).
BLOCKER_DETECTED = "detected"  # a call, not a negative
BLOCKER_NO_CAPTURE = "capture-not-measured"
BLOCKER_DEPTH = "depth-insufficient"


@dataclass(frozen=True)
class SensitivityRecord:
    """Canonical per-virus sensitivity interpretation for one run.

    Every renderer (``sensitivity.tsv``, ``positive_control.json``, ``summary.txt``
    and ``report.html``) reads these fields; none recomputes them (SENS-CORR-03).

    Two molecule counts are kept apart because they differ by the capture term:

    * **required true molecules** (``lod95_molecules``): true viral molecules that
      must be present for a 95 % chance of clearing the gate,
      ``LOD95_MOLECULES * threshold / capture``. For ``depth > 0`` this equals
      ``depth * lod95_per_10k / 1e4``.
    * **expected observed molecules** (``expected_observed_molecules_at_lod95``):
      what those true molecules yield after capture, ``LOD95_MOLECULES *
      threshold``. This is the Poisson mean the 95 % statement is about.

    At capture 0.5 and threshold 1 that is 5.991 required versus 2.996 expected.

    Attributes
    ----------
    virus_name:
        Virus the record describes.
    observed_molecules:
        Molecules attributed to the virus in this run, as reported (for an
        undetected row, the count below the gate when the caller supplied it).
    detection_threshold:
        Sample-level UMI gate that decided the call.
    detected:
        The actual detection decision (the virus cleared the gate). It is passed
        in by the pipeline; it is not inferred from ``observed_molecules``, so a
        nonzero count below the gate is a negative with a nonzero count.
    p_zero_at_lod95:
        P(observing 0 molecules) for a virus sitting exactly at this run's
        LOD95. A negative at this abundance is expected 5 % of the time, so
        ~0.05 by construction; the value is reported to make the arithmetic
        checkable rather than to be interpreted.
    expected_molecules_if_present_at_1_per_10k:
        Expected observed molecules for a virus at 1 UMI per 10k host UMI.
        The natural "is this a real infection or noise" yardstick.
    p_detect_if_present_at_1_per_10k:
        P(reaching the detection threshold) at that abundance.
    lod95_per_10k:
        LOD95 in true viral UMI per 10k host UMI, adjusted by the *effective*
        capture. ``inf`` if depth is 0. Without a measured capture the effective
        capture is 1.0, so this is the depth-only floor, not an estimate.
    lod95_molecules:
        Required **true** molecules at the LOD95 (see above). Not capped at ~3.
    expected_observed_molecules_at_lod95:
        Expected **observed** molecules at the LOD95: ``LOD95_MOLECULES *
        threshold``, independent of capture.
    capture_measured:
        True only when ``capture`` is a valid measurement whose scope covers this
        row. False for rows outside a control's scope and for failed or
        over-recovered controls.
    capture_status:
        Why: ``measured``, ``not-measured``, or a caller-supplied reason such as
        ``control-failed``, ``control-over-recovered``, ``control-out-of-scope``.
    depth_only_sufficient:
        The LOD95 at capture 1.0 is inside the "informative" band: depth alone.
    depth_sufficient:
        The LOD95 *after* the effective capture is inside the "informative" band.
        Legacy column, kept with this definition; it equals ``depth_only_sufficient``
        for any row without a measured capture.
    sensitivity_eligible:
        ``depth_sufficient`` **and** ``capture_measured``: an adequate adjusted
        sensitivity estimate exists for this row. Says nothing about the result.
    informative:
        The informative **negative** (``informative_negative`` column):
        ``sensitivity_eligible`` **and not** ``detected``. A detected positive is
        never an informative negative. This is the covid-run lesson: ample depth
        with no capture term established does not certify an absence.
    negative_blockers:
        Codes explaining why ``informative`` is false (empty when it is true).
    """

    virus_name: str
    observed_molecules: float
    detection_threshold: int
    detected: bool
    p_zero_at_lod95: float
    expected_molecules_if_present_at_1_per_10k: float
    p_detect_if_present_at_1_per_10k: float
    lod95_per_10k: float
    lod95_molecules: float
    expected_observed_molecules_at_lod95: float
    depth_only_sufficient: bool
    depth_sufficient: bool
    capture_measured: bool
    capture_status: str
    sensitivity_eligible: bool
    informative: bool
    capture: float
    lod_interpretation: str
    negative_blockers: tuple[str, ...] = field(default_factory=tuple)
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_row(self) -> dict[str, Any]:
        """Flatten to a TSV row (finite floats, or None where nothing was measured).

        ``capture`` is None unless it was measured: the internal 1.0 on an
        uncertified row is the depth-only floor, and printing it would read as a
        measurement.
        """
        return {
            "virus_name": self.virus_name,
            "observed_molecules": self.observed_molecules,
            "detection_threshold": self.detection_threshold,
            "detected": self.detected,
            "capture": self.capture if self.capture_measured else None,
            "capture_measured": self.capture_measured,
            "capture_status": self.capture_status,
            "lod95_per_10k": self.lod95_per_10k,
            "lod95_molecules": self.lod95_molecules,
            "lod95_expected_observed_molecules": self.expected_observed_molecules_at_lod95,
            "lod_interpretation": self.lod_interpretation,
            "depth_only_sufficient": self.depth_only_sufficient,
            "depth_sufficient": self.depth_sufficient,
            "sensitivity_eligible": self.sensitivity_eligible,
            "informative_negative": self.informative,
            "negative_blockers": ";".join(self.negative_blockers),
            "expected_molecules_at_1_per_10k": (self.expected_molecules_if_present_at_1_per_10k),
            "p_detect_at_1_per_10k": self.p_detect_if_present_at_1_per_10k,
            "p_zero_at_lod95": self.p_zero_at_lod95,
            "notes": "; ".join(self.notes),
        }


SENSITIVITY_COLUMNS: tuple[str, ...] = (
    "virus_name",
    "observed_molecules",
    "detection_threshold",
    "detected",
    "capture",
    "capture_measured",
    "capture_status",
    "lod95_per_10k",
    "lod95_molecules",
    "lod95_expected_observed_molecules",
    "lod_interpretation",
    "depth_only_sufficient",
    "depth_sufficient",
    "sensitivity_eligible",
    "informative_negative",
    "negative_blockers",
    "expected_molecules_at_1_per_10k",
    "p_detect_at_1_per_10k",
    "p_zero_at_lod95",
    "notes",
)


#: Abundance the per-virus yardstick is expressed at, as a fraction.
YARDSTICK_ABUNDANCE = 1e-4  # 1 UMI per 10k host UMI


def sensitivity_record(
    virus_name: str,
    observed_molecules: float,
    depth: float,
    detection_threshold: int = 1,
    capture: float = 1.0,
    capture_measured: bool = False,
    abundance_ceiling: float = DEFAULT_INFORMATIVE_LOD_PER_10K,
    notes: tuple[str, ...] = (),
    detected: bool | None = None,
    capture_status: str | None = None,
) -> SensitivityRecord:
    """Build a :class:`SensitivityRecord` for one virus in one run.

    ``depth`` is total quantified molecules. ``capture`` defaults to 1.0, which
    is the *optimistic* assumption and leaves ``capture_measured`` False; pass
    both when a control or read-level evidence has established a real capture
    term. A record built with the default always carries a note saying so.

    ``detected`` is the pipeline's actual detection decision. Left as None it is
    inferred as ``observed_molecules >= detection_threshold``, which is only
    right for callers that have no better information. ``capture_status``
    records *why* capture is not measured (see :class:`SensitivityRecord`).
    """
    if detection_threshold < 1:
        raise ValueError(f"detection_threshold must be >= 1, got {detection_threshold!r}")
    if capture <= 0:
        raise ValueError(
            f"capture must be > 0, got {capture!r}. A zero capture term means the "
            "reference cannot see the query at all; pass capture=None by leaving "
            "capture_measured=False so the record reports the optimistic bound "
            "plus the reason it is unmeasured."
        )
    if detected is None:
        detected = observed_molecules >= detection_threshold
    if capture_status is None:
        capture_status = CAPTURE_MEASURED if capture_measured else CAPTURE_NOT_MEASURED

    # A capture above 1 means the control recovered more than was planted, so
    # detection is at least as easy as the depth floor; clamp for the LOD rather
    # than letting it tighten the limit below the sampling floor.
    effective_capture = min(1.0, capture)
    lod_10k = lod95_per_10k(depth, capture=effective_capture, threshold=detection_threshold)
    expected_at_yardstick = expected_viral_molecules(depth, YARDSTICK_ABUNDANCE, effective_capture)
    depth_sufficient = lod_10k <= abundance_ceiling
    depth_only_sufficient = (
        lod95_per_10k(depth, capture=1.0, threshold=detection_threshold) <= abundance_ceiling
    )
    eligible = depth_sufficient and capture_measured
    informative = eligible and not detected
    blockers = tuple(
        code
        for code, active in (
            (BLOCKER_DETECTED, detected),
            (BLOCKER_NO_CAPTURE, not capture_measured),
            (BLOCKER_DEPTH, not depth_sufficient),
        )
        if active
    )
    expected_observed = LOD95_MOLECULES * detection_threshold
    all_notes = list(notes)
    if not capture_measured:
        all_notes.append(
            "capture assumed 1.0 (exact reference match) and NOT measured; the "
            "LOD95 shown is a floor — a divergent query needs proportionally more "
            f"molecules ({capture_cliff_text()})"
        )
    elif capture >= 1.0:
        all_notes.append(
            f"measured capture {capture:.4g} >= 1: no capture loss was measurable "
            "at this abundance, so the LOD95 equals the depth floor. This bounds "
            "loss from above only; it does not demonstrate zero divergence"
        )
    if detected:
        all_notes.append(
            "detected: this row is a call, so it is not a negative and "
            "informative_negative is false by definition"
        )
    else:
        if observed_molecules > 0:
            all_notes.append(
                f"{observed_molecules:g} molecule(s) observed but below the detection "
                "gate: not called, and not a true absence either"
            )
        if capture_measured and depth_sufficient:
            tail = (
                "depth is sufficient and capture was measured in scope, so this "
                "negative is certifiable for that scope only"
            )
        elif depth_sufficient:
            tail = (
                "depth is sufficient but capture is unestablished, so this negative "
                "is not certifiable"
            )
        else:
            tail = "depth is insufficient for this row, so this negative is not certifiable"
        all_notes.append(
            ("zero observed: " if observed_molecules == 0 else "not detected: ")
            + "a negative, not an absence; "
            + tail
        )
    return SensitivityRecord(
        virus_name=virus_name,
        observed_molecules=float(observed_molecules),
        detection_threshold=int(detection_threshold),
        detected=bool(detected),
        p_zero_at_lod95=poisson_zero_probability(LOD95_MOLECULES),
        expected_molecules_if_present_at_1_per_10k=expected_at_yardstick,
        p_detect_if_present_at_1_per_10k=poisson_detection_probability(
            expected_at_yardstick, detection_threshold
        ),
        lod95_per_10k=lod_10k,
        lod95_molecules=expected_observed / effective_capture,
        expected_observed_molecules_at_lod95=expected_observed,
        depth_only_sufficient=depth_only_sufficient,
        depth_sufficient=depth_sufficient,
        capture_measured=capture_measured,
        capture_status=capture_status,
        sensitivity_eligible=eligible,
        informative=informative,
        capture=capture,
        lod_interpretation=classify_lod(lod_10k),
        negative_blockers=blockers,
        notes=tuple(all_notes),
    )


_NOT_ESTIMABLE = (
    "Detection limit: NOT ESTIMABLE — this run quantified zero molecules, so "
    "no negative can be interpreted. A 'no virus found' result here is a "
    "quantification failure, not a biological finding."
)

_ASSUMPTIONS = (
    "under the stated assumptions (Poisson sampling of quantified molecules in the "
    "selected count layer, the stated capture, a fixed detection gate)"
)


def _limit_sentences(
    depth: float,
    detection_threshold: int,
    capture: float,
    capture_measured: bool,
    abundance_ceiling: float,
    per_target: bool = False,
) -> tuple[str, str]:
    """The run-level limit and depth-verdict sentences, shared by both statements.

    ``per_target`` words the headline as the depth-only floor that rows without a
    measured capture get, with measured rows qualified separately, so it cannot be
    read as contradicting a per-target capture stated later.
    """
    effective = min(1.0, capture)
    lod_10k = lod95_per_10k(depth, capture=effective, threshold=detection_threshold)
    expected = minimum_molecules_for_detection(0.95, detection_threshold)
    required = expected / effective
    band = classify_lod(lod_10k)
    if per_target:
        capture_text = "depth-only floor, capture ASSUMED 1.0 unless a target below says otherwise"
    elif capture_measured:
        capture_text = f"k-mer capture measured at {capture:.3g}"
    else:
        capture_text = "k-mer capture ASSUMED 1.0, not measured"
    head = (
        f"Detection limit for this run (LOD95, k={DEFAULT_K}, {capture_text}): "
        f"{lod_10k:.4g} viral UMI per 10k host UMI, i.e. about {required:.4g} true "
        f"viral molecules ({expected:.4g} expected observed after capture) at the "
        f"detection threshold. Quantified depth: {depth:,.0f} molecules."
    )
    if lod_10k <= abundance_ceiling:
        verdict = (
            f"Depth is sufficient to resolve a {abundance_ceiling:g}/10k burden "
            f"(band: {band}): a virus present at or above the limit is reported with "
            f"at least 95 % probability {_ASSUMPTIONS}; that is a probability, not a guarantee."
        )
    else:
        verdict = (
            f"Depth is NOT sufficient (band: {band}, limit {lod_10k:.4g}/10k exceeds "
            f"the {abundance_ceiling:g}/10k sufficiency ceiling): a virus below the "
            "limit is missed with >=5% probability. Deepen the library."
        )
    return head, verdict


_NO_CAPTURE_TAIL = (
    "No capture term has been established for this query, so a negative "
    "CANNOT be read as absence regardless of depth: capture falls steeply "
    "with divergence ({cliff}), and the reference panel — "
    "not the sequencing depth — is then the binding limit. Run a positive control "
    "at known abundance, or align the reads to the target directly "
    "(`viralscan evidence`), before reporting a negative."
)


def negative_result_statement(
    depth: float,
    detection_threshold: int = 1,
    capture: float = 1.0,
    capture_measured: bool = False,
    abundance_ceiling: float = DEFAULT_INFORMATIVE_LOD_PER_10K,
    certified_targets: tuple[str, ...] | None = None,
) -> str:
    """Human-readable caveat for a *single* capture term applied run-wide.

    Prefer :func:`sensitivity_statement`, which renders per-target from the
    :class:`SensitivityRecord` list and so cannot describe a mixed-scope run as
    uniformly certifiable. This scalar form remains for callers that have one
    capture value and no records. ``certified_targets`` narrows the "certifiable"
    sentence to exactly those targets; ``None`` keeps the unscoped wording.

    The depth verdict uses the same rule as
    :attr:`SensitivityRecord.depth_sufficient`, so the prose and the
    ``sensitivity.tsv`` column can never disagree about the same run.
    """
    if depth <= 0:
        return _NOT_ESTIMABLE
    head, depth_verdict = _limit_sentences(
        depth, detection_threshold, capture, capture_measured, abundance_ceiling
    )
    if capture_measured and certified_targets is not None:
        tail = (
            "A k-mer capture term was measured for "
            f"{', '.join(certified_targets)} only; negatives for those targets are "
            "certifiable at that capture, and no other virus is covered by it."
        )
    elif capture_measured:
        tail = (
            "A measured capture term exists, so this negative is certifiable at "
            "the stated divergence."
        )
    else:
        tail = _NO_CAPTURE_TAIL.format(cliff=capture_cliff_text())
    return f"{head} {depth_verdict} {tail}"


def sensitivity_statement(
    records: Sequence[SensitivityRecord],
    depth: float,
    detection_threshold: int = 1,
    abundance_ceiling: float = DEFAULT_INFORMATIVE_LOD_PER_10K,
) -> str:
    """Run-level interpretation rendered from the per-target records only.

    The headline limit is the depth-only floor (capture 1.0, not measured), because
    that is what every row without a measured capture gets. Each target whose capture
    *was* measured in scope is then qualified with its own adjusted limit, so a
    mixed-scope panel never reads as uniformly certifiable. Called viruses are
    listed as calls, below-gate nonzero counts are reported as such, and every other
    undetected row is summarised by the reason it is not an informative negative.
    This is the single text renderer behind ``summary.txt`` and ``report.html``.
    """
    if depth <= 0:
        return _NOT_ESTIMABLE
    head, verdict = _limit_sentences(
        depth, detection_threshold, 1.0, False, abundance_ceiling, per_target=True
    )
    parts = [head, verdict]

    called = sorted(r.virus_name for r in records if r.detected)
    if called:
        parts.append(
            f"Called by the detection gate: {', '.join(called)}. A call says nothing "
            "about any virus that was not called."
        )
    informative = sorted((r for r in records if r.informative), key=lambda r: r.virus_name)
    for r in informative:
        parts.append(
            f"Informative negative for {r.virus_name} only: capture measured at "
            f"{r.capture:.3g} on a positive control in its scope, so reporting it needs "
            f"about {r.lod95_molecules:.4g} true molecules "
            f"({r.expected_observed_molecules_at_lod95:.4g} expected observed), "
            f"LOD95 {r.lod95_per_10k:.4g}/10k; a virus at that burden is reported with at "
            f"least 95 % probability under the same assumptions. No other virus, and no untested "
            "genome of this one, is covered."
        )
    other = [r for r in records if not r.detected and not r.informative]
    if other:
        reasons: dict[str, int] = {}
        for r in other:
            key = r.capture_status + (", depth-insufficient" if not r.depth_sufficient else "")
            reasons[key] = reasons.get(key, 0) + 1
        listing = "; ".join(f"{n} x {k}" for k, n in sorted(reasons.items()))
        text = (
            f"{len(other)} undetected virus(es) have no informative negative ({listing}): "
            "their zeros are sampling statements, not evidence of absence."
        )
        if any(not r.capture_measured for r in other):
            text += " " + _NO_CAPTURE_TAIL.format(cliff=capture_cliff_text())
        parts.append(text)
    below = sorted(
        (r for r in records if not r.detected and r.observed_molecules > 0),
        key=lambda r: -r.observed_molecules,
    )
    if below:
        shown = ", ".join(f"{r.virus_name} ({r.observed_molecules:g})" for r in below[:5])
        more = f" and {len(below) - 5} more" if len(below) > 5 else ""
        parts.append(f"Nonzero counts below the detection gate (not called): {shown}{more}.")
    return " ".join(parts)
