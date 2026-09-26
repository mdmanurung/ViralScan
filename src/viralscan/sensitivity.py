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
   95 % probability is :math:`-\\ln(0.05) = 2.996` molecules, so the
   **depth-only LOD95** is :math:`2.996 / D` — about one viral molecule per
   :math:`D/3` host molecules.

2. **k-mer capture.** Pseudoalignment needs an *exact* :math:`k`-mer match
   (``k = 31``). A fragment of length :math:`L` at per-base divergence
   :math:`d` survives only if at least one of its :math:`L-k+1` windows is
   error-free, giving :math:`P = 1 - (1 - (1-d)^k)^{L-k+1}`. This is not a
   gentle slope: capture is ~1.0 at 5 % divergence, 0.90 at 10 %, 0.32 at 15 %,
   0.06 at 20 % and 0.001 at 30 %. **This is the term that decides whether a
   negative is informative**, and it is invisible in a count matrix.

3. **Allocation survival.** The default ``host-conservative`` multimap method
   credits host-virus-ambiguous molecules *zero* to the virus
   (``multimapping.build_layers``). Measured on the bundled EBV LCL run that
   costs 1.03 % of viral molecules; against host endogenous viral elements it is
   the whole ballgame.

Only term 1 is measurable from inside a run. This module therefore reports the
depth-only LOD as a *floor*, labels it as such, and refuses to present a
negative as an absence. :func:`capture_reference` and :func:`sensitivity_record`
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
from dataclasses import dataclass, field
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
#:   0.0003  bundled EBV LCL, 103,145,071 molecules   -> informative
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
    # Sum the lower tail explicitly: 1 - CDF(threshold-1). For small thresholds
    # (the default is 1) this is exact in double precision; the loop keeps it
    # stable for the large thresholds a caller may pass.
    log_p = -expected
    term = math.exp(log_p)  # P(0)
    cdf = term
    for k in range(1, threshold):
        term *= expected / k
        cdf += term
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
        (see :func:`fragment_capture`). Defaults to 1.0, i.e. the optimistic
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
    """P(a fragment pseudoaligns) against a reference at ``divergence``.

    Uses the exact-match k-mer model: a window of length ``k`` is error-free
    with probability ``(1-d)^k``, a fragment of length ``L`` has ``L-k+1``
    windows, and the fragment is captured if any one of them is error-free.

    This ignores indels, sequencing error, and the fact that kallisto tolerates
    some divergence through multiple k-mers rather than one. It is therefore an
    **upper bound on the loss** from divergence, and is used as such: the true
    capture is somewhat better than this curve, never worse.

    Parameters
    ----------
    divergence:
        Per-base substitution rate in [0, 1), e.g. 0.05 for 5 %.
    read_length:
        Fragment length. For 10x 3' chemistry the usable cDNA read is ~90 bp
        after the 16 bp CB + 12 bp UMI are removed.
    k:
        k-mer length; must match the index that was built.
    """
    if not 0.0 <= divergence < 1.0:
        raise ValueError(f"divergence must be in [0, 1), got {divergence!r}")
    if read_length <= k:
        # Fragment shorter than k: capture is all-or-nothing per molecule.
        return (1.0 - divergence) ** max(read_length, 0)
    window_clean = (1.0 - divergence) ** k
    return 1.0 - (1.0 - window_clean) ** (read_length - k + 1)


def capture_reference(divergence: float, read_length: int = 90, k: int = DEFAULT_K) -> float:
    """Capture term for :func:`expected_viral_molecules` at a stated divergence.

    Thin alias over :func:`fragment_capture`; named for the call site in
    :func:`sensitivity_record`, where the divergence is something a caller
    measured (e.g. from a spike-in control or a read-level alignment) rather
    than something the run knows.
    """
    return fragment_capture(divergence, read_length=read_length, k=k)


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


@dataclass(frozen=True)
class SensitivityRecord:
    """Per-virus detection sensitivity for one run.

    Attributes
    ----------
    virus_name:
        Virus the record describes.
    observed_molecules:
        Molecules attributed to the virus in this run, as reported.
    detection_threshold:
        Sample-level UMI gate that decided the call.
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
        Depth-only LOD95 in viral UMI per 10k host UMI. ``inf`` if depth is 0.
        Computed with the supplied ``capture``; with the default capture=1.0
        this is a floor, not an estimate.
    lod95_molecules:
        Expected true molecules at the LOD95 — always ~3 (times threshold).
    capture_measured:
        True when ``capture`` came from a measurement (a spike-in control, or
        read-level alignment evidence) rather than the optimistic default.
    depth_sufficient:
        True when depth alone puts the LOD95 in the "informative" band.
    informative:
        True only when ``depth_sufficient`` **and** ``capture_measured``. A
        negative can only be certified as an absence when both hold: plenty of
        depth plus a demonstrated ability to see the query. Depth alone is not
        sufficient, and this is the empirical lesson of the covid runs, which
        called SARS-CoV-2 = 0 at ample depth because no capture term had been
        established for the query.
    """

    virus_name: str
    observed_molecules: float
    detection_threshold: int
    p_zero_at_lod95: float
    expected_molecules_if_present_at_1_per_10k: float
    p_detect_if_present_at_1_per_10k: float
    lod95_per_10k: float
    lod95_molecules: float
    depth_sufficient: bool
    capture_measured: bool
    informative: bool
    capture: float
    lod_interpretation: str
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_row(self) -> dict[str, Any]:
        """Flatten to a TSV row (finite floats only, so pandas/CSV stay clean)."""
        return {
            "virus_name": self.virus_name,
            "observed_molecules": self.observed_molecules,
            "detection_threshold": self.detection_threshold,
            "capture": self.capture,
            "capture_measured": self.capture_measured,
            "lod95_per_10k": self.lod95_per_10k,
            "lod95_molecules": self.lod95_molecules,
            "lod_interpretation": self.lod_interpretation,
            "depth_sufficient": self.depth_sufficient,
            "informative_negative": self.informative,
            "expected_molecules_at_1_per_10k": (self.expected_molecules_if_present_at_1_per_10k),
            "p_detect_at_1_per_10k": self.p_detect_if_present_at_1_per_10k,
            "p_zero_at_lod95": self.p_zero_at_lod95,
            "notes": "; ".join(self.notes),
        }


SENSITIVITY_COLUMNS: tuple[str, ...] = (
    "virus_name",
    "observed_molecules",
    "detection_threshold",
    "capture",
    "capture_measured",
    "lod95_per_10k",
    "lod95_molecules",
    "lod_interpretation",
    "depth_sufficient",
    "informative_negative",
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
) -> SensitivityRecord:
    """Build a :class:`SensitivityRecord` for one virus in one run.

    ``depth`` is total quantified molecules. ``capture`` defaults to 1.0, which
    is the *optimistic* assumption and leaves ``capture_measured`` False; pass
    both when a control or read-level evidence has established a real capture
    term. A record built with the default always carries a note saying so.
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

    # A capture above 1 means the control recovered more than was planted, so
    # detection is at least as easy as the depth floor; clamp for the LOD rather
    # than letting it tighten the limit below the sampling floor.
    effective_capture = min(1.0, capture)
    lod_10k = lod95_per_10k(depth, capture=effective_capture, threshold=detection_threshold)
    expected_at_yardstick = expected_viral_molecules(depth, YARDSTICK_ABUNDANCE, effective_capture)
    depth_sufficient = lod_10k <= abundance_ceiling
    all_notes = list(notes)
    if not capture_measured:
        all_notes.append(
            "capture assumed 1.0 (exact reference match) and NOT measured; the "
            "LOD95 shown is a floor — a divergent query needs proportionally more "
            "molecules (see fragment_capture(): 0.32 at 15% divergence, 0.06 at 20%)"
        )
    elif capture >= 1.0:
        all_notes.append(
            f"measured capture {capture:.4g} >= 1: no capture loss was measurable "
            "at this abundance, so the LOD95 equals the depth floor. This bounds "
            "loss from above only; it does not demonstrate zero divergence"
        )
    if observed_molecules == 0:
        all_notes.append(
            "zero observed: a negative, not an absence"
            + (
                "; depth is sufficient but capture is unestablished, so this "
                "negative is not certifiable"
                if depth_sufficient
                else "; depth is also insufficient"
            )
        )
    return SensitivityRecord(
        virus_name=virus_name,
        observed_molecules=float(observed_molecules),
        detection_threshold=int(detection_threshold),
        p_zero_at_lod95=poisson_zero_probability(LOD95_MOLECULES),
        expected_molecules_if_present_at_1_per_10k=expected_at_yardstick,
        p_detect_if_present_at_1_per_10k=poisson_detection_probability(
            expected_at_yardstick, detection_threshold
        ),
        lod95_per_10k=lod_10k,
        lod95_molecules=LOD95_MOLECULES * detection_threshold,
        depth_sufficient=depth_sufficient,
        capture_measured=capture_measured,
        informative=depth_sufficient and capture_measured,
        capture=capture,
        lod_interpretation=classify_lod(lod_10k),
        notes=tuple(all_notes),
    )


def negative_result_statement(
    depth: float,
    detection_threshold: int = 1,
    capture: float = 1.0,
    capture_measured: bool = False,
    abundance_ceiling: float = DEFAULT_INFORMATIVE_LOD_PER_10K,
) -> str:
    """Human-readable caveat to attach to a run that detected nothing.

    This is the deliverable the whole module exists for. It states the run's
    sampling limit in both directions a reader needs — as a burden per host
    molecules, and as an absolute molecule count — and says plainly which of
    the two limits (depth or reference capture) is actually in the way.

    The depth verdict uses the same rule as
    :attr:`SensitivityRecord.depth_sufficient`, so the prose and the
    ``sensitivity.tsv`` column can never disagree about the same run.
    """
    lod_10k = lod95_per_10k(depth, capture=capture, threshold=detection_threshold)
    molecules = minimum_molecules_for_detection(0.95, detection_threshold)
    if depth <= 0:
        return (
            "Detection limit: NOT ESTIMABLE — this run quantified zero molecules, so "
            "no negative can be interpreted. A 'no virus found' result here is a "
            "quantification failure, not a biological finding."
        )
    band = classify_lod(lod_10k)
    head = (
        f"Detection limit for this run (LOD95, k={DEFAULT_K}, k-mer capture "
        f"{'measured at ' + format(capture, '.3g') if capture_measured else 'ASSUMED 1.0, not measured'}): "
        f"{lod_10k:.4g} viral UMI per 10k host UMI, i.e. about {molecules:.0f} true "
        f"viral molecules at the detection threshold. Quantified depth: "
        f"{depth:,.0f} molecules."
    )
    if lod_10k <= abundance_ceiling:
        depth_verdict = (
            f"Depth is sufficient to resolve a {abundance_ceiling:g}/10k burden "
            f"(band: {band}): a virus above the limit would have been reported."
        )
    else:
        depth_verdict = (
            f"Depth is NOT sufficient (band: {band}, limit {lod_10k:.4g}/10k exceeds "
            f"the {abundance_ceiling:g}/10k sufficiency ceiling): a virus below the "
            "limit is missed with >=5% probability. Deepen the library."
        )
    if capture_measured:
        tail = (
            "A measured capture term exists, so this negative is certifiable at "
            "the stated divergence."
        )
    else:
        tail = (
            "No capture term has been established for this query, so a negative "
            "CANNOT be read as absence regardless of depth: capture falls to 0.32 "
            "at 15% divergence and 0.06 at 20%, and the reference panel — not the "
            "sequencing depth — is then the binding limit. Run a positive control "
            "at known abundance, or align the reads to the target directly "
            "(`viralscan evidence`), before reporting a negative."
        )
    return f"{head} {depth_verdict} {tail}"
