"""Advise the ``kb count`` parameters ``-x``, ``-w`` and ``--strand`` from the reads (PLAN KBPARAMS-01).

The function :func:`recommend_kb_parameters` combines up to three independent signals:

* the built-in chemistry detector (:mod:`viralscan.chemistry`),
* the ``tchem`` second opinion (:mod:`viralscan.tchem`), and
* optional strand pilot rates (:mod:`viralscan.strand`, ``p_pseudoaligned`` per strand).

A parameter is recommended only when its signals agree. Anything else is reported
with a reason and the run keeps asking for an explicit flag. The function is advisory:
it changes no default. ``run --kb-params auto`` uses it and stops on anything but ``accepted``. Pilot rates must come from the caller, because the
pilot needs a kallisto index.

``status`` per recommendation:

``accepted``    every available signal agrees and the choice cannot lose reads silently
``review``      a value is suggested, but a signal is missing or a limit applies (see ``notes``)
``needs_flag``  no safe value; pass the flag explicitly
``conflict``    two signals disagree; nothing is suggested for that parameter
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from typing import Any, Optional

from viralscan import chemistry, strand, tchem

ACCEPTED = "accepted"
REVIEW = "review"
NEEDS_FLAG = "needs_flag"
CONFLICT = "conflict"


@dataclass
class Param:
    value: Optional[str]
    status: str
    notes: list[str] = field(default_factory=list)


@dataclass
class KbRecommendation:
    technology: Param  # kb -x
    whitelist: Param  # kb -w ("" = kb's own on-list, "None" = no correction)
    strand: Param  # kb --strand
    evidence: dict[str, Any] = field(default_factory=dict)

    @property
    def status(self) -> str:
        """The worst status across the three parameters."""
        order = (ACCEPTED, REVIEW, NEEDS_FLAG, CONFLICT)
        return max(
            (p.status for p in (self.technology, self.whitelist, self.strand)), key=order.index
        )

    def kb_args(self) -> list[str]:
        """``kb count`` arguments for every parameter that has a value (including ``-w ""``)."""
        args: list[str] = []
        for flag, p in (("-x", self.technology), ("-w", self.whitelist), ("--strand", self.strand)):
            if p.value is not None:
                args += [flag, p.value]
        return args

    def as_block(self) -> dict[str, Any]:
        return {"status": self.status, **asdict(self)}


def _technology(
    r1: str,
    r2: Optional[str],
    detect: Callable[..., Any],
    second_opinion: Callable[..., Any],
    tchem_kwargs: dict[str, Any],
    evidence: dict[str, Any],
) -> tuple[Param, Any]:
    det = detect([r1])[0]
    evidence["chemistry"] = det.as_block()
    name = det.chemistry
    if name is None:
        return Param(None, NEEDS_FLAG, [f"built-in detector: {det.reason}"]), None
    op = second_opinion(name, r1, r2, **tchem_kwargs)
    evidence["tchem"] = op.as_block()
    if op.status == "disagree":
        return Param(None, CONFLICT, [f"built-in detector says {name}; {op.reason}"]), op
    if op.status == "agree":
        return Param(name, ACCEPTED, [f"detector and tchem agree: {op.reason}"]), op
    return Param(
        name, REVIEW, [f"only the built-in detector backs {name}: tchem {op.status} ({op.reason})"]
    ), op


def _whitelist(technology: Param) -> Param:
    if technology.value is None:
        return Param(None, technology.status, ["no technology, so no barcode list"])
    value = chemistry.kb_whitelist_arg(technology.value, None)
    note = (
        "kb skips barcode correction (no official on-list)"
        if value == "None"
        else "kb uses its own on-list"
    )
    return Param(value, technology.status, [note])


def _strand(
    op: Any, pilot_rates: Optional[dict[str, float]], tau: float, evidence: dict[str, Any]
) -> Param:
    from_tchem = op.call.strand if op is not None and op.call is not None else None
    from_pilot = None
    if pilot_rates is not None:
        from_pilot = strand.infer_strand(pilot_rates, tau)
        evidence["strand_pilot"] = {"rates": pilot_rates, "tau": tau, "winner": from_pilot}
    notes = []
    if from_tchem:
        notes.append(f"tchem strand: {from_tchem}")
    if from_pilot:
        notes.append(f"pilot winner (tau {tau}): {from_pilot}")
    if from_tchem is None and from_pilot is None:
        return Param(None, NEEDS_FLAG, ["no strand signal: pass --strand or supply pilot rates"])
    if from_tchem and from_pilot and from_tchem != from_pilot:
        # An unstranded pilot against a stranded tchem call keeps every read, so it is the safe fallback.
        return Param("unstranded", CONFLICT, notes + ["signals disagree: fall back to unstranded"])
    value = from_pilot or from_tchem
    assert value is not None
    if from_tchem is None or from_pilot is None:
        return Param(value, REVIEW, notes + ["only one signal available"])
    return Param(value, ACCEPTED, notes + ["tchem and pilot agree"])


def recommend_kb_parameters(
    r1: str,
    r2: Optional[str] = None,
    *,
    pilot_rates: Optional[dict[str, float]] = None,
    tau: float = strand.TAU,
    tchem_mode: str = "auto",
    tchem_bin: Optional[str] = None,
    tchem_whitelists: Optional[str] = None,
    detect: Callable[..., Any] = chemistry.detect,
    second_opinion: Callable[..., Any] = tchem.second_opinion,
) -> KbRecommendation:
    """Recommend ``-x``, ``-w`` and ``--strand`` for one sample; see the module docstring.

    ``tau`` is the existing strand pilot threshold (``strand.TAU``); it is inherited, not new.
    """
    evidence: dict[str, Any] = {}
    kwargs = {"mode": tchem_mode, "binary": tchem_bin, "whitelists": tchem_whitelists}
    technology, op = _technology(r1, r2, detect, second_opinion, kwargs, evidence)
    return KbRecommendation(
        technology=technology,
        whitelist=_whitelist(technology),
        # tchem's strand belongs to the chemistry it named; ignore it when the chemistries conflict
        strand=_strand(None if technology.status == CONFLICT else op, pilot_rates, tau, evidence),
        evidence=evidence,
    )


#: tchem layouts that count as confirmed 3' for the ``--strand`` gate. Drop-seq is
#: treated as 3' (assumption recorded in PLAN ``STRAND-EXPLICIT-01``).
THREE_PRIME_LAYOUTS = ("10x 3'", "Drop-seq")


def is_confirmed_three_prime(tchem_block: Optional[dict[str, Any]]) -> bool:
    """True when the recorded tchem opinion agrees with ``-x`` and names a 3' layout.

    ``-x 10xv2`` names 3' v2 and 5' v1/v2 libraries alike, so the chemistry name cannot
    decide it. Only a tchem match whose layout starts with a 3' prefix counts. A missing,
    unavailable, ambiguous or disagreeing opinion is not a confirmation.
    """
    if not tchem_block or tchem_block.get("status") != "agree":
        return False
    call = tchem_block.get("call") or {}
    layout = call.get("layout") or ""
    return call.get("verdict") == "match" and layout.startswith(THREE_PRIME_LAYOUTS)
