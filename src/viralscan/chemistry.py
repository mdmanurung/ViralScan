"""Chemistry module (PLAN MECH-D, DEF-02): one geometry table, read-based detection.

Every barcode geometry ViralScan uses -- kb's ``-x``, the STARsolo arguments and
the whitelist preflight -- comes from :data:`CHEMISTRIES`. :func:`detect` reads
the chemistry off the first R1 reads (WP1E Q6): which on-list the barcodes match,
and where the template-switch oligo (5') or poly-T (3') starts, which gives the
UMI length. :func:`resolve` then fails on ambiguity or on a mismatch with ``-x``
unless the user forces it.

Measured on the four local libraries (2026-10-03, 100k R1 reads each):

==========================  ======  ==================  ==============================
library                     R1 len  read structure      on-list match
==========================  ======  ==================  ==============================
EBV SRR12682296, 10xv2 3'   26      trimmed             v2 0.973, v3 0.093
HHV-6B SRR20710641, 5'      150     TSO at 26 (0.84)    v2 0.882, v3 0.097
covid x213, GEM-X 5'        28      trimmed             user list 0.672, v4 0.044
HSV-1 SRR8315713, Drop-seq  20      trimmed             none above 0.005
==========================  ======  ==================  ==============================

No bundled list covers GEM-X 5', so such a library is only detected with ``-w``.
"""

from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Optional


@dataclass(frozen=True)
class Chemistry:
    name: str  # kb ``-x``
    cb_len: int
    umi_len: int
    #: ngs_tools' bundled on-list file, or None when the technology has no
    #: official on-list. kb then builds an allowlist from the data, which
    #: ViralScan bypasses (SW-21): cell calling stays ViralScan's own.
    onlist: Optional[str]


CHEMISTRIES: dict[str, Chemistry] = {
    c.name: c
    for c in (
        Chemistry("10xv1", 14, 10, "10x_version1_whitelist.txt.gz"),
        Chemistry("10xv2", 16, 10, "10x_version2_whitelist.txt.gz"),
        Chemistry("10xv3", 16, 12, "10x_version3_whitelist.txt.gz"),
        Chemistry("10xv4", 16, 12, "10x_version4_whitelist.txt.gz"),
        Chemistry("dropseq", 12, 8, None),
    )
}

#: Chemistries :func:`detect` can call from an on-list match. 10xv1 keeps its
#: UMI on a separate read, so it is never inferred from R1.
DETECTABLE = ("10xv2", "10xv3", "10xv4")

SAMPLE_READS = 100_000
MIN_MATCH_RATE = 0.5
#: Fraction of reads that must carry the adapter at one position.
MIN_ADAPTER_RATE = 0.5
TSO = "TTTCTTATATGGG"
POLY_T = "TTTTTTTTTT"
#: Adapter start = 16 bp barcode + 10 or 12 bp UMI.
ADAPTER_STARTS = (26, 28)
#: R1 trimmed to exactly barcode + UMI, by length.
UMI_BY_TRIMMED_LENGTH = {26: 10, 28: 12}
#: The user decision of 2026-10-03: no 10x list matches and R1 is 20 bp means
#: Drop-seq (12 bp barcode + 8 bp UMI).
DROPSEQ_R1_LENGTH = 20


class ChemistryError(ValueError):
    """The library's chemistry is ambiguous or contradicts ``-x``."""


def get(name: str) -> Chemistry:
    try:
        return CHEMISTRIES[name.strip().lower()]
    except KeyError:
        raise ValueError(f"Unknown technology {name!r}") from None


def cb_umi_geometry(technology: str) -> tuple[int, int]:
    """``(cb_len, umi_len)`` for a named chemistry or a kb ``bc:umi:seq`` string.

    Raises ``ValueError`` for anything it cannot resolve -- extracting reads
    with the wrong geometry silently yields nothing, so this fails loudly.
    """
    if technology.strip().lower() in CHEMISTRIES:
        c = get(technology)
        return c.cb_len, c.umi_len
    if ":" in technology:
        try:
            bc, umi, _seq = technology.split(":")[:3]
            _, bcs, bce = (int(x) for x in bc.split(","))
            _, umis, umie = (int(x) for x in umi.split(","))
            return bce - bcs, umie - umis
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Cannot parse barcode geometry from -x {technology!r}: {exc}"
            ) from exc
    raise ValueError(
        f"Unknown technology {technology!r}; add it to CHEMISTRIES or pass an "
        "explicit 'bc:umi:seq' geometry string."
    )


def onlist_path(name: str) -> Optional[str]:
    """Path of the bundled on-list for *name*, or None if kb ships none.

    The lists come with kb-python (via ngs_tools), which is a runtime tool and
    not a pip dependency, so a missing install raises rather than skipping.
    """
    chem = get(name)
    if chem.onlist is None:
        return None
    try:
        from ngs_tools.chemistry import get_chemistry
    except ImportError as exc:
        raise ChemistryError(
            "kb-python's on-lists (ngs_tools) are not importable, so the chemistry "
            "cannot be checked. Install kb-python into this environment, or pass "
            "-x with --force-technology."
        ) from exc
    return get_chemistry(name).whitelist_path


def kb_whitelist_arg(technology: str, whitelist: Optional[str]) -> str:
    """The ``kb count -w`` value: the user's list, ``""`` (kb's on-list) or ``"None"``.

    ``"None"`` makes kb skip barcode correction. That is the SW-21 decision
    (user, 2026-10-03) for a chemistry with no official on-list: kb would
    otherwise build an allowlist from a knee on the data and drop every other
    barcode (Drop-seq HSV-1: 1.8M barcodes to 5,468, -28 % HSV-1 molecules),
    a second, hidden cell caller ahead of ViralScan's own.
    """
    if whitelist:
        return whitelist
    chem = CHEMISTRIES.get(str(technology).strip().lower())
    return "None" if chem is not None and chem.onlist is None else ""


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def sample_r1(path: str, n: int = SAMPLE_READS) -> list[str]:
    seqs: list[str] = []
    with _open(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    return seqs


def match_rates(samples: list[list[str]], lists: dict[str, str]) -> list[dict[str, float]]:
    """Per sample, the fraction of R1 barcodes found in each on-list.

    Streams each list once past every sample's barcode counts, so a 6.8M-entry
    list is never held in memory.
    """
    rates: list[dict[str, float]] = [{} for _ in samples]
    for label, path in lists.items():
        with _open(path) as fh:
            first = fh.readline().strip()
        cb_len = len(first)
        counts = [Counter(s[:cb_len] for s in seqs if len(s) >= cb_len) for seqs in samples]
        hits = [0] * len(samples)
        with _open(path) as fh:
            for line in fh:
                bc = line.strip()
                for i, c in enumerate(counts):
                    hits[i] += c.get(bc, 0)
        for i, seqs in enumerate(samples):
            rates[i][label] = round(hits[i] / len(seqs), 4) if seqs else 0.0
    return rates


def read_structure(seqs: list[str]) -> tuple[int, Optional[str], Optional[int]]:
    """``(modal R1 length, end, umi_len)`` from the R1 sequences.

    *end* is ``"5p"`` (TSO after the UMI), ``"3p"`` (poly-T) or None when R1
    was trimmed to barcode + UMI. The UMI length comes from where the adapter
    starts, else from a trimmed length of 26 or 28.
    """
    r1_len = Counter(len(s) for s in seqs).most_common(1)[0][0] if seqs else 0
    for start in ADAPTER_STARTS:
        for end, motif in (("5p", TSO), ("3p", POLY_T)):
            hits = sum(s[start : start + len(motif)] == motif for s in seqs)
            if seqs and hits / len(seqs) >= MIN_ADAPTER_RATE:
                return r1_len, end, start - 16
    return r1_len, None, UMI_BY_TRIMMED_LENGTH.get(r1_len)


@dataclass
class Detection:
    chemistry: Optional[str]
    basis: str  # "bundled on-list", "user on-list", "R1 length", or "ambiguous"
    reason: str
    r1_length: int
    end: Optional[str]
    umi_len: Optional[int]
    n_sampled: int
    match_rates: dict[str, float] = field(default_factory=dict)

    def as_block(self) -> dict:
        return asdict(self)


def classify(
    seqs: list[str], rates: dict[str, float], user_list: Optional[str] = None
) -> Detection:
    """Call a chemistry from one sample's R1 reads and on-list match rates.

    *rates* is keyed by chemistry name, plus ``"user"`` for the ``-w`` list.
    """
    r1_len, end, umi = read_structure(seqs)

    def out(chem, basis, reason):
        return Detection(chem, basis, reason, r1_len, end, umi, len(seqs), dict(rates))

    if user_list is not None:
        rate = rates.get("user", 0.0)
        if rate < MIN_MATCH_RATE:
            return out(
                None,
                "ambiguous",
                f"only {rate:.1%} of R1 barcodes are in the -w list {user_list}; "
                "it does not belong to this library (F-005)",
            )
        by_umi = {10: "10xv2", 12: "10xv3"}.get(umi or -1)
        if by_umi is None:
            return out(
                None,
                "ambiguous",
                f"the -w list matches ({rate:.1%}) but the UMI length cannot be read "
                f"from R1 (length {r1_len}, no TSO or poly-T at {ADAPTER_STARTS})",
            )
        return out(by_umi, "user on-list", f"-w list matches {rate:.1%}; UMI {umi} bp")

    hits = [c for c in DETECTABLE if rates.get(c, 0.0) >= MIN_MATCH_RATE]
    if len(hits) > 1:
        return out(None, "ambiguous", f"R1 barcodes match several on-lists: {hits}")
    if hits:
        chem = get(hits[0])
        if umi is not None and umi != chem.umi_len:
            return out(
                None,
                "ambiguous",
                f"barcodes match the {chem.name} on-list but R1 implies a {umi} bp "
                f"UMI, not {chem.umi_len}",
            )
        return out(chem.name, "bundled on-list", f"{rates[chem.name]:.1%} on the {chem.name} list")
    if r1_len == DROPSEQ_R1_LENGTH:
        return out("dropseq", "R1 length", "no 10x on-list matches and R1 is 20 bp")
    best = max(rates.items(), key=lambda kv: kv[1], default=("none", 0.0))
    return out(
        None,
        "ambiguous",
        f"no on-list matches (best {best[0]} {best[1]:.1%}) and R1 is {r1_len} bp. "
        "Pass the library's on-list with -w (GEM-X 5' needs Cell Ranger's).",
    )


def detect(r1_paths: list[str], user_list: Optional[str] = None) -> list[Detection]:
    """Detect each sample's chemistry from its first :data:`SAMPLE_READS` R1 reads."""
    lists = {c: p for c in DETECTABLE if (p := onlist_path(c))}
    if user_list:
        lists["user"] = user_list
    samples = [sample_r1(p) for p in r1_paths]
    return [
        classify(seqs, rates, user_list)
        for seqs, rates in zip(samples, match_rates(samples, lists))
    ]


def _same(requested: str, detected: str, basis: str) -> bool:
    """A user on-list fixes only the geometry, so compare geometry there."""
    if basis == "user on-list" or requested.strip().lower() not in CHEMISTRIES:
        return cb_umi_geometry(requested) == cb_umi_geometry(detected)
    return requested.strip().lower() == detected


def resolve(requested: Optional[str], detections: list[Detection], force: bool = False) -> str:
    """The technology to run, per WP1E Q6.

    ``requested`` None means auto. A detection that is ambiguous, that
    disagrees between samples, or that contradicts ``-x`` raises
    :class:`ChemistryError`; ``force`` keeps an explicit ``-x`` regardless.
    """
    if force:
        if not requested:
            raise ChemistryError("--force-technology needs an explicit -x.")
        return requested
    called = {d.chemistry for d in detections}
    problem = None
    if None in called:
        problem = "; ".join(d.reason for d in detections if d.chemistry is None)
    elif len(called) > 1:
        problem = f"samples disagree on the chemistry: {sorted(called)}"
    elif requested and not _same(requested, detections[0].chemistry, detections[0].basis):
        d = detections[0]
        problem = f"-x {requested} but the reads look like {d.chemistry} ({d.reason})"
    if problem:
        raise ChemistryError(
            f"Chemistry check failed: {problem}. Fix -x/-w, or pass -x with "
            "--force-technology to run anyway."
        )
    return requested or detections[0].chemistry
