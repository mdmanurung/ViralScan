"""Library diagnostics: is this single-cell or bulk, 3' or 5', and did the run look right.

Builds on :mod:`viralscan.chemistry` (barcode on-list, read structure) and
:mod:`viralscan.strand` (1M-pair pilot). Three pure helpers plus one entry point
for the ``viralscan check-chemistry`` subcommand (PLAN ``DEF-02`` / ``CHEM-01``):

* :func:`library_kind`: bulk vs single-cell from the barcode repeat rate;
* :func:`infer_end`: 3' vs 5' from R1 structure, else from the strand pilot (F-020);
* :func:`sanity_gate`: post-``kb count`` symptoms of a wrong chemistry or strand.

Every verdict is advice. The only fail-closed paths stay in ``chemistry.resolve``.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Optional

from viralscan import chemistry, strand

#: Barcode geometries to try: (cb_len, umi_len) for 10x (16/10-12) and Drop-seq (12/8).
GEOMETRIES = ((16, 10), (12, 8))
#: Fraction of reads whose barcode prefix occurs on another read. 10x libraries repeat each
#: barcode (EBV 10xv2 0.91, HHV-6B 0.89, tonsil 0.63, Drop-seq at 12 bp 0.68); a bulk R1 is
#: a genomic read and almost never repeats.
BARCODED_MIN_REPEAT = 0.5
BULK_MAX_REPEAT = 0.1
#: Among reads sharing a prefix, the fraction that also share the next UMI-length bases. Real
#: barcode+UMI reads differ there (measured 0.02-0.43); cDNA reads that share a prefix are
#: the same transcript and continue identically (EBV R2 0.94, HHV-6B R2 0.97).
CDNA_MIN_CONTINUATION = 0.7
#: Below this, kb pseudoaligned too little of a host-containing index for the library to be right.
MIN_P_PSEUDOALIGNED = 40.0
FAIL_P_PSEUDOALIGNED = 10.0


def prefix_stats(seqs: list[str], cb_len: int, umi_len: int) -> tuple[float, Optional[float]]:
    """``(repeat rate, continuation rate)`` of the first *cb_len* bases of *seqs*.

    Continuation is None when no prefix repeats.
    """
    long_enough = [s for s in seqs if len(s) >= cb_len + umi_len]
    if not long_enough:
        return 0.0, None
    by_prefix: dict[str, list[str]] = {}
    for s in long_enough:
        by_prefix.setdefault(s[:cb_len], []).append(s[cb_len : cb_len + umi_len])
    repeated = [v for v in by_prefix.values() if len(v) > 1]
    n_rep = sum(len(v) for v in repeated)
    if not n_rep:
        return 0.0, None
    same = sum(c for v in repeated for c in Counter(v).values() if c > 1)
    return n_rep / len(long_enough), same / n_rep


def library_kind(seqs: list[str], on_list_rate: float = 0.0) -> tuple[str, str]:
    """``("single-cell" | "bulk" | "unclear", reason)`` from the first R1 reads.

    ``bulk`` also covers a cDNA mate passed as ``-s1`` (swapped R1/R2).
    """
    r1_len = Counter(len(s) for s in seqs).most_common(1)[0][0] if seqs else 0
    stats = [prefix_stats(seqs, cb, umi) for cb, umi in GEOMETRIES]
    rep, cont = max(stats, key=lambda s: s[0])
    detail = f"barcode repeat {rep:.2f}, continuation {'n/a' if cont is None else f'{cont:.2f}'}, on-list {on_list_rate:.1%}"
    if on_list_rate >= chemistry.MIN_MATCH_RATE:
        return "single-cell", detail
    if rep > BULK_MAX_REPEAT and cont is not None and cont >= CDNA_MIN_CONTINUATION:
        return (
            "bulk",
            f"{detail}: reads sharing a prefix continue identically, so R1 is cDNA (bulk, or -s1/-s2 swapped)",
        )
    if rep >= BARCODED_MIN_REPEAT:
        return "single-cell", detail
    if rep <= BULK_MAX_REPEAT and r1_len >= 40:
        return "bulk", f"{detail}: R1 is {r1_len} bp with no barcode structure"
    return "unclear", f"{detail}, R1 {r1_len} bp"


def infer_end(
    r1_end: Optional[str], pilot_rates: Optional[dict[str, float]]
) -> tuple[Optional[str], str]:
    """``("3p" | "5p" | None, basis)``.

    R1 carrying a TSO or poly-T after the UMI is decisive. A trimmed R1 (barcode + UMI
    only, e.g. tonsil x223 at 28 bp) cannot tell 3' from 5': on-lists are shared. There the
    strand pilot decides (F-020: 5' libraries keep reverse/unstranded >= 0.85, forward
    0.11-0.16; 3' libraries keep forward). A conflict between the two is reported, not resolved.
    """
    from_pilot = None
    if pilot_rates:
        choice = strand.infer_strand(pilot_rates)
        from_pilot = {"reverse": "5p", "forward": "3p"}.get(choice)
    if r1_end and from_pilot and r1_end != from_pilot:
        return None, f"conflict: R1 structure says {r1_end}, strand pilot says {from_pilot}"
    if r1_end:
        return r1_end, "R1 structure (adapter after the UMI)"
    if from_pilot:
        return from_pilot, "strand pilot (R1 is trimmed)"
    return None, "undetermined (R1 trimmed and no strand pilot or an unstranded pilot)"


def sanity_gate(run_info: dict, *, host_in_index: bool = True) -> list[dict[str, str]]:
    """Warnings from a finished ``kb count`` that point at the wrong chemistry or strand.

    *host_in_index* is False after a host-filter (two-step): the reads are already host-free,
    so a low pseudoalignment rate is expected there and is not judged.
    """
    findings: list[dict[str, str]] = []
    p = float(run_info.get("p_pseudoaligned", 0.0))
    if host_in_index and p < MIN_P_PSEUDOALIGNED:
        level = "error" if p < FAIL_P_PSEUDOALIGNED else "warning"
        findings.append(
            {
                "level": level,
                "check": "p_pseudoaligned",
                "message": f"only {p:.1f}% of reads pseudoalign to a host-containing index. "
                "On a 5' library this is the forward-strand default (F-020: 6.5-8.9%); "
                "otherwise suspect the wrong -x/-w (F-005) or the wrong index. "
                "Run `viralscan check-chemistry`.",
            }
        )
    return findings


def sanity_gate_from_dir(kb_dir: str | Path, *, host_in_index: bool = True) -> list[dict[str, str]]:
    """:func:`sanity_gate` on ``<kb_dir>/run_info.json``; a missing file yields no findings."""
    path = Path(kb_dir) / "run_info.json"
    if not path.is_file():
        return []
    return sanity_gate(json.loads(path.read_text()), host_in_index=host_in_index)


def diagnose(
    s1: str,
    s2: Optional[str] = None,
    *,
    whitelist: Optional[str] = None,
    index: Optional[str] = None,
    t2g: Optional[str] = None,
    technology: Optional[str] = None,
    cores: int = 4,
    pilot_reads: int = strand.PILOT_READS,
) -> dict:
    """Everything ``check-chemistry`` reports, as one JSON-ready dict."""
    seqs = chemistry.sample_r1(s1)
    det = chemistry.detect([s1], whitelist)[0]
    on_list = max(det.match_rates.values(), default=0.0)
    kind, kind_reason = library_kind(seqs, on_list)
    report: dict = {
        "library_kind": {"call": kind, "reason": kind_reason},
        "chemistry": det.as_block(),
        "requested_technology": technology,
    }
    pilot = None
    if kind == "single-cell" and index and t2g and s2 and det.chemistry:
        pilot = strand.run_pilot(
            s1, s2, index, t2g, technology or det.chemistry, whitelist, cores, pilot_reads
        )
        report["strand_pilot"] = strand.inference_block(pilot, pilot_reads)
    end, basis = infer_end(det.end, pilot)
    report["end"] = {"call": end, "basis": basis}
    report["advice"] = _advice(report, whitelist)
    return report


def _advice(report: dict, whitelist: Optional[str]) -> list[str]:
    out: list[str] = []
    kind = report["library_kind"]["call"]
    chem = report["chemistry"]
    if kind == "bulk":
        out.append(
            "Bulk library: no cell barcodes. Use `kb count -x BULK` outside ViralScan's per-cell workflow."
        )
        return out
    if kind == "unclear":
        out.append("Neither barcoded nor clearly bulk: check that -s1 is the barcode read (R1).")
    if chem["chemistry"] is None:
        out.append(f"Chemistry not resolved: {chem['reason']}")
    else:
        flag = f"-x {chem['chemistry']}"
        out.append(flag + ("" if whitelist or chem["basis"] != "user on-list" else " -w <on-list>"))
    pilot = report.get("strand_pilot")
    if pilot:
        out.append(f"--strand {pilot['choice']} (pilot rates {pilot['rates']})")
    elif kind == "single-cell":
        out.append(
            "Pass -i/-t/-s2 to add the strand pilot, which also tells 3' from 5' on trimmed R1."
        )
    if report["end"]["call"] is None and "conflict" in report["end"]["basis"]:
        out.append(f"Resolve before running: {report['end']['basis']}.")
    return out
