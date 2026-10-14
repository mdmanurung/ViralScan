"""Q137 read-support flags per molecule, from the competitive alignment (never removes a call).

Gates from the Biomni torque teno read review (2026-10-14), provisional until calibrated. They are FLAGS: every
molecule stays in every table and no detection threshold moves.

A molecule is a corrected (CB, UMI). Its **representative read** is the read with the best viral alignment
(``molecule_verdict.select_representatives``: best viral ``AS``, longest M/=/X, lexical name); a molecule with no
viral alignment is represented by its best host read. Gates 1, 2 and 4 are evaluated on the representative
read's best viral record; gate 3 uses the molecule's best viral and best host ``AS`` over ALL its reads, so it
agrees with ``molecule_verdicts.tsv.gz``.

1. ``aligned_bases >= 60`` (M/=/X only) and ``>= 60 %`` of the trimmed query: the read minus a leading 10x
   TSO/TruSeq reagent and everything from the first >= 15-nt homopolymer on (``anello_align.body_span``).
2. ``dust_unmasked_bases >= 20``: bases of the viral alignment that symmetric DUST (window 64, level 20) leaves
   unmasked, on the whole read. ``aligned_AT_warning`` (not a gate) when max(A, T) of the aligned segment >= 0.70.
3. Competitive: the same minimap2 run scored the read on human and virus; human ``AS`` >= viral ``AS`` is
   ``human_equal_or_better`` (human or ambiguous).
4. ``hotspot_only``: every viral alignment of the representative read lies inside one hotspot of
   ``data/read_support_loci.tsv``. The same file lists the DNAJB6 locus as a human decoy (annotation).

``known_artefact`` per molecule names the pattern that fails it (``human_DNAJB6``: human-equal-or-better with its
best human hit on the decoy locus; ``low_complexity_hotspot``: hotspot-only). The call-level summary reports a
reason only when EVERY molecule of the call fails for it. No genus is blacklisted.

DUST is implemented here (the evidence environment is not guaranteed to carry ``dustmasker``). Its masks equal
``sdust -w 64 -t 20`` on 3153 of 3154 sequences (2,854 round-2 evidence reads plus synthetic repeats); the one
exception is a repeat longer than the window, where sdust extends the mask by 5 nt. ``dustmasker`` differs from
``sdust`` on the same three sequences that it differs from this code (tests/test_read_support.py, skipped
without ``sdust``).

    python -m viralscan.read_support <evidence dir> [<evidence dir> ...] -o <output dir>   # read-only demo
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

from viralscan.anello_align import body_span
from viralscan.molecule_verdict import (
    _CIGAR,
    matched_bases,
    select_representatives,
    side_score,
    verdict_of,
)

MIN_ALIGNED_BASES = 60
MIN_QUERY_FRACTION = 0.60
MIN_DUST_UNMASKED = 20
AT_WARNING_FRACTION = 0.70
DUST_WINDOW = 64
DUST_LEVEL = 20
DUST_METHOD = f"symmetric_dust_w{DUST_WINDOW}_l{DUST_LEVEL}"
LOCI_TSV = Path(__file__).parent / "data" / "read_support_loci.tsv"
_COMP = str.maketrans("ACGTN", "TGCAN")

ROW_FIELDS = (
    "read",
    "cb",
    "umi",
    "count_layer",
    "n_reads",
    "aligned_bases",
    "trimmed_query_bases",
    "query_coverage",
    "dust_unmasked_bases",
    "aligned_max_AT_fraction",
    "virus_as",
    "host_as",
    "verdict",
    "viral_reference",
    "viral_start",
    "viral_end",
    "host_locus",
    "gate_aligned_length",
    "gate_dust",
    "gate_competitive",
    "hotspot_only",
    "host_decoy",
    "aligned_AT_warning",
    "known_artefact",
    "flags",
    "passes_all_gates",
)
SUMMARY_FIELDS = (
    "evidence_dir",
    "n_reads",
    "n_molecules",
    "n_pass_all_gates",
    "fraction_pass_all_gates",
    "n_fail_aligned_length",
    "n_fail_dust",
    "n_human_equal_or_better",
    "n_hotspot_only",
    "n_aligned_AT_warning",
    "n_host_decoy",
    "known_artefact",
    "known_artefact_mixed",
    "dust_method",
)


def dust_mask(seq: str, window: int = DUST_WINDOW, level: int = DUST_LEVEL) -> str:
    """Symmetric DUST (Morgulis et al. 2006): *seq* with low-complexity bases lower-cased.

    A triplet interval of ``l`` triplets scores ``sum c_t (c_t - 1) / 2 / (l - 1)``. It is masked when the score
    exceeds ``level / 10`` and it is *perfect*: every proper subinterval scores strictly lower. Intervals span at
    most *window* bases. Triplets containing N never repeat. Cost is O(len * window); reads are <= 150 nt.
    """
    upper = seq.upper()
    n = len(upper) - 2
    if n < 2:
        return upper
    tri = [t if set(t) <= set("ACGT") else None for t in (upper[i : i + 3] for i in range(n))]
    max_l = window - 2
    score = [[0.0] * (max_l + 1) for _ in range(n)]
    for i in range(n):
        counts: dict[str, int] = {}
        pairs = 0
        for ln in range(1, min(max_l, n - i) + 1):
            t = tri[i + ln - 1]
            if t is not None:
                pairs += counts.get(t, 0)
                counts[t] = counts.get(t, 0) + 1
            score[i][ln] = pairs / (ln - 1) if ln > 1 else 0.0
    sub = [[0.0] * (max_l + 1) for _ in range(n)]
    out = list(upper)
    for ln in range(2, max_l + 1):
        for i in range(n - ln + 1):
            sub[i][ln] = max(
                score[i + 1][ln - 1], score[i][ln - 1], sub[i + 1][ln - 1], sub[i][ln - 1]
            )
            if score[i][ln] > level / 10 and score[i][ln] > sub[i][ln]:
                for p in range(i, i + ln + 2):
                    out[p] = upper[p].lower()
    return "".join(out)


def load_loci(path: str | Path = LOCI_TSV) -> list[dict[str, str]]:
    """Rows of ``data/read_support_loci.tsv`` (kind, reference, start, end, reason, note)."""
    with Path(path).open(newline="") as handle:
        rows = [ln for ln in handle if ln.strip() and not ln.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def _ref_span(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR.findall(cigar) if op in "MDN=X")


def _aligned_positions(cigar: str, reverse: bool) -> tuple[set[int], int]:
    """M/=/X read positions (sequencing orientation) and the full read length, hard clips included."""
    pos = 0
    out: set[int] = set()
    for num, op in _CIGAR.findall(cigar):
        n = int(num)
        if op in "M=X":
            out.update(range(pos, pos + n))
        if op in "MISH=X":
            pos += n
    return ({pos - 1 - p for p in out} if reverse else out), pos


def _sequence(records: list[list[str]], queries: Mapping[str, str] | None) -> str:
    """Read in sequencing orientation: the primary record's SEQ, else the extracted FASTA's."""
    for f in records:
        if not int(f[1]) & 0x900 and f[9] != "*":
            seq = f[9].upper()
            return seq.translate(_COMP)[::-1] if int(f[1]) & 16 else seq
    return (queries or {}).get(records[0][0], "").upper()


def molecule_support_rows(
    sam_lines: Iterable[str],
    queries: Mapping[str, str] | None = None,
    loci: list[dict[str, str]] | None = None,
) -> list[dict[str, object]]:
    """One flag row per corrected (CB, UMI) molecule of a competitive SAM (all reads, secondary included).

    # ponytail: every read's records are held in memory; stream per molecule if a target ever has > 10 M reads.
    """
    loci = load_loci() if loci is None else loci
    hotspots = [r for r in loci if r["kind"] == "hotspot"]
    decoys = [r for r in loci if r["kind"] == "host_decoy"]
    rep, by_read, n_reads = select_representatives(sam_lines)
    members: dict[tuple[str, str], list[str]] = defaultdict(list)
    for name in by_read:
        key = tuple(name.split("_", 2)[:2])
        members[(key[0], key[1])].append(name)
    rows: list[dict[str, object]] = []
    for (cb, umi), rep_name in sorted(rep.items()):
        best: dict[str, int] = {}
        host_hits: list[tuple[int, list[str]]] = []
        for name in members[(cb, umi)]:
            for f in by_read[name]:
                hit = side_score(f)
                if hit is None:
                    continue
                side, score = hit
                best[side] = max(score, best.get(side, score))
                if side == "host":
                    host_hits.append((score, f))
        virus_as, host_as = best.get("virus"), best.get("host")
        verdict = verdict_of(virus_as, host_as)
        records = by_read[rep_name]
        viral = [
            f for f in records if not int(f[1]) & 4 and f[2] != "*" and not f[2].startswith("HOST|")
        ]
        scored = [f for f in viral if side_score(f) is not None]
        hit_rec = min(
            scored,
            key=lambda f: (-(side_score(f) or ("", 0))[1], -matched_bases(f[5]), f[2], int(f[3])),
            default=None,
        )
        seq = _sequence(records, queries)
        positions: set[int] = set()
        if hit_rec is not None:
            positions, length = _aligned_positions(hit_rec[5], bool(int(hit_rec[1]) & 16))
            if seq and length and length != len(seq):
                raise ValueError(
                    f"CIGAR length {length} != read length {len(seq)} for {rep_name!r}"
                )
        aligned = len(positions)
        start, end = body_span(seq) if seq else (0, 0)
        trimmed = end - start
        coverage = min(1.0, aligned / trimmed) if trimmed else 0.0
        masked = dust_mask(seq) if seq else ""
        unmasked = sum(masked[p].isupper() for p in positions) if masked else 0
        at = max(sum(seq[p] == b for p in positions) for b in "AT") / aligned if aligned else 0.0
        ref = hit_rec[2].split("|", 1)[-1] if hit_rec else ""
        v_start = int(hit_rec[3]) if hit_rec else 0
        v_end = v_start + _ref_span(hit_rec[5]) - 1 if hit_rec else 0
        hot = next(
            (
                h
                for h in hotspots
                if viral
                and all(
                    f[2].split("|", 1)[-1] == h["reference"]
                    and int(h["start"]) <= int(f[3])
                    and int(f[3]) + _ref_span(f[5]) - 1 <= int(h["end"])
                    for f in viral
                )
            ),
            None,
        )
        top_host = max((s for s, _ in host_hits), default=None)
        host_best = [f for s, f in host_hits if s == top_host]
        decoy = next(
            (
                d
                for d in decoys
                for f in host_best
                if f[2].split("|", 1)[-1] == d["reference"]
                and int(f[3]) <= int(d["end"])
                and int(f[3]) + _ref_span(f[5]) - 1 >= int(d["start"])
            ),
            None,
        )
        host_locus = f"{host_best[0][2].split('|', 1)[-1]}:{host_best[0][3]}" if host_best else ""
        ok_len = aligned >= MIN_ALIGNED_BASES and coverage >= MIN_QUERY_FRACTION
        ok_dust = unmasked >= MIN_DUST_UNMASKED
        ok_comp = verdict == "virus_best"
        reasons = []
        if decoy and not ok_comp:
            reasons.append(decoy["reason"])
        if hot:
            reasons.append(hot["reason"])
        flags = []
        if not seq:
            flags.append("sequence_unavailable")
        if hit_rec is None:
            flags.append("no_viral_alignment")
        elif not ok_len:
            flags.append("aligned_length_lt60_or_query_lt60pct")
        if hit_rec is not None and not ok_dust:
            flags.append("dust_unmasked_lt20")
        if not ok_comp and hit_rec is not None:
            flags.append("human_equal_or_better")
        if hot:
            flags.append("hotspot_only")
        if at >= AT_WARNING_FRACTION:
            flags.append("aligned_AT_warning")
        rows.append(
            dict(
                read=rep_name,
                cb=cb,
                umi=umi,
                count_layer="umi",
                n_reads=n_reads[(cb, umi)],
                aligned_bases=aligned,
                trimmed_query_bases=trimmed,
                query_coverage=round(coverage, 4),
                dust_unmasked_bases=unmasked,
                aligned_max_AT_fraction=round(at, 4),
                virus_as="" if virus_as is None else virus_as,
                host_as="" if host_as is None else host_as,
                verdict=verdict,
                viral_reference=ref,
                viral_start=v_start or "",
                viral_end=v_end or "",
                host_locus=host_locus,
                gate_aligned_length=ok_len,
                gate_dust=ok_dust,
                gate_competitive=ok_comp,
                hotspot_only=hot is not None,
                host_decoy=decoy is not None,
                aligned_AT_warning=at >= AT_WARNING_FRACTION,
                known_artefact=";".join(reasons),
                flags=";".join(flags),
                passes_all_gates=ok_len and ok_dust and ok_comp and hot is None,
            )
        )
    return rows


def summarise_support(rows: list[dict[str, object]], evidence_dir: str = "") -> dict[str, object]:
    """Call-level summary. ``known_artefact`` lists a reason only when every molecule fails for it."""
    n = len(rows)
    reasons = [set(str(r["known_artefact"]).split(";")) - {""} for r in rows]
    every = sorted(set.intersection(*reasons)) if reasons else []
    union = set().union(*reasons) if reasons else set()
    mixed = sorted(union) if n and all(reasons) and not every and len(union) > 1 else []
    passed = sum(bool(r["passes_all_gates"]) for r in rows)
    return {
        "evidence_dir": evidence_dir,
        "n_reads": sum(int(str(r["n_reads"])) for r in rows),
        "n_molecules": n,
        "n_pass_all_gates": passed,
        "fraction_pass_all_gates": round(passed / n, 4) if n else "",
        "n_fail_aligned_length": sum(not r["gate_aligned_length"] for r in rows),
        "n_fail_dust": sum(not r["gate_dust"] for r in rows),
        "n_human_equal_or_better": sum(not r["gate_competitive"] for r in rows),
        "n_hotspot_only": sum(bool(r["hotspot_only"]) for r in rows),
        "n_aligned_AT_warning": sum(bool(r["aligned_AT_warning"]) for r in rows),
        "n_host_decoy": sum(bool(r["host_decoy"]) for r in rows),
        "known_artefact": ";".join(every),
        "known_artefact_mixed": ";".join(mixed),
        "dust_method": DUST_METHOD,
    }


def _write(path: Path, rows: list[dict[str, object]], fields: Iterable[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_molecule_support(bam: str | Path, out_dir: Path) -> list[dict[str, object]]:
    """Write ``molecule_support.tsv`` and ``read_support_summary.tsv`` for one competitive BAM."""
    with subprocess.Popen(
        ["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdout is not None
        rows = molecule_support_rows(proc.stdout)
    if proc.returncode:
        raise RuntimeError(f"samtools view failed (exit {proc.returncode}) on {bam}")
    out_dir.mkdir(parents=True, exist_ok=True)
    _write(out_dir / "molecule_support.tsv", rows, ROW_FIELDS)
    _write(
        out_dir / "read_support_summary.tsv",
        [summarise_support(rows, out_dir.name)],
        SUMMARY_FIELDS,
    )
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Q137 molecule-support flags for evidence directories (read-only)."
    )
    ap.add_argument("evidence_dir", nargs="+", type=Path)
    ap.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args(argv)
    summaries = []
    for d in args.evidence_dir:
        out = args.output / d.name
        write_molecule_support(d / "competitive_reads.raw.bam", out)
        with (out / "read_support_summary.tsv").open() as handle:
            summaries.append(next(csv.DictReader(handle, delimiter="\t")))
    _write(args.output / "read_support_calls.tsv", summaries, SUMMARY_FIELDS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
