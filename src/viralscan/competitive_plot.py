"""Read-by-read host-versus-virus alignment figure for a ``viralscan evidence`` directory (PLAN ALIGN-PLOT-01).

``plot_competitive_alignment`` answers "what did each read of this call actually align to?". Panel A draws one
lane per read (full read length, original read orientation) with a virus row and a host row, clipped bases grey,
poly-A/T runs hatched and the 10x TSO outlined, and the verdict plus both ``AS`` values at the right. Panel B is
the virus-AS-versus-host-AS scatter and panel C the human loci (1 kb bins) the reads hit.
``cross_call_loci`` pools the alignments of several evidence directories into one TSV and one figure, so a host
locus or viral position that recurs across calls and libraries stands out.

Inputs are read, never modified: ``competitive_reads.raw.bam`` (``samtools view``; supplementary records and
``SA`` tags both count), ``read_verdicts.tsv.gz`` (written by ``scripts/dsr02_validate_call.py``; derived from the
BAM when absent) and ``viral_reads.fasta`` (read sequence, for poly-A/T and TSO and for reads that never aligned).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
import subprocess
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from viralscan.evidence import reference_class, to_int
from viralscan.molecule_verdict import VERDICTS, deduplicate_competitive_alignments, verdict_of

log = logging.getLogger("viralscan")

#: 10x 5' template-switch oligo, the same motif ``scripts/dsr02_verdicts.py`` labels ``repeat_class == "tso"``.
TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
POLY_MIN = 10
MAX_LANES = 60
#: Gate of the Q137 aligned-length rule (read_support.MIN_ALIGNED_BASES), drawn on panel D.
LENGTH_GATE = 60
#: Per-lane text (loci, verdict) only fits while lanes are tall enough to read.
LOCI_TEXT_MAX_LANES = 45
VERDICT_TEXT_MAX_LANES = 60
HOST_BIN = 1000
VIRUS_BIN = 100
FIGURE_STEM = "competitive_alignment"

# Okabe-Ito, colour-blind safe: the verdict colour is the lane-segment colour of the winning side.
COLOURS = {
    "virus": "#D55E00",
    "virus_other": "#F2A97A",
    "host": "#0072B2",
    "host_other": "#8CC3E8",
    "clip": "#E3E3E3",
    "tso": "#009E73",
    "tie": "#CC79A7",
    "unaligned": "#8A8A8A",
}
VERDICT_COLOURS = {
    "virus_best": COLOURS["virus"],
    "host_best": COLOURS["host"],
    "tie": COLOURS["tie"],
    "unaligned": COLOURS["unaligned"],
}
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_COMP = str.maketrans("ACGTN", "TGCAN")


@dataclass
class Segment:
    """One alignment of a read, in the read's original orientation (``q0``/``q1`` are 0-based, half open)."""

    ref: str
    pos: int
    strand: str
    q0: int
    q1: int
    read_length: int
    kind: str  # primary | supplementary | secondary
    score: int | None
    matched: int = 0  # M/=/X bases only

    @property
    def side(self) -> str:
        return "host" if reference_class(self.ref) == "host" else "virus"

    @property
    def label(self) -> str:
        return f"{self.ref.split('|', 1)[-1]}:{self.pos:,}"


@dataclass
class ReadLane:
    name: str
    length: int
    virus_as: int | None
    host_as: int | None
    verdict: str
    repeat_class: str
    counted: str = ""
    segments: list[Segment] = field(default_factory=list)
    sequence: str = ""

    @property
    def library(self) -> str:
        """Run accession from ``<CB>_<UMI>_<n>|<SRR>.<read>``; empty for other read names."""
        return self.name.split("|", 1)[1].split(".")[0] if "|" in self.name else ""

    def viral_aligned(self) -> int:
        """M/=/X bases of the best-scoring viral record (supplementary excluded), the Q137 length; 0 if none."""
        hits = [s for s in self.segments if s.side == "virus" and s.kind != "supplementary"]
        return max(hits, key=lambda s: (s.score or 0, s.matched)).matched if hits else 0

    def poly_runs(self) -> list[tuple[int, int, str]]:
        """``(start, end, base)`` of poly-A/T runs >= ``POLY_MIN`` nt; one-base interruptions are bridged."""
        runs = []
        for base in "AT":
            merged: list[list[int]] = []
            for hit in re.finditer(f"{base}{{5,}}", self.sequence):
                if merged and hit.start() - merged[-1][1] <= 1:
                    merged[-1][1] = hit.end()
                else:
                    merged.append([hit.start(), hit.end()])
            runs += [(a, b, base) for a, b in merged if b - a >= POLY_MIN]
        return sorted(runs)

    def tso_span(self) -> tuple[int, int] | None:
        if self.repeat_class != "tso":
            return None
        for motif in (TSO, TSO.translate(_COMP)[::-1]):
            at = self.sequence.find(motif)
            if at >= 0:
                return at, at + len(motif)
        return None

    def sort_key(self) -> tuple[int, str, int, str]:
        wanted = {"virus_best": "virus", "host_best": "host"}.get(self.verdict)
        primaries = [s for s in self.segments if s.kind == "primary"]
        best = next((s for s in primaries if s.side == wanted), None) or next(iter(primaries), None)
        ref, pos = (best.ref, best.pos) if best else ("~", 0)
        return (
            VERDICTS.index(self.verdict) if self.verdict in VERDICTS else len(VERDICTS),
            ref,
            pos,
            self.name,
        )


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _segment(
    ref: str, pos: int, reverse: bool, cigar: str, kind: str, score: int | None
) -> Segment | None:
    """Aligned query interval of one CIGAR, flipped into the original read orientation for ``-`` strand hits."""
    ops = [(int(n), op) for n, op in _CIGAR.findall(cigar)]
    if not ops:
        return None
    lead = trail = aligned = matched = 0
    seen_aligned = False
    for n, op in ops:
        if op in "MI=X":
            aligned += n
            matched += n if op != "I" else 0
            seen_aligned = True
        elif op in "SH":
            if seen_aligned:
                trail += n
            else:
                lead += n
    length = lead + aligned + trail
    q0, q1 = lead, lead + aligned
    if reverse:
        q0, q1 = length - q1, length - q0
    return Segment(ref, pos, "-" if reverse else "+", q0, q1, length, kind, score, matched)


def _parse_sam(lines: Iterable[str]) -> dict[str, ReadLane]:
    """Per-read segments from SAM text: every mapped record, plus ``SA`` entries no record already carries."""
    reads: dict[str, ReadLane] = {}
    pending_sa: list[tuple[str, str]] = []
    for line in lines:
        if not line or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        lane = reads.setdefault(f[0], ReadLane(f[0], 0, None, None, "unaligned", "none"))
        if not flag & 0x900 and f[9] != "*":
            seq = f[9].upper()
            lane.sequence = seq.translate(_COMP)[::-1] if flag & 16 else seq
            lane.length = len(seq)
        if flag & 4 or f[2] == "*":
            continue
        score = next((int(t[5:]) for t in f[11:] if t.startswith("AS:i:")), None)
        kind = "supplementary" if flag & 0x800 else "secondary" if flag & 0x100 else "primary"
        seg = _segment(f[2], int(f[3]), bool(flag & 16), f[5], kind, score)
        if seg:
            lane.segments.append(seg)
        pending_sa += [(f[0], t[5:]) for t in f[11:] if t.startswith("SA:Z:")]
    for name, sa in pending_sa:
        lane = reads[name]
        for entry in filter(None, sa.split(";")):
            ref, pos, strand, cigar = entry.split(",")[:4]
            if any(
                s.ref == ref and s.pos == int(pos) and s.strand == strand for s in lane.segments
            ):
                continue
            seg = _segment(ref, int(pos), strand == "-", cigar, "supplementary", None)
            if seg:
                lane.segments.append(seg)
    return reads


def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def _opt_int(value: str | None) -> int | None:
    return None if value in (None, "", "None") else to_int(value)


def load_lanes(evidence_dir: str | Path, layer: str = "umi") -> list[ReadLane]:
    """One ``ReadLane`` per read of an evidence directory (verdict table first, BAM-derived if it is absent)."""
    d = Path(evidence_dir)
    if layer not in ("umi", "raw"):
        raise ValueError("layer must be umi or raw")
    bam, verdicts, fasta = (
        d / "competitive_reads.raw.bam",
        d / "read_verdicts.tsv.gz",
        d / "viral_reads.fasta",
    )
    if not bam.is_file() and not verdicts.is_file():
        raise FileNotFoundError(
            f"{d} has neither competitive_reads.raw.bam nor read_verdicts.tsv.gz"
        )
    aligned: dict[str, ReadLane] = {}
    if bam.is_file():
        with subprocess.Popen(
            ["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True
        ) as proc:
            assert proc.stdout is not None
            aligned = _parse_sam(
                deduplicate_competitive_alignments(proc.stdout) if layer == "umi" else proc.stdout
            )
        if proc.returncode:
            raise RuntimeError(f"samtools view failed (exit {proc.returncode}) on {bam}")
    else:
        if layer == "umi":
            raise FileNotFoundError(f"{bam} is required for post-alignment UMI deduplication")
        log.warning("%s has no competitive_reads.raw.bam; lanes show no alignments", d)
    seqs = _read_fasta(fasta) if fasta.is_file() else {}
    lanes: list[ReadLane] = []
    if verdicts.is_file():
        with gzip.open(verdicts, "rt", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            if layer == "umi" and row["read"] not in aligned:
                continue
            lane = aligned.get(
                row["read"], ReadLane(row["read"], 0, None, None, "unaligned", "none")
            )
            lane.length = to_int(row.get("length")) or lane.length
            lane.virus_as, lane.host_as = (
                _opt_int(row.get("virus_as")),
                _opt_int(row.get("host_as")),
            )
            lane.verdict, lane.repeat_class = (
                row.get("verdict") or "unaligned",
                row.get("repeat_class") or "none",
            )
            lane.counted = row.get("counted", "")
            lanes.append(lane)
    else:
        for (
            lane
        ) in aligned.values():  # same rule as molecule_verdict: best non-supplementary AS per side
            for side in ("virus", "host"):
                scores = [
                    s.score
                    for s in lane.segments
                    if s.side == side and s.kind != "supplementary" and s.score is not None
                ]
                setattr(lane, f"{side}_as", max(scores, default=None))
            lane.verdict = verdict_of(lane.virus_as, lane.host_as)
            lanes.append(lane)
    for lane in lanes:
        lane.sequence = seqs.get(lane.name, lane.sequence)
        lane.length = max(
            [lane.length, len(lane.sequence), *(s.read_length for s in lane.segments)]
        )
        if not verdicts.is_file() and TSO in lane.sequence + lane.sequence.translate(_COMP)[::-1]:
            lane.repeat_class = "tso"
    return lanes


def _bin_start(pos: int, size: int) -> int:
    return (pos - 1) // size * size + 1


def locus_reads(lanes: Iterable[ReadLane]) -> dict[tuple[str, str, int], set[str]]:
    """Reads per ``(side, reference, bin start)``: 1 kb bins on the host, 100 bp on viral references."""
    hits: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    for lane in lanes:
        for s in lane.segments:
            hits[
                s.side,
                s.ref.split("|", 1)[-1],
                _bin_start(s.pos, HOST_BIN if s.side == "host" else VIRUS_BIN),
            ].add(lane.name)
    return hits


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _loci_text(lane: ReadLane) -> str:
    parts = []
    for side, tag in (("virus", "V"), ("host", "H")):
        segs = sorted(
            (s for s in lane.segments if s.side == side), key=lambda s: s.kind != "primary"
        )
        if segs:
            parts.append(
                f"{tag} "
                + " ".join(
                    f"{s.label}{s.strand} {s.q1 - s.q0}nt{'' if s.kind == 'primary' else ' ' + s.kind[:4] + '.'}"
                    for s in segs[:3]
                )
            )
    text = "   ".join(parts) or "no alignment"
    return text if len(text) <= 100 else text[:99] + "…"


def _draw_lanes(ax: Any, lanes: list[ReadLane], total: int, unit: str = "molecules") -> None:
    from matplotlib.patches import Rectangle

    n = len(lanes)
    xmax = max((lane.length for lane in lanes), default=100)
    show_loci = n <= LOCI_TEXT_MAX_LANES
    ax.set_xlim(-2, xmax * (2.0 if show_loci else 1.45))
    ax.set_ylim(n - 0.4, -0.6)
    shade = {"primary": 0, "supplementary": 1, "secondary": 2}
    for i, lane in enumerate(lanes):
        ax.add_patch(
            Rectangle((0, i - 0.4), lane.length, 0.8, fc=COLOURS["clip"], ec="#9A9A9A", lw=0.3)
        )
        for s in sorted(lane.segments, key=lambda s: -shade[s.kind]):  # primary last = on top
            colour = COLOURS[s.side if s.kind == "primary" else f"{s.side}_other"]
            y0 = i - 0.4 if s.side == "virus" else i
            ax.add_patch(Rectangle((s.q0, y0), s.q1 - s.q0, 0.4, fc=colour, ec="none"))
        for a, b, base in lane.poly_runs():
            ax.add_patch(
                Rectangle(
                    (a, i - 0.4),
                    b - a,
                    0.8,
                    fill=False,
                    hatch="//" if base == "A" else "\\\\",
                    ec="#333333",
                    lw=0,
                )
            )
        tso = lane.tso_span()
        if tso:
            ax.add_patch(
                Rectangle(
                    (tso[0], i - 0.4), tso[1] - tso[0], 0.8, fill=False, ec=COLOURS["tso"], lw=1.4
                )
            )
        if i and lane.verdict != lanes[i - 1].verdict:
            ax.axhline(i - 0.5, color="#555555", lw=0.6, xmax=0.5 if show_loci else 0.7)
        ax.add_patch(
            Rectangle(
                (xmax * 1.02, i - 0.4),
                xmax * 0.02,
                0.8,
                fc=VERDICT_COLOURS.get(lane.verdict, "#000"),
                ec="none",
            )
        )
        if n <= VERDICT_TEXT_MAX_LANES:
            v, h = (lane.virus_as, lane.host_as)
            ax.text(
                xmax * 1.06,
                i,
                f"{lane.verdict:<10} V={'-' if v is None else v:<4} H={'-' if h is None else h}",
                va="center",
                fontsize=6.5,
                family="monospace",
            )
        if show_loci:
            ax.text(xmax * 1.3, i, _loci_text(lane), va="center", fontsize=6.5, family="monospace")
            ax.text(
                -xmax * 0.012,
                i,
                lane.name.split("|")[0].split("_")[0][:8]
                + "…"
                + lane.name.split("|")[0].split("_")[-1],
                va="center",
                ha="right",
                fontsize=6,
                color="#444444",
            )
    ax.set_yticks([])
    ax.set_xticks(range(0, xmax + 1, 10 if xmax <= 150 else 50))
    ax.tick_params(axis="x", labelsize=7)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_bounds(0, xmax)
    ax.set_xlabel(
        "position along the read (nt, original read orientation); top row virus, bottom row host",
        fontsize=8,
    )
    note = f"{total:,} {unit}"
    ax.set_title(
        f"A   lanes, one per {unit[:-1]}: {note}", loc="left", fontsize=10, fontweight="bold"
    )


def _legend_handles() -> list[object]:
    from matplotlib.patches import Patch

    return [
        Patch(fc=COLOURS["virus"], label="virus, primary"),
        Patch(fc=COLOURS["virus_other"], label="virus, suppl./secondary"),
        Patch(fc=COLOURS["host"], label="host, primary"),
        Patch(fc=COLOURS["host_other"], label="host, suppl./secondary"),
        Patch(fc=COLOURS["clip"], ec="#9A9A9A", label="clipped / unaligned"),
        Patch(fill=False, hatch="//", ec="#333333", label=f"poly-A (>= {POLY_MIN} nt)"),
        Patch(fill=False, hatch="\\\\", ec="#333333", label=f"poly-T (>= {POLY_MIN} nt)"),
        Patch(fill=False, ec=COLOURS["tso"], lw=1.4, label="10x TSO"),
    ]


def _draw_scatter(ax: Any, lanes: list[ReadLane], unit: str = "molecules") -> None:
    import random

    rng = random.Random(0)  # fixed seed: identical points are jittered the same way every run
    top = max([ln.virus_as or 0 for ln in lanes] + [ln.host_as or 0 for ln in lanes] + [60]) * 1.08
    ax.plot([0, top], [0, top], color="#666666", lw=0.8, ls="--", zorder=1)
    ax.axhline(0, color="#BBBBBB", lw=0.6)
    ax.axvline(0, color="#BBBBBB", lw=0.6)
    counts = Counter(ln.verdict for ln in lanes)
    for verdict in VERDICTS:
        pts = [ln for ln in lanes if ln.verdict == verdict and (ln.virus_as or ln.host_as)]
        ax.scatter(
            [(ln.virus_as or 0) + rng.uniform(-1.2, 1.2) for ln in pts],
            [(ln.host_as or 0) + rng.uniform(-1.2, 1.2) for ln in pts],
            s=14,
            alpha=0.5,
            lw=0,
            color=VERDICT_COLOURS[verdict],
            label=f"{verdict} (n={counts[verdict]:,})",
            zorder=2,
        )
    ax.set_xlim(-top * 0.05, top)
    ax.set_ylim(-top * 0.05, top)
    ax.set(xlabel="virus AS (0 = no viral hit)", ylabel="host AS (0 = no host hit)")
    ax.set_title(f"B   score per {unit[:-1]}", loc="left", fontsize=10, fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.set_aspect("equal", adjustable="box")


def _bar_panel(
    ax: Any,
    labels: list[str],
    counts: list[int],
    colour: str,
    title: str,
    empty: str,
    unit: str = "reads",
) -> None:
    from matplotlib.ticker import MaxNLocator

    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    if not counts:
        ax.text(0.5, 0.5, empty, ha="center", va="center", transform=ax.transAxes, color="#555555")
        ax.set_axis_off()
        return
    ax.barh(range(len(counts)), counts, color=colour, height=0.5)
    ax.set_ylim(max(len(counts), 3) - 0.5, -0.5)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_yticks(range(len(counts)), labels, fontsize=7, family="monospace")
    for i, c in enumerate(counts):
        ax.text(c, i, f" {c:,}", va="center", fontsize=7)
    ax.set_xlim(0, max(counts) * 1.2)
    ax.set_xlabel(unit, fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _draw_length_hist(ax: Any, lanes: list[ReadLane], unit: str) -> None:
    """Panel D: viral aligned length per molecule (or read), stacked by verdict, with the 60-base gate."""
    from matplotlib.ticker import MaxNLocator

    ax.set_title(
        f"D   viral aligned length per {unit[:-1]}", loc="left", fontsize=10, fontweight="bold"
    )
    lengths = [ln.viral_aligned() for ln in lanes]
    top = max([*lengths, LENGTH_GATE]) + 5
    bins = list(range(0, top + 5, 5))
    ax.hist(
        [[v for v, ln in zip(lengths, lanes, strict=True) if ln.verdict == k] for k in VERDICTS],
        bins=bins,
        stacked=True,
        color=[VERDICT_COLOURS[k] for k in VERDICTS],
        label=list(VERDICTS),
    )
    ax.axvline(LENGTH_GATE, color="#222222", lw=1, ls="--")
    below = sum(v < LENGTH_GATE for v in lengths)
    ax.text(
        LENGTH_GATE,
        1.0,
        f"gate {LENGTH_GATE} bases: {below:,} of {len(lengths):,} {unit} below ",
        transform=ax.get_xaxis_transform(),
        va="top",
        ha="right",
        fontsize=7,
    )
    ax.set_xlabel("viral aligned bases (M/=/X, no clips; 0 = no viral alignment)", fontsize=8)
    ax.set_ylabel(unit, fontsize=8)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_competitive_alignment(
    evidence_dir: str | Path,
    out_dir: str | Path | None = None,
    max_lanes: int = MAX_LANES,
    layer: str = "umi",
) -> list[str]:
    """Write ``competitive_alignment.png`` and ``.pdf`` for one evidence directory; return the paths.

    *out_dir* defaults to the evidence directory; the evidence files themselves are only read.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    evidence_dir = Path(evidence_dir)
    lanes = sorted(load_lanes(evidence_dir, layer), key=ReadLane.sort_key)
    unit = "molecules" if layer == "umi" else "reads (raw, PCR duplicates included)"
    unit_short = "molecules" if layer == "umi" else "reads"
    with_lanes = (
        0 < len(lanes) <= max_lanes
    )  # the lane block is unreadable above ~60 lanes: panel D shows them
    n = len(lanes) if with_lanes else 0
    lane_h = max(0.04, min(0.32, 8 / n)) if n else 0.32
    a_h = 1.6 + max(n, 3) * lane_h
    fig = plt.figure(figsize=(16, 0.7 + (a_h if with_lanes else 0) + 4.6), layout="constrained")
    rows = [0.7, *([a_h] if with_lanes else []), 4.6]
    grid = fig.add_gridspec(len(rows), 3, height_ratios=rows, width_ratios=[1, 1, 1.2])
    ax_leg = fig.add_subplot(grid[0, :])
    ax_a = fig.add_subplot(grid[1, :]) if with_lanes else None
    ax_b, ax_c, ax_d = (fig.add_subplot(grid[-1, i]) for i in range(3))
    ax_leg.axis("off")
    ax_leg.legend(handles=_legend_handles(), loc="center", ncol=8, frameon=False, fontsize=8)
    fig.suptitle(
        f"Competitive alignment, human vs virus, counted as {unit}: {evidence_dir.name}",
        fontsize=12,
        fontweight="bold",
    )
    loci = Counter(
        {key: len(names) for key, names in locus_reads(lanes).items() if key[0] == "host"}
    )
    if lanes:
        if ax_a is not None:
            _draw_lanes(ax_a, lanes, len(lanes), unit_short)
        _draw_scatter(ax_b, lanes, unit_short)
        _draw_length_hist(ax_d, lanes, unit_short)
    else:
        for ax, title in (
            (ax_b, f"B   score per {unit_short[:-1]}"),
            (ax_d, f"D   viral aligned length per {unit_short[:-1]}"),
        ):
            ax.text(
                0.5,
                0.5,
                f"no {unit_short} in this evidence directory",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="#555555",
            )
            ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
            ax.set_axis_off()
    top_loci = sorted(loci.items(), key=lambda kv: (-kv[1], kv[0][1], kv[0][2]))[:10]
    _bar_panel(
        ax_c,
        [f"{ref}:{start:,}" for (_, ref, start), _ in top_loci],
        [c for _, c in top_loci],
        COLOURS["host"],
        f"C   human loci hit ({HOST_BIN // 1000} kb bins; {len(loci)} loci, top {len(top_loci)})",
        "no host alignments",
        unit_short,
    )
    out = Path(out_dir) if out_dir is not None else evidence_dir
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "pdf"):
        path = out / f"{FIGURE_STEM}.{ext}"
        with plt.rc_context({"hatch.linewidth": 0.5}):
            fig.savefig(path, dpi=160)
        paths.append(str(path))
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------
# Cross-call view
# ---------------------------------------------------------------------------

CROSS_CALL_FIELDS = (
    "kind",
    "reference",
    "bin_start",
    "bin_end",
    "n_reads",
    "count_layer",
    "n_molecules",
    "n_calls",
    "n_libraries",
    "libraries",
    "calls",
)


def cross_call_loci(
    evidence_dirs: Iterable[str | Path],
    out_prefix: str | Path,
    top: int = 20,
    layer: str = "umi",
) -> list[dict[str, Any]]:
    """Pool alignments over several evidence directories; write ``<prefix>.tsv`` and ``<prefix>.png``.

    One TSV row per host locus (1 kb bin) or viral position (100 bp bin): distinct reads, calls (evidence
    directories) and libraries (run accession in the read name) hitting it, recurrent loci first. The figure shows
    the top ``top`` host loci and viral positions as bars of reads, stacked by library.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dirs = [Path(d) for d in evidence_dirs]
    reads: dict[tuple[str, str, int], set[str | tuple[str, ...]]] = defaultdict(set)
    calls: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    per_lib: dict[tuple[str, str, int], Counter[str]] = defaultdict(Counter)
    for d in dirs:
        lanes = load_lanes(d, layer)
        library_of = {lane.name: lane.library or d.name for lane in lanes}
        identities = {
            lane.name: (library_of[lane.name], *lane.name.split("_")[:2])
            if layer == "umi"
            else lane.name
            for lane in lanes
        }
        for key, names in locus_reads(lanes).items():
            calls[key].add(d.name)
            for name in names:
                identity = identities[name]
                if identity not in reads[key]:
                    reads[key].add(identity)
                    per_lib[key][library_of[name]] += 1
    rows: list[dict[str, Any]] = []
    for key in reads:
        kind, ref, start = key
        size = HOST_BIN if kind == "host" else VIRUS_BIN
        rows.append(
            {
                "kind": kind,
                "reference": ref,
                "bin_start": start,
                "bin_end": start + size - 1,
                "n_reads": len(reads[key]),
                "count_layer": layer,
                "n_molecules": len(reads[key]) if layer == "umi" else "",
                "n_calls": len(calls[key]),
                "n_libraries": len(per_lib[key]),
                "libraries": ";".join(f"{k}={v}" for k, v in sorted(per_lib[key].items())),
                "calls": ";".join(sorted(calls[key])),
            }
        )
    rows.sort(
        key=lambda r: (
            r["kind"] != "host",
            -int(r["n_libraries"]),
            -int(r["n_reads"]),
            str(r["reference"]),
            int(r["bin_start"]),
        )
    )
    prefix = Path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with Path(f"{prefix}.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CROSS_CALL_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    libraries = sorted({lib for c in per_lib.values() for lib in c})
    palette = [
        "#0072B2",
        "#E69F00",
        "#009E73",
        "#CC79A7",
        "#56B4E9",
        "#D55E00",
        "#F0E442",
        "#999999",
    ]
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, 0.9 + 0.3 * max(1, min(top, len(rows))) + 1.8),
        layout="constrained",
        squeeze=False,
    )
    for ax, kind, title in (
        (axes[0][0], "host", f"Host loci ({HOST_BIN // 1000} kb bins)"),
        (axes[0][1], "virus", f"Viral positions ({VIRUS_BIN} bp bins)"),
    ):
        sel = [r for r in rows if r["kind"] == kind][:top]
        ax.set_title(
            f"{title}: {len(dirs)} calls, {len(libraries)} libraries",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        if not sel:
            ax.text(
                0.5,
                0.5,
                f"no {kind} alignments",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="#555555",
            )
            ax.set_axis_off()
            continue
        left = [0] * len(sel)
        for j, lib in enumerate(libraries):
            vals = [
                per_lib[(kind, str(r["reference"]), int(r["bin_start"]))].get(lib, 0) for r in sel
            ]
            if any(vals):
                ax.barh(
                    range(len(sel)),
                    vals,
                    left=left,
                    color=palette[j % len(palette)],
                    height=0.7,
                    label=lib,
                )
                left = [a + b for a, b in zip(left, vals)]
        ax.set_yticks(
            range(len(sel)),
            [f"{r['reference']}:{int(r['bin_start']):,}" for r in sel],
            fontsize=7,
            family="monospace",
        )
        ax.invert_yaxis()
        for i, r in enumerate(sel):
            ax.text(
                left[i],
                i,
                f" {r['n_reads']} {'molecules' if layer == 'umi' else 'raw reads'}, {r['n_libraries']} lib, {r['n_calls']} call",
                va="center",
                fontsize=7,
            )
        ax.set_xlim(0, max(left) * 1.45)
        ax.set_xlabel(
            f"distinct {'molecules' if layer == 'umi' else 'raw reads'} (stacked by library)",
            fontsize=8,
        )
        ax.legend(frameon=False, fontsize=7, title="library", title_fontsize=7, loc="lower right")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(
        "Loci hit across calls (recurrent loci sort first)", fontsize=12, fontweight="bold"
    )
    fig.savefig(Path(f"{prefix}.png"), dpi=160)
    plt.close(fig)
    return rows


def run_evidence_plot(args: argparse.Namespace) -> None:
    """``viralscan evidence-plot``: one figure per evidence directory, optionally the cross-call view."""
    if getattr(args, "coverage_map", False):
        from viralscan.evidence import coverage_map_figure

        if args.cross_call:
            raise ValueError("--cross-call cannot be combined with --coverage-map")
        for directory in map(Path, args.evidence_dir):
            out = Path(args.output) / directory.name if args.output else None
            coverage_map_figure(directory, args.target, args.annotation, out)
        return
    if getattr(args, "target", None) or getattr(args, "annotation", None):
        raise ValueError("--target and --annotation require --coverage-map")
    for d in map(Path, args.evidence_dir):
        out = Path(args.output) / d.name if args.output else None
        for path in plot_competitive_alignment(d, out, args.max_lanes, args.layer):
            log.info("Wrote %s", path)
    if args.cross_call:
        rows = cross_call_loci(args.evidence_dir, args.cross_call, layer=args.layer)
        log.info("Cross-call loci: %d rows -> %s.tsv/.png", len(rows), args.cross_call)
