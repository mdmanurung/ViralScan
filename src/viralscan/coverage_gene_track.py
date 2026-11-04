"""Gene track of the ``evidence-plot --coverage-map`` figure (PLAN EVID-05).

Pure layout logic, no matplotlib: GTF rows become gene models, the most covered models get
labels, and labels are packed into rows so none overlap. ``evidence.coverage_map_figure``
draws the result. The rules come from the round 2026-11-02 prototype (CDS-only genes reach
the gene row's UTR ends, repeat copies split, long spliced transcripts left out).
"""

import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

#: Transcripts or genes longer than this are left out and counted. EBV latency transcripts
#: span most of the genome and would hide every other box; 20 kb keeps all lytic genes.
MAX_FEATURE_BP = 20_000
#: CDS blocks of one gene further apart than this are separate copies (genes annotated in
#: both terminal repeats of HHV-6 or HSV), not one spliced gene.
COPY_GAP = 10_000
#: Genes labelled on the figure. Eight labels stay readable on a 10-inch panel even for the
#: crowded combined EBV GTF; every drawn gene stays in the JSON sidecar.
TOP_N_LABELS = 8
#: Label rows per strand lane. Labels that need a fourth row are dropped and counted.
MAX_LABEL_ROWS = 3
#: Width of one label character as a fraction of the genome axis (7 pt text on 10 inches,
#: with the axis about 0.85 of the figure). Slightly generous, so labels never touch.
LABEL_CHAR_FRACTION = 0.0085
#: Empty space kept between two labels in one row, as a fraction of the genome.
LABEL_PAD_FRACTION = 0.005


@dataclass(frozen=True)
class GtfRow:
    """One GTF row of the target contig (1-based, inclusive)."""

    feature: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]


@dataclass(frozen=True)
class GeneModel:
    """One drawn gene: blocks (exons, else CDS, else the gene span) of one copy."""

    name: str
    strand: str
    blocks: tuple[tuple[int, int], ...]
    size_bp: int

    @property
    def start(self) -> int:
        return self.blocks[0][0]

    @property
    def end(self) -> int:
        return max(end for _, end in self.blocks)


@dataclass(frozen=True)
class PlacedLabel:
    """A label with its row (0 is nearest the gene box) and clamped centre."""

    index: int
    row: int
    centre: float
    low: float
    high: float


def read_gtf_rows(path: str | Path, seqids: Iterable[str]) -> list[GtfRow]:
    """GTF rows of the given seqids with a feature type the gene track uses."""
    wanted = set(seqids)
    rows: list[GtfRow] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9 or fields[0] not in wanted:
                continue
            if fields[2] in {"gene", "transcript", "exon", "CDS"}:
                attributes = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
                rows.append(
                    GtfRow(fields[2], int(fields[3]), int(fields[4]), fields[6], attributes)
                )
    return rows


def _label(attributes: dict[str, str], feature: str) -> str:
    """Name, else ``gene``, else the last token of the gene ID (HUM_HERP6B_U14 is U14)."""
    name = attributes.get("gene_name") or attributes.get("gene")
    if name:
        return name
    gene_id = attributes.get("gene_id")
    return gene_id.split("_")[-1] if gene_id else feature


def build_gene_models(
    rows: Iterable[GtfRow],
    length: int,
    max_bp: int | None = None,
    copy_gap: int | None = None,
) -> tuple[list[GeneModel], int]:
    """Gene models and the number of models left out for being longer than ``max_bp``.

    A gene with transcript rows is drawn per transcript; its gene row is skipped. Blocks come
    from exons, else CDS, else the gene or transcript span. CDS-only genes with blocks more
    than ``copy_gap`` apart split into copies, and each copy reaches out to the ends of the
    gene row that contains it (its UTRs). A row with start > end wraps the origin and is
    split at the reference end; its size counts both parts.
    """
    max_bp = MAX_FEATURE_BP if max_bp is None else max_bp
    copy_gap = COPY_GAP if copy_gap is None else copy_gap
    materialised = list(rows)
    with_transcripts = {
        r.attributes.get("gene_id") for r in materialised if r.feature == "transcript"
    }
    groups: dict[str, dict[str, list[tuple[int, int]]]] = {}
    meta: dict[str, tuple[str, str]] = {}
    wrapping: set[str] = set()
    for number, row in enumerate(materialised):
        gene_id = row.attributes.get("gene_id")
        if row.feature == "gene" and gene_id in with_transcripts:
            continue
        per_transcript = gene_id in with_transcripts and row.feature != "gene"
        key = row.attributes.get("transcript_id") if per_transcript else gene_id
        key = key or f"{row.feature}:{number}"
        meta.setdefault(key, (_label(row.attributes, row.feature), row.strand))
        group = groups.setdefault(key, {"exon": [], "CDS": [], "span": []})
        kind = "span" if row.feature in {"gene", "transcript"} else row.feature
        if row.start > row.end:
            wrapping.add(key)
            group[kind] += [(1, row.end), (row.start, length)]
        else:
            group[kind].append((row.start, row.end))
    models: list[GeneModel] = []
    n_long = 0
    for key, group in groups.items():
        name, strand = meta[key]
        blocks = sorted(group["exon"] or group["CDS"] or group["span"])
        if not blocks:
            continue
        copies = [[blocks[0]]]
        for block in blocks[1:]:
            if (
                not group["exon"]
                and key not in wrapping
                and block[0] - copies[-1][-1][1] > copy_gap
            ):
                copies.append([block])
            else:
                copies[-1].append(block)
        for copy in copies:
            if not group["exon"] and key not in wrapping:
                for s0, s1 in group["span"]:
                    if s0 <= copy[0][0] and copy[-1][1] <= s1 and s1 - s0 + 1 < copy_gap:
                        copy[0] = (s0, copy[0][1])
                        copy[-1] = (copy[-1][0], s1)
            if key in wrapping:
                size, copy = sum(e - s + 1 for s, e in copy), copy[::-1]  # origin side last
            else:
                size = max(e for _, e in copy) - copy[0][0] + 1
            if size > max_bp:
                n_long += 1
                continue
            models.append(GeneModel(name, strand, tuple(copy), size))
    return sorted(models, key=lambda m: (m.start, m.end, m.name)), n_long


def mean_block_depth(depth: Any, bin_size: int, blocks: Iterable[tuple[int, int]]) -> float:
    """Mean of binned depth over the bases of the blocks (1-based, inclusive)."""
    total = bases = 0.0
    for start, end in blocks:
        for position_bin in range((start - 1) // bin_size, (end - 1) // bin_size + 1):
            low = max(start, position_bin * bin_size + 1)
            high = min(end, (position_bin + 1) * bin_size)
            total += depth[position_bin] * (high - low + 1)
            bases += high - low + 1
    return total / bases if bases else 0.0


def select_labelled(coverage: Sequence[float], top_n: int = TOP_N_LABELS) -> list[int]:
    """Indices of the ``top_n`` most covered genes, best first; uncovered genes are never named."""
    ranked = sorted((i for i, c in enumerate(coverage) if c > 0), key=lambda i: (-coverage[i], i))
    return ranked[:top_n]


def place_labels(
    wanted: Sequence[tuple[int, float, float]],
    length: int,
    max_rows: int = MAX_LABEL_ROWS,
    pad: float | None = None,
) -> tuple[list[PlacedLabel], list[int]]:
    """Greedy, collision-free label rows; ``wanted`` is (index, centre, width) in priority order.

    Each label takes the first row where it keeps ``pad`` bp clear of every label already in
    that row; its centre is clamped so it stays inside 1..length. A label that fits in none
    of the ``max_rows`` rows is returned in the second list.
    """
    gap = LABEL_PAD_FRACTION * length if pad is None else pad
    rows: list[list[tuple[float, float]]] = [[] for _ in range(max_rows)]
    placed: list[PlacedLabel] = []
    unplaced: list[int] = []
    for index, centre, width in wanted:
        half = min(width, length) / 2
        centre = min(max(centre, 1 + half), length - half)
        low, high = centre - half, centre + half
        row = next(
            (
                r
                for r, taken in enumerate(rows)
                if all(high + gap <= lo or hi + gap <= low for lo, hi in taken)
            ),
            None,
        )
        if row is None:
            unplaced.append(index)
            continue
        rows[row].append((low, high))
        placed.append(PlacedLabel(index, row, centre, low, high))
    return placed, unplaced


def label_width(name: str, length: int) -> float:
    """Width of a label in reference bp."""
    return (len(name) + 1) * LABEL_CHAR_FRACTION * length
