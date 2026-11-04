"""Gene-relative strand of viral reads, calibrated per library (PLAN ``EVID-03``).

``alignment_qc.tsv`` reports ``reverse_strand_fraction``: the share of reads aligned to the reverse strand of the
reference. That number says nothing about genes. A 10x 3' library puts R2 reads on the transcript's own strand, a
10x 5' library puts them on the opposite strand, so "antisense" is only meaningful against what the same library
does to genes that are known to be unambiguous. This module writes ``strand_qc.tsv`` and leaves every other
evidence output untouched.

Unit (``docs/evidence_unit.md``): the member reads of ``counted_virus_best`` molecules, found by joining each BAM
read to ``read_lineage.tsv.gz`` on ``read_number``. Never the representative read. A read's alignment is its best
viral record (highest ``AS``, secondary records included; ties go to the smallest reference and position). An
aligned block is a stretch of the alignment between ``N`` gaps, so spliced reads give several blocks. Each block is
classified by the strand of the annotated features it overlaps:

* ``sense``: every overlapped feature is on the read's strand;
* ``antisense``: every overlapped feature is on the opposite strand;
* ``ambiguous_overlap``: overlapped features on both strands (or with no strand);
* ``outside_annotation``: no overlapped feature.

Viral features are the ``gene`` rows of the panel GTF (``exon`` then ``CDS`` rows when it has no gene rows).

Calibration. The expected orientation comes from the library's host reads: member reads of ``counted_host_best``
molecules, best host record, blocks that overlap exons of exactly one gene on one strand. If at least
:data:`MIN_HOST_BLOCKS` such blocks exist and at least :data:`MIN_ORIENTATION_FRACTION` of them agree on one
orientation, that orientation is the expected one (``sense`` for a 3' library, ``antisense`` for a 5' one);
the viral blocks are then reported as matching it or not. These host reads are only those that the kallisto
index sent to the target virus and that aligned better to the host, so the calibration set is small and can be
too small to use.

``status`` is ``assessed`` or ``not_assessed`` (with ``reason``), or ``not_evaluable`` when no counted
virus-best molecule has a viral alignment. ``not_assessed`` reasons: ``chemistry_unknown`` (the run's
technology is not one of :data:`viralscan.chemistry.CHEMISTRIES`; nothing is computed),
``verdicts_not_produced_in_competitor_mode`` (nothing is computed), ``host_gtf_missing``,
``too_few_host_blocks``, ``host_orientation_mixed``. For the last three the viral class counts are still written,
the expected-orientation columns are blank.
"""

from __future__ import annotations

import csv
import gzip
import re
import subprocess
from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

from viralscan.chemistry import CHEMISTRIES
from viralscan.molecule_verdict import SKIP_FLAGS

#: Fewest unambiguous host blocks a library needs before its orientation is trusted.
MIN_HOST_BLOCKS = 200
#: Share of host blocks that must agree for a library to count as having one orientation.
MIN_ORIENTATION_FRACTION = 0.8
CLASSES = ("sense", "antisense", "ambiguous_overlap", "outside_annotation")
ALL_VIRAL = "ALL_VIRAL"
VIRAL_PREFIXES = ("VIRUS|", "TARGET|", "SAMEGROUP|", "RELATED|", "DECOY|")
CATEGORIES = (
    "counted_virus_best",
    "counted_host_best",
    "counted_tie",
    "counted_unaligned",
    "not_counted",
)
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_GENE_ID = re.compile(r'gene_id\s+"([^"]*)"')

FIELDS = (
    "reference",
    "status",
    "reason",
    "technology",
    "expected_orientation",
    "host_blocks_used",
    "host_blocks_sense",
    "host_blocks_antisense",
    "host_blocks_ambiguous",
    "host_blocks_outside",
    "host_sense_fraction",
    *(f"blocks_{c}" for c in CLASSES),
    *(f"reads_{c}" for c in CLASSES),
    *(f"molecules_{c}" for c in CLASSES),
    "blocks_expected_orientation",
    "blocks_unexpected_orientation",
    "fraction_expected_orientation",
    *(f"n_molecules_{c}" for c in CATEGORIES),
    *(f"n_reads_{c}" for c in CATEGORIES),
    "n_reads_virus_best_without_viral_alignment",
    "verdict_ties_resolvable",
)


class FeatureIndex:
    """Annotated intervals of one contig, searchable by overlap (1-based inclusive coordinates)."""

    def __init__(self, features: Iterable[tuple[int, int, str, str]]) -> None:
        rows = sorted(features)
        self.starts = [r[0] for r in rows]
        self.ends = [r[1] for r in rows]
        self.strands = [r[2] for r in rows]
        self.genes = [r[3] for r in rows]
        self.prefix_max_end: list[int] = []
        top = 0
        for end in self.ends:
            top = max(top, end)
            self.prefix_max_end.append(top)

    def overlapping(self, start: int, end: int) -> list[tuple[str, str]]:
        """``(strand, gene)`` of every feature overlapping ``[start, end]``."""
        hits: list[tuple[str, str]] = []
        i = bisect_right(self.starts, end) - 1
        while i >= 0 and self.prefix_max_end[i] >= start:
            if self.ends[i] >= start:
                hits.append((self.strands[i], self.genes[i]))
            i -= 1
        return hits


def load_features(
    gtf_paths: Iterable[str | Path], kinds: Sequence[str], seqids: set[str] | None = None
) -> dict[str, FeatureIndex]:
    """Index GTF rows by seqid, using the first of *kinds* that a seqid has rows for.

    *seqids* limits what is read (a host GTF has millions of rows); ``None`` reads all.
    """
    by_kind: dict[str, dict[str, list[tuple[int, int, str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for path in gtf_paths:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9 or f[2] not in kinds or (seqids is not None and f[0] not in seqids):
                    continue
                match = _GENE_ID.search(f[8])
                by_kind[f[0]][f[2]].append(
                    (int(f[3]), int(f[4]), f[6], match.group(1) if match else f"{f[0]}:{f[3]}")
                )
    index: dict[str, FeatureIndex] = {}
    for seqid, kinds_found in by_kind.items():
        kind = next(k for k in kinds if k in kinds_found)
        index[seqid] = FeatureIndex(kinds_found[kind])
    return index


def aligned_blocks(pos: int, cigar: str) -> list[tuple[int, int]]:
    """1-based inclusive reference blocks of an alignment; an ``N`` gap starts a new block."""
    blocks: list[tuple[int, int]] = []
    start, cursor = pos, pos
    for n, op in ((int(n), op) for n, op in _CIGAR.findall(cigar)):
        if op in "M=XD":
            cursor += n
        elif op == "N":
            if cursor > start:
                blocks.append((start, cursor - 1))
            cursor += n
            start = cursor
    if cursor > start:
        blocks.append((start, cursor - 1))
    return blocks


def classify_block(
    reverse: bool, hits: Sequence[tuple[str, str]], single_gene: bool = False
) -> str:
    """Class of one block from the features it overlaps; *single_gene* also requires one gene."""
    if not hits:
        return "outside_annotation"
    strands = {s for s, _ in hits}
    if len(strands) != 1 or strands == {"."} or (single_gene and len({g for _, g in hits}) > 1):
        return "ambiguous_overlap"
    return "sense" if strands == {"-" if reverse else "+"} else "antisense"


class _Hit(NamedTuple):
    """One alignment record; tuple order is the tie-break order (best ``AS`` first)."""

    neg_score: int
    ref: str
    pos: int
    flag: int
    cigar: str

    def blocks(self) -> list[tuple[int, int]]:
        return aligned_blocks(self.pos, self.cigar)


def _best_records(
    sam_lines: Iterable[str], wanted: Mapping[int, tuple[str, str]], side: str
) -> dict[int, _Hit]:
    """Best record per wanted read on *side* (``host`` / ``virus``).

    Every record of a wanted read is a candidate except unmapped and supplementary ones.
    """
    best: dict[int, _Hit] = {}
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or int(f[1]) & SKIP_FLAGS or f[2] == "*":
            continue
        number = _number(f[0])
        if number not in wanted or f[2].startswith("HOST|") != (side == "host"):
            continue
        score = next((int(t[5:]) for t in f[11:] if t.startswith("AS:i:")), None)
        if score is None:
            continue
        hit = _Hit(-score, f[2], int(f[3]), int(f[1]), f[5])
        if number not in best or hit < best[number]:
            best[number] = hit
    return best


def _number(qname: str) -> int:
    return int(qname.split("_", 2)[2].split("|", 1)[0])


def _blank(**extra: object) -> dict[str, object]:
    return {field: "" for field in FIELDS} | extra


def _fraction(a: int, b: int) -> float | str:
    return a / b if b else ""


def strand_qc_rows(
    sam_lines: Iterable[str],
    verdicts: Iterable[Mapping[str, str]],
    lineage: Iterable[Mapping[str, str]],
    *,
    technology: str | None,
    host_features: Mapping[str, FeatureIndex] | None,
    viral_features: Mapping[str, FeatureIndex],
) -> list[dict[str, object]]:
    """Rows of ``strand_qc.tsv`` from competitive SAM lines, verdict rows and lineage rows.

    *host_features* is None when no host GTF was given. Feature keys are bare GTF seqids; BAM references
    lose their ``HOST|`` / ``VIRUS|``-style prefix before lookup.
    """
    sam = list(sam_lines)
    category: dict[tuple[str, str], str] = {}
    for vrow in verdicts:
        category[vrow["cb"], vrow["ub"]] = (
            f"counted_{vrow['verdict']}" if vrow["counted"] == "True" else "not_counted"
        )
    members: dict[int, tuple[str, str]] = {}
    n_reads: Counter[str] = Counter()
    for lrow in lineage:
        key = (lrow["cb"], lrow["ub"])
        members[int(lrow["read_number"])] = key
        n_reads[category.get(key, "not_counted")] += 1
    n_molecules = Counter(category.values())
    counts = {
        **{f"n_molecules_{c}": n_molecules[c] for c in CATEGORIES},
        **{f"n_reads_{c}": n_reads[c] for c in CATEGORIES},
    }
    base = _blank(technology=technology or "", verdict_ties_resolvable=True, **counts)
    if (technology or "").lower() not in CHEMISTRIES:
        return [
            base | {"reference": ALL_VIRAL, "status": "not_assessed", "reason": "chemistry_unknown"}
        ]

    def reads_of(cat: str) -> dict[int, tuple[str, str]]:
        return {n: k for n, k in members.items() if category.get(k) == cat}

    virus_reads, host_reads = reads_of("counted_virus_best"), reads_of("counted_host_best")
    viral_best = _best_records(sam, virus_reads, "virus")
    host_best = _best_records(sam, host_reads, "host") if host_features is not None else {}
    base["n_reads_virus_best_without_viral_alignment"] = len(virus_reads) - len(viral_best)
    if not viral_best:
        return [
            base
            | {
                "reference": ALL_VIRAL,
                "status": "not_evaluable",
                "reason": "no_counted_virus_best_alignment",
            }
        ]

    host: Counter[str] = Counter()
    for hit in host_best.values():
        idx = (host_features or {}).get(hit.ref.removeprefix("HOST|"))
        for start, end in hit.blocks():
            host[
                classify_block(
                    bool(hit.flag & 0x10),
                    idx.overlapping(start, end) if idx else [],
                    single_gene=True,
                )
            ] += 1
    used = host["sense"] + host["antisense"]
    sense_fraction = _fraction(host["sense"], used)
    reason, expected = "", ""
    if host_features is None:
        reason = "host_gtf_missing"
    elif used < MIN_HOST_BLOCKS:
        reason = "too_few_host_blocks"
    elif max(host["sense"], host["antisense"]) / used < MIN_ORIENTATION_FRACTION:
        reason = "host_orientation_mixed"
    else:
        expected = "sense" if host["sense"] > host["antisense"] else "antisense"
    calibration = {
        "host_blocks_used": used,
        "host_blocks_sense": host["sense"],
        "host_blocks_antisense": host["antisense"],
        "host_blocks_ambiguous": host["ambiguous_overlap"],
        "host_blocks_outside": host["outside_annotation"],
        "host_sense_fraction": sense_fraction,
        "expected_orientation": expected,
    }

    blocks: Counter[tuple[str, str]] = Counter()
    read_sets: dict[tuple[str, str], set[int]] = defaultdict(set)
    mol_sets: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
    for number, hit in viral_best.items():
        bare = hit.ref.split("|", 1)[1] if hit.ref.startswith(VIRAL_PREFIXES) else hit.ref
        idx = viral_features.get(hit.ref) or viral_features.get(bare)
        for start, end in hit.blocks():
            cls = classify_block(bool(hit.flag & 0x10), idx.overlapping(start, end) if idx else [])
            for ref in (hit.ref, ALL_VIRAL):
                blocks[ref, cls] += 1
                read_sets[ref, cls].add(number)
                mol_sets[ref, cls].add(virus_reads[number])

    rows: list[dict[str, object]] = []
    for ref in sorted({r for r, _ in blocks}, key=lambda r: (r == ALL_VIRAL, r)):
        row = base | calibration | {"reference": ref}
        for cls in CLASSES:
            row[f"blocks_{cls}"] = blocks[ref, cls]
            row[f"reads_{cls}"] = len(read_sets[ref, cls])
            row[f"molecules_{cls}"] = len(mol_sets[ref, cls])
        row["status"] = "assessed" if expected else "not_assessed"
        row["reason"] = reason
        if expected:
            other = "antisense" if expected == "sense" else "sense"
            row["blocks_expected_orientation"] = blocks[ref, expected]
            row["blocks_unexpected_orientation"] = blocks[ref, other]
            row["fraction_expected_orientation"] = _fraction(
                blocks[ref, expected], blocks[ref, expected] + blocks[ref, other]
            )
        rows.append(row)
    return rows


def write_strand_qc(
    bam: str | Path,
    out_dir: Path,
    *,
    technology: str | None,
    host_gtf: str | Path | None,
    viral_gtfs: Iterable[str | Path],
    competitor_mode: bool = False,
) -> list[dict[str, object]]:
    """Write ``strand_qc.tsv`` in *out_dir* from the competitive BAM and the verdict and lineage tables."""
    if competitor_mode:
        rows = [
            _blank(reference=ALL_VIRAL, status="not_assessed", technology=technology or "",
                   reason="verdicts_not_produced_in_competitor_mode")
        ]  # fmt: skip
    else:
        with subprocess.Popen(
            ["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True
        ) as proc:
            assert proc.stdout is not None
            sam = list(proc.stdout)
        if proc.returncode:
            raise RuntimeError(f"samtools view failed (exit {proc.returncode}) on {bam}")
        references = {line.split("\t", 3)[2] for line in sam if not line.startswith("@")}
        host_ids = {r.removeprefix("HOST|") for r in references if r.startswith("HOST|")}
        host_features = (
            load_features([host_gtf], ("exon",), host_ids) if host_gtf is not None else None
        )
        viral_features = load_features(viral_gtfs, ("gene", "exon", "CDS"))
        with (
            gzip.open(out_dir / "molecule_verdicts.tsv.gz", "rt", newline="") as v,
            gzip.open(out_dir / "read_lineage.tsv.gz", "rt", newline="") as lin,
        ):
            rows = strand_qc_rows(
                sam,
                csv.DictReader(v, delimiter="\t"),
                csv.DictReader(lin, delimiter="\t"),
                technology=technology,
                host_features=host_features,
                viral_features=viral_features,
            )
    with (out_dir / "strand_qc.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows
