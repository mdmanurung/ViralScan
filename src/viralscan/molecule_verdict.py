"""Per-molecule host-versus-virus verdicts from the competitive alignment (PLAN VERDICT-01).

``alignment_qc`` counts each read's primary alignment, so a read that scores the same on host and virus is
assigned to one of them at random, and ``host_competitive_fraction`` is one number for the whole run. A call
needs a per-molecule answer instead: for each (cell, UMI), did the best viral alignment beat the best host
alignment? ``minimap2 -ax sr`` keeps secondary alignments within 80 % of the primary score, so both sides are
in the raw competitive BAM and no new alignment pass is needed.

A molecule is ``virus_best`` when its best viral ``AS`` exceeds its best host ``AS``, ``host_best`` for the
reverse, ``tie`` when equal, and ``unaligned`` when neither aligned. Only one side present means the other
scored below the secondary cutoff, so that side wins.

Lineage rows with ``assigned_weight > 0`` are the reads the call counted (host/virus-ambiguous reads get 0
under ``host-conservative``); ``counted`` marks molecules with at least one such read.
"""

from __future__ import annotations

import csv
import gzip
import subprocess
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

#: SAM flags skipped: unmapped (0x4) and supplementary (0x800). Secondary (0x100) is kept on purpose.
SKIP_FLAGS = 0x804

VERDICTS = ("virus_best", "host_best", "tie", "unaligned")


def _score(fields: list[str]) -> int | None:
    for tag in fields[11:]:
        if tag.startswith("AS:i:"):
            return int(tag[5:])
    return None


def _read_number(qname: str) -> int | None:
    """Read number from an extracted-read name ``<CB>_<UMI>_<n>|<read id>``, else None."""
    parts = qname.split("_")
    if len(parts) < 3:
        return None
    try:
        return int(parts[2].split("|", 1)[0])
    except ValueError:
        return None


def best_scores(sam_lines: Iterable[str]) -> dict[str, dict[str, int]]:
    """Best ``AS`` per read name and side (``host`` / ``virus``), secondary alignments included."""
    scores: dict[str, dict[str, int]] = {}
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11 or int(fields[1]) & SKIP_FLAGS or fields[2] == "*":
            continue
        score = _score(fields)
        if score is None:
            continue
        side = "host" if fields[2].startswith("HOST|") else "virus"
        best = scores.setdefault(fields[0], {})
        best[side] = max(score, best.get(side, score))
    return scores


@dataclass
class _Molecule:
    reads: int = 0
    counted: bool = False
    virus_as: int | None = None
    host_as: int | None = None


def molecule_rows(sam_lines: Iterable[str], counted_reads: set[int]) -> Iterator[dict[str, object]]:
    """One row per (cell, UMI) with its best host and virus score and a verdict."""
    molecules: dict[tuple[str, str], _Molecule] = {}
    for qname, best in best_scores(sam_lines).items():
        parts = qname.split("_")
        if len(parts) < 3 or not parts[0] or not parts[1]:
            continue
        mol = molecules.setdefault((parts[0], parts[1]), _Molecule())
        mol.reads += 1
        mol.counted = mol.counted or _read_number(qname) in counted_reads
        if "virus" in best:
            mol.virus_as = (
                best["virus"] if mol.virus_as is None else max(mol.virus_as, best["virus"])
            )
        if "host" in best:
            mol.host_as = best["host"] if mol.host_as is None else max(mol.host_as, best["host"])
    for (cb, ub), mol in molecules.items():
        virus_as, host_as = mol.virus_as, mol.host_as
        if virus_as is None and host_as is None:
            verdict = "unaligned"
        elif host_as is None or (virus_as is not None and virus_as > host_as):
            verdict = "virus_best"
        elif virus_as is None or host_as > virus_as:
            verdict = "host_best"
        else:
            verdict = "tie"
        yield {
            "cb": cb,
            "ub": ub,
            "reads": mol.reads,
            "counted": mol.counted,
            "virus_as": virus_as,
            "host_as": host_as,
            "verdict": verdict,
        }


def counted_read_numbers(lineage_tsv_gz: str | Path) -> set[int]:
    """Read numbers of lineage rows the call counted (``assigned_weight > 0``)."""
    with gzip.open(lineage_tsv_gz, "rt", newline="") as handle:
        return {
            int(row["read_number"])
            for row in csv.DictReader(handle, delimiter="\t")
            if row["assigned_weight"] not in ("", "None") and float(row["assigned_weight"]) > 0
        }


def summarise(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    """Molecule counts per scope (``all`` / ``counted``) and verdict, every verdict listed."""
    counts: Counter[tuple[str, str]] = Counter()
    for row in rows:
        counts["all", str(row["verdict"])] += 1
        if row["counted"]:
            counts["counted", str(row["verdict"])] += 1
    return [
        {"scope": scope, "verdict": verdict, "molecules": counts[scope, verdict]}
        for scope in ("all", "counted")
        for verdict in VERDICTS
    ]


def write_molecule_verdicts(
    sam_lines: Iterable[str], lineage_tsv_gz: str | Path, out_dir: Path
) -> list[dict[str, object]]:
    """Write ``molecule_verdicts.tsv.gz`` and ``molecule_verdict_summary.tsv``; return the summary."""
    rows = list(molecule_rows(sam_lines, counted_read_numbers(lineage_tsv_gz)))
    fields = ["cb", "ub", "reads", "counted", "virus_as", "host_as", "verdict"]
    with gzip.open(out_dir / "molecule_verdicts.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = summarise(rows)
    with (out_dir / "molecule_verdict_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["scope", "verdict", "molecules"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary)
    return summary


def write_from_bam(
    bam: str | Path, lineage_tsv_gz: str | Path, out_dir: Path
) -> list[dict[str, object]]:
    """``write_molecule_verdicts`` on a competitive BAM, streaming ``samtools view`` (KSHV has tens of millions of reads)."""
    with subprocess.Popen(
        ["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdout is not None
        summary = write_molecule_verdicts(proc.stdout, lineage_tsv_gz, out_dir)
    if proc.returncode:
        raise RuntimeError(f"samtools view failed (exit {proc.returncode}) on {bam}")
    return summary
