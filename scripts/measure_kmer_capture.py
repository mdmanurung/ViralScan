#!/usr/bin/env python3
"""Measure how much of a viral population's k-mer space a reference panel captures.

PLAN `CAT-11`. This is the measurement behind `REF-01`'s decision to make the
expanded anellovirus panel the default, promoted from a throwaway script to a
reproducible tool.

Why this and not "count the genomes": a panel can hold thousands of genomes and
still be blind to a real isolate, because pseudoalignment matches exact k-mers.
The bundled 20-genome TTV panel held 20 genomes and shared **zero** 31-mers with
85.8 % of the 2,042 real anellovirus genomes. Genome count is not coverage.

Three numbers, in increasing order of honesty:

``coverage``
    Fraction of a target genome's distinct k-mers that also occur in the panel.
    With the target *in* the panel this is 1.0 by construction, which is why it
    is reported but never used as the gate.

``leave_one_out``
    The same fraction with the target's own genome removed from the panel. This
    is the estimate for a strain nobody has sequenced yet, and it is the number
    that matters when asking "would we see a new isolate?".

``p_fragment``
    Fraction of length-``--read-length`` windows of the target that contain at
    least one panel k-mer — i.e. P(a read pseudoaligns), leave-one-out. A panel
    can have low ``leave_one_out`` coverage and still capture every read, if the
    shared k-mers are spread evenly; this is the number that predicts detection.

Usage
-----
    python scripts/measure_kmer_capture.py \\
        --panel ref/panel.fa --population all_genomes.fa \\
        --groups src/viralscan/data/anellovirus_accessions.tsv \\
        --group-key accession --group-value genus \\
        --out results/kmer_capture

Self-check (no arguments needed):

    python scripts/measure_kmer_capture.py --self-check
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from viralscan.sensitivity import DEFAULT_K  # noqa: E402

DEFAULT_READ_LENGTH = 90
_ACGT = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    """The lexicographically smaller of a k-mer and its reverse complement.

    kallisto indexes canonical k-mers, so a panel genome deposited in the
    opposite orientation to a target still matches. Measuring forward-strand
    k-mers only therefore understates what the tool would actually detect.
    """
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(accession, sequence)``; the accession is the first header token."""
    name: Optional[str] = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def kmers(sequence: str, k: int, *, strand: str = "canonical") -> set[str]:
    """Distinct k-mers over the ACGT alphabet.

    Ambiguity codes are skipped rather than expanded: an ``N`` cannot match a
    concrete k-mer in an index either, so counting it would overstate coverage.

    ``strand="canonical"`` (the default) collapses each k-mer with its reverse
    complement, which is what kallisto's index does. ``strand="forward"``
    reproduces the pre-fix behaviour and is kept only so the superseded numbers
    stay regenerable.
    """
    seq = sequence.upper()
    fold = strand == "canonical"
    out: set[str] = set()
    for i in range(len(seq) - k + 1):
        window = seq[i : i + k]
        if _ACGT.issuperset(window):
            out.add(canonical(window) if fold else window)
    return out


def fragment_hit_fraction(
    sequence: str,
    panel: set[str],
    k: int,
    read_length: int,
    *,
    exclude: Optional[set[str]] = None,
    strand: str = "canonical",
) -> float:
    """Fraction of length-``read_length`` windows holding ≥ 1 panel k-mer.

    Computed by walking the genome once and counting, for each fragment window,
    whether any of its k-mer start positions is a panel hit. A prefix sum over
    the per-position hit flags keeps this linear rather than quadratic.

    ``exclude`` holds the k-mers that only the genome under test contributes to
    the panel, so leave-one-out is evaluated without materialising a separate
    panel set per genome (which would be quadratic in the panel size).
    """
    seq = sequence.upper()
    n_kmer_positions = len(seq) - k + 1
    if n_kmer_positions <= 0:
        return 0.0
    fold = strand == "canonical"
    hits = [0] * n_kmer_positions
    for i in range(n_kmer_positions):
        window = seq[i : i + k]
        if not _ACGT.issuperset(window):
            continue
        if fold:
            window = canonical(window)
        if window in panel and not (exclude and window in exclude):
            hits[i] = 1
    prefix = [0] * (n_kmer_positions + 1)
    for i, flag in enumerate(hits):
        prefix[i + 1] = prefix[i] + flag
    # A fragment starting at s spans k-mer start positions s .. s+read_length-k.
    kmers_per_fragment = read_length - k + 1
    if kmers_per_fragment <= 0:
        return 0.0
    n_fragments = n_kmer_positions - kmers_per_fragment + 1
    if n_fragments <= 0:
        # Genome shorter than one fragment: treat the whole genome as one window.
        return 1.0 if prefix[-1] > 0 else 0.0
    captured = sum(1 for s in range(n_fragments) if prefix[s + kmers_per_fragment] - prefix[s] > 0)
    return captured / n_fragments


def load_groups(path: Path, key_col: str, value_col: str) -> dict[str, str]:
    """Load an accession → group label map from a TSV (e.g. accession → genus)."""
    groups: dict[str, str] = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row.get(key_col) or "").strip()
            value = (row.get(value_col) or "").strip()
            if not key:
                continue
            groups[key] = value or "unclassified"
            groups.setdefault(key.split(".")[0], value or "unclassified")
    return groups


def measure(
    panel_path: Path,
    population_path: Path,
    *,
    k: int = DEFAULT_K,
    read_length: int = DEFAULT_READ_LENGTH,
    groups: Optional[dict[str, str]] = None,
    strand: str = "canonical",
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Return ``(per_genome_rows, summary)``."""
    panel_kmers_by_acc: dict[str, set[str]] = {}
    for acc, seq in read_fasta(panel_path):
        panel_kmers_by_acc[acc] = kmers(seq, k, strand=strand)
    panel_all: set[str] = set()
    for value in panel_kmers_by_acc.values():
        panel_all |= value
    # Count how many panel genomes contribute each k-mer, so leave-one-out can
    # subtract a genome's private k-mers without rebuilding the union each time.
    multiplicity: dict[str, int] = {}
    for value in panel_kmers_by_acc.values():
        for km in value:
            multiplicity[km] = multiplicity.get(km, 0) + 1

    panel_bare = {a.split(".")[0]: a for a in panel_kmers_by_acc}
    rows: list[dict[str, object]] = []
    for acc, seq in read_fasta(population_path):
        target = kmers(seq, k, strand=strand)
        if not target:
            continue
        in_panel_key = acc if acc in panel_kmers_by_acc else panel_bare.get(acc.split(".")[0])
        shared_kmers = target & panel_all
        coverage = len(shared_kmers) / len(target)
        if in_panel_key is not None:
            # Leave-one-out without rebuilding the panel: the only k-mers that
            # disappear when this genome is removed are the ones it alone
            # contributes. Rebuilding the union per genome would be quadratic in
            # the panel size and does not finish on a 2,042-genome panel.
            own = panel_kmers_by_acc[in_panel_key]
            private = {km for km in own if multiplicity.get(km, 0) == 1}
        else:
            private = set()
        loo_shared = len(shared_kmers - private) if private else len(shared_kmers)
        rows.append(
            {
                "accession": acc,
                "group": (groups or {}).get(acc, (groups or {}).get(acc.split(".")[0], "")),
                "in_panel": int(in_panel_key is not None),
                "length": len(seq),
                "n_kmers": len(target),
                "coverage": round(coverage, 6),
                "leave_one_out": round(loo_shared / len(target), 6),
                "p_fragment": round(
                    fragment_hit_fraction(
                        seq, panel_all, k, read_length, exclude=private, strand=strand
                    ),
                    6,
                ),
            }
        )

    def _mean(values: list[float]) -> float:
        return sum(values) / len(values) if values else 0.0

    def _median(values: list[float]) -> float:
        if not values:
            return 0.0
        ordered = sorted(values)
        mid = len(ordered) // 2
        if len(ordered) % 2:
            return ordered[mid]
        return (ordered[mid - 1] + ordered[mid]) / 2

    cov = [float(r["coverage"]) for r in rows]
    loo = [float(r["leave_one_out"]) for r in rows]
    pf = [float(r["p_fragment"]) for r in rows]
    summary: dict[str, object] = {
        "k": k,
        "strand": strand,
        "read_length": read_length,
        "panel_genomes": len(panel_kmers_by_acc),
        "panel_kmers": len(panel_all),
        "population_genomes": len(rows),
        "median_coverage": round(_median(cov), 6),
        "median_leave_one_out": round(_median(loo), 6),
        "mean_leave_one_out": round(_mean(loo), 6),
        "zero_coverage_fraction": round(sum(1 for v in loo if v == 0.0) / len(loo), 6)
        if loo
        else 0.0,
        "median_p_fragment": round(_median(pf), 6),
        "mean_p_fragment": round(_mean(pf), 6),
    }
    if groups:
        by_group: dict[str, list[dict[str, object]]] = {}
        for row in rows:
            by_group.setdefault(str(row["group"]) or "unclassified", []).append(row)
        summary["groups"] = {
            name: {
                "genomes": len(members),
                "median_leave_one_out": round(
                    _median([float(m["leave_one_out"]) for m in members]), 6
                ),
                "median_p_fragment": round(_median([float(m["p_fragment"]) for m in members]), 6),
                "zero_coverage_fraction": round(
                    sum(1 for m in members if float(m["leave_one_out"]) == 0.0) / len(members),
                    6,
                ),
            }
            for name, members in sorted(by_group.items())
        }
    return rows, summary


def self_check() -> int:
    """Assert the arithmetic on sequences whose answers are known by hand."""
    import random
    import tempfile

    rng = random.Random(0)
    genome = "".join(rng.choice("ACGT") for _ in range(2000))
    other = "".join(rng.choice("ACGT") for _ in range(2000))
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        panel = tmp_path / "panel.fa"
        pop = tmp_path / "pop.fa"
        panel.write_text(f">A.1 self\n{genome}\n>B.1 other\n{other}\n")
        pop.write_text(f">A.1 self\n{genome}\n")

        # A genome against a panel containing itself: coverage 1.0 by construction.
        rows, _ = measure(panel, pop)
        assert rows[0]["coverage"] == 1.0, rows[0]
        # Leave-one-out against an unrelated panel member: ~0 shared k-mers.
        assert rows[0]["leave_one_out"] < 0.01, rows[0]
        assert rows[0]["p_fragment"] < 0.01, rows[0]

        # A genome present twice is still covered leave-one-out (the duplicate
        # keeps every k-mer's multiplicity above 1).
        panel.write_text(f">A.1 self\n{genome}\n>C.1 copy\n{genome}\n")
        rows, _ = measure(panel, pop)
        assert rows[0]["leave_one_out"] == 1.0, rows[0]
        assert rows[0]["p_fragment"] == 1.0, rows[0]

        # Positive control from REF-01: a genome against itself alone is 100 %.
        panel.write_text(f">A.1 self\n{genome}\n")
        rows, _ = measure(panel, pop)
        assert rows[0]["coverage"] == 1.0

        # A panel holding only the reverse complement of the target must still
        # capture it, because kallisto indexes canonical k-mers. Measuring the
        # forward strand alone reported these genomes as undetectable; 350 of the
        # 2,042 panel anelloviruses were affected and 55 of them by >0.25
        # fragment capture.
        panel.write_text(f">R.1 revcomp\n{reverse_complement(genome)}\n")
        rows, _ = measure(panel, pop)
        assert rows[0]["coverage"] == 1.0, rows[0]
        assert rows[0]["p_fragment"] == 1.0, rows[0]
        rows, _ = measure(panel, pop, strand="forward")
        assert rows[0]["coverage"] < 0.01, rows[0]
    print("self-check passed")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panel", type=Path, help="Reference panel FASTA.")
    parser.add_argument("--population", type=Path, help="Genomes to measure against the panel.")
    parser.add_argument("--groups", type=Path, help="TSV mapping accessions to a group label.")
    parser.add_argument("--group-key", default="accession")
    parser.add_argument("--group-value", default="genus")
    parser.add_argument("--k", type=int, default=DEFAULT_K)
    parser.add_argument(
        "--strand",
        choices=("canonical", "forward"),
        default="canonical",
        help=(
            "canonical (default) folds each k-mer with its reverse complement, "
            "matching kallisto's index; forward reproduces the superseded "
            "single-strand numbers."
        ),
    )
    parser.add_argument("--read-length", type=int, default=DEFAULT_READ_LENGTH)
    parser.add_argument("--out", type=Path, help="Output prefix for .tsv and .json.")
    parser.add_argument("--self-check", action="store_true", help="Run arithmetic self-tests.")
    args = parser.parse_args(argv)

    if args.self_check:
        return self_check()
    if not args.panel or not args.population:
        parser.error("--panel and --population are required (or pass --self-check)")

    groups = load_groups(args.groups, args.group_key, args.group_value) if args.groups else None
    rows, summary = measure(
        args.panel,
        args.population,
        k=args.k,
        read_length=args.read_length,
        groups=groups,
        strand=args.strand,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        tsv_path = args.out.with_suffix(".tsv")
        with open(tsv_path, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]) if rows else ["accession"], delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)
        args.out.with_suffix(".json").write_text(json.dumps(summary, indent=2, sort_keys=True))
        print(f"wrote {tsv_path} and {args.out.with_suffix('.json')}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
