#!/usr/bin/env python3
"""DSR-02: one verdict per `viralscan evidence` directory.

    python scripts/dsr02_verdicts.py <evidence_root> > verdicts.tsv

Rules, first match wins (raw alignment_qc.tsv rows, reads summed by reference_class):
  no_support     no read aligned to host or virus
  host_best      no read aligned to a viral reference, or host_homology flagged
  low_complexity viral reads, but low_complexity flagged or complex_body_fraction < 0.5
  viral_best     otherwise
The flags alone are not enough: when no read aligns to the virus (or BLAST masks all reads) the
evidence flags read `not_flagged` although every read is host (GSM5725695, F-028), so the verdict
counts alignments itself.
`complex_body_fraction` is the share of reads WITH a templated body (`anello_align.is_complex_body`), so a LOW
value means low complexity. It was read the other way until VERDICT-01 (2026-10-07), which labelled EBV
(1.52 M clean viral reads) `low_complexity`. `tso_reads` is a diagnostic only: reads carrying the 10x TSO, the reads
STAR's 0.9 match filter lets through the host subtraction (F-028).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
COLS = ["evidence_dir", "target_reads", "virus_reads", "host_reads", "host_homology", "low_complexity",
        "sibling_or_host_ambiguity", "tso_reads", "virus_best_molecules", "host_best_molecules",
        "tie_molecules", "verdict"]

#: PROVISIONAL (VERDICT-01): a call needs at least this share of its counted molecules to score higher on the
#: virus than on the host. Frozen, or replaced, at the calibration gate on known positives and artefacts.
MIN_VIRUS_BEST_FRACTION = 0.5


def _tsv(path: Path):
    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def read_weighted_complex_fraction(virus_rows) -> float | None:
    """Reads-weighted share of viral reads with a templated body; None when no row reports one."""
    rows = [(int(r["reads"]), float(r["complex_body_fraction"])) for r in virus_rows if r["complex_body_fraction"]]
    total = sum(n for n, _ in rows)
    return sum(n * f for n, f in rows) / total if total else None


def molecule_counts(evidence_dir: Path) -> dict[str, int] | None:
    """Counted molecules per verdict from ``molecule_verdict_summary.tsv``; None for older evidence."""
    path = evidence_dir / "molecule_verdict_summary.tsv"
    if not path.is_file():
        return None
    return {r["verdict"]: int(r["molecules"]) for r in _tsv(path) if r["scope"] == "counted"}


def verdict(evidence_dir: Path) -> dict:
    qc = [r for r in _tsv(evidence_dir / "alignment_qc.tsv") if r["count_layer"] == "raw"]
    virus = [r for r in qc if r["reference_class"] == "virus"]
    virus_reads = sum(int(r["reads"]) for r in virus)
    host_reads = sum(int(r["reads"]) for r in qc if r["reference_class"] == "host")
    flags = {r["flag"]: r["status"] for r in _tsv(evidence_dir / "interpretation_flags.tsv")}
    complex_body = read_weighted_complex_fraction(virus)
    seqs = [ln.strip() for ln in (evidence_dir / "viral_reads.fasta").open() if not ln.startswith(">")]
    mols = molecule_counts(evidence_dir)
    if mols is not None:
        # Per molecule, not the run-level host_homology flag: that one is a single number for every call.
        aligned = mols["virus_best"] + mols["host_best"] + mols["tie"]
        host_side = not aligned or mols["virus_best"] / aligned < MIN_VIRUS_BEST_FRACTION
    else:
        host_side = not virus_reads or flags.get("host_homology") == "flagged"
    if (not virus_reads and not host_reads) or (mols is not None and not aligned):
        v = "no_support"
    elif host_side:
        v = "host_best"
    elif flags.get("low_complexity") == "flagged" or (complex_body is not None and complex_body < 0.5):
        v = "low_complexity"
    else:
        v = "viral_best"
    return {"evidence_dir": evidence_dir.name, "target_reads": len(seqs), "virus_reads": virus_reads,
            "host_reads": host_reads, "host_homology": flags.get("host_homology", ""),
            "low_complexity": flags.get("low_complexity", ""),
            "sibling_or_host_ambiguity": flags.get("sibling_or_host_ambiguity", ""),
            "tso_reads": sum(TSO in s for s in seqs),
            "virus_best_molecules": "" if mols is None else mols["virus_best"],
            "host_best_molecules": "" if mols is None else mols["host_best"],
            "tie_molecules": "" if mols is None else mols["tie"], "verdict": v}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("evidence_root", type=Path)
    args = ap.parse_args(argv)
    out = csv.DictWriter(sys.stdout, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    out.writeheader()
    for d in sorted(p for p in args.evidence_root.iterdir() if (p / "interpretation_flags.tsv").is_file()):
        out.writerow(verdict(d))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
