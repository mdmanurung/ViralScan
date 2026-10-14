#!/usr/bin/env python3
"""DSR-02: one verdict per `viralscan evidence` directory.

    python scripts/dsr02_verdicts.py <evidence_root> [--reference-cells <round>/reference_cells] > verdicts.tsv

Rules, first match wins (UMI-deduplicated alignment_qc.tsv rows, molecules summed by reference_class; raw-layer read counts are
reported as `virus_reads`/`host_reads` and decide nothing):
  no_support     no read aligned to host or virus
  host_best      no read aligned to a viral reference, or host_homology flagged
  low_complexity viral reads, but low_complexity flagged or complex_body_fraction < 0.5
  viral_best     otherwise
`molecule_scope` says which molecules decided: `restricted` (counted molecules of the frozen reference cells,
CELLS-01), `all_cells` (every cell the evidence run saw), `no_molecule_table` / `no_frozen_cells` (a restriction was
asked for and could not be applied; the verdict then rests on all cells and must not be read as a cell-level one).
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
import gzip
import sys
from collections import Counter
from pathlib import Path

TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
COLS = ["evidence_dir", "target_reads", "virus_reads", "host_reads", "virus_molecules", "host_molecules", "host_homology", "low_complexity",
        "sibling_or_host_ambiguity", "tso_reads", "virus_best_molecules", "host_best_molecules",
        "tie_molecules", "molecule_scope", "verdict"]
VERDICTS = ("virus_best", "host_best", "tie", "unaligned")

#: PROVISIONAL (VERDICT-01): a call needs at least this share of its counted molecules to score higher on the
#: virus than on the host. Frozen, or replaced, at the calibration gate on known positives and artefacts.
MIN_VIRUS_BEST_FRACTION = 0.5
#: PROVISIONAL, same status: below this reads-weighted share of templated bodies the viral reads are low complexity.
MIN_COMPLEX_BODY_FRACTION = 0.5


def _tsv(path: Path):
    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def read_weighted_complex_fraction(virus_rows) -> float | None:
    """Count-weighted share with a templated body; the rows are the deduplicated layer, so molecules. None when no row reports one."""
    rows = [(int(r["reads"]), float(r["complex_body_fraction"])) for r in virus_rows if r["complex_body_fraction"]]
    total = sum(n for n, _ in rows)
    return sum(n * f for n, f in rows) / total if total else None


def molecule_counts(evidence_dir: Path, cells: set[str] | None = None) -> dict[str, int] | None:
    """Counted molecules per verdict; None when the evidence has no molecule verdicts.

    Without ``cells`` this is ``molecule_verdict_summary.tsv`` (every cell the evidence run saw). With a cell
    set (the frozen reference cells, CELLS-01) only counted molecules whose barcode is in it are tallied, from
    ``molecule_verdicts.tsv.gz``; None when that table is absent, so a caller never mistakes an unrestricted
    tally for a restricted one.
    """
    if cells is None:
        path = evidence_dir / "molecule_verdict_summary.tsv"
        if not path.is_file():
            return None
        return {r["verdict"]: int(r["molecules"]) for r in _tsv(path) if r["scope"] == "counted"}
    table = evidence_dir / "molecule_verdicts.tsv.gz"
    if not table.is_file():
        return None
    counts: Counter[str] = Counter()
    with gzip.open(table, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["counted"] == "True" and r["cb"] in cells:
                counts[r["verdict"]] += 1
    return {v: counts[v] for v in VERDICTS}


def frozen_cells(reference_cells_dir: Path, evidence_name: str) -> set[str] | None:
    """Frozen reference cells (`<ds>__<sample>.tsv`, scripts/dsr06_reference_cells.py) of an evidence dir.

    Evidence directories are `<dataset>__<sample>__<arm>__<virus slug>`; None when no list is frozen.
    """
    path = reference_cells_dir / ("__".join(evidence_name.split("__")[:2]) + ".tsv")
    if not path.is_file():
        return None
    return {ln.strip() for ln in path.read_text().splitlines()[1:] if ln.strip()}


def verdict(evidence_dir: Path, cells: set[str] | None = None) -> dict:
    all_qc = list(_tsv(evidence_dir / "alignment_qc.tsv"))
    qc = [r for r in all_qc if r["count_layer"] == "deduplicated"]
    if all_qc and not qc:
        raise ValueError(f"{evidence_dir}: missing UMI-deduplicated QC; regenerate evidence")
    virus = [r for r in qc if r["reference_class"] == "virus"]
    # Decisions use the UMI-deduplicated layer: one representative per corrected (CB, UMI, reference), so these
    # are molecules. `virus_reads` / `host_reads` below are the raw layer, labelled reads, and decide nothing.
    virus_reads = sum(int(r["reads"]) for r in virus)
    host_reads = sum(int(r["reads"]) for r in qc if r["reference_class"] == "host")
    raw = [r for r in all_qc if r["count_layer"] == "raw"]
    raw_virus = sum(int(r["reads"]) for r in raw if r["reference_class"] == "virus") if raw else ""
    raw_host = sum(int(r["reads"]) for r in raw if r["reference_class"] == "host") if raw else ""
    flags = {r["flag"]: r["status"] for r in _tsv(evidence_dir / "interpretation_flags.tsv")}
    complex_body = read_weighted_complex_fraction(virus)
    seqs = [ln.strip() for ln in (evidence_dir / "viral_reads.fasta").open() if not ln.startswith(">")]
    mols = molecule_counts(evidence_dir, cells)
    # "restricted" = counted molecules of the frozen reference cells; "all_cells" = every cell the evidence run
    # saw (no frozen list given); "no_molecule_table" = a restriction was asked for but cannot be applied.
    scope = (
        "all_cells" if cells is None and mols is not None
        else "restricted" if mols is not None
        else "no_molecule_table" if cells is not None
        else ""
    )
    if mols is None and cells is not None:
        mols = molecule_counts(evidence_dir)  # unrestricted fallback, labelled by `molecule_scope`
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
    elif flags.get("low_complexity") == "flagged" or (complex_body is not None and complex_body < MIN_COMPLEX_BODY_FRACTION):
        v = "low_complexity"
    else:
        v = "viral_best"
    return {"evidence_dir": evidence_dir.name, "target_reads": len(seqs), "virus_reads": raw_virus,
            "host_reads": raw_host, "virus_molecules": virus_reads, "host_molecules": host_reads, "host_homology": flags.get("host_homology", ""),
            "low_complexity": flags.get("low_complexity", ""),
            "sibling_or_host_ambiguity": flags.get("sibling_or_host_ambiguity", ""),
            "tso_reads": sum(TSO in s for s in seqs),
            "virus_best_molecules": "" if mols is None else mols["virus_best"],
            "host_best_molecules": "" if mols is None else mols["host_best"],
            "tie_molecules": "" if mols is None else mols["tie"], "molecule_scope": scope, "verdict": v}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("evidence_root", type=Path)
    ap.add_argument("--reference-cells", type=Path, help="<round>/reference_cells: count only molecules of "
                    "the frozen reference cells (needs molecule_verdicts.tsv.gz in each evidence dir)")
    args = ap.parse_args(argv)
    out = csv.DictWriter(sys.stdout, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    out.writeheader()
    for d in sorted(p for p in args.evidence_root.iterdir() if (p / "interpretation_flags.tsv").is_file()):
        cells = frozen_cells(args.reference_cells, d.name) if args.reference_cells else None
        row = verdict(d, cells)
        if args.reference_cells and cells is None:
            row["molecule_scope"] = "no_frozen_cells"  # said so, not silently an all-cells tally
        out.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
