#!/usr/bin/env python3
"""DSR-02: read-validate ONE viral call, for any dataset, from arguments (no hand-edited copies).

Chain: the run's equivalence classes -> bus records of the call -> exact reads (UMI exact, CB as in the run)
-> competitive minimap2 against GRCh38 + panel -> per-read verdict (virus_best / host_best / tie / unaligned,
plus a repeat class). Steps 1-4 are `viralscan evidence` (which also writes molecule verdicts and alignment QC);
this script adds the per-read table with the repeat class, which is what showed F-025 (an HHV-6B "call" that
was human telomere repeats). Dataset-specific values are arguments or table rows, never code.

    # one call
    python scripts/dsr02_validate_call.py --run-dir <sample run dir> --virus "<label or gene id>" \\
        --panel-fasta <viral fasta of that run's index> --host-fasta <GRCh38 fasta> --out <evidence dir>
    # row N of calls.tsv (scripts/dsr02_enumerate_calls.py output), e.g. a SLURM array task
    python scripts/dsr02_validate_call.py --calls calls.tsv --row $SLURM_ARRAY_TASK_ID \\
        --panel-fasta <fasta> [--panel-for DATASET=<fasta> ...] --host-fasta <fasta>
    # recompute only the per-read table from an existing evidence directory
    python scripts/dsr02_validate_call.py --postprocess-only --out <evidence dir>

Outputs next to the `viralscan evidence` files: `read_verdicts.tsv.gz` (one row per read),
`read_verdict_summary.tsv` (verdict x repeat class) and `code_sha.txt` (the code that produced the evidence).
Then `scripts/dsr02_verdicts.py <evidence root>` gives the per-call verdict.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:  # run from a checkout without `pip install -e .`
    sys.path.insert(0, str(REPO / "src"))

from viralscan.evidence import _is_low_complexity  # noqa: E402
from viralscan.molecule_verdict import (  # noqa: E402
    counted_read_numbers,
    side_score,
    verdict_of,
)
from viralscan.run_safety import software_identity  # noqa: E402

READ_COLS = [
    "read",
    "cb",
    "umi",
    "read_number",
    "length",
    "virus_as",
    "host_as",
    "verdict",
    "repeat_class",
    "counted",
]
#: PROVISIONAL label rule (annotation only, it gates no verdict): a read is `telomeric` when at least this
#: share of its bases lie in TTAGGG / CCCTAA hexamers (either strand, any frame) and it has two or more.
TELOMERIC_MIN_FRACTION = 0.5
REPEAT_CLASSES = ("tso", "telomeric", "low_complexity", "none")
_COMP = str.maketrans("ACGT", "TGCA")
#: Evidence-directory files of ~GB size that the TSVs and BAMs make redundant.
BULKY = ("competitive_host_target.fasta", "blast/competitive_db*")


def _tso() -> str:
    spec = importlib.util.spec_from_file_location(
        "dsr02_verdicts", REPO / "scripts" / "dsr02_verdicts.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return str(mod.TSO)


TSO = _tso()
TSO_RC = TSO.translate(_COMP)[::-1]


def repeat_class(seq: str) -> str:
    """tso > telomeric > low_complexity > none; strand-agnostic (BAM stores reverse-strand reads reversed)."""
    seq = seq.upper()
    if TSO in seq or TSO_RC in seq:
        return "tso"
    hexamers = seq.count("TTAGGG") + seq.count("CCCTAA")
    if seq and hexamers >= 2 and 6 * hexamers / len(seq) >= TELOMERIC_MIN_FRACTION:
        return "telomeric"
    return "low_complexity" if _is_low_complexity(seq) else "none"


def read_verdict_rows(
    sam_lines: Iterable[str], counted_reads: set[int] | None = None
) -> Iterator[dict[str, object]]:
    """One row per read: best host and virus ``AS``, verdict, repeat class. Streams SAM text.

    Reads without a usable alignment on either side are `unaligned`. Sequence comes from the
    primary record (secondary records carry ``*``).
    """
    # ponytail: per-read dicts in memory (the molecule table does the same); chunk if a target ever has > 50 M reads
    best: dict[str, dict[str, int]] = {}
    seqs: dict[str, str] = {}
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        if not int(f[1]) & 0x900 and f[9] != "*":
            seqs[f[0]] = f[9]
        hit = side_score(f)
        if hit is not None:
            side, score = hit
            cur = best.setdefault(f[0], {})
            cur[side] = max(score, cur.get(side, score))
    for name, seq in seqs.items():
        parts = name.split("_")
        cb, umi = (parts[0], parts[1]) if len(parts) >= 3 else ("", "")
        try:
            number: int | None = int(parts[2].split("|", 1)[0]) if len(parts) >= 3 else None
        except ValueError:
            number = None
        scores = best.get(name, {})
        yield {
            "read": name,
            "cb": cb,
            "umi": umi,
            "read_number": "" if number is None else number,
            "length": len(seq),
            "virus_as": scores.get("virus", ""),
            "host_as": scores.get("host", ""),
            "verdict": verdict_of(scores.get("virus"), scores.get("host")),
            "repeat_class": repeat_class(seq),
            "counted": "" if counted_reads is None else int(number in counted_reads),
        }


def summarise(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    counts: Counter[tuple[str, str]] = Counter()
    for r in rows:
        counts[str(r["verdict"]), str(r["repeat_class"])] += 1
    return [
        {"verdict": v, "repeat_class": c, "reads": n} for (v, c), n in sorted(counts.items())
    ]


def write_read_verdicts(
    sam_lines: Iterable[str], out_dir: Path, lineage: Path | None = None
) -> list[dict[str, object]]:
    """Write ``read_verdicts.tsv.gz`` and ``read_verdict_summary.tsv`` into ``out_dir``."""
    counted = counted_read_numbers(lineage) if lineage and lineage.is_file() else None
    rows = list(read_verdict_rows(sam_lines, counted))
    with gzip.open(out_dir / "read_verdicts.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=READ_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    summary = summarise(rows)
    with (out_dir / "read_verdict_summary.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["verdict", "repeat_class", "reads"], delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(summary)
    return summary


def postprocess(out_dir: Path, *, keep_intermediates: bool = False) -> list[dict[str, object]]:
    """Per-read verdicts from the competitive BAM, the code stamp, and the bulky-file cleanup."""
    bam = out_dir / "competitive_reads.raw.bam"
    if not bam.is_file():
        raise SystemExit(f"{bam} not found: run `viralscan evidence --viral-fasta --host-fasta` first")
    with subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        summary = write_read_verdicts(proc.stdout, out_dir, out_dir / "read_lineage.tsv.gz")
    if proc.returncode:
        raise RuntimeError(f"samtools view failed (exit {proc.returncode}) on {bam}")
    ident = software_identity()
    stamp = ident.get("git_commit") or f"build_sha256:{ident.get('build_sha256', 'unknown')}"
    (out_dir / "code_sha.txt").write_text(f"{stamp}\n")
    if not keep_intermediates:
        for pattern in BULKY:
            for path in out_dir.glob(pattern):
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
    return summary


def slug(virus: str) -> str:
    """Same rule as `scripts/dsr05_dossier.py`, which looks the directory up by this name."""
    return re.sub(r"[^A-Za-z0-9]+", "_", virus).strip("_")


def call_row(calls_tsv: Path, row: int) -> dict[str, str]:
    """Data row ``row`` (0-based) of a `dsr02_enumerate_calls.py` table."""
    with calls_tsv.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not 0 <= row < len(rows):
        raise SystemExit(f"--row {row} is outside {calls_tsv} ({len(rows)} calls)")
    return rows[row]


def evidence_command(
    run_dir: str, virus: str, out: str, panel_fasta: str, host_fasta: str, cores: int = 8
) -> list[str]:
    """`viralscan evidence` argv for one call (list form, never a shell string)."""
    return [
        sys.executable,
        "-m",
        "viralscan.menu",
        "evidence",
        "--run-dir",
        run_dir,
        "-o",
        out,
        "--virus",
        virus,
        "--viral-fasta",
        panel_fasta,
        "--host-fasta",
        host_fasta,
        "--blast",
        "--cell-tags",
        "--cores",
        str(cores),
    ]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--calls", type=Path, help="calls.tsv from dsr02_enumerate_calls.py (with --row)")
    src.add_argument("--run-dir", help="a finished ViralScan sample run directory (with --virus)")
    ap.add_argument("--row", type=int, default=0, help="0-based data row of --calls (array task id)")
    ap.add_argument("--virus", help="exact label, accession or gene id (with --run-dir)")
    ap.add_argument("--panel-fasta", help="viral FASTA of the run's index (default for every dataset)")
    ap.add_argument(
        "--panel-for",
        action="append",
        default=[],
        metavar="DATASET=FASTA",
        help="panel FASTA for one dataset of --calls, overriding --panel-fasta (repeatable)",
    )
    ap.add_argument("--host-fasta", help="host genome FASTA (GRCh38 primary assembly)")
    ap.add_argument("--out", type=Path, help="evidence directory (single call, or --postprocess-only)")
    ap.add_argument("--out-root", type=Path, help="with --calls: <root>/<dataset>__<sample>__<arm>__<virus slug>")
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--force", action="store_true", help="redo a finished evidence directory")
    ap.add_argument("--keep-intermediates", action="store_true", help="keep the ~GB competitive FASTA and BLAST db")
    ap.add_argument("--postprocess-only", action="store_true", help="only (re)write the per-read verdicts of --out")
    ap.add_argument("--dry-run", action="store_true", help="print the `viralscan evidence` command and stop")
    args = ap.parse_args(argv)

    if args.postprocess_only:
        if args.out is None:
            ap.error("--postprocess-only needs --out")
        postprocess(args.out, keep_intermediates=True)
        return 0
    if args.calls is not None:
        call = call_row(args.calls, args.row)
        run_dir, virus = call["run_dir"], call["virus"]
        panels = dict(p.split("=", 1) for p in args.panel_for)
        panel = panels.get(call["dataset"], args.panel_fasta)
        root = args.out_root or args.calls.parent
        out = root / f"{call['dataset']}__{call['sample']}__{call['arm']}__{slug(virus)}"
    elif args.run_dir:
        if not args.virus or args.out is None:
            ap.error("--run-dir needs --virus and --out")
        run_dir, virus, panel, out = args.run_dir, args.virus, args.panel_fasta, args.out
    else:
        ap.error("give --calls (with --row) or --run-dir")
    if not panel or not args.host_fasta:
        ap.error("--host-fasta and a panel FASTA (--panel-fasta, or --panel-for for this dataset) are required")
    cmd = evidence_command(run_dir, virus, str(out), panel, args.host_fasta, args.cores)
    if args.dry_run:
        print(shlex.join(cmd))
        return 0
    if (out / "interpretation_flags.tsv").is_file() and not args.force:
        print(f"already done: {out}", file=sys.stderr)
    else:
        subprocess.run(cmd, check=True)
    postprocess(out, keep_intermediates=args.keep_intermediates)
    print(f"DONE {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
