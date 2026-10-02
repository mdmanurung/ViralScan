"""Opt-in ``--strand auto``: infer read strandedness from a small pilot (PLAN DEF-02, F-020)."""

from __future__ import annotations

import gzip
import itertools
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Optional

PILOT_READS = 1_000_000
#: A strand wins if it keeps >= TAU of the reads the unstranded run pseudoaligns.
TAU = 0.8
STRANDS = ("forward", "reverse", "unstranded")


def infer_strand(rates: dict[str, float], tau: float = TAU) -> str:
    """Pick a strand from pilot p_pseudoaligned rates keyed by STRANDS.

    reverse/forward win if rate/unstranded >= tau; the larger ratio wins and a tie
    falls to ``unstranded`` (fail toward not discarding reads).
    """
    unstr = rates["unstranded"]
    if unstr <= 0:
        return "unstranded"
    ratios = {s: rates[s] / unstr for s in ("forward", "reverse")}
    passing = {s: r for s, r in ratios.items() if r >= tau}
    if not passing:
        return "unstranded"
    best = max(passing.values())
    winners = [s for s, r in passing.items() if r == best]
    return winners[0] if len(winners) == 1 else "unstranded"


def head_fastq(src: str, dest: Path, n_reads: int) -> None:
    """Copy the first ``n_reads`` records of a (gzipped or plain) FASTQ to ``dest``."""
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt") as fin, open(dest, "w") as fout:
        fout.writelines(itertools.islice(fin, 4 * n_reads))


def _kb_pilot(
    strand: str,
    r1: Path,
    r2: Path,
    out: Path,
    index: str,
    t2g: str,
    technology: str,
    whitelist: Optional[str],
    cores: int,
) -> float:
    cmd = [
        "kb", "count", "-i", index, "-g", t2g, "-x", technology, "-t", str(cores),
        "--overwrite", "-o", str(out), "--strand", strand,
    ]  # fmt: skip
    if whitelist:
        cmd += ["-w", whitelist]
    subprocess.run(cmd + [str(r1), str(r2)], check=True, capture_output=True)
    return float(json.loads((out / "run_info.json").read_text())["p_pseudoaligned"])


def run_pilot(
    s1: str,
    s2: str,
    index: str,
    t2g: str,
    technology: str,
    whitelist: Optional[str],
    cores: int,
    n_reads: int = PILOT_READS,
    runner: Callable[..., float] = _kb_pilot,
) -> dict[str, float]:
    """Run the three explicit-strand pilots on the first ``n_reads`` pairs."""
    with tempfile.TemporaryDirectory(prefix="viralscan_strand_") as tmp:
        tmp_dir = Path(tmp)
        r1, r2 = tmp_dir / "R1.fastq", tmp_dir / "R2.fastq"
        head_fastq(s1, r1, n_reads)
        head_fastq(s2, r2, n_reads)
        return {
            s: runner(s, r1, r2, tmp_dir / s, index, t2g, technology, whitelist, cores)
            for s in STRANDS
        }


def inference_block(rates: dict[str, float], n_reads: int = PILOT_READS) -> dict:
    """The per-sample ``strand_inference`` manifest record."""
    return {
        "choice": infer_strand(rates),
        "rates": {s: rates[s] for s in STRANDS},
        "pilot_reads": n_reads,
        "tau": TAU,
    }
