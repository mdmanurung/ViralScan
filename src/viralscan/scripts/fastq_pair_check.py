"""
Fail the run when R1 and R2 hold different numbers of reads.

kallisto pairs reads by position and exits 0 when R2 is shorter, so a truncated
R2 silently drops reads. This rule counts both raw inputs (no dependency on
``kb_count``) and writes ``log/fastq_pair_check.json`` with the pair count.
An old run resumed after its FASTQs were deleted records ``"pairs": null``.
"""

import gzip
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any


def _count_reads(path: str | Path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        n_lines = sum(1 for _ in fh)  # constant memory
    if n_lines % 4:
        raise ValueError(f"{path} has {n_lines} lines, not a multiple of 4: truncated FASTQ")
    return n_lines // 4


def check_pairs(r1: str | Path, r2: str | Path) -> int:
    n1, n2 = _count_reads(r1), _count_reads(r2)
    if n1 != n2:
        raise ValueError(f"R1 has {n1} reads, R2 has {n2} ({r1} vs {r2})")
    return n1


def run(r1: str | Path, r2: str | Path, done_file: str | Path) -> None:
    """Entry point: count both raw FASTQs and write the pair-count receipt."""
    result: dict[str, object] = {"pairs": None, "skipped": "raw FASTQs missing"}
    if Path(r1).exists() and Path(r2).exists():
        result = {"pairs": check_pairs(r1, r2)}
    else:
        print(f"WARNING fastq_pair_check skipped: {r1} or {r2} no longer exists.")
    Path(done_file).write_text(json.dumps(result) + "\n")


if TYPE_CHECKING:
    snakemake: Any  # injected by the Snakemake runner; declared for the type checker only

if "snakemake" in globals():
    run(
        snakemake.config["sample1"],  # noqa: F821
        snakemake.config["sample2"],  # noqa: F821
        snakemake.output[0],  # noqa: F821
    )
