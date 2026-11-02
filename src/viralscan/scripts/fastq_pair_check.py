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


if "snakemake" in globals():
    smk = snakemake  # type: ignore[name-defined]  # noqa: F821
    r1, r2 = smk.config["sample1"], smk.config["sample2"]
    result: dict[str, object] = {"pairs": None, "skipped": "raw FASTQs missing"}
    if Path(r1).exists() and Path(r2).exists():
        result = {"pairs": check_pairs(r1, r2)}
    else:
        print(f"WARNING fastq_pair_check skipped: {r1} or {r2} no longer exists.")
    Path(smk.output[0]).write_text(json.dumps(result) + "\n")
