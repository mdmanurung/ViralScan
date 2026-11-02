"""
Fail the run when R1 and R2 hold different numbers of reads.

kallisto pairs reads by position and exits 0 when R2 is shorter, so a truncated
R2 silently drops reads. This rule counts both raw inputs (no dependency on
``kb_count``) and writes ``log/fastq_pair_check.json`` with the pair count,
which ``multimap`` reads into ``count_audit.tsv`` as ``fastq_read_pairs``.
"""

import gzip
import json
from pathlib import Path


def _count_reads(path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        n_lines = sum(1 for _ in fh)  # constant memory
    if n_lines % 4:
        raise ValueError(f"{path} has {n_lines} lines, not a multiple of 4: truncated FASTQ")
    return n_lines // 4


def check_pairs(r1, r2) -> int:
    n1, n2 = _count_reads(r1), _count_reads(r2)
    if n1 != n2:
        raise ValueError(f"R1 has {n1} reads, R2 has {n2} ({r1} vs {r2})")
    return n1


if "snakemake" in globals():
    pairs = check_pairs(snakemake.config["sample1"], snakemake.config["sample2"])  # noqa: F821
    Path(snakemake.output[0]).write_text(json.dumps({"pairs": pairs}) + "\n")  # noqa: F821
