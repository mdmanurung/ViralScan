"""One SAM layer: CIGAR parsing, tag reads and the ``samtools view`` stream (PLAN SAM-01)."""

import re
import subprocess
from collections.abc import Iterator

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    """``[(length, op), ...]`` of one CIGAR string; empty for ``*``."""
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def tag_value(fields: list[str], prefix: str) -> str | None:
    """Text after the first optional tag that starts with *prefix* (e.g. ``"NH:i:"``)."""
    for tag in fields[11:]:
        if tag.startswith(prefix):
            return tag[len(prefix) :]
    return None


def int_tag(fields: list[str], name: str, *, strict: bool = False) -> int | None:
    """Integer value of optional tag *name* (e.g. ``"AS"``), or None if absent.

    A value that is not an integer gives None, or ``ValueError`` when *strict*
    (use strict where a malformed score must stop the run, not drop a record).
    """
    value = tag_value(fields, f"{name}:i:")
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        if strict:
            raise
        return None


def view_lines(
    bam: str | object, *opts: str, region: str | None = None, samtools: str = "samtools"
) -> Iterator[str]:
    """SAM text lines from ``samtools view`` in constant memory.

    Raises ``subprocess.CalledProcessError`` when samtools exits non-zero.
    """
    command = [samtools, "view", *opts, str(bam), *([region] if region else [])]
    with subprocess.Popen(
        command, stdout=subprocess.PIPE, text=True, encoding="utf-8", bufsize=1 << 20
    ) as proc:
        assert proc.stdout is not None
        yield from proc.stdout
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, proc.args)


def view_fields(
    bam: str | object, *opts: str, region: str | None = None, samtools: str = "samtools"
) -> Iterator[list[str]]:
    """Tab-split fields of each record from :func:`view_lines`."""
    for line in view_lines(bam, *opts, region=region, samtools=samtools):
        yield line.rstrip("\n").split("\t")
