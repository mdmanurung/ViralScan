"""Shared FASTQ test-fixture helper.

Imported by test_host_filter and test_evidence (and any future test that needs
to synthesise small FASTQ files in a tmp_path).  Keeping it here avoids a
pytest-fixture dependency so call sites stay as plain function calls.
"""

from __future__ import annotations

import gzip
from pathlib import Path


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    """Write *records* (name, seq) as a minimal FASTQ file at *path*.

    If *path* ends in ``.gz`` the file is gzip-compressed.  Quality scores
    are all ``I`` (Phred 40) so they satisfy FASTQ parsers.
    """
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
