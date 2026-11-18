"""The one FASTA reader (PLAN ``SEQIO-01``).

A record is ``(identifier, sequence)``: the identifier is the first whitespace-delimited token of the
``>`` line, the sequence is the joined, stripped lines below it. Plain and ``.gz`` files read alike.
Blank lines are skipped, a line that is not a header and comes before the first header is ignored, and a
``>`` line with no identifier raises ``FastaError``.

What a caller does with a repeated identifier is its own rule, so ``iter_fasta`` yields every record and
``read_fasta`` takes ``unique=`` (``False``: the last record wins; ``True``: raise).
"""

from collections.abc import Iterator
from pathlib import Path

from viralscan.utils import open_maybe_gzip


class FastaError(ValueError):
    """A FASTA file the reader cannot make records of."""


def iter_fasta(path: str | Path, *, upper: bool = False) -> Iterator[tuple[str, str]]:
    """Yield ``(identifier, sequence)`` for each record of a FASTA, in file order."""
    name: str | None = None
    chunks: list[str] = []
    with open_maybe_gzip(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper() if upper else "".join(chunks)
                tokens = line[1:].split()
                if not tokens:
                    raise FastaError(f"FASTA header without an identifier in {path}")
                name, chunks = tokens[0], []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper() if upper else "".join(chunks)


def read_fasta(path: str | Path, *, upper: bool = False, unique: bool = False) -> dict[str, str]:
    """``{identifier: sequence}``; with ``unique=True`` a repeated identifier raises ``FastaError``."""
    records: dict[str, str] = {}
    for name, sequence in iter_fasta(path, upper=upper):
        if unique and name in records:
            raise FastaError(f"Duplicate FASTA ID {name!r} in {path}")
        records[name] = sequence
    return records
