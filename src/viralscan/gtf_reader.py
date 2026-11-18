"""The one GTF reader (PLAN ``SEQIO-01``): data rows, attributes, gene IDs.

Two attribute rules exist in the package and give different answers on real inputs, so both are kept
and named; the PLAN row ``SEQIO-01`` lists the disagreement for the maintainer to settle:

- ``parse_attributes`` reads every ``key "value"`` pair into a dict (``\\w+`` key, one or more spaces,
  quoted value, empty allowed; the last of a repeated key wins; unquoted values are not seen).
- ``gene_id_of`` searches the text for ``gene_id "value"`` (one space, non-empty value; the first match
  wins; ``ref_gene_id "x"`` also matches).

``iter_gtf_fields`` skips comments, blank lines and short rows. A row with more than nine columns (a
trailing tab) is kept by default and dropped with ``exact=True``, because the callers split on this.
"""

import re
from collections.abc import Iterable, Iterator
from pathlib import Path

from viralscan.utils import open_maybe_gzip

_PAIR = re.compile(r'(\w+)\s+"([^"]*)"')
_GENE_ID = re.compile(r'gene_id "([^"]+)"')


def iter_gtf_fields(lines: Iterable[str], *, exact: bool = False) -> Iterator[list[str]]:
    """The tab-split columns of each data line (not ``#``, at least nine columns; exactly nine if ``exact``)."""
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 if exact else len(fields) < 9:
            continue
        yield fields


def parse_attributes(text: str) -> dict[str, str]:
    """``{key: value}`` of the quoted ``key "value"`` pairs in a GTF attribute column."""
    return dict(_PAIR.findall(text))


def gene_id_of(text: str) -> str | None:
    """The first ``gene_id "value"`` in a GTF attribute column, or ``None``."""
    match = _GENE_ID.search(text)
    return match.group(1) if match else None


def gene_ids(lines: Iterable[str]) -> set[str]:
    """Every ``gene_id`` (``gene_id_of`` rule) of the data lines of a GTF."""
    found = (gene_id_of(fields[8]) for fields in iter_gtf_fields(lines))
    return {gene_id for gene_id in found if gene_id}


def gtf_gene_ids(path: str | Path) -> set[str]:
    """Every ``gene_id`` of a GTF file (plain or ``.gz``)."""
    with open_maybe_gzip(path, encoding="utf-8", errors="replace") as handle:
        return gene_ids(handle)
