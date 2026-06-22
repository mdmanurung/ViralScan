"""Anellovirus reference accession table and name-map helpers.

Public API
----------
load_accession_table(path=None) -> list[Row]
    Load the packaged anellovirus_accessions.tsv.  Each Row is a dict with
    keys ``accession``, ``virus_name``, ``genus``, ``family``, ``source``.

anello_name_map() -> dict[str, str]
    Build an accession → genus-level display name map from the packaged TSV.
    Each accession (bare and versioned) maps to its Anelloviridae genus name
    (e.g. ``"Alphatorquevirus"``, ``"Betatorquevirus"``).  Unclassified entries
    fall back to ``"Anelloviridae"``.  The map is compatible with
    :func:`viralscan.virus_grouping.virus_name_for_gene`: exact-match handles
    bare accessions; the existing boundary-aware prefix rule handles
    ``{accession}_geneN`` variants.

merged_name_map() -> dict[str, str]
    Return ``{**VIRUS_NAME_MAP, **anello_name_map()}``.  Use this when calling
    :func:`viralscan.virus_grouping.virus_name_for_gene` or
    :func:`viralscan.virus_grouping.group_genes_by_virus` so that both legacy
    TTV-prefix gene IDs *and* new accession-keyed gene IDs resolve correctly.

The packaged TSV lives at ``src/viralscan/data/anellovirus_accessions.tsv`` and
is built from:
  * ViralScan's 20 bundled RefSeq anellovirus GTFs (``viralscan-refseq`` source).
  * ~2,022 CD-HIT representative genomes from the Clareau lab's
    ``clareaulab/anellovirus_reference`` (GitHub, 2025) derived from NCBI Virus
    complete human Anelloviridae sequences (``clareaulab`` source).

Sources:
  Lareau et al., clareaulab/anellovirus_reference,
  https://github.com/clareaulab/anellovirus_reference
"""

from __future__ import annotations

import csv
import importlib.resources
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from viralscan.constants import VIRUS_NAME_MAP

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_TSV_FILENAME = "anellovirus_accessions.tsv"
_PACKAGE_DATA = "viralscan.data"

# Genus values that indicate "no genus known" — fall back to family
_UNCLASSIFIED = {"unclassified", ""}


def _default_tsv_path() -> Path:
    """Return the path to the packaged TSV via importlib.resources."""
    try:
        # Python ≥ 3.9: importlib.resources.files() is the stable API
        ref = importlib.resources.files(_PACKAGE_DATA).joinpath(_TSV_FILENAME)
        return Path(str(ref))
    except AttributeError:
        # Python 3.8 fallback
        import importlib.resources as pkg_resources

        with pkg_resources.path(_PACKAGE_DATA, _TSV_FILENAME) as p:
            return Path(p)


def _strip_version(accession: str) -> str:
    """Return accession without trailing .N version suffix."""
    return re.sub(r"\.\d+$", "", accession)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

Row = dict[str, Any]


def load_accession_table(path: str | Path | None = None) -> list[Row]:
    """Load the anellovirus accession table.

    Parameters
    ----------
    path:
        Explicit path to a TSV file.  Defaults to the packaged
        ``anellovirus_accessions.tsv``.

    Returns
    -------
    List of dicts with keys ``accession``, ``virus_name``, ``genus``,
    ``family``, ``source``.
    """
    tsv_path = Path(path) if path is not None else _default_tsv_path()
    rows: list[Row] = []
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


@lru_cache(maxsize=1)
def anello_name_map() -> dict[str, str]:
    """Return an accession → genus display name map for Anelloviridae.

    Both the versioned accession (e.g. ``"PQ438050.1"``) and the bare
    accession without version (``"PQ438050"``) are added as keys so that
    :func:`~viralscan.virus_grouping.virus_name_for_gene` matches regardless
    of whether the gene ID carries a version suffix.

    The boundary-aware prefix rule in ``virus_grouping`` also handles
    ``"{acc}_geneN"`` variants without any extra entries here, because
    ``virus_name_for_gene("PQ438050.1_gene1")`` first checks whether
    ``"PQ438050.1_gene1"`` starts with ``"PQ438050.1"`` and the next
    character is ``_`` → it does.
    """
    rows = load_accession_table()
    name_map: dict[str, str] = {}
    for row in rows:
        acc = row["accession"].strip()
        genus = row["genus"].strip()
        display = genus if genus.lower() not in _UNCLASSIFIED else "Anelloviridae"
        # Add versioned accession (e.g. "PQ438050.1")
        name_map[acc] = display
        # Add bare accession without version (e.g. "PQ438050")
        bare = _strip_version(acc)
        if bare != acc:
            name_map.setdefault(bare, display)
    return name_map


def merged_name_map() -> dict[str, str]:
    """Return VIRUS_NAME_MAP merged with anello_name_map().

    Accession-level keys (from ``anello_name_map``) take priority over the
    existing prefix-based keys so that accession-keyed gene IDs (emitted by
    ``build_reference._genome_as_transcript_gtf``) are resolved correctly.
    Legacy TTV-prefix gene IDs remain functional via the inherited ``"TTV"``
    entry in :data:`~viralscan.constants.VIRUS_NAME_MAP`.

    Returns
    -------
    A new dict (does not mutate ``VIRUS_NAME_MAP``).
    """
    merged: dict[str, str] = dict(VIRUS_NAME_MAP)
    merged.update(anello_name_map())
    return merged
