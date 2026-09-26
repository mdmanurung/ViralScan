"""Anellovirus reference accession table and name-map helpers.

Public API
----------
load_accession_table(path=None) -> list[Row]
    Load the packaged anellovirus_accessions.tsv.  Each Row is a dict with
    keys ``accession``, ``virus_name``, ``genus``, ``family``, ``source``.

load_gene_table(path=None) -> list[GeneRow]
    Load the packaged anellovirus_genes.tsv — one row per real NCBI CDS
    feature, generated from cached GenBank flatfiles by
    ``extras/build_anellovirus_genes.py``.

gtf_text_for(accessions) -> str
    Emit a merged viral GTF carrying **real gene structure** for every
    accession covered by the packaged gene table, and a whole-genome
    single-exon placeholder for any accession that is not.

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

Why the gene table exists
-------------------------
Before it, every one of the panel's 2,042 genomes was reduced to a single
placeholder gene ``{accession}_gene1`` spanning the whole genome.  A
counting-matrix column built that way is a per-genome *competition* bucket, not
a measurement: because a whole-genome transcript shares sequence with every
other genome, reads cross-map in proportion to conservation, and the single
most-conserved genome absorbs the panel's whole viral signal.  Measured on a
COVID scRNA-seq run, 1,167,103 of 1,169,272 anellovirus UMI (99.8 %) landed in
``MW455439.1_gene1`` and were reported as "Alphatorquevirus".  The packaged
gene table replaces the placeholder with the CDS features NCBI actually
annotates, so a hit is attributable to a locus rather than to a conservation
rank.

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
from collections.abc import Iterable
from functools import lru_cache
from pathlib import Path
from typing import Any

from viralscan.constants import VIRUS_NAME_MAP

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_TSV_FILENAME = "anellovirus_accessions.tsv"
_GENE_TSV_FILENAME = "anellovirus_genes.tsv"
_PACKAGE_DATA = "viralscan.data"

# Genus values that indicate "no genus known" — fall back to family
_UNCLASSIFIED = {"unclassified", ""}

#: Columns of ``anellovirus_genes.tsv`` that hold integers.
_GENE_INT_COLUMNS = ("n_exons", "genome_length")

#: Columns of ``anellovirus_genes.tsv`` that hold booleans rendered as text.
_GENE_BOOL_COLUMNS = ("origin_spanning",)


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


def _default_gene_tsv_path() -> Path:
    """Return the path to the packaged gene catalogue via importlib.resources."""
    try:
        ref = importlib.resources.files(_PACKAGE_DATA).joinpath(_GENE_TSV_FILENAME)
        return Path(str(ref))
    except AttributeError:
        import importlib.resources as pkg_resources

        with pkg_resources.path(_PACKAGE_DATA, _GENE_TSV_FILENAME) as p:
            return Path(p)


def _strip_version(accession: str) -> str:
    """Return accession without trailing .N version suffix."""
    return re.sub(r"\.\d+$", "", accession)


def _coerce_gene_row(row: dict[str, str]) -> GeneRow:
    coerced: dict[str, Any] = dict(row)
    for column in _GENE_INT_COLUMNS:
        value = str(coerced.get(column, "")).strip()
        coerced[column] = int(value) if value else 0
    for column in _GENE_BOOL_COLUMNS:
        coerced[column] = str(coerced.get(column, "")).strip().lower() == "true"
    return coerced  # type: ignore[return-value]


def _exon_blocks(exons: str) -> list[tuple[int, int]]:
    """Parse the catalogue's ``start:end,start:end`` exon field."""
    blocks: list[tuple[int, int]] = []
    for piece in exons.split(","):
        piece = piece.strip()
        if not piece:
            continue
        start, _, end = piece.partition(":")
        blocks.append((int(start), int(end)))
    return blocks


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

Row = dict[str, Any]
GeneRow = dict[str, Any]


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


def load_gene_table(path: str | Path | None = None) -> list[GeneRow]:
    """Load the packaged per-genome gene catalogue.

    One row per real NCBI CDS feature, with keys ``accession``, ``gene_id``,
    ``transcript_id``, ``gene``, ``product``, ``n_exons``, ``exons``,
    ``strand``, ``genome_length``, ``topology``, ``origin_spanning``,
    ``source_genotype``, ``source_isolate``, ``source_strain``.  ``n_exons`` and
    ``genome_length`` are coerced to ``int`` and ``origin_spanning`` to
    ``bool``.

    Returns an empty list when the catalogue is not packaged, so a partial
    install degrades to whole-genome placeholders rather than failing.
    """
    tsv_path = Path(path) if path is not None else _default_gene_tsv_path()
    if not tsv_path.is_file():
        return []
    rows: list[GeneRow] = []
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(_coerce_gene_row(dict(row)))
    return rows


def annotated_accessions(path: str | Path | None = None) -> set[str]:
    """Accessions the packaged gene catalogue resolves to at least one real gene."""
    return {str(row["accession"]).strip() for row in load_gene_table(path)}


def gtf_text_for(
    accessions: Iterable[str],
    fasta_texts: dict[str, str] | None = None,
    gene_path: str | Path | None = None,
) -> str:
    """Return a merged viral GTF with real gene structure where NCBI has it.

    For every accession present in the packaged gene catalogue the real exon
    structure is emitted, one ``exon`` row per CDS interval in transcript
    order.  Accessions the catalogue does not cover fall back to the whole-genome
    single-exon placeholder that :func:`viralscan.scripts.build_reference.
    _genome_as_transcript_gtf` produces, so the merged GTF stays complete and
    ``kb ref`` never sees an unannotated sequence.

    Parameters
    ----------
    accessions:
        Accessions to annotate, in the order they should appear in the GTF.
    fasta_texts:
        Optional ``{accession: fasta_text}`` used only to size the placeholder
        exon for an uncovered accession.  Missing entries get a 1 bp exon, which
        ``kb ref`` tolerates but which keeps an unannotated genome detectable.
    gene_path:
        Override the gene-catalogue TSV (tests).
    """
    catalogue = load_gene_table(gene_path)
    by_accession: dict[str, list[GeneRow]] = {}
    for row in catalogue:
        by_accession.setdefault(str(row["accession"]).strip(), []).append(row)

    lines: list[str] = []
    covered: set[str] = set()
    for accession in accessions:
        accession = accession.strip()
        if not accession:
            continue
        rows = by_accession.get(accession)
        if not rows:
            continue
        covered.add(accession)
        for row in rows:
            gene_id = str(row["gene_id"])
            transcript_id = str(row["transcript_id"])
            product = str(row["product"]).replace('"', "'")
            symbol = str(row["gene_symbol"]).replace('"', "'")
            label = symbol or product or str(row["gene"])
            attrs = (
                f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; '
                f'gene_name "{label}"; gene_biotype "protein_coding"; '
                f'product "{product}"; n_exons "{int(row["n_exons"])}";'
            )
            if symbol:
                attrs += f' gene "{symbol}";'
            for start, end in _exon_blocks(str(row["exons"])):
                lines.append(
                    f"{accession}\tNCBI\texon\t{start}\t{end}\t.\t{str(row['strand'])}\t0\t{attrs}"
                )

    from viralscan.scripts.build_reference import _genome_as_transcript_gtf

    for accession in accessions:
        accession = accession.strip()
        if not accession or accession in covered:
            continue
        text = (fasta_texts or {}).get(accession, "")
        lines.extend(_genome_as_transcript_gtf(text, accession).splitlines())
    return "\n".join(lines) + ("\n" if lines else "")


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


def candidate_gene_ids(segments: int = 1) -> set[str]:
    """Gene IDs the expanded anellovirus panel can contribute to a reference.

    Returns three sets' worth of IDs, in one call:

    1. Every real gene in the packaged gene catalogue
       (:func:`load_gene_table`) — the names a reference built with
       :func:`gtf_text_for` actually emits.
    2. ``{accession}_geneN`` for N in ``1..segments``, the whole-genome
       placeholder convention, for references built from older GTFs.
    3. ``{accession}_txN`` companions, which the placeholder convention also
       writes.

    Over-inclusion is harmless and over-counting is the safe direction:
    ``detect_genes`` only reports IDs that are actually columns of the count
    matrix, so a name that does not exist is simply never seen. Under-inclusion
    is the bug this exists to prevent — 2,022 of 2,042 anellovirus genomes,
    measured as 91 % of the panel, were invisible to detection when the panel
    contributed only ``_gene1`` names.

    Parameters
    ----------
    segments:
        How many ``_geneN`` indices to emit per accession. Anelloviridae genomes
        are overwhelmingly single-segment, so 1 is correct for almost all of
        them; raise it to cover multi-segment records.
    """
    if segments < 1:
        raise ValueError(f"segments must be >= 1, got {segments!r}")
    gene_ids: set[str] = {str(row["gene_id"]) for row in load_gene_table()}
    for row in load_accession_table():
        accession = row["accession"].strip()
        if not accession:
            continue
        for n in range(1, segments + 1):
            gene_ids.add(f"{accession}_gene{n}")
            gene_ids.add(f"{accession}_tx{n}")
    return gene_ids


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
