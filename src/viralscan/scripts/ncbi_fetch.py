"""
Download a FASTA + GTF reference from NCBI given one or more nucleotide
accession numbers (e.g. RefSeq IDs like ``NC_002021.3``) and return paths
suitable as input to ``kb ref``.

The implementation deliberately avoids heavy third-party deps (no Biopython):
it uses NCBI E-utilities ``efetch`` over plain HTTP via ``requests`` and
includes a minimal GenBank-flatfile parser that extracts CDS and ``misc_RNA``
features and emits a GTF annotation. This is sufficient for ``kb ref`` /
``gffread``, which only needs ``gene_id`` / ``transcript_id`` attributes on exon
records.

Downloads are cached under ``~/.cache/viralscan/ncbi/<accession>/`` so that
re-running the workflow does not re-hit NCBI.  Each accession directory holds
``<accession>.fasta``, ``<accession>.gb`` (the raw GenBank flatfile) and
``<accession>.gtf``, each with a ``.sha256`` sidecar.  Retaining the flatfile is
what makes the annotation re-derivable offline: the gene catalogue shipped in
``src/viralscan/data/anellovirus_genes.tsv`` is generated from it, so a re-run
of :func:`fetch_genbank` after an NCBI re-annotation never re-hits the network.

Gene identifiers are **genome-scoped**
--------------------------------------
``gene_id`` is ``<accession>_<ncbi token>``, where the token is the CDS's
``/locus_tag``, then ``/gene``, then ``/protein_id``, then ``cds<N>``.  The
NCBI symbols themselves are *not* unique — across the packaged
Anelloviridae panel ``ORF1`` is the ``/product`` of 150 different genomes and
``/gene="orf1"`` of 63 more — so emitting them bare collapses thousands of
distinct genomes onto a handful of counting-matrix columns.  The bare symbol
is preserved verbatim in the ``gene_name`` attribute.

A non-coding (``misc_RNA``) feature follows the same scheme with ``/protein_id``
replaced by ``/product`` and an ``rna<N>`` ordinal, so the scheme is never
special-cased per virus.  This is what makes EBV's EBER1/EBER2 countable at
all: they are ``misc_RNA``, they are the highest-abundance latent EBV
transcripts, and no protein identifier exists for them.

Circular topology
-----------------
Anelloviridae are circular single-stranded DNA, but NCBI records annotate
features in a *linear* representation: an origin-spanning feature is written
as an explicit ``join(<high>..<end>,<low>..<low>)`` whose interval order is the
5'→3' transcript order.  The parser therefore never re-sorts intervals; it only
reverses them for minus-strand features, and tags an origin-spanning join with
``origin_spanning="true"`` so the wrap is auditable in the shipped catalogue.
"""

from __future__ import annotations

import os
import re
import time
import warnings
from collections.abc import Callable, Iterable, Iterator
from pathlib import Path

import requests

from viralscan.run_safety import sha256_file

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ACCESSION_RE = re.compile(r"^[A-Za-z]{1,3}_?\d+(\.\d+)?$")
DEFAULT_CACHE_DIR = (
    Path(os.environ.get("VIRALSCAN_CACHE", Path.home() / ".cache" / "viralscan")) / "ncbi"
)

LOCUS_LENGTH_RE = re.compile(r"(\d+)\s+(?:bp|aa)\b")
UNSAFE_ID_CHARS = re.compile(r"[^A-Za-z0-9_.-]")

GENE_NAME_QUALIFIERS = ("locus_tag", "gene", "protein_id")
NONCODING_GENE_NAME_QUALIFIERS = ("locus_tag", "gene", "product")
NONCODING_FEATURE_KEY = "misc_RNA"


class NCBIFetchError(RuntimeError):
    """Raised when an NCBI download or parse fails."""


def _validate_accession(accession: str) -> str:
    acc = accession.strip()
    if not ACCESSION_RE.match(acc):
        raise NCBIFetchError(
            f"Invalid NCBI accession {acc!r}. Expected e.g. 'NC_002021.3' or 'KX020937.1'."
        )
    return acc


#: Minimum seconds between efetch requests. NCBI's published limit is 3
#: requests/second without an API key and 10/second with one (NBK25497);
#: proactive spacing keeps multi-accession bursts under the limit instead of
#: relying on the post-failure 429 backoff alone.
_RATE_INTERVAL_S = 1.0 / 3.0
_RATE_INTERVAL_WITH_KEY_S = 1.0 / 10.0
_last_request_at = 0.0


def _throttle(api_key: str | None) -> None:
    global _last_request_at
    interval = _RATE_INTERVAL_WITH_KEY_S if api_key else _RATE_INTERVAL_S
    wait = interval - (time.monotonic() - _last_request_at)
    if wait > 0:
        time.sleep(wait)
    _last_request_at = time.monotonic()


def _efetch(accession: str, rettype: str, email: str | None, api_key: str | None) -> str:
    """Call NCBI efetch and return the response body as text.

    Retries with exponential backoff on transient errors (429, 5xx, network).
    """
    params: dict[str, str] = {
        "db": "nuccore",
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
    }
    if email:
        params["email"] = email
        params["tool"] = "ViralScan"
    if api_key:
        params["api_key"] = api_key

    last_err: Exception | None = None
    for attempt in range(4):
        _throttle(api_key)
        try:
            resp = requests.get(EUTILS_BASE, params=params, timeout=60)
        except requests.RequestException as exc:
            last_err = exc
        else:
            if resp.status_code == 200 and resp.text.strip():
                return str(resp.text)
            if resp.status_code in (429, 500, 502, 503, 504):
                last_err = NCBIFetchError(
                    f"NCBI efetch returned {resp.status_code} for {accession} ({rettype})"
                )
            else:
                raise NCBIFetchError(
                    f"NCBI efetch failed for {accession} ({rettype}): "
                    f"HTTP {resp.status_code} — {resp.text[:200]}"
                )
        time.sleep(2**attempt)
    raise NCBIFetchError(f"NCBI efetch gave up on {accession} ({rettype}): {last_err}")


def _parse_location(loc: str) -> list[tuple[int, int, str]]:
    """Parse a GenBank feature location string into (start, end, strand) tuples.

    Handles simple, complement, and join forms. 1-based inclusive coordinates
    are preserved for GTF output.
    """
    s = loc.strip()
    strand = "+"
    if s.startswith("complement("):
        strand = "-"
        s = s[len("complement(") : -1]
    if s.startswith("join("):
        s = s[len("join(") : -1]
    parts: list[tuple[int, int, str]] = []
    for piece in s.split(","):
        piece = piece.strip().lstrip("<").lstrip(">")
        m = re.match(r"^[<>]?(\d+)\.\.[<>]?(\d+)$", piece)
        if not m:
            m2 = re.match(r"^[<>]?(\d+)$", piece)
            if not m2:
                # Nested join(complement(...)), order(...), remote segments and
                # ^ between-base locations land here: the feature is emitted
                # with those exons missing, so say so loudly.
                warnings.warn(
                    f"_parse_location: dropping unparseable interval {piece!r} "
                    f"from location {loc!r}; the GTF feature will lack it",
                    stacklevel=2,
                )
                continue
            start = end = int(m2.group(1))
        else:
            start, end = int(m.group(1)), int(m.group(2))
        parts.append((start, end, strand))
    return parts


def _locus_fields(genbank_text: str) -> dict[str, str | int]:
    """Return ``{"version", "genome_length", "topology", "molecule"}`` for a record.

    The LOCUS line is authoritative for length and for whether the submitter
    declared the sequence circular.  Most Anelloviridae records declare
    ``linear`` even though the genome is circular in vivo, so ``topology`` is
    recorded as a submitter declaration and is *not* used to decide how
    coordinates are interpreted.
    """
    fields: dict[str, str | int] = {
        "version": "",
        "genome_length": 0,
        "topology": "",
        "molecule": "",
    }
    for line in genbank_text.splitlines():
        if line.startswith("VERSION"):
            tokens = line.split()
            if len(tokens) >= 2:
                fields["version"] = tokens[1]
        elif line.startswith("LOCUS"):
            tokens = line.split()
            if len(tokens) > 1 and not fields["version"]:
                fields["version"] = tokens[1]
            match = LOCUS_LENGTH_RE.search(line)
            if match:
                fields["genome_length"] = int(match.group(1))
            fields["topology"] = "circular" if "circular" in line else "linear"
            for token in tokens:
                if token.endswith("DNA") or token.endswith("RNA"):
                    fields["molecule"] = token
    return fields


def iter_features(genbank_text: str) -> Iterator[tuple[str, str, dict[str, str]]]:
    """Yield ``(key, location, qualifiers)`` for every feature in a flatfile.

    Only single-line qualifier values are returned; multi-line values such as
    ``/translation`` keep their first physical line.  That is sufficient for the
    qualifiers this module consumes (``/gene``, ``/locus_tag``, ``/product``,
    ``/protein_id``, ``/isolate``, ``/strain``), all of which are short.
    """
    lines = genbank_text.splitlines()
    in_features = False
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("FEATURES"):
            in_features = True
            index += 1
            continue
        if line.startswith("ORIGIN") or line.startswith("//"):
            in_features = False
        if not in_features or not line.startswith("     ") or not line[5:6].strip():
            index += 1
            continue
        key = line[5:21].strip()
        if not key:
            index += 1
            continue
        location = line[21:].strip()
        cursor = index + 1
        while (
            cursor < len(lines)
            and lines[cursor].startswith(" " * 21)
            and lines[cursor][21:22] != "/"
        ):
            location += lines[cursor][21:].strip()
            cursor += 1
        qualifiers: dict[str, str] = {}
        while cursor < len(lines) and lines[cursor].startswith(" " * 21):
            text = lines[cursor][21:].strip()
            if text.startswith("/"):
                name, _, value = text[1:].partition("=")
                qualifiers.setdefault(name, value.strip().strip('"'))
            cursor += 1
        yield key, location, qualifiers
        index = cursor


def _source_qualifiers(genbank_text: str) -> dict[str, str]:
    """Return the ``source`` feature's qualifiers (isolate, strain, host, …)."""
    for key, _location, qualifiers in iter_features(genbank_text):
        if key == "source":
            return qualifiers
    return {}


def _transcript_order(parts: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    """Return *parts* in 5'→3' transcript order.

    GenBank writes a feature's intervals in transcript order already, so a
    plus-strand feature is left untouched — including an origin-spanning
    ``join(<high>..<end>,<low>..<low>)``, whose first interval is the
    high-coordinate one.  A minus-strand feature reads the reverse, so its
    interval list is reversed.  Sorting is never applied: a coordinate sort
    would silently reorder an origin-spanning join into a scrambled transcript.
    """
    if parts and parts[0][2] == "-":
        return list(reversed(parts))
    return list(parts)


def _origin_spans(parts: list[tuple[int, int, str]], genome_length: int) -> bool:
    """True when a multi-interval feature crosses the origin of a circular genome.

    The test is structural rather than arithmetic on the genome length: in
    transcript order the first interval must start *after* the last interval
    ends, which is only possible if the transcript ran off the end of the
    linear representation and back to its start.  An ordinary internal join is
    monotonically increasing and cannot satisfy this.

    ``genome_length`` is unused by the test but is part of the signature so that
    callers pass the record's declared length alongside the exons; it is what
    :mod:`extras.build_anellovirus_genes` records for its out-of-range audit.
    """
    del genome_length
    if len(parts) < 2:
        return False
    return parts[0][0] > parts[-1][1]


def _gene_token(qualifiers: dict[str, str], index: int) -> str:
    """Shortest NCBI identifier available for a CDS, sanitised for a gene ID.

    ``/locus_tag`` is preferred because it is the submitter's stable name for the
    locus, then ``/gene`` (the GenBank symbol, e.g. ``orf2/5``), then
    ``/protein_id``, which is globally unique across the whole of NCBI and was
    observed to be reused by zero of the 206 sampled Anelloviridae accessions.
    A CDS with none of the three falls back to its 1-based ordinal.
    """
    for name in GENE_NAME_QUALIFIERS:
        value = qualifiers.get(name, "").strip()
        if value:
            return UNSAFE_ID_CHARS.sub("_", value)
    return f"cds{index}"


def _noncoding_token(qualifiers: dict[str, str], index: int) -> str:
    """Gene-ID token for a non-coding feature, sanitised exactly as a CDS token is.

    Same precedence as :func:`_gene_token` minus ``/protein_id`` — a
    ``misc_RNA`` has no protein — plus ``/product``, which is what carries the
    transcript's name: EBV's EBER features have no ``/locus_tag`` and no
    ``/gene``, only ``/product="EBER-1 (pol III transcript)"`` (RefSeq
    ``NC_007605.1``) or ``/locus_tag="HHV4tp2_gs01"`` + ``/product="EBER-1"``
    (``NC_009334.1``).  Without ``/product`` in the precedence both EBERs would
    fall through to an opaque ordinal.

    The ordinal fallback is ``rna<N>`` rather than ``cds<N>`` so that a
    non-coding feature with no identifier at all is still distinguishable from
    a CDS with none.
    """
    for name in NONCODING_GENE_NAME_QUALIFIERS:
        value = qualifiers.get(name, "").strip()
        if value:
            return UNSAFE_ID_CHARS.sub("_", value)
    return f"rna{index}"


def _panel_gene_ids(
    cds_features: list[tuple[str, dict[str, str]]],
    accession: str,
    token_for: Callable[[dict[str, str], int], str] = _gene_token,
    seen: dict[str, int] | None = None,
) -> list[tuple[str, str, str, str, dict[str, str]]]:
    """Return ``(gene_id, transcript_id, label, location, qualifiers)`` per CDS.

    ``gene_id`` is ``<accession>_<token>``; ``token`` is the bare identifier the
    gene ID is built from, so ``gene_id`` is always decomposable; ``transcript_id``
    is the CDS's ``/protein_id`` where one exists (as before) and
    ``<gene_id>_t<N>`` otherwise.

    Duplicate tokens within one genome — RefSeq annotates the EBNA-1 locus twice,
    and some Anelloviridae records repeat a ``/gene`` across features — are
    disambiguated with a ``_dup<N>`` suffix starting at 2.  ``_dup`` rather than a
    bare ordinal because a bare ordinal is indistinguishable from a real token
    (``orf12``), whereas no NCBI symbol contains ``_dup``.

    ``token_for`` chooses which identifier the token is built from:
    :func:`_gene_token` for CDS, :func:`_noncoding_token` for non-coding features.
    ``<N>`` is the 1-based ordinal *within the group passed in*, so passing CDS
    and non-coding features as two separate groups leaves every CDS ordinal — and
    with it every ``cds<N>`` fallback and ``_t<N>`` transcript ID — unchanged.

    ``seen`` is the cross-group token counter.  Sharing it across both groups is
    what keeps a non-coding feature from minting a gene ID a CDS already owns:
    HHV-8 annotates the gp79 glycoprotein CDS and the T0.7 transcript — one
    ``misc_RNA`` — with the same ``/locus_tag="HHV8GK18_gp79"``, and HTLV-2
    annotates a CDS and a genome-scale ``misc_RNA`` both as ``HTLV2gs1``, so
    de-duplicating each group in isolation would emit two rows under one
    ``gene_id`` and merge a protein and a transcript into a single
    counting-matrix column.  Omitting it leaves a CDS-only caller exactly as
    before.
    """
    out: list[tuple[str, str, str, str, dict[str, str]]] = []
    seen = {} if seen is None else seen
    for index, (location, qualifiers) in enumerate(cds_features, 1):
        token = token_for(qualifiers, index)
        seen[token] = seen.get(token, 0) + 1
        gene_id = token if seen[token] == 1 else f"{token}_dup{seen[token]}"
        transcript_id = qualifiers.get("protein_id", "").strip() or f"{gene_id}_t{index}"
        out.append((f"{accession}_{gene_id}", transcript_id, token, location, qualifiers))
    return out


def _display_label(token: str, qualifiers: dict[str, str]) -> str:
    """Human-readable gene name: ``/gene``, then ``/locus_tag``, then ``/product``.

    The previous writer used the same precedence, so a CDS annotated with only a
    product string reports ``ORF1`` rather than a synthetic ordinal.
    """
    return (
        qualifiers.get("gene", "").strip()
        or qualifiers.get("locus_tag", "").strip()
        or qualifiers.get("product", "").strip()
        or token
    )


def _gtf_attributes(
    gene_id: str,
    transcript_id: str,
    token: str,
    qualifiers: dict[str, str],
    exons: list[tuple[int, int, str]],
    genome_length: int,
    biotype: str = "protein_coding",
) -> str:
    parts = [
        f'gene_id "{gene_id}"',
        f'transcript_id "{transcript_id}"',
        f'gene_name "{_display_label(token, qualifiers)}"',
        f'gene_biotype "{biotype}"',
    ]
    for qualifier in ("product", "locus_tag", "gene", "protein_id", "note"):
        value = qualifiers.get(qualifier, "").strip().replace('"', "'")
        if value:
            parts.append(f'{qualifier} "{value}"')
    parts.append(f'n_exons "{len(exons)}"')
    if _origin_spans(exons, genome_length):
        parts.append('origin_spanning "true"')
    # GTF spec: every attribute is terminated by a semicolon. Matches the
    # repo's other writers (build_reference.py, _whole_genome_gtf_from_fasta).
    return "; ".join(parts) + ";"


def _genbank_to_gtf(genbank_text: str, accession: str) -> str:
    """Minimal GenBank → GTF converter.

    Extracts CDS features and, because the most abundant viral transcripts are
    routinely non-coding, ``misc_RNA`` features, and emits one ``exon`` line per
    location interval with ``gene_id`` and ``transcript_id`` attributes. This is
    the minimum that ``kb ref``/``gffread`` need to extract transcript sequences.

    ``misc_RNA`` is the key EBV annotates EBER1/EBER2 under — RefSeq
    ``NC_007605.1`` 6629..6795 and 6956..7128, ``NC_009334.1`` 6634..6800 and
    6961..7133 — and those are the highest-abundance latent EBV transcripts, so
    a CDS-only reference cannot see latent infection at all.  ``ncRNA``/miRNA
    features are deliberately *not* emitted: EBV's are 18–24 nt, shorter than
    kallisto's ``k=31``, so they contain no 31-mer and would be dead targets.

    CDS rows are emitted first in flatfile order and the non-coding rows after
    them, so every CDS gene ID, transcript ID, row order and byte position is
    identical to a CDS-only build.  That is what keeps the 1,995-accession
    Anelloviridae catalogue and every already-cached ``<accession>.gtf``
    untouched; it also matches how the vendor EBV reference is assembled, with
    the two EBER records appended after the 94 CDS records.

    A single-exon non-coding feature is emitted as the one exon NCBI annotates;
    no second exon is fabricated.  Single-exon targets are not at risk from
    ``kb ref``: the vendor EBV reference already carries both EBERs as
    single-exon ``exon`` records, 2,446 of the 2,515 packaged Anelloviridae
    genes are single-exon, and the single-exon whole-genome placeholder
    ``MW455439.1_gene1`` took 1,167,103 UMI in a COVID run.

    Two consequences of admitting the whole ``misc_RNA`` class are worth
    knowing, both measured across the 2,249 cached flatfiles (6 records carry
    any ``misc_RNA`` at all; the 1,995-accession Anelloviridae panel carries
    none, so its GTFs and the packaged catalogue are unaffected):

    * ``k`` is the length floor, not the exon count.  HHV-8 ``NC_009333.1``
      annotates 10 miRNAs as ``misc_RNA`` at 20–23 nt, which contain no 31-mer,
      so ``kallisto index`` skips them; they stay in the GTF so the annotation
      records that NCBI calls them, rather than hiding a filter in the writer.
    * Two records annotate a genome-scale ``misc_RNA``: hepatitis A
      ``NC_001489.1`` 1..7478 is the whole 7,478-nt genome, and HTLV-2
      ``NC_001488.1`` 316..8751 is 94 % of its 8,952 nt.  Emitted into a
      ``viral_gene`` index these are exactly the whole-genome pseudo-transcripts
      that collapsed 99.8 % of anellovirus UMI into one bucket, so resolving
      them belongs to the two-index-tier split in
      ``docs/plans/2026-09-27-viral-reference-panel-expansion.md`` §5.1, not to
      a span threshold invented here.

    Gene IDs are genome-scoped — see the module docstring.  Raises
    :class:`NCBIFetchError` when the record carries neither a CDS nor a
    ``misc_RNA`` feature, which is the caller's signal to fall back to
    :func:`_whole_genome_gtf_from_fasta`.
    """
    locus = _locus_fields(genbank_text)
    seqid_field = str(locus["version"]) or accession
    genome_length = int(locus["genome_length"])
    cds_features = [
        (location, qualifiers)
        for key, location, qualifiers in iter_features(genbank_text)
        if key == "CDS"
    ]
    noncoding_features = [
        (location, qualifiers)
        for key, location, qualifiers in iter_features(genbank_text)
        if key == NONCODING_FEATURE_KEY
    ]
    if not cds_features and not noncoding_features:
        raise NCBIFetchError(
            f"No CDS or {NONCODING_FEATURE_KEY} features found in GenBank record for "
            f"{accession}; cannot build a GTF for kb ref."
        )

    out: list[str] = []
    seen_tokens: dict[str, int] = {}
    for group, token_for, biotype in (
        (cds_features, _gene_token, "protein_coding"),
        (noncoding_features, _noncoding_token, NONCODING_FEATURE_KEY),
    ):
        for gene_id, transcript_id, token, location, qualifiers in _panel_gene_ids(
            group, accession, token_for=token_for, seen=seen_tokens
        ):
            exons = _transcript_order(_parse_location(location))
            if not exons:
                continue
            attrs = _gtf_attributes(
                gene_id, transcript_id, token, qualifiers, exons, genome_length, biotype
            )
            for start, end, strand in exons:
                out.append(f"{seqid_field}\tNCBI\texon\t{start}\t{end}\t.\t{strand}\t0\t{attrs}")

    if not out:
        raise NCBIFetchError(
            f"CDS and {NONCODING_FEATURE_KEY} features in GenBank record for {accession} "
            "carry no parseable location; cannot build a GTF for kb ref."
        )
    return "\n".join(out) + "\n"


def _whole_genome_gtf_from_fasta(fasta_text: str, accession: str) -> str:
    """Fallback GTF for sequences with no CDS annotations.

    Treats each FASTA record as a single gene/transcript/exon spanning the full
    sequence length.  Mirrors the ``_genome_as_transcript_gtf`` convention used
    by ``build_reference.py`` for anellovirus whole-genome references.
    """
    lines: list[str] = []
    seq_idx = 0
    current_header = ""
    current_length = 0

    def _flush() -> None:
        nonlocal seq_idx
        if not current_header:
            return
        seq_idx += 1
        gene_id = f"{accession}_gene{seq_idx}"
        tx_id = f"{accession}_tx{seq_idx}"
        attrs = (
            f'gene_id "{gene_id}"; transcript_id "{tx_id}"; '
            f'gene_name "{accession}"; gene_biotype "whole_genome";'
        )
        seqname = current_header.split()[0]
        end = current_length if current_length > 0 else 1
        for feat in ("gene", "transcript", "exon"):
            lines.append(f"{seqname}\tViralScan\t{feat}\t1\t{end}\t.\t+\t.\t{attrs}")

    for raw in fasta_text.splitlines():
        row = raw.strip()
        if not row:
            continue
        if row.startswith(">"):
            _flush()
            current_header = row[1:]
            current_length = 0
        else:
            current_length += len(row)
    _flush()
    # Trailing newline matters: mode-3 multi-accession runs concatenate
    # per-accession GTF chunks with ``"".join(...)``, so a chunk without a
    # final newline glues its last line onto the next chunk's first line.
    return "\n".join(lines) + "\n"


# CAT-37: bump whenever ``_genbank_to_gtf`` / ``_whole_genome_gtf_from_fasta`` change what
# they emit, so cached GTFs from an older generator are rebuilt instead of reused.
GTF_FORMAT_VERSION = 3


def _gtf_stamp(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".fmt")


def _gtf_cache_valid(path: Path) -> bool:
    """``_cache_valid`` plus a match on the generator format-version stamp."""
    stamp = _gtf_stamp(path)
    return (
        _cache_valid(path)
        and stamp.exists()
        and stamp.read_text().strip() == str(GTF_FORMAT_VERSION)
    )


def _cache_valid(path: Path) -> bool:
    """Return True iff *path* exists, is non-empty, and its .sha256 sidecar matches.

    A missing or mismatched sidecar causes this function to return False, which
    triggers a fresh download in ``_fetch_one()``.  This guards against silent
    reuse of truncated files left by interrupted downloads.

    Audit §3.2: audits/2026-05-08-full-pipeline.md
    """
    if not path.exists() or path.stat().st_size == 0:
        return False
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if not sidecar.exists():
        return False
    return sidecar.read_text().strip() == sha256_file(path)


def _write_cached(path: Path, content: str) -> None:
    """Write *content* to *path* and create/update the companion .sha256 sidecar."""
    path.write_text(content)
    sidecar = path.with_suffix(path.suffix + ".sha256")
    sidecar.write_text(sha256_file(path))


def genbank_cache_path(accession: str, cache_dir: str | os.PathLike[str] | None = None) -> Path:
    """Return the cache path of an accession's raw GenBank flatfile."""
    acc = _validate_accession(accession)
    root = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    return root / acc / f"{acc}.gb"


def fetch_genbank(
    accession: str,
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
) -> tuple[Path, str]:
    """Return ``(cached_flatfile_path, genbank_text)`` for one accession.

    Cache-first: an existing flatfile with a matching ``.sha256`` sidecar is
    reused verbatim, so a re-run of the gene-catalogue generator costs no NCBI
    requests.  A missing or mismatched sidecar triggers a fresh ``efetch`` via
    the same retry/rate-limit plumbing as :func:`_efetch`.

    Keeping the flatfile rather than only its derived GTF is what makes the
    annotation re-derivable: ``extras/build_anellovirus_genes.py`` reads these
    files directly, and gene IDs, exon blocks and product strings can all be
    recomputed after a code change without touching NCBI.

    Raises :class:`NCBIFetchError` when the payload is not a GenBank record.
    """
    acc = _validate_accession(accession)
    path = genbank_cache_path(acc, cache_dir)
    if _cache_valid(path):
        return path, path.read_text()

    resolved_email = email or os.environ.get("NCBI_EMAIL")
    if not resolved_email:
        raise NCBIFetchError("NCBI requires an email address. Pass --ncbi-email or set NCBI_EMAIL.")
    resolved_key = api_key or os.environ.get("NCBI_API_KEY")

    genbank_text = _efetch(acc, "gb", resolved_email, resolved_key)
    if "FEATURES" not in genbank_text:
        raise NCBIFetchError(f"Unexpected GenBank payload for {acc}: {genbank_text[:120]!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_cached(path, genbank_text)
    return path, genbank_text


def catalogue_rows(accession: str, genbank_text: str) -> list[dict[str, object]]:
    """Return one row per CDS feature of *accession*, for the shipped gene TSV.

    Columns mirror ``src/viralscan/data/anellovirus_genes.tsv``:
    ``accession``, ``gene_id``, ``transcript_id``, ``gene`` (the token the gene ID
    is built from), ``gene_symbol`` (the verbatim ``/gene``), ``product``,
    ``n_exons``, ``exons`` (``start:end,start:end`` in transcript order),
    ``strand``, ``genome_length``, ``topology``, ``origin_spanning``,
    ``source_genotype``, ``source_isolate``, ``source_strain``.

    ``source_genotype`` is NCBI's own ``/genotype`` qualifier on the ``source``
    feature, copied verbatim.  It is **not** a derived or inferred genogroup.
    Measured over all 2,042 panel accessions: no record carries a ``/genogroup``
    qualifier, and exactly two carry ``/genotype`` — ``NC_014081.1`` (``"6"``,
    whose ``/organism`` is ``Torque teno virus 3``) and ``NC_014094.1``
    (``"28"``, organism ``Torque teno virus 6``).  Those two disagree with their
    own organism number, which is the concrete reason genogroup must not be
    inferred from ``/organism``: the species and genogroup namespaces are
    separate and neither is derivable from the other.  The other 2,040 ship
    empty rather than being back-filled from free text.

    Returns ``[]`` for a record with no CDS feature; the caller is expected to
    record the accession as unannotated rather than to invent a gene.

    Non-coding features are deliberately excluded here even though
    :func:`_genbank_to_gtf` emits them: this table is the packaged *gene*
    catalogue whose rows are consumed as protein-coding loci
    (``viralscan.anellovirus`` → ``gtf_text_for``), and no accession in the
    Anelloviridae panel carries a ``misc_RNA`` feature, so including them would
    change no shipped row while widening that contract.
    """
    locus = _locus_fields(genbank_text)
    genome_length = int(locus["genome_length"])
    topology = str(locus["topology"])
    source = _source_qualifiers(genbank_text)
    cds_features = [
        (location, qualifiers)
        for key, location, qualifiers in iter_features(genbank_text)
        if key == "CDS"
    ]
    rows: list[dict[str, object]] = []
    for gene_id, transcript_id, token, location, qualifiers in _panel_gene_ids(
        cds_features, accession
    ):
        exons = _transcript_order(_parse_location(location))
        if not exons:
            continue
        rows.append(
            {
                "accession": accession,
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "gene": token,
                "gene_symbol": qualifiers.get("gene", "").strip(),
                "product": qualifiers.get("product", "").strip(),
                "n_exons": len(exons),
                "exons": ",".join(f"{start}:{end}" for start, end, _ in exons),
                "strand": exons[0][2],
                "genome_length": genome_length,
                "topology": topology,
                "origin_spanning": str(_origin_spans(exons, genome_length)).lower(),
                "source_genotype": source.get("genotype", "").strip(),
                "source_isolate": source.get("isolate", "").strip(),
                "source_strain": source.get("strain", "").strip(),
            }
        )
    return rows


def _fetch_one(
    accession: str,
    cache_dir: Path,
    email: str | None,
    api_key: str | None,
) -> tuple[Path, Path]:
    """Return (fasta_path, gtf_path) for a single accession, using the cache.

    Cache hits are validated with a SHA-256 sidecar file.  If the sidecar is
    missing or the hash does not match (e.g. interrupted previous download),
    the file is re-fetched from NCBI.

    The raw GenBank flatfile is cached alongside the two derived artefacts as
    ``<accession>.gb``, so the annotation stays re-derivable offline: deleting
    only ``<accession>.gtf`` rebuilds it from the retained flatfile with no
    network access.

    Audit §3.2: audits/2026-05-08-full-pipeline.md
    """
    acc = _validate_accession(accession)
    acc_dir = cache_dir / acc
    acc_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = acc_dir / f"{acc}.fasta"
    gtf_path = acc_dir / f"{acc}.gtf"

    if not _cache_valid(fasta_path):
        fasta_text = _efetch(acc, "fasta", email, api_key)
        if not fasta_text.startswith(">"):
            raise NCBIFetchError(f"Unexpected FASTA payload for {acc}: {fasta_text[:120]!r}")
        _write_cached(fasta_path, fasta_text)

    if not _gtf_cache_valid(gtf_path):
        # fetch_genbank reuses the retained .gb when present, else refetches.
        _gb_path, genbank_text = fetch_genbank(acc, email, api_key, cache_dir)
        try:
            gtf_content = _genbank_to_gtf(genbank_text, acc)
        except NCBIFetchError:
            # No CDS annotations — fall back to a whole-genome single-exon GTF
            gtf_content = _whole_genome_gtf_from_fasta(fasta_path.read_text(), acc)
        _write_cached(gtf_path, gtf_content)
        _gtf_stamp(gtf_path).write_text(str(GTF_FORMAT_VERSION))

    return fasta_path, gtf_path


def fetch_reference(
    accessions: Iterable[str],
    out_dir: str | os.PathLike[str],
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
) -> tuple[Path, Path]:
    """Download FASTA + GTF for one or more NCBI nucleotide accessions.

    Multiple accessions are concatenated into a single merged FASTA and a
    single merged GTF, mirroring ViralScan's existing comma-separated
    ``--fasta``/``--gtf`` semantics.

    Parameters
    ----------
    accessions:
        Iterable of nucleotide accession strings (RefSeq or GenBank).
    out_dir:
        Directory where merged ``reference.fasta`` and ``reference.gtf``
        will be written.
    email:
        Contact email; required by NCBI's E-utilities terms of service.
        Falls back to the ``NCBI_EMAIL`` environment variable.
    api_key:
        Optional NCBI API key for higher rate limits. Falls back to
        ``NCBI_API_KEY`` environment variable.
    cache_dir:
        Override the per-accession cache directory.

    Returns
    -------
    (fasta_path, gtf_path) as :class:`pathlib.Path` objects, ready for
    ``kb ref``.
    """
    accessions = [a for a in accessions if a]
    if not accessions:
        raise NCBIFetchError("At least one NCBI accession must be provided.")

    email = email or os.environ.get("NCBI_EMAIL")
    if not email:
        raise NCBIFetchError("NCBI requires an email address. Pass --ncbi-email or set NCBI_EMAIL.")
    api_key = api_key or os.environ.get("NCBI_API_KEY")

    cache = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    merged_fasta = out / "reference.fasta"
    merged_gtf = out / "reference.gtf"

    fasta_chunks: list[str] = []
    gtf_chunks: list[str] = []
    for acc in accessions:
        fasta_path, gtf_path = _fetch_one(acc, cache, email, api_key)
        fasta_chunks.append(fasta_path.read_text())
        gtf_chunks.append(gtf_path.read_text())

    merged_fasta.write_text("".join(fasta_chunks))
    merged_gtf.write_text("".join(gtf_chunks))

    if merged_fasta.stat().st_size == 0 or merged_gtf.stat().st_size == 0:
        raise NCBIFetchError("Merged reference files are empty after download.")

    return merged_fasta, merged_gtf
