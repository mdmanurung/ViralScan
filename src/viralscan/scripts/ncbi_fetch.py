"""
Download a FASTA + GTF reference from NCBI given one or more nucleotide
accession numbers (e.g. RefSeq IDs like ``NC_002021.3``) and return paths
suitable as input to ``kb ref``.

The implementation deliberately avoids heavy third-party deps (no Biopython):
it uses NCBI E-utilities ``efetch`` over plain HTTP via ``requests`` and
includes a minimal GenBank-flatfile parser that extracts CDS features and
emits a GTF annotation. This is sufficient for ``kb ref`` / ``gffread``,
which only needs ``gene_id`` / ``transcript_id`` attributes on exon records.

Downloads are cached under ``~/.cache/viralscan/ncbi/<accession>/`` so that
re-running the workflow does not re-hit NCBI.
"""

from __future__ import annotations

import hashlib
import os
import re
import time
from pathlib import Path
from typing import Iterable

import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ACCESSION_RE = re.compile(r"^[A-Za-z]{1,3}_?\d+(\.\d+)?$")
DEFAULT_CACHE_DIR = Path(os.environ.get("VIRALSCAN_CACHE", Path.home() / ".cache" / "viralscan")) / "ncbi"


class NCBIFetchError(RuntimeError):
    """Raised when an NCBI download or parse fails."""


def _validate_accession(accession: str) -> str:
    acc = accession.strip()
    if not ACCESSION_RE.match(acc):
        raise NCBIFetchError(
            f"Invalid NCBI accession {acc!r}. Expected e.g. 'NC_002021.3' or 'KX020937.1'."
        )
    return acc


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
        try:
            resp = requests.get(EUTILS_BASE, params=params, timeout=60)
        except requests.RequestException as exc:
            last_err = exc
        else:
            if resp.status_code == 200 and resp.text.strip():
                return resp.text
            if resp.status_code in (429, 500, 502, 503, 504):
                last_err = NCBIFetchError(
                    f"NCBI efetch returned {resp.status_code} for {accession} ({rettype})"
                )
            else:
                raise NCBIFetchError(
                    f"NCBI efetch failed for {accession} ({rettype}): "
                    f"HTTP {resp.status_code} — {resp.text[:200]}"
                )
        time.sleep(2 ** attempt)
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
        s = s[len("complement("):-1]
    if s.startswith("join("):
        s = s[len("join("):-1]
    parts: list[tuple[int, int, str]] = []
    for piece in s.split(","):
        piece = piece.strip().lstrip("<").lstrip(">")
        m = re.match(r"^[<>]?(\d+)\.\.[<>]?(\d+)$", piece)
        if not m:
            m2 = re.match(r"^[<>]?(\d+)$", piece)
            if not m2:
                continue
            start = end = int(m2.group(1))
        else:
            start, end = int(m.group(1)), int(m.group(2))
        parts.append((start, end, strand))
    return parts


def _genbank_to_gtf(genbank_text: str, accession: str) -> str:
    """Minimal GenBank → GTF converter.

    Extracts CDS features and emits one ``exon`` line per location interval
    with ``gene_id`` and ``transcript_id`` attributes. This is the minimum
    that ``kb ref``/``gffread`` need to extract transcript sequences.
    """
    lines = genbank_text.splitlines()
    seqid_field = accession
    in_features = False
    out: list[str] = []
    cds_count = 0

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("VERSION"):
            tokens = line.split()
            if len(tokens) >= 2:
                seqid_field = tokens[1]
        if line.startswith("FEATURES"):
            in_features = True
            i += 1
            continue
        if line.startswith("ORIGIN") or line.startswith("//"):
            in_features = False

        if in_features and line.startswith("     CDS "):
            loc = line[21:].strip()
            j = i + 1
            while j < len(lines) and lines[j].startswith(" " * 21) and not lines[j][21:22] == "/":
                loc += lines[j][21:].strip()
                j += 1

            gene_name = None
            product = None
            protein_id = None
            while j < len(lines) and lines[j].startswith(" " * 21):
                qline = lines[j].strip()
                if qline.startswith("/gene="):
                    gene_name = qline.split("=", 1)[1].strip().strip('"')
                elif qline.startswith("/product="):
                    product = qline.split("=", 1)[1].strip().strip('"')
                elif qline.startswith("/protein_id="):
                    protein_id = qline.split("=", 1)[1].strip().strip('"')
                j += 1

            cds_count += 1
            gene_id = gene_name or protein_id or f"{accession}_cds{cds_count}"
            transcript_id = protein_id or f"{gene_id}_t{cds_count}"
            label = (gene_name or product or gene_id).replace('"', "'")

            for start, end, strand in _parse_location(loc):
                attrs = (
                    f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; '
                    f'gene_name "{label}";'
                )
                out.append(
                    f"{seqid_field}\tNCBI\texon\t{start}\t{end}\t.\t{strand}\t0\t{attrs}"
                )
            i = j
            continue
        i += 1

    if not out:
        raise NCBIFetchError(
            f"No CDS features found in GenBank record for {accession}; "
            "cannot build a GTF for kb ref."
        )
    return "\n".join(out) + "\n"


def _checksum(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch_one(
    accession: str,
    cache_dir: Path,
    email: str | None,
    api_key: str | None,
) -> tuple[Path, Path]:
    """Return (fasta_path, gtf_path) for a single accession, using the cache."""
    acc = _validate_accession(accession)
    acc_dir = cache_dir / acc
    acc_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = acc_dir / f"{acc}.fasta"
    gtf_path = acc_dir / f"{acc}.gtf"

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        fasta_text = _efetch(acc, "fasta", email, api_key)
        if not fasta_text.startswith(">"):
            raise NCBIFetchError(f"Unexpected FASTA payload for {acc}: {fasta_text[:120]!r}")
        fasta_path.write_text(fasta_text)

    if not gtf_path.exists() or gtf_path.stat().st_size == 0:
        genbank_text = _efetch(acc, "gb", email, api_key)
        if "FEATURES" not in genbank_text:
            raise NCBIFetchError(f"Unexpected GenBank payload for {acc}: {genbank_text[:120]!r}")
        gtf_path.write_text(_genbank_to_gtf(genbank_text, acc))

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
        raise NCBIFetchError(
            "NCBI requires an email address. Pass --ncbi-email or set NCBI_EMAIL."
        )
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
