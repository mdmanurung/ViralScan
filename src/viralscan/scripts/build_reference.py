"""build_reference.py — build a combined host + virus kallisto reference.

Public API
----------
build_combined_reference(
    host_species, virus_accessions, out_dir,
    email=None, api_key=None, cache_dir=None, run_kb_ref=True,
) -> dict

    Downloads host cDNA FASTA + GTF from Ensembl and viral FASTA from NCBI
    (via ncbi_fetch.fetch_reference), concatenates both, and optionally runs
    ``kb ref`` to produce a kallisto index + t2g mapping.

fetch_host_cdna(species, out_dir, cache_dir=None) -> (fasta_path, gtf_path)

    Download Ensembl cDNA FASTA + GTF for a supported host species.
    Results are cached under ~/.cache/viralscan/ensembl/<species>/.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from viralscan.constants import ENSEMBL_SPECIES
from viralscan.run_safety import sha256_file
from viralscan.sensitivity import DEFAULT_K
from viralscan.validation import require_schema_valid

log = logging.getLogger("viralscan")

# ---------------------------------------------------------------------------
# Ensembl HTTPS-FTP mirror helpers
# ---------------------------------------------------------------------------

_ENSEMBL_VERSION_URL = "https://ftp.ensembl.org/pub/VERSION"
_ENSEMBL_FALLBACK_RELEASE = 116
# current_fasta / current_gtf symlinks are unreliable; use release-pinned paths.
_ENSEMBL_FTP = "https://ftp.ensembl.org/pub/release-{release}/fasta/{species}/cdna/"
_ENSEMBL_GTF = "https://ftp.ensembl.org/pub/release-{release}/gtf/{species}/"


def _ensembl_release() -> int:
    """Return the current Ensembl release number, falling back to a hardcoded value."""
    try:
        with urllib.request.urlopen(_ENSEMBL_VERSION_URL, timeout=15) as resp:  # noqa: S310
            return int(resp.read().strip())
    except Exception as exc:
        log.warning(
            "Could not fetch Ensembl release from %s (%s); defaulting to %d",
            _ENSEMBL_VERSION_URL,
            exc,
            _ENSEMBL_FALLBACK_RELEASE,
        )
        return _ENSEMBL_FALLBACK_RELEASE


def _ensembl_species_key(species: str) -> str:
    """Normalise user-supplied species name to ENSEMBL_SPECIES key."""
    key = species.strip().lower().replace(" ", "_")
    if key in ENSEMBL_SPECIES:
        return key
    # Try looking up by Ensembl name (e.g. 'homo_sapiens')
    for short, (ens, _) in ENSEMBL_SPECIES.items():
        if key == ens:
            return short
    supported = ", ".join(sorted(ENSEMBL_SPECIES))
    raise ValueError(f"Unknown host species {species!r}. Supported values: {supported}")


def _download(url: str, dest: Path, timeout: int = 120, retries: int = 3) -> Path:
    """Download *url* to *dest* with simple retry logic."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            log.info("Downloading %s", url)
            with urllib.request.urlopen(url, timeout=timeout) as resp, open(dest, "wb") as fh:  # noqa: S310
                shutil.copyfileobj(resp, fh)
            return dest
        except Exception as exc:
            if attempt < retries - 1:
                wait = 2**attempt
                log.warning("Download error (%s); retrying in %ds …", exc, wait)
                time.sleep(wait)
            else:
                raise RuntimeError(f"Failed to download {url}: {exc}") from exc
    return dest  # unreachable


def _fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence: list[str] = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence).upper()))
                identifier = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if identifier is not None:
        records.append((identifier, "".join(sequence).upper()))
    return records


#: Maximum homopolymer run tolerated inside a reference k-mer.  kallisto's own
#: build-time guard clips poly-A tails longer than 10, so 11 is deliberately one
#: step looser than the tool and one step tighter than a typical poly-A tail.
LOW_COMPLEXITY_MAX_RUN = 11

#: Minimum number of distinct bases in a reference k-mer.
LOW_COMPLEXITY_MIN_BASES = 3

#: A perfect tandem repeat of a unit this short or shorter is low-complexity.
LOW_COMPLEXITY_MAX_TANDEM = 5


def low_complexity_kmer_counts(
    sequence: str,
    k: int = DEFAULT_K,
    max_run: int = LOW_COMPLEXITY_MAX_RUN,
    min_bases: int = LOW_COMPLEXITY_MIN_BASES,
    max_tandem: int = LOW_COMPLEXITY_MAX_TANDEM,
) -> dict[str, int]:
    """Break the k-mers of *sequence* down by why they are low-complexity.

    ``pure_homopolymer`` is the class that actually causes harm, so it is counted
    separately.  Measured on the deployed panel: 170 pure 31-mers across 9
    records, and the k-mers that a 10x poly-A tail matched were literally
    ``A``*31 and its near neighbours.  A panel can carry only a handful of them
    and still manufacture >1 % of R2 reads, because poly-A reads are not rare —
    which is why the gate is an absolute count and not a fraction.
    """
    counts = {"pure_homopolymer": 0, "long_run": 0, "few_bases": 0, "tandem": 0, "total": 0}
    seq = sequence.upper()
    for start in range(len(seq) - k + 1):
        window = seq[start : start + k]
        if set(window) - set("ACGT"):
            continue
        counts["total"] += 1
        distinct = len(set(window))
        if distinct == 1:
            counts["pure_homopolymer"] += 1
            continue
        longest = run = 1
        for i in range(1, k):
            run = run + 1 if window[i] == window[i - 1] else 1
            longest = max(longest, run)
        if longest > max_run:
            counts["long_run"] += 1
            continue
        if distinct < min_bases:
            counts["few_bases"] += 1
            continue
        if any(
            k % unit == 0 and window == window[:unit] * (k // unit)
            for unit in range(1, max_tandem + 1)
        ):
            counts["tandem"] += 1
    return counts


def low_complexity_kmer_fraction(
    sequence: str,
    k: int = DEFAULT_K,
    max_run: int = LOW_COMPLEXITY_MAX_RUN,
    min_bases: int = LOW_COMPLEXITY_MIN_BASES,
    max_tandem: int = LOW_COMPLEXITY_MAX_TANDEM,
) -> tuple[int, int]:
    """Return ``(n_low_complexity, n_total)`` k-mers in *sequence*.

    A k-mer is low-complexity when it has a homopolymer run longer than *max_run*,
    fewer than *min_bases* distinct bases, or is a perfect tandem repeat of a unit
    of at most *max_tandem* bases.  Windows containing a non-ACGT character are
    skipped, because kallisto replaces them and they cannot reach the index as
    written.

    This is deliberately a *k-mer space* property, not an N-masking property.  A
    panel can be almost entirely unmasked and still contribute homopolymer
    k-mers, and those k-mers are what match the long poly-A/poly-T tails that
    dominate 10x R2 reads.
    """
    counts = low_complexity_kmer_counts(sequence, k, max_run, min_bases, max_tandem)
    return (
        counts["pure_homopolymer"] + counts["long_run"] + counts["few_bases"] + counts["tandem"],
        counts["total"],
    )


def low_complexity_report(
    fasta: Path,
    k: int = DEFAULT_K,
    max_run: int = LOW_COMPLEXITY_MAX_RUN,
) -> dict[str, tuple[int, int]]:
    """Map each record identifier to its ``(low_complexity_kmers, total_kmers)``."""
    return {
        identifier: low_complexity_kmer_fraction(sequence, k=k, max_run=max_run)
        for identifier, sequence in _fasta_records(fasta)
    }


def validate_reference_records(
    fasta: Path,
    *,
    max_low_complexity_fraction: float | None = None,
    max_pure_homopolymer_kmers: int | None = None,
    k: int = DEFAULT_K,
    max_run: int = LOW_COMPLEXITY_MAX_RUN,
) -> list[tuple[str, str]]:
    """Fail before indexing on empty, duplicate-ID, or duplicate-sequence records.

    Two optional low-complexity gates, both off by default so existing callers are
    unaffected:

    *max_low_complexity_fraction*
        Per-record ceiling on the fraction of k-mers that are low-complexity.

    *max_pure_homopolymer_kmers*
        Per-record ceiling on k-mers that are a single base repeated. This is the
        class a 10x poly-A or poly-T tail matches exactly, so it is gated as an
        absolute count rather than a fraction: 170 such k-mers across 9 records
        are enough to manufacture >1 % of R2 reads, because the reads are not
        rare even though the k-mers are.

    Pass ``0`` to require a fully masked panel.
    """
    records = _fasta_records(fasta)
    if not records:
        raise ValueError(f"Reference FASTA has no sequences: {fasta}")
    seen_ids: set[str] = set()
    seen_sequences: dict[str, str] = {}
    for identifier, sequence in records:
        if not sequence:
            raise ValueError(f"Reference sequence {identifier!r} is empty.")
        if identifier in seen_ids:
            raise ValueError(f"Duplicate FASTA identifier before index construction: {identifier}")
        seen_ids.add(identifier)
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest in seen_sequences:
            raise ValueError(
                f"Exact duplicate sequences before index construction: "
                f"{seen_sequences[digest]} and {identifier}"
            )
        seen_sequences[digest] = identifier
    if max_pure_homopolymer_kmers is not None or max_low_complexity_fraction is not None:
        homopolymer_offenders: list[tuple[str, int]] = []
        fraction_offenders: list[tuple[str, int, int]] = []
        for identifier, sequence in records:
            counts = low_complexity_kmer_counts(sequence, k=k, max_run=max_run)
            if (
                max_pure_homopolymer_kmers is not None
                and counts["pure_homopolymer"] > max_pure_homopolymer_kmers
            ):
                homopolymer_offenders.append((identifier, counts["pure_homopolymer"]))
            if max_low_complexity_fraction is not None and counts["total"]:
                low = (
                    counts["pure_homopolymer"]
                    + counts["long_run"]
                    + counts["few_bases"]
                    + counts["tandem"]
                )
                if low / counts["total"] > max_low_complexity_fraction:
                    fraction_offenders.append((identifier, low, counts["total"]))
        problems = []
        if homopolymer_offenders:
            homopolymer_offenders.sort(key=lambda row: -row[1])
            shown = ", ".join(f"{i} ({n})" for i, n in homopolymer_offenders[:10])
            more = (
                ""
                if len(homopolymer_offenders) <= 10
                else f" (+{len(homopolymer_offenders) - 10} more)"
            )
            problems.append(
                f"{len(homopolymer_offenders)} record(s) contain pure-homopolymer "
                f"{k}-mers above the limit of {max_pure_homopolymer_kmers}: {shown}{more}"
            )
        if fraction_offenders:
            fraction_offenders.sort(key=lambda row: -(row[1] / row[2]))
            shown = ", ".join(f"{i} ({a}/{b})" for i, a, b in fraction_offenders[:10])
            more = (
                "" if len(fraction_offenders) <= 10 else f" (+{len(fraction_offenders) - 10} more)"
            )
            problems.append(
                f"{len(fraction_offenders)} record(s) exceed the low-complexity fraction "
                f"limit of {max_low_complexity_fraction}: {shown}{more}"
            )
        if problems:
            raise ValueError(
                "Reference k-mers are low-complexity and will match the poly-A and poly-T "
                "tails that dominate 10x reads. " + "; ".join(problems) + ". Mask the panel "
                "with `dustmasker` (windows 64 and 30, merged, masked to N) before indexing, "
                "or raise the limits deliberately. A kallisto D-list cannot fix this: it "
                "filters host-homologous k-mers, not self-similarity inside a viral contig."
            )
    return records


def write_reference_manifest(
    fasta: Path,
    output: Path,
    *,
    profile: str,
    host_species: str,
    viral_identifiers: set[str],
    annotations: Optional[dict[str, dict[str, object]]] = None,
    genome_dlist: Optional[Path] = None,
) -> Path:
    """Write machine-readable per-sequence provenance for a frozen reference."""
    created_at = datetime.now(timezone.utc).isoformat()
    sequences = []
    for identifier, sequence in validate_reference_records(fasta):
        is_viral = identifier in viral_identifiers
        counts = [sequence.count(base) for base in "ACGT"]
        fractions = [count / len(sequence) for count in counts if count]
        if not fractions:
            fractions = [1.0]
        entropy = -sum(fraction * math.log2(fraction) for fraction in fractions)
        low_kmers, total_kmers = low_complexity_kmer_fraction(sequence)
        record: dict[str, object] = {
            "accession_version": identifier,
            "taxonomy": "virus" if is_viral else host_species,
            "source": "NCBI nucleotide" if is_viral else "Ensembl cDNA",
            "source_snapshot": "retrieved build input",
            "retrieved_at": created_at,
            "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "length": len(sequence),
            "source_licence": "source database terms apply",
            "cluster": None,
            "representative_status": "input" if is_viral else "host_transcript",
            "inclusion_rationale": (
                "requested viral accession" if is_viral else "competitive host transcriptome"
            ),
            "low_complexity_max_base_fraction": max(fractions),
            "low_complexity_entropy": entropy,
            "low_complexity_flag": max(fractions) >= 0.80 or entropy < 1.20,
            "low_complexity_kmers": low_kmers,
            "low_complexity_kmer_fraction": (
                round(low_kmers / total_kmers, 6) if total_kmers else 0.0
            ),
        }
        if annotations and identifier in annotations:
            record.update(annotations[identifier])
        elif is_viral:
            record["host_homology_status"] = "not_assessed_no_host_genome"
        sequences.append(record)
    manifest = {
        "schema_version": "3.0.0",
        "profile": profile,
        "created_at": created_at,
        "host_species": host_species,
        "fasta_sha256": sha256_file(fasta),
        "genome_dlist": (
            {
                "path": str(genome_dlist.resolve()),
                "sha256": sha256_file(genome_dlist),
                "purpose": "mask host-genomic k-mers shared with viral sequences",
            }
            if genome_dlist
            else None
        ),
        "sequences": sequences,
    }
    require_schema_valid(manifest, "reference_manifest.schema.json", output)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return output


def _parse_host_homology_paf(
    paf_text: str, viral_lengths: dict[str, int]
) -> dict[str, dict[str, object]]:
    """Reduce raw minimap2 PAF alignments to maximum per-query host homology."""
    annotations = {
        identifier: {
            "host_homology_status": "measured",
            "host_homology_max_identity": 0.0,
            "host_homology_max_query_coverage": 0.0,
            "host_homology_max_aligned_bases": 0,
            "host_homology_best_target": "",
        }
        for identifier in viral_lengths
    }
    for line in paf_text.splitlines():
        fields = line.split("\t")
        if len(fields) < 12 or fields[0] not in annotations:
            continue
        query, query_length, query_start, query_end = (
            fields[0],
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
        )
        matches, block_length = int(fields[9]), int(fields[10])
        identity = matches / block_length if block_length else 0.0
        query_coverage = (query_end - query_start) / query_length if query_length else 0.0
        current = annotations[query]
        if block_length > int(current["host_homology_max_aligned_bases"]):
            current.update(
                host_homology_max_identity=identity,
                host_homology_max_query_coverage=query_coverage,
                host_homology_max_aligned_bases=block_length,
                host_homology_best_target=fields[5],
            )
    return annotations


def measure_host_homology(
    viral_fasta: Path, host_genome: Path, output_tsv: Path
) -> dict[str, dict[str, object]]:
    """Measure viral-sequence homology to the full host genome and retain raw metrics."""
    minimap2 = shutil.which("minimap2")
    if minimap2 is None:
        raise RuntimeError(
            "--genome-dlist requires minimap2 to annotate host-genome homology. "
            "Install the full ViralScan environment."
        )
    viral_lengths = {
        identifier: len(sequence) for identifier, sequence in _fasta_records(viral_fasta)
    }
    proc = subprocess.run(  # noqa: S603
        [minimap2, "-x", "asm10", str(host_genome), str(viral_fasta)],
        check=True,
        capture_output=True,
        text=True,
    )
    annotations = _parse_host_homology_paf(proc.stdout, viral_lengths)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession_version", *next(iter(annotations.values()), {}).keys()]
    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for identifier, values in sorted(annotations.items()):
            writer.writerow({"accession_version": identifier, **values})
    return annotations


def _list_ensembl_files(species_name: str, url_base: str, retries: int = 3) -> list[str]:
    """Scrape the Ensembl HTTP index page and return file-name links."""
    import html.parser

    class _Parser(html.parser.HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.links: list[str] = []

        def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
            if tag == "a":
                for k, v in attrs:
                    if k == "href" and v and not v.startswith("?") and not v.startswith("/"):
                        self.links.append(v)

    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url_base, timeout=30) as resp:  # noqa: S310
                html_bytes = resp.read()
            break
        except Exception as exc:
            if attempt < retries - 1:
                wait = 2**attempt
                log.warning(
                    "Could not list Ensembl directory %s (%s); retrying in %ds …",
                    url_base,
                    exc,
                    wait,
                )
                time.sleep(wait)
            else:
                raise RuntimeError(f"Could not list Ensembl directory {url_base}: {exc}") from exc

    parser = _Parser()
    parser.feed(html_bytes.decode("utf-8", errors="replace"))
    return parser.links


def fetch_host_cdna(
    species: str,
    out_dir: os.PathLike[str] | str,
    cache_dir: Optional[os.PathLike[str] | str] = None,
) -> tuple[Path, Path]:
    """Download Ensembl cDNA FASTA (gzipped) and GTF for *species*.

    Parameters
    ----------
    species:
        Short species name, e.g. ``"human"``, ``"mouse"``.  Run
        ``viralscan build-ref --list-species`` to see all supported names.
    out_dir:
        Directory where downloaded files will be *copied* (symlinked from cache).
    cache_dir:
        Root of the download cache.  Defaults to ``~/.cache/viralscan/ensembl``.

    Returns
    -------
    (fasta_path, gtf_path): paths to the local gzipped FASTA and GTF.
    """
    key = _ensembl_species_key(species)
    ens_name, assembly = ENSEMBL_SPECIES[key]

    if cache_dir is None:
        cache_dir = Path.home() / ".cache" / "viralscan" / "ensembl" / key
    else:
        cache_dir = Path(cache_dir) / key
    cache_dir.mkdir(parents=True, exist_ok=True)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── cDNA FASTA ──────────────────────────────────────────────────────────
    release = _ensembl_release()
    log.info("Using Ensembl release %d", release)
    cdna_base = _ENSEMBL_FTP.format(release=release, species=ens_name)
    cdna_links = _list_ensembl_files(ens_name, cdna_base)
    cdna_files = [f for f in cdna_links if re.search(r"\.cdna\.all\.fa\.gz$", f)]
    if not cdna_files:
        raise RuntimeError(
            f"Could not find a cdna.all.fa.gz file at {cdna_base}. "
            "Ensembl may have reorganised their FTP layout."
        )
    cdna_filename = cdna_files[0]
    cdna_cache = cache_dir / cdna_filename
    if not cdna_cache.exists():
        _download(cdna_base + cdna_filename, cdna_cache)
    else:
        log.info("Using cached cDNA FASTA: %s", cdna_cache)

    cdna_out = out_dir / cdna_filename
    if not cdna_out.exists():
        shutil.copy2(cdna_cache, cdna_out)

    # ── GTF ─────────────────────────────────────────────────────────────────
    gtf_base = _ENSEMBL_GTF.format(release=release, species=ens_name)
    gtf_links = _list_ensembl_files(ens_name, gtf_base)
    # We want the toplevel (not abinitio, not chr patch_hapl_scaff, not README)
    gtf_files = [
        f
        for f in gtf_links
        if re.search(r"\.\d+\.gtf\.gz$", f) and "abinitio" not in f and "chr_patch" not in f
    ]
    if not gtf_files:
        raise RuntimeError(f"Could not find a release-numbered .gtf.gz at {gtf_base}.")
    gtf_filename = gtf_files[0]
    gtf_cache = cache_dir / gtf_filename
    if not gtf_cache.exists():
        _download(gtf_base + gtf_filename, gtf_cache)
    else:
        log.info("Using cached GTF: %s", gtf_cache)

    gtf_out = out_dir / gtf_filename
    if not gtf_out.exists():
        shutil.copy2(gtf_cache, gtf_out)

    return cdna_out, gtf_out


# ---------------------------------------------------------------------------
# Viral GTF helper (port of extras/Viral_GTF_maker.py)
# ---------------------------------------------------------------------------


def _genome_as_transcript_gtf(fasta_text: str, accession: str) -> str:
    """Convert a whole-genome FASTA to a minimal GTF.

    Each sequence in *fasta_text* is represented as one gene, one transcript,
    and one exon spanning the entire sequence.  The feature biotype is
    ``whole_genome``, matching the convention in ``extras/Viral_GTF_maker.py``.

    Parameters
    ----------
    fasta_text:
        Plain-text (not gzip) FASTA content for a single viral genome.
    accession:
        NCBI accession used to construct stable gene/transcript IDs.

    Returns
    -------
    GTF lines as a single string (no trailing newline).
    """
    lines: list[str] = []
    current_header: str = ""
    current_length: int = 0
    seq_idx: int = 0

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
        for feature in ("gene", "transcript", "exon"):
            lines.append(f"{seqname}\tViralScan\t{feature}\t1\t{end}\t.\t+\t.\t{attrs}")

    for raw in fasta_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            _flush()
            current_header = line[1:]
            current_length = 0
        else:
            current_length += len(line)

    _flush()
    return "\n".join(lines)


def host_cdna_as_gtf(
    host_fasta_gz: os.PathLike[str] | str, out_path: os.PathLike[str] | str
) -> int:
    """Write a cDNA-level GTF from an Ensembl cDNA FASTA (seqname = transcript ID).

    ``kb ref`` extracts cDNA by matching each GTF seqname against a FASTA sequence
    header.  Ensembl ships a *cDNA* FASTA (headers are ENST transcript IDs) but its
    companion GTF is *chromosomal* (seqnames ``1``, ``2``, ``X`` …).  Handing that
    pair to ``kb ref`` makes it hang forever at "Splitting genome" because no
    chromosomal seqname matches a cDNA header.  Emitting one gene/transcript/exon
    per cDNA record — seqname = transcript ID, ``gene_id`` = the ``gene:ENSG…``
    field, coordinates ``1..length`` — makes the GTF consistent with the FASTA.

    Ensembl cDNA header example::

        >ENST00000632684.1 cdna chromosome:GRCh38:… gene:ENSG00000273663.1 gene_biotype:…

    Parameters
    ----------
    host_fasta_gz:
        Path to the gzip-compressed Ensembl cDNA FASTA.
    out_path:
        Destination path for the generated (plain-text) GTF.

    Returns
    -------
    Number of transcript records written.
    """
    n = 0
    current_id: Optional[str] = None
    current_gene = ""
    current_len = 0

    with gzip.open(host_fasta_gz, "rt") as fasta, open(out_path, "w") as out:

        def _flush() -> None:
            nonlocal n
            if current_id is None:
                return
            attrs = f'gene_id "{current_gene}"; transcript_id "{current_id}";'
            for feature in ("gene", "transcript", "exon"):
                out.write(
                    f"{current_id}\tEnsembl_cDNA\t{feature}\t1\t{current_len}\t.\t+\t.\t{attrs}\n"
                )
            n += 1

        for raw in fasta:
            raw = raw.rstrip()
            if raw.startswith(">"):
                _flush()
                parts = raw[1:].split()
                current_id = parts[0]
                current_gene = next(
                    (p[len("gene:") :] for p in parts if p.startswith("gene:")), parts[0]
                )
                current_len = 0
            else:
                current_len += len(raw)

        _flush()

    return n


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------


def build_combined_reference(
    host_species: str,
    virus_accessions: list[str],
    out_dir: os.PathLike[str] | str,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
    cache_dir: Optional[os.PathLike[str] | str] = None,
    run_kb_ref: bool = True,
    include_anellovirus: bool = False,
    allow_partial_panel: bool = False,
    profile: str = "curated",
    genome_dlist: Optional[os.PathLike[str] | str] = None,
) -> dict[str, Optional[Path]]:
    """Build a combined host + virus kallisto reference.

    Steps
    -----
    1. Download Ensembl cDNA FASTA + GTF for *host_species*.
    2. Download NCBI FASTA for each accession in *virus_accessions*
       (via :func:`viralscan.scripts.ncbi_fetch.fetch_reference`).
    2b. If *include_anellovirus*, fetch the 2,042 packaged anellovirus accessions
        (skip-and-log on individual failures; abort only if >50% fail).
    3. Synthesise a ``whole_genome`` GTF for each viral sequence.
    4. Concatenate host cDNA FASTA + all viral FASTAs → ``combined.fa``
       (gzip-encoded; the viral sequences are plain-text, appended after
       decompression of the host FASTA).
    5. Concatenate host GTF + viral GTFs → ``combined.gtf``.
    6. If *run_kb_ref* is ``True`` and ``kb`` is on ``$PATH``:
       ``kb ref -i index.idx -g t2g.txt -f1 cdna.fa combined.fa combined.gtf``

    Parameters
    ----------
    host_species:
        Short species name, e.g. ``"human"``.
    virus_accessions:
        List of NCBI accession numbers (e.g. ``["NC_045512.2"]``).
    out_dir:
        Destination directory for all output files.
    email:
        E-mail address for NCBI E-utilities (recommended, avoids throttling).
    api_key:
        NCBI API key for higher request rate.
    cache_dir:
        Cache root; defaults to ``~/.cache/viralscan``.
    run_kb_ref:
        Whether to run ``kb ref`` after concatenating files.
    include_anellovirus:
        When ``True``, union the full packaged anellovirus accession
        table into the reference.  Accessions already in *virus_accessions* are
        de-duplicated so they are not fetched twice.  Use ``--no-anellovirus``
        (via :func:`build_ref_main`) to skip. The v3 default is ``False``.

    Returns
    -------
    dict with keys: ``fasta``, ``gtf``, ``index`` (None if *run_kb_ref* is
    False), ``t2g`` (None if *run_kb_ref* is False).
    """
    # Late import to avoid circular dependency at module load time.
    from viralscan.scripts.ncbi_fetch import fetch_reference as _ncbi_fetch

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genome_dlist_path = Path(genome_dlist) if genome_dlist else None
    if genome_dlist_path and not genome_dlist_path.is_file():
        raise ValueError(f"Genome D-list FASTA does not exist: {genome_dlist_path}")

    ncbi_cache = Path(cache_dir) / "ncbi" if cache_dir else None

    log.info("Step 1/5  Fetching host cDNA for '%s' …", host_species)
    # NB: the chromosomal GTF returned here is intentionally NOT used for the combined
    # GTF (see Step 5) — it is kept only for provenance. The combined GTF is generated
    # from the cDNA FASTA headers so seqnames match.
    host_fasta_gz, _host_gtf_gz = fetch_host_cdna(host_species, out_dir / "host", cache_dir)

    log.info("Step 2/5  Fetching %d viral accessions from NCBI …", len(virus_accessions))
    viral_fasta_path, viral_gtf_path = _ncbi_fetch(
        virus_accessions,
        out_dir=out_dir / "viral",
        email=email,
        api_key=api_key,
        cache_dir=ncbi_cache,
    )

    if include_anellovirus:
        from viralscan.anellovirus import load_accession_table as _load_anello_table
        from viralscan.scripts.ncbi_fetch import (
            DEFAULT_CACHE_DIR as _NCBI_DEFAULT_CACHE,
        )
        from viralscan.scripts.ncbi_fetch import (
            NCBIFetchError as _NCBIFetchError,
        )
        from viralscan.scripts.ncbi_fetch import (
            _fetch_one,
        )

        anello_rows = _load_anello_table()
        all_anello_accs = {row["accession"].strip() for row in anello_rows}
        new_anello = sorted(all_anello_accs - set(virus_accessions))
        anello_cache = ncbi_cache if ncbi_cache is not None else _NCBI_DEFAULT_CACHE

        log.info(
            "Step 2b/5  Including %d anellovirus accessions (use --no-anellovirus to skip).",
            len(new_anello),
        )

        anello_failures: list[str] = []
        with open(viral_fasta_path, "a") as _anello_fh:
            for i, acc in enumerate(new_anello, 1):
                if i % 250 == 0:
                    log.info("  Anellovirus fetch progress: %d / %d", i, len(new_anello))
                try:
                    fasta_p, _ = _fetch_one(acc, anello_cache, email, api_key)
                    text = fasta_p.read_text()
                    if text and not text.endswith("\n"):
                        text += "\n"
                    _anello_fh.write(text)
                except _NCBIFetchError as exc:
                    anello_failures.append(f"{acc}: {exc}")

        if anello_failures:
            missing_report = out_dir / "missing_accessions.tsv"
            missing_report.write_text(
                "accession\terror\n"
                + "\n".join(
                    failure.split(": ", 1)[0] + "\t" + failure.split(": ", 1)[-1]
                    for failure in anello_failures
                )
                + "\n",
                encoding="utf-8",
            )
            fail_frac = len(anello_failures) / len(new_anello) if new_anello else 0.0
            log.warning(
                "Anellovirus fetch: %d / %d accessions failed (%.0f%%).",
                len(anello_failures),
                len(new_anello),
                fail_frac * 100,
            )
            if not allow_partial_panel:
                raise RuntimeError(
                    f"Anellovirus fetch failed for {len(anello_failures)}/{len(new_anello)} "
                    f"accessions. See {missing_report}. Re-run after fixing retrieval, or "
                    "explicitly use --allow-partial-panel."
                )
            log.warning(
                "Continuing with %d successfully-fetched anellovirus accessions.",
                len(new_anello) - len(anello_failures),
            )

    log.info("Step 3/5  Building viral GTF …")
    with open(viral_fasta_path) as fh:
        viral_fasta_text = fh.read()

    # Build per-accession GTF blocks using accession-specific FASTA
    # (ncbi_fetch returns a concatenated FASTA; we split on accession headers)
    anello_accessions: set[str] = set()
    if include_anellovirus:
        from viralscan.anellovirus import load_accession_table

        anello_accessions = {row["accession"].strip() for row in load_accession_table()}

    # CAT-01: reuse the real CDS structure NCBI already returned. Before this,
    # the fetched GTF was discarded and every non-anellovirus accession became a
    # single whole-genome gene, which silently disabled gene programmes.
    real_gtf_blocks: dict[str, list[str]] = {}
    try:
        with open(viral_gtf_path) as gtf_fh:
            real_gtf_blocks = index_gtf_by_seqname(gtf_fh.read())
    except OSError as exc:
        log.warning("Could not read the fetched viral GTF (%s); using placeholders.", exc)

    viral_gtf_lines: list[str] = []
    annotation_sources: Counter[str] = Counter()
    placeholder_accessions: list[str] = []
    current_acc = None
    current_lines: list[str] = []

    def _flush_block() -> None:
        if not current_acc or not current_lines:
            return
        block_gtf, source = viral_gtf_block(
            "\n".join(current_lines),
            current_acc,
            anello_accessions=anello_accessions,
            real_gtf_blocks=real_gtf_blocks,
        )
        annotation_sources[source] += 1
        if source == "placeholder":
            placeholder_accessions.append(current_acc)
        if block_gtf:
            viral_gtf_lines.append(block_gtf)

    for raw in viral_fasta_text.splitlines():
        line = raw.strip()
        if line.startswith(">"):
            _flush_block()
            # Extract accession from header (first token, strip ">")
            header_token = line[1:].split()[0]
            # Keep the versioned accession (e.g. "NC_045512.2") so it matches
            # the GTF gene_id and anello_name_map keys.
            current_acc = header_token
            current_lines = [line]
        else:
            current_lines.append(line)

    _flush_block()

    total_blocks = sum(annotation_sources.values())
    log.info(
        "Viral annotation: %d real CDS from the anellovirus catalogue, %d real CDS from NCBI, "
        "%d whole-genome placeholders (of %d accessions).",
        annotation_sources["catalogue"],
        annotation_sources["ncbi"],
        annotation_sources["placeholder"],
        total_blocks,
    )
    if placeholder_accessions:
        # A placeholder is a competition bucket, not a measurement: gene
        # programmes cannot resolve a marker against a whole-genome feature.
        log.warning(
            "%d accession(s) have no CDS annotation and fall back to a whole-genome "
            "feature; gene programmes will not resolve for them: %s%s",
            len(placeholder_accessions),
            ", ".join(placeholder_accessions[:10]),
            " …" if len(placeholder_accessions) > 10 else "",
        )

    our_viral_gtf = out_dir / "viral" / "viral_whole_genome.gtf"
    our_viral_gtf.parent.mkdir(parents=True, exist_ok=True)
    with open(our_viral_gtf, "w") as fh:
        fh.write("\n".join(viral_gtf_lines))
        if viral_gtf_lines:
            fh.write("\n")

    log.info("Step 4/5  Concatenating FASTA …")
    combined_fasta = out_dir / "combined.fa"
    with open(combined_fasta, "wb") as out_fh:
        # Decompress host cDNA gzip into combined
        with gzip.open(host_fasta_gz, "rb") as gz_fh:
            shutil.copyfileobj(gz_fh, out_fh)
        # Append viral FASTA (plain text from NCBI fetch)
        with open(viral_fasta_path, "rb") as vf:
            shutil.copyfileobj(vf, out_fh)

    log.info("Step 5/5  Concatenating GTF …")
    # The Ensembl companion GTF (host_gtf_gz) is *chromosomal* (seqnames 1/2/X) and
    # does NOT match the cDNA FASTA headers (ENST…), which would make kb ref hang at
    # "Splitting genome". Generate a cDNA-level host GTF from the FASTA instead.
    host_cdna_gtf = out_dir / "host" / "host_cdna.gtf"
    host_cdna_gtf.parent.mkdir(parents=True, exist_ok=True)
    n_host_tx = host_cdna_as_gtf(host_fasta_gz, host_cdna_gtf)
    log.info("  Host cDNA GTF: %s (%d transcripts)", host_cdna_gtf, n_host_tx)

    combined_gtf = out_dir / "combined.gtf"
    with open(combined_gtf, "wb") as out_fh:
        with open(host_cdna_gtf, "rb") as host_fh:
            shutil.copyfileobj(host_fh, out_fh)
        # Append viral GTF
        with open(our_viral_gtf, "rb") as vf:
            shutil.copyfileobj(vf, out_fh)

    log.info("Combined FASTA: %s", combined_fasta)
    log.info("Combined GTF:   %s", combined_gtf)

    viral_identifiers = {
        line[1:].split()[0] for line in viral_fasta_text.splitlines() if line.startswith(">")
    }
    manifest_profile = "anellovirus-expanded" if include_anellovirus else profile
    homology_annotations = (
        measure_host_homology(
            viral_fasta_path,
            genome_dlist_path,
            out_dir / "host_homology_annotations.tsv",
        )
        if genome_dlist_path
        else None
    )
    manifest_path = write_reference_manifest(
        combined_fasta,
        out_dir / "reference_manifest.json",
        profile=manifest_profile,
        host_species=host_species,
        viral_identifiers=viral_identifiers,
        annotations=homology_annotations,
        genome_dlist=genome_dlist_path,
    )

    index_path: Optional[Path] = None
    t2g_path: Optional[Path] = None

    if run_kb_ref:
        kb_bin = shutil.which("kb")
        if kb_bin is None:
            raise RuntimeError(
                "'kb' not found on PATH but index construction was requested. "
                "Install the full ViralScan environment or pass --no-kb-ref explicitly."
            )
        else:
            index_path = out_dir / "index.idx"
            t2g_path = out_dir / "t2g.txt"
            cdna_fa = out_dir / "cdna.fa"  # kb ref -f1 output
            cmd = [
                kb_bin,
                "ref",
                "-i",
                str(index_path),
                "-g",
                str(t2g_path),
                "-f1",
                str(cdna_fa),
            ]
            if genome_dlist_path:
                cmd.extend(["--d-list", str(genome_dlist_path)])
            cmd.extend([str(combined_fasta), str(combined_gtf)])
            log.info("Running: %s", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True)  # noqa: S603
                log.info("kb ref complete. Index: %s", index_path)
            except subprocess.CalledProcessError as exc:
                log.error(
                    "kb ref failed (exit %d); combined files are still available.", exc.returncode
                )
                raise  # propagate — caller decides whether to abort

    return {
        "fasta": combined_fasta,
        "gtf": combined_gtf,
        "index": index_path,
        "t2g": t2g_path,
        "manifest": manifest_path,
    }


# ---------------------------------------------------------------------------
# Anellovirus-specific reference builder  (B.1 – B.4)
# ---------------------------------------------------------------------------


def _run_dustmasker(fasta_in: Path, fasta_out: Path) -> bool:
    """Run dustmasker to hard-mask low-complexity regions.

    Uses window=64, level=30 (clareaulab parameters). Returns True if masking
    ran successfully, False if dustmasker is not on PATH (masked → unmasked
    copy is written to *fasta_out* in the False case via the caller).
    """
    binary = shutil.which("dustmasker")
    if binary is None:
        log.warning(
            "dustmasker not found on PATH — skipping hard-masking. "
            "Install NCBI BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs"
            "&DOC_TYPE=Download) to enable this step."
        )
        return False
    cmd = [
        binary,
        "-in",
        str(fasta_in),
        "-out",
        str(fasta_out),
        "-outfmt",
        "fasta",
        "-window",
        "64",
        "-level",
        "30",
    ]
    log.info("Running: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)  # noqa: S603
    except subprocess.CalledProcessError as exc:
        log.error(
            "dustmasker failed (exit %d); continuing without hard-masking.",
            exc.returncode,
        )
        return False
    log.info("Hard-masking complete: %s", fasta_out)
    return True


def _run_cdhit_est(fasta_in: Path, fasta_out: Path, identity: float = 0.95) -> bool:
    """Run cd-hit-est to cluster near-identical sequences.

    Returns True if clustering ran, False if cd-hit-est is not on PATH.
    """
    binary = shutil.which("cd-hit-est")
    if binary is None:
        log.warning(
            "cd-hit-est not found on PATH — skipping clustering. "
            "Install CD-HIT (https://github.com/weizhongli/cdhit) to enable."
        )
        return False
    word_size = 8 if identity >= 0.9 else (7 if identity >= 0.88 else 6)
    cmd = [
        binary,
        "-i",
        str(fasta_in),
        "-o",
        str(fasta_out),
        "-c",
        str(identity),
        "-n",
        str(word_size),
        "-M",
        "8000",
        "-T",
        "0",
        "-d",
        "0",  # keep full sequence name
    ]
    log.info("Running: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)  # noqa: S603
    except subprocess.CalledProcessError as exc:
        log.error(
            "cd-hit-est failed (exit %d); continuing without clustering.",
            exc.returncode,
        )
        return False
    log.info("Clustering complete: %s", fasta_out)
    return True


def index_gtf_by_seqname(gtf_text: str) -> dict[str, list[str]]:
    """Group GTF lines by their seqname (column 1).

    ``ncbi_fetch.fetch_reference`` returns a *merged* GTF carrying the real CDS
    structure for every accession it fetched, with genome-scoped gene IDs like
    ``NC_001526.4_HpV16gp3``. Splitting it per accession lets the combined-
    reference builder reuse that annotation instead of discarding it (PLAN
    `CAT-01`). Comment and blank lines are dropped.
    """
    blocks: dict[str, list[str]] = {}
    for raw in gtf_text.splitlines():
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        seqname = line.split("\t", 1)[0].strip()
        if not seqname:
            continue
        blocks.setdefault(seqname, []).append(line)
    return blocks


def viral_gtf_block(
    fasta_text: str,
    accession: str,
    *,
    anello_accessions: Optional[set[str]] = None,
    real_gtf_blocks: Optional[dict[str, list[str]]] = None,
) -> tuple[str, str]:
    """Return ``(gtf_text, source)`` for one viral accession, best annotation first.

    Order, most informative first (PLAN `CAT-01`):

    1. ``catalogue`` — the packaged anellovirus CDS catalogue.
    2. ``ncbi`` — the real GTF NCBI already returned for this accession.
    3. ``placeholder`` — one whole-genome gene, used **only** when neither of the
       above covers the record.

    Before `CAT-01` every non-anellovirus accession took branch 3, so a natively
    built index carried one ``{accession}_gene1`` bucket per genome and gene
    programmes had nothing to resolve against.
    """
    if anello_accessions and accession in anello_accessions:
        block = _catalogued_anello_gtf(fasta_text, accession)
        if block.strip():
            return block, "catalogue"
    if real_gtf_blocks:
        lines = real_gtf_blocks.get(accession)
        if lines is None:
            # Match across a version mismatch in either direction: the FASTA
            # header and the GTF seqname do not always agree on the suffix.
            bare = accession.split(".")[0]
            for key, value in real_gtf_blocks.items():
                if key.split(".")[0] == bare:
                    lines = value
                    break
        if lines:
            return "\n".join(lines), "ncbi"
    return _genome_as_transcript_gtf(fasta_text, accession), "placeholder"


def _catalogued_anello_gtf(fasta_text: str, accession: str) -> str:
    """Real CDS structure for one anellovirus accession, else the placeholder.

    A thin wrapper over :func:`viralscan.anellovirus.gtf_text_for` so the
    combined-reference path can annotate one accession at a time while streaming
    a merged FASTA.
    """
    from viralscan.anellovirus import gtf_text_for

    return gtf_text_for([accession], fasta_texts={accession: fasta_text})


def _gtf_from_merged_fasta(fasta_path: Path, gtf_path: Path) -> None:
    """Split *fasta_path* by accession header and emit whole-genome GTF to *gtf_path*.

    This is the *fallback* annotation.  When the packaged
    ``anellovirus_genes.tsv`` covers the accessions, prefer
    :func:`viralscan.anellovirus.gtf_text_for`, which emits the real NCBI CDS
    structure; a whole-genome single-exon gene is a competition bucket rather
    than a measurement, because a whole-genome transcript shares sequence with
    every other genome in the panel and reads cross-map in proportion to
    conservation.
    """
    gtf_blocks: list[str] = []
    current_acc: Optional[str] = None
    current_lines: list[str] = []

    def _flush() -> None:
        if current_acc and current_lines:
            block = _genome_as_transcript_gtf("\n".join(current_lines), current_acc)
            if block:
                gtf_blocks.append(block)

    with open(fasta_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                current_acc = line[1:].split()[0]
                current_lines = [line]
            else:
                current_lines.append(line)
    _flush()

    with open(gtf_path, "w") as fh:
        fh.write("\n".join(gtf_blocks))
        if gtf_blocks:
            fh.write("\n")


def _anellovirus_gtf(fasta_path: Path, gtf_path: Path) -> tuple[int, int]:
    """Write an anellovirus GTF, preferring the packaged real-gene catalogue.

    Returns ``(annotated_accessions, placeholder_accessions)``.  Any accession the
    catalogue does not cover still gets a whole-genome placeholder, because
    ``kb ref`` silently drops a sequence that has no GTF row and the genome would
    then be neither quantified nor detectable.
    """
    from viralscan.anellovirus import gtf_text_for, load_gene_table

    order: list[str] = []
    texts: dict[str, str] = {}
    current_acc: Optional[str] = None
    current_lines: list[str] = []

    def _flush() -> None:
        if current_acc and current_lines:
            order.append(current_acc)
            texts[current_acc] = "\n".join(current_lines)

    with open(fasta_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                current_acc = line[1:].split()[0]
                current_lines = [line]
            else:
                current_lines.append(line)
    _flush()

    catalogue = {str(row["accession"]).strip() for row in load_gene_table()}
    annotated = {acc for acc in order if acc in catalogue}
    with open(gtf_path, "w") as fh:
        fh.write(gtf_text_for(order, fasta_texts=texts))
    return len(annotated), len(order) - len(annotated)


def build_anellovirus_reference(
    out_dir: os.PathLike[str] | str,
    accessions: Optional[list[str]] = None,
    mask: bool = True,
    cluster: bool = False,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
    cache_dir: Optional[os.PathLike[str] | str] = None,
    run_kb_ref: bool = True,
    fasta_path: Optional[Path] = None,
    genome_dlist: Optional[os.PathLike[str] | str] = None,
) -> dict[str, Optional[Path]]:
    """Build a kallisto-ready Anelloviridae reference.

    When *fasta_path* is ``None`` (the default), sequences are downloaded from
    NCBI for all accessions in the packaged ``anellovirus_accessions.tsv`` (or
    the explicit *accessions* subset).  When *fasta_path* is supplied (e.g.
    the bundled FASTA from ``viralscan data fetch``), the NCBI download step is
    skipped entirely and the provided FASTA is used as the starting point.

    After obtaining the FASTA the pipeline optionally hard-masks low-complexity
    regions with ``dustmasker``, optionally clusters with ``cd-hit-est``,
    regenerates per-genome GTFs, and optionally runs ``kb ref`` to produce a
    kallisto index.

    Parameters
    ----------
    out_dir:
        Destination directory for output files.
    accessions:
        Explicit list of NCBI accession numbers (only used when *fasta_path*
        is ``None``).  Defaults to all ~2,042 accessions in the packaged TSV.
    mask:
        Hard-mask low-complexity regions with ``dustmasker -window 64
        -level 30``.  A requested mask step fails if ``dustmasker`` is absent.
    cluster:
        Cluster near-identical sequences with ``cd-hit-est -c 0.95``.
        Off by default because the packaged table already uses CD-HIT
        representatives.  A requested clustering step fails if ``cd-hit-est``
        is absent.
    email:
        E-mail address for NCBI E-utilities (only used when *fasta_path* is
        ``None``; required per NCBI policy).
    api_key:
        NCBI API key for higher request rates (only used when *fasta_path* is
        ``None``).
    cache_dir:
        Cache root; defaults to ``~/.cache/viralscan``.
    run_kb_ref:
        Build a kallisto index + t2g via ``kb ref``.  A requested index step
        fails if ``kb`` is absent from PATH.
    fasta_path:
        Pre-built merged FASTA to use instead of downloading from NCBI.  When
        provided the NCBI fetch step (Step 1) is skipped.

    Returns
    -------
    dict with keys ``fasta``, ``gtf``, ``index`` (None if not built), ``t2g``
    (None if not built).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genome_dlist_path = Path(genome_dlist) if genome_dlist else None
    if genome_dlist_path and not genome_dlist_path.is_file():
        raise ValueError(f"Genome D-list FASTA does not exist: {genome_dlist_path}")

    if fasta_path is not None:
        log.info("Step 1/4  Using provided FASTA, skipping NCBI download: %s", fasta_path)
        working_fasta = Path(fasta_path)
    else:
        from viralscan.anellovirus import load_accession_table
        from viralscan.scripts.ncbi_fetch import fetch_reference as _ncbi_fetch

        ncbi_out = out_dir / "ncbi"
        ncbi_cache = Path(cache_dir) / "ncbi" if cache_dir else None

        if accessions is None:
            rows = load_accession_table()
            accessions = [r["accession"] for r in rows]
            log.info(
                "Loaded %d accessions from packaged anellovirus_accessions.tsv", len(accessions)
            )

        log.info("Step 1/4  Fetching %d Anelloviridae accessions from NCBI …", len(accessions))
        merged_fasta, _ncbi_gtf = _ncbi_fetch(
            accessions,
            out_dir=ncbi_out,
            email=email,
            api_key=api_key,
            cache_dir=ncbi_cache,
        )
        working_fasta = merged_fasta

    # Start with the raw/provided FASTA; optionally mask then cluster.

    if mask:
        log.info("Step 2/4  Hard-masking with dustmasker …")
        masked_fasta = out_dir / "anellovirus.masked.fa"
        ran = _run_dustmasker(working_fasta, masked_fasta)
        if ran:
            working_fasta = masked_fasta
        else:
            raise RuntimeError(
                "Anellovirus masking was requested but dustmasker was unavailable or failed. "
                "Install BLAST+ or pass --no-mask explicitly."
            )
    else:
        log.info("Step 2/4  Masking disabled — skipping dustmasker.")

    # A requested mask step failing is not the only way an unmasked panel reaches
    # the index: `--no-mask`, or a build path that never called dustmasker at all.
    # Verify the property that actually matters — whether any *k-mer* is
    # low-complexity — rather than trusting that a mask ran.
    max_low_complexity_fraction = 0.0 if mask else 0.05
    max_pure_homopolymer_kmers = 0 if mask else 2
    report = low_complexity_report(working_fasta)
    overall_low = sum(low for low, _ in report.values())
    overall_total = sum(total for _, total in report.values())
    if report:
        worst_id, (worst_low, worst_total) = max(
            report.items(), key=lambda row: (row[1][0] / row[1][1]) if row[1][1] else 0.0
        )
        pure = sum(
            low_complexity_kmer_counts(sequence)["pure_homopolymer"]
            for _, sequence in _fasta_records(working_fasta)
        )
        log.info(
            "Low-complexity k-mer scan: %s/%s (%.4f%%) low-complexity, %d pure-homopolymer, "
            "across %d records; worst record %s at %d/%s",
            f"{overall_low:,}",
            f"{overall_total:,}",
            100 * overall_low / overall_total if overall_total else 0.0,
            pure,
            len(report),
            worst_id,
            worst_low,
            worst_total,
        )
    # Structural problems (empty, duplicate ID, duplicate sequence) keep their own
    # ValueError; only the low-complexity verdict is reported as a gate failure.
    validate_reference_records(working_fasta)
    try:
        validate_reference_records(
            working_fasta,
            max_low_complexity_fraction=max_low_complexity_fraction,
            max_pure_homopolymer_kmers=max_pure_homopolymer_kmers,
        )
    except ValueError as exc:
        raise RuntimeError(
            "Anellovirus panel failed the low-complexity k-mer gate "
            f"(pure-homopolymer limit {max_pure_homopolymer_kmers}, fraction limit "
            f"{max_low_complexity_fraction:.2f} per record): {exc}"
        ) from exc
    log.info(
        "Low-complexity k-mer gate passed (pure-homopolymer <= %d, fraction <= %.2f)",
        max_pure_homopolymer_kmers,
        max_low_complexity_fraction,
    )

    if cluster:
        log.info("Step 3/4  Clustering with cd-hit-est …")
        clustered_fasta = out_dir / "anellovirus.clustered.fa"
        ran = _run_cdhit_est(working_fasta, clustered_fasta)
        if ran:
            working_fasta = clustered_fasta
        else:
            raise RuntimeError(
                "Anellovirus clustering was requested but cd-hit-est was unavailable or failed. "
                "Install CD-HIT or omit --cluster explicitly."
            )
    else:
        log.info("Step 3/4  Clustering disabled — skipping cd-hit-est.")

    # Write the final FASTA to the canonical output name.
    final_fasta = out_dir / "anellovirus.fa"
    if working_fasta != final_fasta:
        shutil.copy2(working_fasta, final_fasta)
    records = validate_reference_records(final_fasta)

    log.info("Step 4/4  Building viral GTF …")
    final_gtf = out_dir / "anellovirus.gtf"
    n_annotated, n_placeholder = _anellovirus_gtf(final_fasta, final_gtf)
    log.info(
        "  %d/%d genomes carry real NCBI CDS structure; %d fall back to a "
        "whole-genome placeholder (record has no CDS feature).",
        n_annotated,
        n_annotated + n_placeholder,
        n_placeholder,
    )

    log.info("Anellovirus FASTA: %s", final_fasta)
    log.info("Anellovirus GTF:   %s", final_gtf)

    homology_annotations = (
        measure_host_homology(
            final_fasta,
            genome_dlist_path,
            out_dir / "host_homology_annotations.tsv",
        )
        if genome_dlist_path
        else None
    )
    manifest_path = write_reference_manifest(
        final_fasta,
        out_dir / "reference_manifest.json",
        profile="anellovirus-representative" if cluster else "anellovirus-expanded",
        host_species="none",
        viral_identifiers={identifier for identifier, _sequence in records},
        annotations=homology_annotations,
        genome_dlist=genome_dlist_path,
    )

    index_path: Optional[Path] = None
    t2g_path: Optional[Path] = None

    if run_kb_ref:
        kb_bin = shutil.which("kb")
        if kb_bin is None:
            raise RuntimeError(
                "'kb' not found on PATH but index construction was requested. "
                "Install the full ViralScan environment or pass --no-kb-ref explicitly."
            )
        else:
            index_path = out_dir / "index.idx"
            t2g_path = out_dir / "t2g.txt"
            cdna_fa = out_dir / "cdna.fa"
            cmd = [
                kb_bin,
                "ref",
                "-i",
                str(index_path),
                "-g",
                str(t2g_path),
                "-f1",
                str(cdna_fa),
            ]
            if genome_dlist_path:
                cmd.extend(["--d-list", str(genome_dlist_path)])
            cmd.extend([str(final_fasta), str(final_gtf)])
            log.info("Running: %s", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True)  # noqa: S603
                log.info("kb ref complete. Index: %s", index_path)
            except subprocess.CalledProcessError as exc:
                log.error(
                    "kb ref failed (exit %d); FASTA and GTF are still available.", exc.returncode
                )
                raise  # propagate — caller decides whether to abort

    return {
        "fasta": final_fasta,
        "gtf": final_gtf,
        "index": index_path,
        "t2g": t2g_path,
        "manifest": manifest_path,
    }


# ---------------------------------------------------------------------------
# CLI entry point (called from menu.py build-ref subcommand)
# ---------------------------------------------------------------------------


def build_ref_main(args: argparse.Namespace) -> None:
    """Orchestrator called by ``viralscan build-ref``."""
    from viralscan.utils import configure_logging

    configure_logging(
        verbose=bool(getattr(args, "verbose", False)),
        quiet=bool(getattr(args, "quiet", False)),
    )

    if getattr(args, "list_species", False):
        print("Supported host species:")
        for key, (ens, asm) in sorted(ENSEMBL_SPECIES.items()):
            print(f"  {key:<16} ({ens}, {asm})")
        sys.exit(0)

    # Fail before any download if the requested index cannot be constructed.
    if not getattr(args, "no_kb_ref", False) and shutil.which("kb") is None:
        log.error(
            "'kb' is not on PATH but index construction was requested. Install the full "
            "ViralScan environment or pass --no-kb-ref explicitly."
        )
        sys.exit(2)
    genome_dlist = getattr(args, "genome_dlist", None)
    if genome_dlist and not Path(genome_dlist).is_file():
        log.error("--genome-dlist does not exist or is not a file: %s", genome_dlist)
        sys.exit(2)
    if genome_dlist and shutil.which("minimap2") is None:
        log.error("--genome-dlist requires minimap2 for host-homology annotation.")
        sys.exit(2)

    reference_panel = getattr(args, "reference_panel", None)
    if reference_panel == "anellovirus":
        bundled_fasta: Optional[Path] = None
        from viralscan.data_fetch import (
            ViralScanDataError as _DataError,
        )
        from viralscan.data_fetch import (
            bundled_anellovirus_fasta,
        )

        try:
            bundled_fasta = bundled_anellovirus_fasta(getattr(args, "cache_dir", None))
            log.info("Using bundled anellovirus FASTA from Zenodo cache: %s", bundled_fasta)
        except _DataError as exc:
            log.info("Bundled FASTA not available (%s); falling back to NCBI download.", exc)
        try:
            result = build_anellovirus_reference(
                out_dir=args.output,
                accessions=getattr(args, "virus_accessions", None),
                mask=not getattr(args, "no_mask", False),
                cluster=getattr(args, "cluster", False),
                email=getattr(args, "ncbi_email", None),
                api_key=getattr(args, "ncbi_api_key", None),
                cache_dir=getattr(args, "cache_dir", None),
                run_kb_ref=not getattr(args, "no_kb_ref", False),
                fasta_path=bundled_fasta,
                genome_dlist=genome_dlist,
            )
        except (subprocess.CalledProcessError, RuntimeError, ValueError) as exc:
            # Error already logged by the builder.
            log.error("Anellovirus reference build failed: %s", exc)
            sys.exit(1)
        print("\nAnellovirus reference build complete.")
        print(f"  FASTA          : {result['fasta']}")
        print(f"  GTF            : {result['gtf']}")
        if result["index"]:
            print(f"  kallisto index : {result['index']}")
            print(f"  t2g mapping    : {result['t2g']}")
        else:
            print("  kallisto index : not built (run 'kb ref' manually if needed)")
        return

    if not args.host:
        log.error(
            "--host is required (e.g. --host human). "
            "Use --reference-panel anellovirus for an anellovirus-only reference."
        )
        sys.exit(1)

    if not args.virus_accessions:
        log.error("--virus-accessions is required")
        sys.exit(1)

    include_anello = getattr(args, "anellovirus", False)
    if include_anello:
        log.info(
            "Anellovirus accessions will be included in the combined reference (explicit opt-in)."
        )

    try:
        result = build_combined_reference(
            host_species=args.host,
            virus_accessions=args.virus_accessions,
            out_dir=args.output,
            email=getattr(args, "ncbi_email", None),
            api_key=getattr(args, "ncbi_api_key", None),
            cache_dir=getattr(args, "cache_dir", None),
            run_kb_ref=not getattr(args, "no_kb_ref", False),
            include_anellovirus=include_anello,
            allow_partial_panel=getattr(args, "allow_partial_panel", False),
            profile=getattr(args, "profile", "curated"),
            genome_dlist=getattr(args, "genome_dlist", None),
        )
    except (subprocess.CalledProcessError, RuntimeError, ValueError) as exc:
        # Error already logged by the builder.
        log.error("Reference build failed: %s", exc)
        sys.exit(1)

    print("\nReference build complete.")
    print(f"  Combined FASTA : {result['fasta']}")
    print(f"  Combined GTF   : {result['gtf']}")
    if result["index"]:
        print(f"  kallisto index : {result['index']}")
        print(f"  t2g mapping    : {result['t2g']}")
    else:
        print("  kallisto index : not built (run 'kb ref' manually if needed)")
