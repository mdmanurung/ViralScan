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
import gzip
import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Optional

from viralscan.constants import ENSEMBL_SPECIES

log = logging.getLogger("viralscan")

# ---------------------------------------------------------------------------
# Ensembl HTTPS-FTP mirror helpers
# ---------------------------------------------------------------------------

_ENSEMBL_FTP = "https://ftp.ensembl.org/pub/current_fasta/{species}/cdna/"
_ENSEMBL_GTF = "https://ftp.ensembl.org/pub/current_gtf/{species}/"


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


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _list_ensembl_files(species_name: str, url_base: str) -> list[str]:
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

    try:
        with urllib.request.urlopen(url_base, timeout=30) as resp:  # noqa: S310
            html_bytes = resp.read()
    except Exception as exc:
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
    cdna_base = _ENSEMBL_FTP.format(species=ens_name)
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
    gtf_base = _ENSEMBL_GTF.format(species=ens_name)
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
) -> dict[str, Optional[Path]]:
    """Build a combined host + virus kallisto reference.

    Steps
    -----
    1. Download Ensembl cDNA FASTA + GTF for *host_species*.
    2. Download NCBI FASTA for each accession in *virus_accessions*
       (via :func:`viralscan.scripts.ncbi_fetch.fetch_reference`).
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

    Returns
    -------
    dict with keys: ``fasta``, ``gtf``, ``index`` (None if *run_kb_ref* is
    False), ``t2g`` (None if *run_kb_ref* is False).
    """
    # Late import to avoid circular dependency at module load time.
    from viralscan.scripts.ncbi_fetch import fetch_reference as _ncbi_fetch

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ncbi_cache = Path(cache_dir) / "ncbi" if cache_dir else None

    log.info("Step 1/5  Fetching host cDNA for '%s' …", host_species)
    host_fasta_gz, host_gtf_gz = fetch_host_cdna(host_species, out_dir / "host", cache_dir)

    log.info("Step 2/5  Fetching %d viral accessions from NCBI …", len(virus_accessions))
    viral_fasta_path, _viral_gtf_path = _ncbi_fetch(
        virus_accessions,
        out_dir=out_dir / "viral",
        email=email,
        api_key=api_key,
        cache_dir=ncbi_cache,
    )

    log.info("Step 3/5  Building whole-genome viral GTF …")
    # ncbi_fetch already writes a GTF, but we regenerate from our helper to
    # ensure consistent gene_biotype = "whole_genome" formatting.
    with open(viral_fasta_path) as fh:
        viral_fasta_text = fh.read()

    # Build per-accession GTF blocks using accession-specific FASTA
    # (ncbi_fetch returns a concatenated FASTA; we split on accession headers)
    viral_gtf_lines: list[str] = []
    current_acc = None
    current_lines: list[str] = []

    for raw in viral_fasta_text.splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if current_acc and current_lines:
                block_gtf = _genome_as_transcript_gtf("\n".join(current_lines), current_acc)
                if block_gtf:
                    viral_gtf_lines.append(block_gtf)
            # Extract accession from header (first token, strip ">")
            header_token = line[1:].split()[0]
            # Keep the versioned accession (e.g. "NC_045512.2") so it matches
            # the GTF gene_id and anello_name_map keys.
            current_acc = header_token
            current_lines = [line]
        else:
            current_lines.append(line)

    if current_acc and current_lines:
        block_gtf = _genome_as_transcript_gtf("\n".join(current_lines), current_acc)
        if block_gtf:
            viral_gtf_lines.append(block_gtf)

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
    combined_gtf = out_dir / "combined.gtf"
    with open(combined_gtf, "wb") as out_fh:
        # Decompress host GTF gzip into combined
        with gzip.open(host_gtf_gz, "rb") as gz_fh:
            shutil.copyfileobj(gz_fh, out_fh)
        # Append viral GTF
        with open(our_viral_gtf, "rb") as vf:
            shutil.copyfileobj(vf, out_fh)

    log.info("Combined FASTA: %s", combined_fasta)
    log.info("Combined GTF:   %s", combined_gtf)

    index_path: Optional[Path] = None
    t2g_path: Optional[Path] = None

    if run_kb_ref:
        kb_bin = shutil.which("kb")
        if kb_bin is None:
            log.warning(
                "'kb' not found on PATH; skipping kb ref. "
                "Install kb-python and re-run with the same output directory."
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
                str(combined_fasta),
                str(combined_gtf),
            ]
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
        "-in", str(fasta_in),
        "-out", str(fasta_out),
        "-outfmt", "fasta",
        "-window", "64",
        "-level", "30",
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
        "-i", str(fasta_in),
        "-o", str(fasta_out),
        "-c", str(identity),
        "-n", str(word_size),
        "-M", "8000",
        "-T", "0",
        "-d", "0",  # keep full sequence name
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


def _gtf_from_merged_fasta(fasta_path: Path, gtf_path: Path) -> None:
    """Split *fasta_path* by accession header and emit whole-genome GTF to *gtf_path*."""
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


def build_anellovirus_reference(
    out_dir: os.PathLike[str] | str,
    accessions: Optional[list[str]] = None,
    mask: bool = True,
    cluster: bool = False,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
    cache_dir: Optional[os.PathLike[str] | str] = None,
    run_kb_ref: bool = True,
) -> dict[str, Optional[Path]]:
    """Build a kallisto-ready Anelloviridae reference from NCBI.

    Downloads sequences for all accessions in the packaged
    ``anellovirus_accessions.tsv`` (or a user-supplied subset), regenerates
    per-genome GTFs using ``_genome_as_transcript_gtf``, optionally hard-masks
    low-complexity regions with ``dustmasker``, optionally clusters
    near-identical sequences with ``cd-hit-est``, and optionally runs
    ``kb ref`` to produce a kallisto index.

    Parameters
    ----------
    out_dir:
        Destination directory for output files.
    accessions:
        Explicit list of NCBI accession numbers.  Defaults to all ~2,042
        accessions in the packaged ``anellovirus_accessions.tsv``.
    mask:
        Hard-mask low-complexity regions with ``dustmasker -window 64
        -level 30`` (stolen from clareaulab/anellovirus_reference).
        Silently skipped when ``dustmasker`` is not on PATH.
    cluster:
        Cluster near-identical sequences with ``cd-hit-est -c 0.95``.
        Off by default because the packaged table already uses CD-HIT
        representatives.  Silently skipped when ``cd-hit-est`` is not
        on PATH.
    email:
        E-mail address for NCBI E-utilities (required; see NCBI policy).
    api_key:
        NCBI API key for higher request rates.
    cache_dir:
        Cache root; defaults to ``~/.cache/viralscan``.
    run_kb_ref:
        Build a kallisto index + t2g via ``kb ref`` (skipped if ``kb`` is
        absent from PATH).

    Returns
    -------
    dict with keys ``fasta``, ``gtf``, ``index`` (None if not built), ``t2g``
    (None if not built).
    """
    from viralscan.anellovirus import load_accession_table
    from viralscan.scripts.ncbi_fetch import fetch_reference as _ncbi_fetch

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ncbi_out = out_dir / "ncbi"
    ncbi_cache = Path(cache_dir) / "ncbi" if cache_dir else None

    if accessions is None:
        rows = load_accession_table()
        accessions = [r["accession"] for r in rows]
        log.info("Loaded %d accessions from packaged anellovirus_accessions.tsv", len(accessions))

    log.info("Step 1/4  Fetching %d Anelloviridae accessions from NCBI …", len(accessions))
    merged_fasta, _ncbi_gtf = _ncbi_fetch(
        accessions,
        out_dir=ncbi_out,
        email=email,
        api_key=api_key,
        cache_dir=ncbi_cache,
    )

    # Start with the raw merged FASTA; optionally mask then cluster.
    working_fasta = merged_fasta

    if mask:
        log.info("Step 2/4  Hard-masking with dustmasker …")
        masked_fasta = out_dir / "anellovirus.masked.fa"
        ran = _run_dustmasker(working_fasta, masked_fasta)
        if ran:
            working_fasta = masked_fasta
        else:
            log.info("Masking skipped — continuing with unmasked FASTA.")
    else:
        log.info("Step 2/4  Masking disabled — skipping dustmasker.")

    if cluster:
        log.info("Step 3/4  Clustering with cd-hit-est …")
        clustered_fasta = out_dir / "anellovirus.clustered.fa"
        ran = _run_cdhit_est(working_fasta, clustered_fasta)
        if ran:
            working_fasta = clustered_fasta
        else:
            log.info("Clustering skipped — continuing with unclustered FASTA.")
    else:
        log.info("Step 3/4  Clustering disabled — skipping cd-hit-est.")

    # Write the final FASTA to the canonical output name.
    final_fasta = out_dir / "anellovirus.fa"
    if working_fasta != final_fasta:
        shutil.copy2(working_fasta, final_fasta)

    log.info("Step 4/4  Building whole-genome GTF …")
    final_gtf = out_dir / "anellovirus.gtf"
    _gtf_from_merged_fasta(final_fasta, final_gtf)

    log.info("Anellovirus FASTA: %s", final_fasta)
    log.info("Anellovirus GTF:   %s", final_gtf)

    index_path: Optional[Path] = None
    t2g_path: Optional[Path] = None

    if run_kb_ref:
        kb_bin = shutil.which("kb")
        if kb_bin is None:
            log.warning(
                "'kb' not found on PATH; skipping kb ref. "
                "Install kb-python and re-run with the same output directory."
            )
        else:
            index_path = out_dir / "index.idx"
            t2g_path = out_dir / "t2g.txt"
            cdna_fa = out_dir / "cdna.fa"
            cmd = [
                kb_bin,
                "ref",
                "-i", str(index_path),
                "-g", str(t2g_path),
                "-f1", str(cdna_fa),
                str(final_fasta),
                str(final_gtf),
            ]
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

    if getattr(args, "anellovirus", False):
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
            )
        except subprocess.CalledProcessError:
            # Error already logged by the builder.
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
        log.error("--host is required (e.g. --host human). Use --anellovirus for anellovirus-only.")
        sys.exit(1)

    if not args.virus_accessions:
        log.error("--virus-accessions is required")
        sys.exit(1)

    try:
        result = build_combined_reference(
            host_species=args.host,
            virus_accessions=args.virus_accessions,
            out_dir=args.output,
            email=getattr(args, "ncbi_email", None),
            api_key=getattr(args, "ncbi_api_key", None),
            cache_dir=getattr(args, "cache_dir", None),
            run_kb_ref=not getattr(args, "no_kb_ref", False),
        )
    except subprocess.CalledProcessError:
        # Error already logged by the builder.
        sys.exit(1)

    print("\nReference build complete.")
    print(f"  Combined FASTA : {result['fasta']}")
    print(f"  Combined GTF   : {result['gtf']}")
    if result["index"]:
        print(f"  kallisto index : {result['index']}")
        print(f"  t2g mapping    : {result['t2g']}")
    else:
        print("  kallisto index : not built (run 'kb ref' manually if needed)")
