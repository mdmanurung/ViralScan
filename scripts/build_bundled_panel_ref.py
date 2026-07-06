#!/usr/bin/env python3
"""
Build the ViralScan bundled-panel reference index — one-time setup.

Builds a combined HOST + VIRUS kallisto index for joint quantification.
Human cDNA is downloaded from Ensembl via fetch_host_cdna; viral FASTAs are
fetched from NCBI via ncbi_fetch._fetch_one.  Both host and viral transcripts
are quantified together, enabling coexpression analysis between host genes
(ENST*) and viral genes (ADENO_*, EPSTEIN_*, etc.).

Usage:
    NCBI_EMAIL=you@example.com PYTHONPATH=/path/to/ViralScan/src \\
        python scripts/build_bundled_panel_ref.py \\
        --out /exports/para-lipg-hpc/mdmanurung/viralscan_bulk_gse128078/ref

SLURM (activate test_viralscan conda env before sbatch):
    sbatch --job-name=build_panel_ref --cpus-per-task=2 --mem=16G --time=08:00:00 \\
        --wrap "NCBI_EMAIL=... PYTHONPATH=$PWD/src python $PWD/scripts/build_bundled_panel_ref.py --out $WORKDIR/ref"

Pre-conditions verified at startup:
  1. Bundled GTFs contain exon features (kb ref extracts cDNA from exon rows;
     missing exons yield an empty/partial index that wc -l alone won't catch).
  2. Every GTF seqname is a valid NCBI accession (flag non-matching names early).
  3. Every GTF seqname has a matching FASTA header after download (a suppressed
     or missing NCBI accession would otherwise let kb ref silently drop viruses).
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def _find_repo_root() -> Path:
    """Walk up from this script to find the repo root (contains src/viralscan/)."""
    here = Path(__file__).resolve().parent
    for candidate in [here.parent, here]:
        if (candidate / "src" / "viralscan").is_dir():
            return candidate
    sys.exit(
        "ERROR: cannot locate ViralScan repo root (src/viralscan/ not found). "
        "Run from inside the ViralScan repo."
    )


_ACC_RE = re.compile(r"^[A-Za-z]{1,3}_?\d+(\.\d+)?$")


def _extract_gtf_seqnames(gtf_files: list[Path]) -> tuple[set[str], set[str]]:
    """Return (all_seqnames, feature_types) from the bundled GTF corpus."""
    all_seqnames: set[str] = set()
    feature_types: set[str] = set()
    for gtf in gtf_files:
        for line in gtf.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            all_seqnames.add(cols[0])
            feature_types.add(cols[2])
    return all_seqnames, feature_types


def _fasta_seq_ids(fasta_paths: list[Path]) -> set[str]:
    """Collect all sequence IDs (text before first space on '>' lines) from FASTAs."""
    ids: set[str] = set()
    for fp in fasta_paths:
        for line in fp.read_text().splitlines():
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def main() -> None:
    repo = _find_repo_root()
    sys.path.insert(0, str(repo / "src"))

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--out", required=True, type=Path,
        help="Output directory for panel.idx, panel.t2g, cdna.fa, combined.*",
    )
    p.add_argument(
        "--host-species", default="human",
        help="Ensembl host species for the combined host+viral reference (default: human). "
             "Host transcripts (ENST*) are quantified alongside viral ones for coexpression analysis.",
    )
    p.add_argument(
        "--ncbi-email", default=os.environ.get("NCBI_EMAIL"),
        help="Contact email for NCBI E-utilities (or set NCBI_EMAIL env var)",
    )
    p.add_argument(
        "--ncbi-api-key", default=os.environ.get("NCBI_API_KEY"),
        help="NCBI API key for higher rate limits (optional)",
    )
    p.add_argument(
        "--cache-dir", type=Path, default=None,
        help="Override the NCBI fetch cache directory",
    )
    args = p.parse_args()

    if not args.ncbi_email:
        sys.exit("ERROR: NCBI requires an email. Pass --ncbi-email or set NCBI_EMAIL.")

    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    from viralscan.scripts.ncbi_fetch import (  # noqa: E402
        DEFAULT_CACHE_DIR, _fetch_one, NCBIFetchError,
    )
    from viralscan.scripts.build_reference import (  # noqa: E402
        fetch_host_cdna,
        host_cdna_as_gtf,
        _genome_as_transcript_gtf,
    )
    from viralscan.anellovirus import load_accession_table as _load_anello_table  # noqa: E402

    cache_dir: Path = args.cache_dir or DEFAULT_CACHE_DIR

    # ── 1. Download human (host) reference from Ensembl ──────────────────────
    print(f"Step 1/7  Downloading {args.host_species} cDNA + GTF from Ensembl …")
    host_cache = Path.home() / ".cache" / "viralscan" / "ensembl" / args.host_species
    host_fasta_gz, host_gtf_gz = fetch_host_cdna(args.host_species, out / "host", host_cache)
    print(f"  cDNA  : {host_fasta_gz}")
    print(f"  GTF   : {host_gtf_gz}")

    # ── 2. Discover bundled GTFs ──────────────────────────────────────────────
    print("Step 2/7  Scanning bundled viral GTFs …")
    gtf_dir = repo / "src" / "viralscan" / "data"
    gtf_files = sorted(gtf_dir.glob("*.gtf"))
    if not gtf_files:
        sys.exit(f"ERROR: no *.gtf files found in {gtf_dir}")
    print(f"  Found {len(gtf_files)} bundled GTFs in {gtf_dir}")

    # ── 3. Pre-checks ─────────────────────────────────────────────────────────
    print("Step 3/7  Pre-checks …")
    all_seqnames, feature_types = _extract_gtf_seqnames(gtf_files)
    if "exon" not in feature_types:
        sys.exit(
            f"ERROR: bundled GTFs contain no 'exon' features (found: {sorted(feature_types)}). "
            "kb ref extracts cDNA from exon rows — the index would be empty."
        )
    print(f"  'exon' present. All feature types: {sorted(feature_types)}")

    non_accession = sorted(s for s in all_seqnames if not _ACC_RE.match(s))
    if non_accession:
        sys.exit(
            f"ERROR: {len(non_accession)} GTF seqname(s) do not look like NCBI accessions:\n"
            + "\n".join(f"  {s}" for s in non_accession)
        )
    accessions = sorted(all_seqnames)
    print(f"  {len(accessions)} unique NCBI accessions in GTF seqnames")

    # ── 4. Download viral FASTAs from NCBI ───────────────────────────────────
    print(f"Step 4/7  Downloading {len(accessions)} viral FASTAs from NCBI …")
    fasta_paths: list[Path] = []
    errors: list[str] = []
    for i, acc in enumerate(accessions, 1):
        print(f"  [{i:3d}/{len(accessions)}] {acc} … ", end="", flush=True)
        try:
            fasta_path, _ = _fetch_one(acc, cache_dir, args.ncbi_email, args.ncbi_api_key)
            fasta_paths.append(fasta_path)
            print("OK")
        except NCBIFetchError as exc:
            errors.append(f"{acc}: {exc}")
            print(f"FAILED: {exc}")

    if errors:
        sys.exit(
            f"\nERROR: {len(errors)} accession(s) failed to download "
            "(fix and re-run — successful downloads are cached):\n"
            + "\n".join(f"  {e}" for e in errors)
        )

    # ── 4b. Fetch anellovirus panel (clareaulab accessions) ──────────────────
    print("Step 4b/7  Fetching anellovirus panel (clareaulab accessions) …")
    anello_rows = _load_anello_table()
    anello_accs = sorted(
        row["accession"].strip()
        for row in anello_rows
        if row.get("source", "").strip() == "clareaulab"
    )
    print(f"  {len(anello_accs)} clareaulab anellovirus accessions to fetch")

    anello_fasta_texts: list[str] = []
    anello_gtf_texts: list[str] = []
    anello_errors: list[str] = []

    for i, acc in enumerate(anello_accs, 1):
        if i % 250 == 0:
            print(f"  Anellovirus fetch progress: {i} / {len(anello_accs)}", flush=True)
        try:
            fasta_path, _ = _fetch_one(acc, cache_dir, args.ncbi_email, args.ncbi_api_key)
            text = fasta_path.read_text()
            if text and not text.endswith("\n"):
                text += "\n"
            anello_fasta_texts.append(text)
            anello_gtf_texts.append(_genome_as_transcript_gtf(text, acc))
        except NCBIFetchError as exc:
            anello_errors.append(f"{acc}: {exc}")

    if anello_errors:
        fail_frac = len(anello_errors) / len(anello_accs) if anello_accs else 0.0
        print(
            f"  WARNING: {len(anello_errors)} / {len(anello_accs)} anellovirus accessions "
            f"failed to download ({fail_frac:.0%}). "
            "Skipping failures — re-run to retry (successful downloads are cached).",
            flush=True,
        )
        if fail_frac > 0.5:
            sys.exit(
                f"ERROR: >50% of anellovirus accessions failed "
                f"({len(anello_errors)}/{len(anello_accs)}). "
                "Check NCBI connectivity and re-run (cached downloads will be reused)."
            )

    print(
        f"  {len(anello_fasta_texts)} anellovirus FASTAs fetched, "
        f"{len(anello_gtf_texts)} whole-genome GTFs generated",
        flush=True,
    )

    # ── 5. Seqname coverage check ─────────────────────────────────────────────
    print("Step 5/7  Seqname coverage check …")
    fasta_ids = _fasta_seq_ids(fasta_paths)
    missing = set(accessions) - fasta_ids
    if missing:
        sys.exit(
            f"ERROR: {len(missing)} GTF seqname(s) have no matching FASTA record "
            "(kb ref would silently drop them):\n"
            + "\n".join(f"  {s}" for s in sorted(missing))
        )
    print(f"  OK: all {len(accessions)} GTF seqnames have FASTA records")

    # ── 6. Concatenate: human (gzip) + viral FASTAs → combined.fa ────────────
    print("Step 6/7  Concatenating references …")
    combined_fa = out / "combined.fa"
    combined_gtf = out / "combined.gtf"

    with open(combined_fa, "wb") as fh:
        # Human cDNA first (gzip-encoded from Ensembl)
        with gzip.open(host_fasta_gz, "rb") as gz:
            shutil.copyfileobj(gz, fh)
        # Curated viral FASTAs (plain text from NCBI cache)
        for fp in fasta_paths:
            data = fp.read_bytes()
            fh.write(data)
            if not data.endswith(b"\n"):
                fh.write(b"\n")
        # Anellovirus FASTAs (plain text, fetched in Step 4b)
        for text in anello_fasta_texts:
            fh.write(text.encode())
    print(f"  combined.fa  → {combined_fa}")

    # The Ensembl companion GTF (host_gtf_gz) is *chromosomal* (seqnames 1/2/X) and does
    # NOT match the cDNA FASTA headers (ENST…) — handing that pair to kb ref makes it hang
    # forever at "Splitting genome". Generate a cDNA-level host GTF from the FASTA instead.
    host_cdna_gtf = out / "host" / "host_cdna.gtf"
    host_cdna_gtf.parent.mkdir(parents=True, exist_ok=True)
    n_host_tx = host_cdna_as_gtf(host_fasta_gz, host_cdna_gtf)
    print(f"  host cDNA GTF → {host_cdna_gtf} ({n_host_tx} transcripts)")

    with open(combined_gtf, "wb") as fh:
        # Host cDNA-level GTF first (seqname = ENST, matches the cDNA FASTA)
        with open(host_cdna_gtf, "rb") as host_fh:
            shutil.copyfileobj(host_fh, fh)
        # Bundled viral GTFs (plain text, curated gene_ids preserved)
        for gtf in gtf_files:
            data = gtf.read_bytes()
            fh.write(data)
            if not data.endswith(b"\n"):
                fh.write(b"\n")
        # Anellovirus GTFs (synthesized whole-genome GTFs; gene_ids = {acc}_geneN)
        for gtf_text in anello_gtf_texts:
            encoded = gtf_text.encode()
            fh.write(encoded)
            if not encoded.endswith(b"\n"):
                fh.write(b"\n")
    print(f"  combined.gtf → {combined_gtf}")

    # ── 7. Run kb ref ─────────────────────────────────────────────────────────
    print("Step 7/7  Running kb ref …")
    kb_bin = shutil.which("kb")
    if kb_bin is None:
        sys.exit(
            "ERROR: 'kb' not found on PATH. "
            "Activate the test_viralscan conda env before running this script."
        )

    panel_idx = out / "panel.idx"
    panel_t2g = out / "panel.t2g"
    cdna_fa = out / "cdna.fa"

    cmd = [
        kb_bin, "ref",
        "-i", str(panel_idx),
        "-g", str(panel_t2g),
        "-f1", str(cdna_fa),
        "--overwrite",
        str(combined_fa),
        str(combined_gtf),
    ]
    print(f"  {' '.join(cmd)}")
    subprocess.run(cmd, check=True)  # noqa: S603

    # ── Verify output ─────────────────────────────────────────────────────────
    for path in [panel_idx, panel_t2g, cdna_fa]:
        if not path.exists() or path.stat().st_size == 0:
            sys.exit(f"ERROR: expected output missing or empty: {path}")

    t2g_lines = sum(1 for _ in panel_t2g.open())
    human_lines = sum(1 for ln in panel_t2g.open() if ln.startswith("ENST"))
    viral_lines = t2g_lines - human_lines
    print(f"\nIndex build complete:")
    print(f"  panel.idx : {panel_idx}  ({panel_idx.stat().st_size // (1024 * 1024)} MB)")
    print(f"  panel.t2g : {panel_t2g}  ({t2g_lines:,} total  |  {human_lines:,} human ENST*  |  {viral_lines:,} viral)")
    print(f"  cdna.fa   : {cdna_fa}")

    epstein = sum(1 for ln in panel_t2g.open() if "EPSTEIN" in ln)
    print(f"  Spot check: {epstein} EPSTEIN (EBV) entries in panel.t2g", end="")
    if epstein == 0:
        print("  ← WARNING: expected >0; check EBV GTF/FASTA inclusion")
    else:
        print("  ✓")

    anello_t2g = sum(1 for ln in panel_t2g.open() if "_gene" in ln and not ln.startswith("ENST"))
    print(f"  Spot check: {anello_t2g} anellovirus-style (_geneN) entries in panel.t2g", end="")
    if anello_fasta_texts and anello_t2g == 0:
        print("  ← WARNING: expected >0; check Step 4b anellovirus fetch and GTF append")
    else:
        print("  ✓")

    print(
        "\nNext step: run format probe on one sample to confirm kb count -x BULK output layout.\n"
        "The output h5ad/mtx will contain both ENST* (host) and viral gene_ids —\n"
        "bulk_viral_summarize.py emits a viral summary TSV; the h5ad is the full\n"
        "host+viral matrix for coexpression analysis.\n"
        "See the plan notes in scripts/bulk_viral_summarize.py."
    )


if __name__ == "__main__":
    main()
