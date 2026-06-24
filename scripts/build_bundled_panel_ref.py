#!/usr/bin/env python3
"""
Build the ViralScan bundled-panel reference index — one-time setup.

Downloads FASTA for each accession in src/viralscan/data/*.gtf from NCBI
(reusing ncbi_fetch for retry + SHA-256-validated caching), concatenates with
the curated bundled GTFs (preserving ADENO_*/EPSTEIN_* gene_ids), then runs
kb ref to produce panel.idx + panel.t2g.

Usage:
    NCBI_EMAIL=you@example.com PYTHONPATH=/path/to/ViralScan/src \\
        python scripts/build_bundled_panel_ref.py \\
        --out /exports/para-lipg-hpc/mdmanurung/viralscan_bulk_gse128078/ref

SLURM (activate test_viralscan conda env before sbatch):
    sbatch --job-name=build_panel_ref --cpus-per-task=2 --mem=8G --time=06:00:00 \\
        --wrap "NCBI_EMAIL=... PYTHONPATH=$PWD/src python $PWD/scripts/build_bundled_panel_ref.py --out $WORKDIR/ref"

Pre-conditions verified at startup:
  1. Bundled GTFs contain exon features (kb ref extracts cDNA from exon rows;
     missing exons yield an empty/partial index that wc -l alone won't catch).
  2. Every GTF seqname has a matching FASTA header after download (a suppressed
     or missing NCBI accession would otherwise let kb ref silently drop viruses).
"""

from __future__ import annotations

import argparse
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
    # Ensure PYTHONPATH includes src/ so viralscan.scripts.ncbi_fetch is importable
    sys.path.insert(0, str(repo / "src"))

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--out", required=True, type=Path,
        help="Output directory for panel.idx, panel.t2g, cdna.fa, combined.*",
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

    # Import after sys.path is set
    from viralscan.scripts.ncbi_fetch import (  # noqa: E402
        DEFAULT_CACHE_DIR, _fetch_one, NCBIFetchError,
    )

    cache_dir: Path = args.cache_dir or DEFAULT_CACHE_DIR

    # ── 1. Discover bundled GTFs ──────────────────────────────────────────────
    gtf_dir = repo / "src" / "viralscan" / "data"
    gtf_files = sorted(gtf_dir.glob("*.gtf"))
    if not gtf_files:
        sys.exit(f"ERROR: no *.gtf files found in {gtf_dir}")
    print(f"Found {len(gtf_files)} bundled GTFs in {gtf_dir}")

    # ── 2. Pre-check: confirm exon features ──────────────────────────────────
    all_seqnames, feature_types = _extract_gtf_seqnames(gtf_files)
    if "exon" not in feature_types:
        sys.exit(
            f"ERROR: bundled GTFs contain no 'exon' features (found: {sorted(feature_types)}). "
            "kb ref extracts cDNA from exon rows — the index would be empty."
        )
    print(f"Pre-check OK: 'exon' present. All feature types: {sorted(feature_types)}")

    # Extract accession-format seqnames (skip any non-accession metadata lines)
    accessions = sorted(s for s in all_seqnames if _ACC_RE.match(s))
    print(f"Extracted {len(accessions)} unique NCBI accessions from GTF seqnames")

    # ── 3. Download FASTAs from NCBI ─────────────────────────────────────────
    fasta_paths: list[Path] = []
    errors: list[str] = []
    for i, acc in enumerate(accessions, 1):
        print(f"[{i:3d}/{len(accessions)}] {acc} … ", end="", flush=True)
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

    # ── 4. Seqname coverage check ─────────────────────────────────────────────
    fasta_ids = _fasta_seq_ids(fasta_paths)
    missing = {s for s in all_seqnames if _ACC_RE.match(s)} - fasta_ids
    if missing:
        sys.exit(
            f"ERROR: {len(missing)} GTF seqname(s) have no matching FASTA record "
            "(kb ref would silently drop them):\n"
            + "\n".join(f"  {s}" for s in sorted(missing))
        )
    print(f"Coverage check OK: all {len(accessions)} GTF seqnames have FASTA records")

    # ── 5. Concatenate FASTAs + bundled GTFs ─────────────────────────────────
    combined_fa = out / "combined.fa"
    combined_gtf = out / "combined.gtf"

    print(f"\nConcatenating {len(fasta_paths)} FASTAs → {combined_fa}")
    with open(combined_fa, "w") as fh:
        for fp in fasta_paths:
            text = fp.read_text()
            fh.write(text)
            if not text.endswith("\n"):
                fh.write("\n")

    print(f"Concatenating {len(gtf_files)} bundled GTFs → {combined_gtf}")
    with open(combined_gtf, "w") as fh:
        for gtf in gtf_files:
            text = gtf.read_text()
            fh.write(text)
            if not text.endswith("\n"):
                fh.write("\n")

    # ── 6. Run kb ref ─────────────────────────────────────────────────────────
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
    print(f"\nRunning: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)  # noqa: S603

    # ── 7. Verify output ──────────────────────────────────────────────────────
    for path in [panel_idx, panel_t2g, cdna_fa]:
        if not path.exists() or path.stat().st_size == 0:
            sys.exit(f"ERROR: expected output missing or empty: {path}")

    t2g_lines = sum(1 for _ in panel_t2g.open())
    print(f"\nIndex build complete:")
    print(f"  panel.idx : {panel_idx}  ({panel_idx.stat().st_size // (1024 * 1024)} MB)")
    print(f"  panel.t2g : {panel_t2g}  ({t2g_lines:,} transcript-gene mappings)")
    print(f"  cdna.fa   : {cdna_fa}")

    epstein = sum(1 for ln in panel_t2g.open() if "EPSTEIN" in ln)
    print(f"  Spot check: {epstein} EPSTEIN (EBV) entries in panel.t2g", end="")
    if epstein == 0:
        print("  ← WARNING: expected >0; check EBV GTF/FASTA inclusion")
    else:
        print("  ✓")

    print(
        "\nNext step: run format probe on one sample to confirm kb count -x BULK output layout,\n"
        "then adjust bulk_viral_summarize.py to match the actual files.\n"
        "See the plan notes in scripts/bulk_viral_summarize.py."
    )


if __name__ == "__main__":
    main()
