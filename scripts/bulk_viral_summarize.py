#!/usr/bin/env python3
"""
Aggregate bulk kb count outputs for GSE128078 into a per-virus summary TSV.

Scans $OUTDIR/<SRR>/ for kb count results, sums counts per gene_id, maps
gene_ids to virus names via viralscan.virus_grouping, normalizes by
n_pseudoaligned, and emits bulk_viral_summary.tsv.

IMPORTANT — run the format probe first:
    ls -R $OUTDIR/<any_SRR>/counts_unfiltered/
    confirm whether adata.h5ad or cells_x_genes.mtx (+ .genes.txt) exists
    before running this script on the full cohort.

Usage:
    PYTHONPATH=/path/to/ViralScan/src python scripts/bulk_viral_summarize.py \\
        --out-dir /exports/para-lipg-hpc/mdmanurung/viralscan_bulk_gse128078/out \\
        --summary /exports/para-lipg-hpc/mdmanurung/viralscan_bulk_gse128078/bulk_viral_summary.tsv

Scientific caveats (printed to stderr and written to the TSV header):
  - Viral-only index (no host decoy): human k-mers can match viral sequences,
    inflating absolute counts. Treat top hits with suspicion, especially those
    concentrated in low-complexity regions.
  - Counts are read-level pseudoalignment (no UMI deduplication).
  - Cross-virus multimapping is handled by kallisto's EC mechanism, not
    ViralScan's single-cell EM correction.
  - Whole-blood bulk RNA-seq has very low viral fractions; most samples will be
    near-zero. Use RPM (reads per million pseudoaligned) for comparison.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def _load_counts(sample_dir: Path) -> tuple[list[str], list[float]]:
    """Return (gene_ids, counts) from a kb count output directory.

    Tries h5ad first (adata.var_names + X matrix), falls back to
    cells_x_genes.mtx + cells_x_genes.genes.txt.  Both are produced by
    kb count -x BULK; which exists depends on the kb-python version and flags.

    Verify the actual output layout with `ls -R <sample>/counts_unfiltered/`
    before trusting this function on a new kb version.
    """
    counts_dir = sample_dir / "counts_unfiltered"
    if not counts_dir.is_dir():
        raise FileNotFoundError(f"counts_unfiltered/ not found in {sample_dir}")

    # Try h5ad (produced by kb count --h5ad)
    h5ad = counts_dir / "adata.h5ad"
    if h5ad.exists():
        try:
            import anndata  # type: ignore[import]
            import numpy as np
            adata = anndata.read_h5ad(h5ad)
            gene_ids = list(adata.var_names)
            # BULK mode collapses cells to 1 row; sum across that single row
            mat = adata.X
            if hasattr(mat, "toarray"):
                arr = mat.toarray()
            else:
                arr = np.asarray(mat)
            counts = list(arr.sum(axis=0))
            return gene_ids, counts
        except Exception as exc:  # noqa: BLE001
            print(f"  WARNING: h5ad load failed ({exc}); falling back to mtx", file=sys.stderr)

    # Fallback: cells_x_genes.mtx + genes.txt (format-stable across kb versions)
    mtx_path = counts_dir / "cells_x_genes.mtx"
    genes_path = counts_dir / "cells_x_genes.genes.txt"
    if mtx_path.exists() and genes_path.exists():
        from scipy.io import mmread  # type: ignore[import]
        import numpy as np
        gene_ids = [ln.strip().split("\t")[0] for ln in genes_path.read_text().splitlines() if ln.strip()]
        mat = mmread(mtx_path)
        arr = np.asarray(mat.toarray())
        # cells_x_genes.mtx is barcodes × genes (same orientation as anndata.X)
        counts = list(arr.sum(axis=0))
        return gene_ids, counts

    raise FileNotFoundError(
        f"No count files found in {counts_dir}.\n"
        "Expected adata.h5ad or cells_x_genes.mtx + cells_x_genes.genes.txt.\n"
        "Run: ls -R <sample>/counts_unfiltered/ to inspect the actual layout."
    )


def _load_run_info(sample_dir: Path) -> dict:
    run_info = sample_dir / "run_info.json"
    if not run_info.exists():
        return {}
    with open(run_info) as fh:
        return json.load(fh)


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--out-dir", required=True, type=Path,
        help="Directory containing <SRR>/ subdirs (output of bulk_viral_scan.sh)",
    )
    p.add_argument(
        "--summary", required=True, type=Path,
        help="Output TSV path for the per-sample per-virus summary",
    )
    p.add_argument(
        "--pythonpath", type=str, default=None,
        help="Prepend to sys.path for viralscan import (alternative to PYTHONPATH env var)",
    )
    args = p.parse_args()

    if args.pythonpath:
        sys.path.insert(0, args.pythonpath)

    try:
        from viralscan.virus_grouping import virus_name_for_gene  # noqa: E402
    except ImportError:
        sys.exit(
            "ERROR: cannot import viralscan. Set PYTHONPATH=<repo>/src or pass --pythonpath."
        )

    out_dir = args.out_dir
    if not out_dir.is_dir():
        sys.exit(f"ERROR: --out-dir {out_dir} does not exist")

    # Discover sample directories (each is a SRR accession subdir)
    sample_dirs = sorted(d for d in out_dir.iterdir() if d.is_dir())
    if not sample_dirs:
        sys.exit(f"ERROR: no subdirectories found in {out_dir}")
    print(f"Found {len(sample_dirs)} sample directories in {out_dir}", file=sys.stderr)

    caveats = [
        "# CAVEATS",
        "# viral-only index (no host decoy): absolute counts suspect due to host k-mer leakage",
        "# read-level counts only (no UMI deduplication)",
        "# cross-virus multimapping uncorrected (no single-cell EM)",
        "# RPM = reads per million pseudoaligned; use for cross-sample comparison",
        "#",
    ]

    # First pass: collect all virus names across all samples to build unified columns
    all_virus_names: set[str] = set()
    sample_data: list[dict] = []

    for sd in sample_dirs:
        srr = sd.name
        print(f"  Loading {srr} … ", end="", file=sys.stderr)
        try:
            gene_ids, counts = _load_counts(sd)
        except FileNotFoundError as exc:
            print(f"SKIP ({exc})", file=sys.stderr)
            continue

        run_info = _load_run_info(sd)
        n_pseudo = int(run_info.get("n_pseudoaligned", 0) or 0)

        # Aggregate counts by virus name
        virus_counts: dict[str, float] = {}
        for gene_id, count in zip(gene_ids, counts):
            virus = virus_name_for_gene(gene_id)
            virus_counts[virus] = virus_counts.get(virus, 0.0) + float(count)
        all_virus_names.update(virus_counts)

        total_viral = sum(virus_counts.values())
        sample_data.append({
            "srr": srr,
            "n_pseudoaligned": n_pseudo,
            "total_viral_reads": total_viral,
            "virus_counts": virus_counts,
        })
        print(
            f"OK (n_pseudo={n_pseudo:,}, total_viral={int(total_viral)})",
            file=sys.stderr,
        )

    if not sample_data:
        sys.exit("ERROR: no samples loaded successfully")

    # Sort virus columns alphabetically for reproducibility
    virus_cols = sorted(all_virus_names)

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with open(args.summary, "w", newline="") as fh:
        for line in caveats:
            fh.write(line + "\n")
        writer = csv.writer(fh, delimiter="\t")
        header = ["srr", "n_pseudoaligned", "total_viral_reads", "total_viral_rpm"] + virus_cols
        writer.writerow(header)
        for sd in sample_data:
            n_pseudo = sd["n_pseudoaligned"]
            total = sd["total_viral_reads"]
            rpm_scale = (1_000_000 / n_pseudo) if n_pseudo > 0 else 0.0
            total_rpm = round(total * rpm_scale, 4)
            row = [
                sd["srr"],
                n_pseudo,
                int(total),
                total_rpm,
            ] + [round(sd["virus_counts"].get(v, 0.0) * rpm_scale, 4) for v in virus_cols]
            writer.writerow(row)

    print(f"\nWrote {len(sample_data)} samples × {len(virus_cols)} viruses → {args.summary}", file=sys.stderr)
    print("Sanity checks:", file=sys.stderr)
    for sd in sample_data[:3]:
        n = sd["n_pseudoaligned"]
        t = sd["total_viral_reads"]
        if n > 0 and t > n:
            print(f"  WARNING: {sd['srr']}: total_viral ({int(t)}) > n_pseudoaligned ({n}) — check index", file=sys.stderr)
        else:
            top = sorted(sd["virus_counts"].items(), key=lambda x: x[1], reverse=True)[:3]
            top_str = ", ".join(f"{v}:{int(c)}" for v, c in top) or "(none)"
            print(f"  {sd['srr']}: n_pseudo={n:,}  top_hits=[{top_str}]", file=sys.stderr)


if __name__ == "__main__":
    main()
