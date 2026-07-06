#!/usr/bin/env python3
"""
Compare STARsolo EBV quantification vs ViralScan kallisto output.

Reads the STARsolo GeneFull filtered matrix, counts cells with ≥1 and ≥10
EBV UMIs (features whose gene_id starts with EPSTEIN_), and optionally
cross-tabulates against a ViralScan viral_summary TSV.

Usage:
    python scripts/compare_starsolo_viralscan.py \\
        --starsolo-dir starsolo_p22_6/starsolo_ebv/Solo.out \\
        --out starsolo_p22_6/comparison_starsolo_vs_viralscan.tsv \\
        [--viralscan-summary results/.../viral_summary_EBV.tsv]
"""
import argparse
import csv
import gzip
import os
import sys

import numpy as np
import scipy.io
import scipy.sparse


def _open(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_starsolo_matrix(solo_out_dir, feature_type="GeneFull", filtered=True):
    """Return (csr_matrix[genes x cells], barcodes, gene_ids, gene_names)."""
    subdir = "filtered" if filtered else "raw"
    mdir = os.path.join(solo_out_dir, feature_type, subdir)

    # barcodes
    bc_path = _find(mdir, "barcodes.tsv.gz", "barcodes.tsv")
    with _open(bc_path) as fh:
        barcodes = [line.strip() for line in fh]

    # features: gene_id \t gene_name \t feature_type
    feat_path = _find(mdir, "features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    gene_ids, gene_names = [], []
    with _open(feat_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            gene_ids.append(parts[0])
            gene_names.append(parts[1] if len(parts) > 1 else parts[0])

    # matrix (Market Exchange format, genes x cells)
    mtx_path = _find(mdir, "matrix.mtx.gz", "matrix.mtx")
    if mtx_path.endswith(".gz"):
        import io
        with gzip.open(mtx_path, "rb") as gz:
            mat = scipy.io.mmread(io.BytesIO(gz.read()))
    else:
        mat = scipy.io.mmread(mtx_path)
    mat = scipy.sparse.csr_matrix(mat)

    return mat, barcodes, gene_ids, gene_names


def _find(directory, *names):
    for name in names:
        p = os.path.join(directory, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"None of {names} found in {directory}. "
        "Did STARsolo finish successfully?"
    )


def ebv_per_cell(mat, gene_ids):
    """Sum UMIs for EBV genes (gene_id starts with EPSTEIN_) per cell."""
    mask = np.array([g.startswith("EPSTEIN_") for g in gene_ids])
    if not mask.any():
        print("WARNING: no genes with gene_id starting with 'EPSTEIN_' found in features. "
              "Check that the EBV GTF was concatenated into the combined GTF.",
              file=sys.stderr)
        return np.zeros(mat.shape[1], dtype=np.float64)
    ebv_mat = mat[mask, :]  # genes x cells
    return np.asarray(ebv_mat.sum(axis=0)).flatten()


def parse_viralscan_summary(path):
    """Read a ViralScan viral_summary TSV; return dict with key metrics."""
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if "EBV" in row.get("virus", "") or "Epstein" in row.get("virus", ""):
                result = row
                break
    if not result and os.path.exists(path):
        # try first non-header row regardless of virus column
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                result = row
                break
    return result


def main():
    ap = argparse.ArgumentParser(description="Compare STARsolo vs ViralScan EBV detection")
    ap.add_argument("--starsolo-dir", required=True,
                    help="Path to the Solo.out directory produced by STARsolo")
    ap.add_argument("--viralscan-summary", default=None,
                    help="Optional ViralScan viral_summary TSV for comparison")
    ap.add_argument("--out", required=True,
                    help="Output TSV path for comparison table")
    ap.add_argument("--feature-type", default="GeneFull",
                    choices=["Gene", "GeneFull"],
                    help="STARsolo feature type (default: GeneFull)")
    args = ap.parse_args()

    print(f"Loading STARsolo {args.feature_type}/filtered matrix ...")
    mat, barcodes, gene_ids, gene_names = load_starsolo_matrix(
        args.starsolo_dir, feature_type=args.feature_type
    )
    n_cells = mat.shape[1]
    n_genes = mat.shape[0]
    n_ebv_genes = sum(1 for g in gene_ids if g.startswith("EPSTEIN_"))
    print(f"  Cells: {n_cells:,}   Genes: {n_genes:,}   EBV genes: {n_ebv_genes}")

    upc = ebv_per_cell(mat, gene_ids)

    n1 = int(np.sum(upc >= 1))
    n10 = int(np.sum(upc >= 10))
    pct1 = 100.0 * n1 / n_cells if n_cells else 0.0
    pct10 = 100.0 * n10 / n_cells if n_cells else 0.0

    print(f"\nSTARsolo EBV detection (GeneFull, filtered cells):")
    print(f"  Total cells  : {n_cells:,}")
    print(f"  EBV ≥1 UMI  : {n1:,}  ({pct1:.2f}%)")
    print(f"  EBV ≥10 UMI : {n10:,}  ({pct10:.2f}%)")

    # ViralScan reference numbers from 1M-read dry-run (hard-coded from prior session)
    VS_REFERENCE = {
        "tool": "ViralScan (1M read subsample)",
        "n_cells": "~8,523 (estimated at 1M reads)",
        "n_ebv_1umi": "~285",
        "pct_ebv_1umi": "3.34",
        "n_ebv_10umi": "285",
        "pct_ebv_10umi": "3.34",
        "note": "1M-read subsample; super-expressors (SE) defined as ≥10 UMI",
    }

    vs_parsed = {}
    if args.viralscan_summary:
        vs_parsed = parse_viralscan_summary(args.viralscan_summary)
        print(f"\nViralScan summary from: {args.viralscan_summary}")
        if vs_parsed:
            print(f"  {vs_parsed}")

    rows = [
        {
            "tool": "STARsolo (full depth, GeneFull filtered)",
            "n_cells": n_cells,
            "n_ebv_1umi": n1,
            "pct_ebv_1umi": f"{pct1:.2f}",
            "n_ebv_10umi": n10,
            "pct_ebv_10umi": f"{pct10:.2f}",
            "note": "CellRanger2 knee filter; no 10x whitelist (permissive)",
        },
        VS_REFERENCE,
    ]

    if vs_parsed:
        rows.append({
            "tool": "ViralScan (from summary file)",
            **{k: vs_parsed.get(k, "") for k in
               ["n_cells", "n_ebv_1umi", "pct_ebv_1umi", "n_ebv_10umi", "pct_ebv_10umi"]},
            "note": args.viralscan_summary,
        })

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fieldnames = ["tool", "n_cells", "n_ebv_1umi", "pct_ebv_1umi",
                  "n_ebv_10umi", "pct_ebv_10umi", "note"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nComparison table written to: {args.out}")


if __name__ == "__main__":
    main()
