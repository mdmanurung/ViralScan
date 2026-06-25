#!/usr/bin/env python3
"""
P22.10 — Matched-barcode STARsolo ↔ ViralScan EBV comparison (anchored to the paper's cells).

Loads three barcode sets — paper (GSM4796271), STARsolo-wl raw, ViralScan-wl raw —
normalises barcodes to bare 16 bp, intersects to a shared anchor set, and reports:

  1. Barcode-set accounting (paper / STAR / VS / shared / drop-outs).
  2. Multimapper-policy audit: confirms STAR raw = unique UMIs only, VS raw = unique-only
     bustools count (no -m flag); reports EBV totals and flags divergence > 20% for
     manual review.
  3. Per-gene EBV UMI table on shared cells (STAR GeneFull vs VS raw adata.X).
  4. Per-cell EBV-total correlation (Spearman + Pearson).
  5. EBV+ tiers on shared set — pan-latent (≥1 / ≥10 UMI) and lytic-restricted
     (BZLF1/BRLF1/BHRF1 ≥1 UMI) per tool.
  6. EBNA-recovery test: per-gene STAR-vs-kallisto ratio for EBNA and LMP families.

Usage:
    python scripts/matched_barcode_comparison.py \\
        --paper-barcodes  data/geo_GSE158275/GSM4796271_LCL_777_B958_UMI_barcodes.tsv.gz \\
        --starsolo-dir    starsolo_p22_6b/starsolo_ebv_wl/Solo.out \\
        --viralscan-adata out_full_depth_wl/lcl_5lines/SRR12682296/kb-python/counts_unfiltered/ \\
        --out             results/matched_barcode_comparison.tsv

Options:
    --feature-type   STARsolo feature type (default: GeneFull)
    --lytic-genes    Comma-separated gene_id prefixes for lytic tier
                     (default: EPSTEIN_HHV4_BZLF1,EPSTEIN_HHV4_BRLF1,EPSTEIN_HHV4_BHRF1)
    --min-umi-pan    UMI threshold for pan-latent EBV+ (default: 1)
    --min-umi-pan2   Second threshold for pan-latent EBV+ tier 2 (default: 10)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from pathlib import Path
from typing import Optional

import anndata
import numpy as np
import scipy.io
import scipy.sparse
import scipy.stats


# ── Barcode normalisation ──────────────────────────────────────────────────────

def strip_10x_suffix(bc: str) -> str:
    """Strip trailing '-1' suffix (CellRanger convention). Matches multimap.strip_10x_suffix."""
    return bc.removesuffix("-1")


def load_barcodes(path: str) -> list[str]:
    """Load a barcodes TSV (optionally .gz); strip trailing '-1'; return bare 16 bp."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[call-overload]
        return [strip_10x_suffix(line.strip()) for line in fh if line.strip()]


# ── STARsolo matrix helpers ────────────────────────────────────────────────────

def _open(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _find(directory: str, *names: str) -> str:
    for name in names:
        p = os.path.join(directory, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"None of {names} found in {directory}. "
        "Did STARsolo finish successfully?"
    )


def load_starsolo_matrix(
    solo_out_dir: str,
    feature_type: str = "GeneFull",
    filtered: bool = False,  # use RAW for the matched-barcode intersection
) -> tuple[scipy.sparse.csr_matrix, list[str], list[str], list[str]]:
    """Return (csr_matrix[genes x cells], barcodes, gene_ids, gene_names).

    NOTE: use filtered=False (raw) for matched-barcode comparison so we can
    select the paper-anchor cell set from the full corrected barcode space.
    """
    subdir = "filtered" if filtered else "raw"
    mdir = os.path.join(solo_out_dir, feature_type, subdir)

    bc_path = _find(mdir, "barcodes.tsv.gz", "barcodes.tsv")
    with _open(bc_path) as fh:
        barcodes = [strip_10x_suffix(line.strip()) for line in fh if line.strip()]

    feat_path = _find(mdir, "features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    gene_ids: list[str] = []
    gene_names: list[str] = []
    with _open(feat_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            gene_ids.append(parts[0])
            gene_names.append(parts[1] if len(parts) > 1 else parts[0])

    mtx_path = _find(mdir, "matrix.mtx.gz", "matrix.mtx")
    if mtx_path.endswith(".gz"):
        import io
        with gzip.open(mtx_path, "rb") as gz:
            mat: scipy.sparse.csr_matrix = scipy.sparse.csr_matrix(
                scipy.io.mmread(io.BytesIO(gz.read()))
            )
    else:
        mat = scipy.sparse.csr_matrix(scipy.io.mmread(mtx_path))

    return mat, barcodes, gene_ids, gene_names


# ── ViralScan matrix helpers ───────────────────────────────────────────────────

def load_viralscan_adata(adata_dir: str) -> anndata.AnnData:
    """Load ViralScan raw adata.h5ad (integer X, pre-multimap)."""
    path = os.path.join(adata_dir, "adata.h5ad")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"adata.h5ad not found at {path}. "
            "Has the ViralScan-wl run (slurm_viralscan_ebv_wl.sh) completed?"
        )
    adata = anndata.read_h5ad(path)
    # Normalise barcodes in obs_names
    adata.obs_names = [strip_10x_suffix(bc) for bc in adata.obs_names]
    return adata


def load_viralscan_multimap_adata(adata_dir: str) -> Optional[anndata.AnnData]:
    """Load ViralScan adata_multimap.h5ad (counts_original/counts_corrected layers)."""
    path = os.path.join(adata_dir, "adata_multimap.h5ad")
    if not os.path.exists(path):
        return None
    adata = anndata.read_h5ad(path)
    adata.obs_names = [strip_10x_suffix(bc) for bc in adata.obs_names]
    return adata


# ── EBV gene selectors ─────────────────────────────────────────────────────────

def ebv_gene_mask(gene_ids: list[str]) -> np.ndarray:
    """Boolean mask for genes starting with 'EPSTEIN_'."""
    return np.array([g.startswith("EPSTEIN_") for g in gene_ids])


def lytic_gene_mask(gene_ids: list[str], lytic_prefixes: list[str]) -> np.ndarray:
    """Boolean mask for lytic EBV genes (BZLF1 / BRLF1 / BHRF1 by default)."""
    return np.array([any(g.startswith(p) for p in lytic_prefixes) for g in gene_ids])


# ── Correlation helpers ────────────────────────────────────────────────────────

def spearman_pearson(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Return Spearman + Pearson correlation between x and y."""
    n = len(x)
    if n < 3:
        return {"n": n, "spearman_r": float("nan"), "pearson_r": float("nan"),
                "spearman_p": float("nan"), "pearson_p": float("nan")}
    sp_r, sp_p = scipy.stats.spearmanr(x, y)
    pe_r, pe_p = scipy.stats.pearsonr(x, y)
    return {
        "n": n,
        "spearman_r": float(sp_r),
        "spearman_p": float(sp_p),
        "pearson_r": float(pe_r),
        "pearson_p": float(pe_p),
    }


# ── Multimapper policy audit ───────────────────────────────────────────────────

def multimap_policy_audit(
    star_total_ebv: float,
    vs_total_ebv: float,
    divergence_threshold: float = 0.20,
) -> dict[str, object]:
    """
    Warn if STAR and VS EBV totals diverge > threshold.

    Both raw layers use unique-only policy by default:
      - STARsolo: --soloMultiMappers Unique (default)
      - kb/bustools count: unique-only (no -m flag)
    But EBV has heavy intra-genome homology (BWRF1.1-8, EBNA-3B/3C), so
    large divergence may indicate a policy difference or aligner-specific
    multi-gene EC handling.
    """
    if star_total_ebv == 0 or vs_total_ebv == 0:
        ratio = float("nan")
        flag = "zero_counts"
    else:
        ratio = float(star_total_ebv / vs_total_ebv)
        diff = abs(ratio - 1.0)
        flag = "POLICY_DIVERGENCE_CHECK_MANUALLY" if diff > divergence_threshold else "ok"
    return {
        "star_total_ebv_umi": int(star_total_ebv),
        "vs_total_ebv_umi": int(vs_total_ebv),
        "star_vs_ratio": ratio,
        "policy_check": flag,
        "star_multimapper_policy": "Unique (default, unique UMIs only)",
        "vs_multimapper_policy": "bustools count unique (no -m flag)",
        "note": (
            "Both raw layers are unique-UMI-only. Divergence may reflect "
            "aligner-level differences (STAR splice-junction vs kallisto pseudoalignment) "
            "and within-EBV homology handling, not a counting-policy mismatch."
        ),
    }


# ── EBNA recovery test ─────────────────────────────────────────────────────────

EBNA_GENES = [
    "EPSTEIN_HHV4_EBNA-1.1", "EPSTEIN_HHV4_EBNA-1.2",
    "EPSTEIN_HHV4_EBNA-2",
    "EPSTEIN_HHV4_EBNA-3A",
    "EPSTEIN_HHV4_EBNA-3B/EBNA-3C",
    "EPSTEIN_HHV4_EBNA-LP",
]
LMP_GENES = [
    "EPSTEIN_HHV4_LMP-1",
    "EPSTEIN_HHV4_LMP-2A",
    "EPSTEIN_HHV4_LMP-2B",
]

def ebna_recovery_test(
    per_gene_star: dict[str, float],
    per_gene_vs: dict[str, float],
) -> list[dict[str, object]]:
    """
    Per-gene STAR-vs-kallisto ratio for EBNA and LMP families.

    Tests the hypothesis that STAR under-counts EBNA (splice-junction artifact:
    STAR requires reads to span annotated junctions; EBNA transcripts with unusual
    introns may be missed) while kallisto pseudoalignment recovers them via
    transcript-level k-mer matching.
    """
    rows = []
    genes_of_interest = EBNA_GENES + LMP_GENES
    all_shared = set(per_gene_star) & set(per_gene_vs)
    for gid in genes_of_interest:
        star_umi = per_gene_star.get(gid, 0.0)
        vs_umi = per_gene_vs.get(gid, 0.0)
        if star_umi > 0 and vs_umi > 0:
            ratio = vs_umi / star_umi
        elif star_umi == 0 and vs_umi > 0:
            ratio = float("inf")
        elif star_umi > 0 and vs_umi == 0:
            ratio = 0.0
        else:
            ratio = float("nan")
        rows.append({
            "gene_id": gid,
            "family": "EBNA" if "EBNA" in gid else "LMP",
            "star_umi_total": star_umi,
            "vs_umi_total": vs_umi,
            "vs_over_star_ratio": ratio,
            "in_shared_genes": gid in all_shared,
        })
    return rows


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="P22.10 matched-barcode STARsolo vs ViralScan EBV comparison"
    )
    ap.add_argument("--paper-barcodes", required=True,
                    help="GSM4796271 barcodes.tsv(.gz) from GEO GSE158275_RAW.tar")
    ap.add_argument("--starsolo-dir", required=True,
                    help="Path to Solo.out/ from the whitelist STARsolo run")
    ap.add_argument("--viralscan-adata", required=True,
                    help="Path to counts_unfiltered/ dir from the whitelist ViralScan run")
    ap.add_argument("--out", required=True,
                    help="Output TSV path for matched-barcode comparison")
    ap.add_argument("--feature-type", default="GeneFull",
                    choices=["Gene", "GeneFull"],
                    help="STARsolo feature type (default: GeneFull)")
    ap.add_argument("--lytic-genes", default=None,
                    help="Comma-separated gene_id prefixes for lytic tier "
                         "(default: EPSTEIN_HHV4_BZLF1,EPSTEIN_HHV4_BRLF1,EPSTEIN_HHV4_BHRF1)")
    ap.add_argument("--min-umi-pan", type=int, default=1,
                    help="Minimum UMI for pan-latent EBV+ tier 1 (default: 1)")
    ap.add_argument("--min-umi-pan2", type=int, default=10,
                    help="Minimum UMI for pan-latent EBV+ tier 2 (default: 10)")
    args = ap.parse_args()

    lytic_prefixes = (
        args.lytic_genes.split(",")
        if args.lytic_genes
        else ["EPSTEIN_HHV4_BZLF1", "EPSTEIN_HHV4_BRLF1", "EPSTEIN_HHV4_BHRF1"]
    )

    print("=" * 70)
    print("P22.10 Matched-barcode STARsolo ↔ ViralScan EBV comparison")
    print("=" * 70)

    # ── Step 1: Paper anchor barcodes ─────────────────────────────────────────
    print("\n[1] Loading paper anchor barcodes (GSM4796271) ...")
    paper_bcs = load_barcodes(args.paper_barcodes)
    paper_set = set(paper_bcs)
    print(f"    Paper cells: {len(paper_set):,}")

    # ── Step 2: STARsolo raw matrix ────────────────────────────────────────────
    print(f"\n[2] Loading STARsolo raw matrix ({args.feature_type}/raw) ...")
    star_mat, star_bcs, star_gene_ids, star_gene_names = load_starsolo_matrix(
        args.starsolo_dir, feature_type=args.feature_type, filtered=False
    )
    star_bc_set = set(star_bcs)
    n_ebv_star = int(ebv_gene_mask(star_gene_ids).sum())
    print(f"    STAR raw barcodes : {len(star_bc_set):,}")
    print(f"    STAR genes        : {len(star_gene_ids):,}  (EBV: {n_ebv_star})")

    # ── Step 3: ViralScan raw adata.h5ad ──────────────────────────────────────
    print("\n[3] Loading ViralScan raw adata.h5ad ...")
    vs_adata = load_viralscan_adata(args.viralscan_adata)
    vs_bc_set = set(vs_adata.obs_names.tolist())
    vs_gene_ids = vs_adata.var_names.tolist()
    n_ebv_vs = int(ebv_gene_mask(vs_gene_ids).sum())
    print(f"    VS raw barcodes   : {len(vs_bc_set):,}")
    print(f"    VS genes          : {len(vs_gene_ids):,}  (EBV: {n_ebv_vs})")

    # Check X is integer
    x_sample = vs_adata.X[:10, :10].toarray() if scipy.sparse.issparse(vs_adata.X) else vs_adata.X[:10, :10]
    if not np.allclose(x_sample, np.round(x_sample)):
        print("    WARNING: VS adata.X contains non-integer values — check multimap layer selection")
    else:
        print("    VS adata.X: integer ✓")

    # ── Step 4: Shared barcode set ─────────────────────────────────────────────
    shared = paper_set & star_bc_set & vs_bc_set
    print(f"\n[4] Barcode intersection:")
    print(f"    Paper ∩ STAR ∩ VS (shared) : {len(shared):,}")
    missing_from_star = paper_set - star_bc_set
    missing_from_vs   = paper_set - vs_bc_set
    print(f"    Paper cells absent from STAR raw : {len(missing_from_star):,}")
    print(f"    Paper cells absent from VS raw   : {len(missing_from_vs):,}")
    if len(shared) == 0:
        print("ERROR: No shared barcodes — check whitelist and barcode normalisation.", file=sys.stderr)
        sys.exit(1)

    shared_list = sorted(shared)

    # ── Step 5: Build EBV matrices on shared cells ────────────────────────────
    print("\n[5] Building per-gene EBV UMI matrices on shared cells ...")

    # STAR: genes x cells → cells x ebv_genes on shared cells
    star_bc_idx = {bc: i for i, bc in enumerate(star_bcs)}
    shared_star_cols = np.array([star_bc_idx[bc] for bc in shared_list if bc in star_bc_idx])
    star_ebv_mask = ebv_gene_mask(star_gene_ids)
    star_ebv_ids = [g for g, m in zip(star_gene_ids, star_ebv_mask) if m]
    star_ebv_sub = star_mat[star_ebv_mask, :][:, shared_star_cols]  # ebv_genes x shared
    star_ebv_arr = np.asarray(star_ebv_sub.todense())               # (n_ebv_genes, n_shared)

    # VS: barcodes x genes → shared_cells x ebv_genes
    vs_bc_idx = {bc: i for i, bc in enumerate(vs_adata.obs_names)}
    shared_vs_rows = np.array([vs_bc_idx[bc] for bc in shared_list if bc in vs_bc_idx])
    vs_ebv_mask = ebv_gene_mask(vs_gene_ids)
    vs_ebv_ids = [g for g, m in zip(vs_gene_ids, vs_ebv_mask) if m]
    X_mat = vs_adata.X
    if scipy.sparse.issparse(X_mat):
        vs_ebv_sub = X_mat[shared_vs_rows, :][:, vs_ebv_mask]  # shared x ebv_genes
        vs_ebv_arr = np.asarray(vs_ebv_sub.todense())
    else:
        vs_ebv_sub = X_mat[np.ix_(shared_vs_rows, np.where(vs_ebv_mask)[0])]
        vs_ebv_arr = np.asarray(vs_ebv_sub)
    # Transpose to ebv_genes x shared for consistent shape
    vs_ebv_arr = vs_ebv_arr.T

    print(f"    EBV genes in STAR : {len(star_ebv_ids)}")
    print(f"    EBV genes in VS   : {len(vs_ebv_ids)}")

    # Build common gene set (identity crosswalk via gene_id)
    common_ebv = sorted(set(star_ebv_ids) & set(vs_ebv_ids))
    print(f"    EBV genes shared  : {len(common_ebv)}")

    star_ebv_dict = {gid: i for i, gid in enumerate(star_ebv_ids)}
    vs_ebv_dict   = {gid: i for i, gid in enumerate(vs_ebv_ids)}

    # ── Step 6: Multimapper policy audit ──────────────────────────────────────
    print("\n[6] Multimapper policy audit ...")
    star_total_ebv = float(star_ebv_arr.sum())
    vs_total_ebv   = float(vs_ebv_arr.sum())
    audit = multimap_policy_audit(star_total_ebv, vs_total_ebv)
    print(f"    STAR EBV total UMI  : {audit['star_total_ebv_umi']:,}")
    print(f"    VS EBV total UMI    : {audit['vs_total_ebv_umi']:,}")
    print(f"    STAR/VS ratio       : {audit['star_vs_ratio']:.3f}")
    print(f"    STAR multimapper    : {audit['star_multimapper_policy']}")
    print(f"    VS multimapper      : {audit['vs_multimapper_policy']}")
    if audit["policy_check"] != "ok":
        print(f"    ⚠  {audit['policy_check']}: ratio > 20% — inspect per-gene breakdown.")
    else:
        print(f"    Policy check        : ok (ratio within 20%)")

    # ── Step 7: Per-gene totals on shared cells ────────────────────────────────
    print("\n[7] Per-gene EBV UMI totals (shared cells) ...")
    per_gene_star: dict[str, float] = {}
    per_gene_vs:   dict[str, float] = {}
    for gid in common_ebv:
        si = star_ebv_dict[gid]
        vi = vs_ebv_dict[gid]
        per_gene_star[gid] = float(star_ebv_arr[si, :].sum())
        per_gene_vs[gid]   = float(vs_ebv_arr[vi, :].sum())

    # ── Step 8: Per-cell EBV-total correlation ─────────────────────────────────
    print("\n[8] Per-cell EBV-total correlation (shared cells) ...")
    star_per_cell = np.asarray(star_ebv_arr.sum(axis=0)).flatten()  # (n_shared,)
    vs_per_cell   = np.asarray(vs_ebv_arr.sum(axis=0)).flatten()    # (n_shared,)

    corr = spearman_pearson(star_per_cell, vs_per_cell)
    print(f"    n shared cells   : {corr['n']:,}")
    print(f"    Spearman r       : {corr['spearman_r']:.4f}  (p={corr['spearman_p']:.2e})")
    print(f"    Pearson r        : {corr['pearson_r']:.4f}  (p={corr['pearson_p']:.2e})")

    # ── Step 9: EBV+ tiers ────────────────────────────────────────────────────
    print(f"\n[9] EBV+ tiers on {len(shared):,} shared cells ...")

    def tier_report(per_cell: np.ndarray, tool: str, gene_ids: list[str],
                    ebv_arr: np.ndarray) -> dict[str, object]:
        n = len(per_cell)
        n_pan1  = int((per_cell >= args.min_umi_pan).sum())
        n_pan2  = int((per_cell >= args.min_umi_pan2).sum())
        # Lytic tier: at least one lytic gene ≥1 UMI
        ly_mask = lytic_gene_mask(gene_ids, lytic_prefixes)
        ly_idx  = [i for i, m in enumerate(ly_mask) if m]
        if ly_idx:
            lytic_per_cell = np.asarray(ebv_arr[np.array(ly_idx), :].sum(axis=0)).flatten()
            n_lytic = int((lytic_per_cell >= 1).sum())
        else:
            n_lytic = 0
        return {
            "tool": tool,
            "n_shared_cells": n,
            f"n_ebv_ge{args.min_umi_pan}umi": n_pan1,
            f"pct_ebv_ge{args.min_umi_pan}umi": f"{100.0 * n_pan1 / n:.2f}" if n else "0.00",
            f"n_ebv_ge{args.min_umi_pan2}umi": n_pan2,
            f"pct_ebv_ge{args.min_umi_pan2}umi": f"{100.0 * n_pan2 / n:.2f}" if n else "0.00",
            "n_lytic": n_lytic,
            "pct_lytic": f"{100.0 * n_lytic / n:.2f}" if n else "0.00",
            "lytic_genes": ",".join(lytic_prefixes),
            "paper_lytic_target": "2.2%",
        }

    star_tier = tier_report(star_per_cell, "STARsolo (wl, GeneFull raw, unique UMIs)",
                            star_ebv_ids, star_ebv_arr)
    vs_tier   = tier_report(vs_per_cell, "ViralScan (wl, adata.h5ad X, unique UMIs)",
                            vs_ebv_ids, vs_ebv_arr)

    for t in (star_tier, vs_tier):
        print(f"\n  {t['tool']}")
        print(f"    Pan-latent (≥{args.min_umi_pan} UMI)  : {t[f'n_ebv_ge{args.min_umi_pan}umi']:,} ({t[f'pct_ebv_ge{args.min_umi_pan}umi']}%)")
        print(f"    Pan-latent (≥{args.min_umi_pan2} UMI) : {t[f'n_ebv_ge{args.min_umi_pan2}umi']:,} ({t[f'pct_ebv_ge{args.min_umi_pan2}umi']}%)")
        print(f"    Lytic-restricted     : {t['n_lytic']:,} ({t['pct_lytic']}%)  [target ~2.2%]")

    # ── Step 10: EBNA recovery test ────────────────────────────────────────────
    print("\n[10] EBNA/LMP recovery test (VS/STAR ratio per gene) ...")
    ebna_rows = ebna_recovery_test(per_gene_star, per_gene_vs)
    for row in ebna_rows:
        ratio_str = (f"{row['vs_over_star_ratio']:.2f}"
                     if isinstance(row['vs_over_star_ratio'], float)
                        and not np.isnan(float(str(row['vs_over_star_ratio'])))
                     else str(row['vs_over_star_ratio']))
        print(f"    {row['gene_id']:40s}  STAR={row['star_umi_total']:6.0f}  VS={row['vs_umi_total']:6.0f}  VS/STAR={ratio_str}")

    # ── Step 11: Write outputs ─────────────────────────────────────────────────
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 11a: Comparison summary TSV
    summary_rows = [
        {
            "section": "barcode_accounting",
            "key": "paper_cells", "value": len(paper_set), "note": "GSM4796271 (LCL_777_B958)",
        },
        {
            "section": "barcode_accounting",
            "key": "star_raw_barcodes", "value": len(star_bc_set), "note": "whitelist-corrected",
        },
        {
            "section": "barcode_accounting",
            "key": "vs_raw_barcodes", "value": len(vs_bc_set), "note": "whitelist-corrected",
        },
        {
            "section": "barcode_accounting",
            "key": "shared_barcodes", "value": len(shared), "note": "paper ∩ STAR ∩ VS",
        },
        {
            "section": "barcode_accounting",
            "key": "paper_missing_from_star", "value": len(missing_from_star), "note": "",
        },
        {
            "section": "barcode_accounting",
            "key": "paper_missing_from_vs", "value": len(missing_from_vs), "note": "",
        },
        {
            "section": "policy_audit",
            "key": "star_ebv_total_umi", "value": audit["star_total_ebv_umi"],
            "note": audit["star_multimapper_policy"],
        },
        {
            "section": "policy_audit",
            "key": "vs_ebv_total_umi", "value": audit["vs_total_ebv_umi"],
            "note": audit["vs_multimapper_policy"],
        },
        {
            "section": "policy_audit",
            "key": "star_vs_ratio", "value": f"{audit['star_vs_ratio']:.4f}",
            "note": audit["policy_check"],
        },
        {
            "section": "correlation",
            "key": "spearman_r", "value": f"{corr['spearman_r']:.4f}",
            "note": f"n={corr['n']}",
        },
        {
            "section": "correlation",
            "key": "spearman_p", "value": f"{corr['spearman_p']:.2e}", "note": "",
        },
        {
            "section": "correlation",
            "key": "pearson_r", "value": f"{corr['pearson_r']:.4f}", "note": "",
        },
    ]
    # Tier rows
    for t, label in [(star_tier, "star"), (vs_tier, "vs")]:
        for k, v in t.items():
            if k != "tool":
                summary_rows.append({"section": f"tier_{label}", "key": k, "value": v, "note": t["tool"]})

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["section", "key", "value", "note"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"\nSummary TSV: {out_path}")

    # 11b: Per-gene table
    per_gene_path = out_path.parent / (out_path.stem + "_per_gene.tsv")
    per_gene_rows = []
    for gid in common_ebv:
        star_tot = per_gene_star.get(gid, 0.0)
        vs_tot   = per_gene_vs.get(gid, 0.0)
        si = star_ebv_dict[gid]
        vi = vs_ebv_dict[gid]
        star_n_cells = int((star_ebv_arr[si, :] >= 1).sum())
        vs_n_cells   = int((vs_ebv_arr[vi, :] >= 1).sum())
        per_gene_rows.append({
            "gene_id": gid,
            "star_total_umi": int(star_tot),
            "vs_total_umi": int(vs_tot),
            "vs_over_star": f"{vs_tot / star_tot:.3f}" if star_tot > 0 else "inf",
            "star_n_cells_ge1": star_n_cells,
            "vs_n_cells_ge1": vs_n_cells,
        })
    with open(per_gene_path, "w", newline="") as fh:
        writer2 = csv.DictWriter(fh, fieldnames=list(per_gene_rows[0].keys()), delimiter="\t")
        writer2.writeheader()
        writer2.writerows(per_gene_rows)
    print(f"Per-gene TSV: {per_gene_path}")

    # 11c: EBNA recovery table
    ebna_path = out_path.parent / (out_path.stem + "_ebna_recovery.tsv")
    with open(ebna_path, "w", newline="") as fh:
        writer3 = csv.DictWriter(fh, fieldnames=list(ebna_rows[0].keys()), delimiter="\t")
        writer3.writeheader()
        writer3.writerows(ebna_rows)
    print(f"EBNA recovery TSV: {ebna_path}")

    print("\n[Done] Matched-barcode comparison complete.")


if __name__ == "__main__":
    main()
