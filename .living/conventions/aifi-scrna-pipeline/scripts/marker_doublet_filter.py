#!/usr/bin/env python3
"""
Marker-based doublet filtering within cell classes.

Removes clusters expressing markers of the wrong lineage, catching
cross-lineage doublets that scrublet missed.

Usage:
    python marker_doublet_filter.py \
        --input /path/to/clustered.h5ad \
        --output /path/to/filtered.h5ad \
        --l2_key AIFI_L2_prediction \
        --cluster_key leiden_1.5
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# ============================================================
# Default AIFI marker-based doublet thresholds
# Format: {L2_type: [(reason, gene, threshold, metric), ...]}
# metric: "frac" = fraction of cells detected
#          "mean" = mean expression value
# ============================================================
DEFAULT_THRESHOLDS = {
    "ASDC": [
        ("B cell doublet", "MS4A1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
        ("T cell doublet", "CD3E", 0.4, "frac"),
    ],
    "CD14 monocyte": [
        ("B cell doublet", "MS4A1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
        ("T cell doublet", "CD3E", 0.1, "frac"),
        ("Erythrocyte doublet", "HBB", 1.0, "mean"),
    ],
    "CD16 monocyte": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
        ("T cell doublet", "CD3E", 0.2, "frac"),
    ],
    "CD56bright NK cell": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.2, "frac"),
        ("Myeloid doublet", "FCN1", 0.4, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
        ("T cell doublet", "CD3D", 0.4, "frac"),
    ],
    "CD56dim NK cell": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.7, "mean"),
        ("Myeloid doublet", "FCN1", 0.4, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
        ("T cell doublet", "IL7R", 0.4, "frac"),
    ],
    "Memory CD4 T cell": [
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
    ],
    "Memory CD8 T cell": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.7, "mean"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
    ],
    "Naive CD4 T cell": [
        ("B cell doublet", "MS4A1", 0.2, "frac"),
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
    ],
    "Naive CD8 T cell": [
        ("B cell doublet", "MS4A1", 0.2, "frac"),
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
    ],
    "Naive B cell": [
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
        ("T cell doublet", "CD3D", 0.4, "frac"),
    ],
    "Memory B cell": [
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
        ("T cell doublet", "CD3D", 0.4, "frac"),
    ],
    "Effector B cell": [
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
        ("T cell doublet", "CD3D", 0.2, "frac"),
    ],
    "Treg": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
    ],
    "MAIT": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
    ],
    "gdT": [
        ("Erythrocyte doublet", "HBB", 0.2, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.4, "frac"),
    ],
    "pDC": [
        ("B cell doublet", "MS4A1", 0.2, "frac"),
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Myeloid doublet", "FCN1", 0.2, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
        ("T cell doublet", "CD3D", 0.4, "frac"),
    ],
    "Intermediate monocyte": [
        ("B cell doublet", "MS4A1", 0.4, "frac"),
        ("Erythrocyte doublet", "HBB", 0.4, "frac"),
        ("Platelet doublet", "PPBP", 0.2, "frac"),
        ("T cell doublet", "CD3E", 0.2, "frac"),
    ],
}


def get_gene_expression(adata, gene):
    """Extract gene expression as dense array."""
    if gene not in adata.var_names:
        return None
    expr = adata[:, gene].X
    if sparse.issparse(expr):
        return np.asarray(expr.toarray()).flatten()
    return np.asarray(expr).flatten()


def compute_cluster_stats(adata, cluster_key, gene):
    """Compute per-cluster fraction detected and mean expression."""
    expr = get_gene_expression(adata, gene)
    if expr is None:
        return pd.DataFrame()

    df = pd.DataFrame({
        'cluster': adata.obs[cluster_key].values,
        'expr': expr,
        'detected': expr > 0,
    })
    stats = df.groupby('cluster').agg(
        frac=('detected', 'mean'),
        mean=('expr', 'mean'),
        n_cells=('expr', 'count'),
    )
    return stats


def filter_population(adata, population, cluster_key, thresholds):
    """Apply marker-based doublet filtering to one L2 population."""
    audit = []
    clusters_to_remove = {}

    for reason, gene, threshold, metric in thresholds:
        stats = compute_cluster_stats(adata, cluster_key, gene)
        if stats.empty:
            continue

        flagged = stats[stats[metric] > threshold].index.tolist()
        for cl in flagged:
            if cl not in clusters_to_remove:
                clusters_to_remove[cl] = []
            clusters_to_remove[cl].append(reason)

        if flagged:
            n_cells = stats.loc[flagged, 'n_cells'].sum()
            audit.append({
                'population': population,
                'reason': reason,
                'gene': gene,
                'threshold': threshold,
                'metric': metric,
                'n_clusters': len(flagged),
                'n_cells': int(n_cells),
            })

    # Remove flagged clusters
    all_flagged = set(clusters_to_remove.keys())
    mask = ~adata.obs[cluster_key].isin(all_flagged)
    n_removed = (~mask).sum()

    return adata[mask].copy(), audit, n_removed


def run_marker_doublet_filter(adata, l2_key='AIFI_L2_prediction',
                                cluster_key='leiden_1.5',
                                thresholds=None,
                                min_genes=750, max_low_frac=0.3):
    """
    Run marker-based doublet filtering across all L2 populations.

    Also applies low-quality cluster removal.
    """
    if thresholds is None:
        thresholds = DEFAULT_THRESHOLDS

    all_audit = []
    total_removed = 0
    barcodes_to_keep = []

    populations = adata.obs[l2_key].unique()
    for pop in sorted(populations):
        pop_mask = adata.obs[l2_key] == pop
        adata_pop = adata[pop_mask].copy()

        if pop in thresholds:
            adata_pop, audit, n_removed = filter_population(
                adata_pop, pop, cluster_key, thresholds[pop]
            )
            all_audit.extend(audit)
            total_removed += n_removed
            if n_removed > 0:
                print(f"  {pop}: removed {n_removed:,} marker doublet cells")

        # Low-quality cluster removal (skip for Erythrocyte/Platelet)
        if pop not in ['Erythrocyte', 'Platelet']:
            low_gene_frac = adata_pop.obs.groupby(cluster_key).apply(
                lambda x: (x['n_genes_by_counts'] < min_genes).mean()
            )
            low_quality = low_gene_frac[low_gene_frac > max_low_frac].index
            if len(low_quality) > 0:
                n_lq = adata_pop.obs[cluster_key].isin(low_quality).sum()
                adata_pop = adata_pop[
                    ~adata_pop.obs[cluster_key].isin(low_quality)
                ].copy()
                total_removed += n_lq
                all_audit.append({
                    'population': pop,
                    'reason': 'Low gene detection',
                    'gene': 'n_genes',
                    'threshold': min_genes,
                    'metric': f'>{max_low_frac:.0%} below threshold',
                    'n_clusters': len(low_quality),
                    'n_cells': int(n_lq),
                })
                print(f"  {pop}: removed {n_lq:,} low-quality cells")

        barcodes_to_keep.extend(adata_pop.obs.index.tolist())

    # Filter original adata
    adata_filtered = adata[adata.obs.index.isin(barcodes_to_keep)].copy()

    print(f"\n  Total removed: {total_removed:,}")
    print(f"  Cells retained: {adata_filtered.n_obs:,}")

    return adata_filtered, pd.DataFrame(all_audit)


def main():
    parser = argparse.ArgumentParser(
        description='Marker-based doublet filtering'
    )
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--l2_key', default='AIFI_L2_prediction')
    parser.add_argument('--cluster_key', default='leiden_1.5')
    parser.add_argument('--thresholds_json', default=None,
                        help='Custom thresholds JSON file')
    parser.add_argument('--min_genes', type=int, default=750)
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    thresholds = DEFAULT_THRESHOLDS
    if args.thresholds_json:
        with open(args.thresholds_json) as f:
            thresholds = json.load(f)

    print("Running marker-based doublet filtering...")
    adata_filtered, audit_df = run_marker_doublet_filter(
        adata,
        l2_key=args.l2_key,
        cluster_key=args.cluster_key,
        thresholds=thresholds,
        min_genes=args.min_genes,
    )

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    adata_filtered.write_h5ad(args.output)

    audit_path = args.output.replace('.h5ad', '_doublet_audit.csv')
    audit_df.to_csv(audit_path, index=False)
    print(f"Saved to {args.output}")
    print(f"Audit trail: {audit_path}")


if __name__ == '__main__':
    main()
