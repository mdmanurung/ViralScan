#!/usr/bin/env python3
"""
AIFI-standard scanpy clustering pipeline.

Normalize → HVG → Scale → PCA → Harmony → Neighbors → UMAP → Leiden

Usage:
    python aifi_cluster_pipeline.py \
        --input /path/to/qc_filtered.h5ad \
        --output /path/to/clustered.h5ad \
        --harmony_key cohort \
        --resolutions 1.0 1.5 2.0
"""

import argparse
import os

import numpy as np
import scanpy as sc
import scanpy.external as sce


def aifi_cluster_pipeline(adata, harmony_key='cohort',
                           n_neighbors=50, n_pcs=30,
                           resolutions=None,
                           exclude_gene_prefixes=None,
                           layer='counts'):
    """
    AIFI-standard scanpy processing pipeline.

    Parameters
    ----------
    adata : AnnData
        Input data. Raw counts should be in the specified layer.
    harmony_key : str
        obs column for Harmony integration (e.g., 'cohort', 'batch_id')
    n_neighbors : int
        Number of neighbors (AIFI default: 50)
    n_pcs : int
        PCs for neighbor computation (AIFI default: 30)
    resolutions : list of float
        Leiden clustering resolutions
    exclude_gene_prefixes : list of str or None
        Gene prefixes to exclude from HVG (e.g., ['IGH', 'IGL', 'IGK'] for B cells)
    layer : str
        Layer containing raw counts
    """
    if resolutions is None:
        resolutions = [1.0, 1.5, 2.0]

    # 1. Start from raw counts
    if layer and layer in adata.layers:
        adata.X = adata.layers[layer].copy()
    print(f"  Starting with {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # 2. Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 3. HVG selection
    sc.pp.highly_variable_genes(adata)
    n_hvg = adata.var['highly_variable'].sum()

    # Optional: exclude gene families
    if exclude_gene_prefixes:
        for prefix in exclude_gene_prefixes:
            mask = adata.var_names.str.startswith(prefix)
            n_excluded = (mask & adata.var['highly_variable']).sum()
            adata.var.loc[mask, 'highly_variable'] = False
            if n_excluded > 0:
                print(f"  Excluded {n_excluded} {prefix}* genes from HVGs")
        n_hvg_after = adata.var['highly_variable'].sum()
        print(f"  HVGs: {n_hvg} → {n_hvg_after} (after exclusions)")
    else:
        print(f"  HVGs: {n_hvg}")

    # 4. Scale HVG data
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg)

    # 5. PCA
    sc.tl.pca(adata_hvg, svd_solver='arpack')
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

    # Transfer PCs to full object for later use
    n_pcs_computed = adata_hvg.varm['PCs'].shape[1]
    adata.varm['PCs'] = np.zeros((adata.n_vars, n_pcs_computed))
    hvg_idx = np.where(adata.var['highly_variable'])[0]
    adata.varm['PCs'][hvg_idx] = adata_hvg.varm['PCs']
    print(f"  PCA: {n_pcs_computed} components computed")

    # 6. Harmony integration
    if harmony_key and harmony_key in adata.obs.columns:
        n_groups = adata.obs[harmony_key].nunique()
        if n_groups > 1:
            sce.pp.harmony_integrate(adata, key=harmony_key)
            use_rep = 'X_pca_harmony'
            print(f"  Harmony: integrated across {n_groups} groups "
                  f"(key: {harmony_key})")
        else:
            use_rep = 'X_pca'
            print(f"  Harmony: skipped (only 1 group for {harmony_key})")
    else:
        use_rep = 'X_pca'
        print(f"  Harmony: skipped (key '{harmony_key}' not found)")

    # 7. Neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors,
                    use_rep=use_rep, n_pcs=n_pcs)
    print(f"  Neighbors: n={n_neighbors}, n_pcs={n_pcs}, rep={use_rep}")

    # 8. UMAP
    sc.tl.umap(adata)

    # 9. Leiden at multiple resolutions
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
        n_clusters = adata.obs[f'leiden_{res}'].nunique()
        print(f"  Leiden (res={res}): {n_clusters} clusters")

    # 10. Compute marker genes for the middle resolution
    mid_res = resolutions[len(resolutions) // 2]
    mid_key = f'leiden_{mid_res}'
    sc.tl.rank_genes_groups(adata, groupby=mid_key, method='wilcoxon')
    print(f"  Marker genes computed for {mid_key}")

    return adata


def main():
    parser = argparse.ArgumentParser(
        description='AIFI-standard clustering pipeline'
    )
    parser.add_argument('--input', required=True, help='Input h5ad')
    parser.add_argument('--output', required=True, help='Output h5ad')
    parser.add_argument('--harmony_key', default='cohort',
                        help='Metadata column for Harmony (default: cohort)')
    parser.add_argument('--n_neighbors', type=int, default=50)
    parser.add_argument('--n_pcs', type=int, default=30)
    parser.add_argument('--resolutions', type=float, nargs='+',
                        default=[1.0, 1.5, 2.0])
    parser.add_argument('--exclude_prefixes', nargs='*', default=None,
                        help='Gene prefixes to exclude from HVGs')
    parser.add_argument('--layer', default='counts',
                        help='Layer with raw counts')
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    print(f"Running AIFI clustering pipeline...")
    adata = aifi_cluster_pipeline(
        adata,
        harmony_key=args.harmony_key,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        resolutions=args.resolutions,
        exclude_gene_prefixes=args.exclude_prefixes,
        layer=args.layer,
    )

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"\nSaved to {args.output}")


if __name__ == '__main__':
    main()
