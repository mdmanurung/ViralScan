#!/usr/bin/env python3
"""
AIFI-standard serial QC filtering pipeline.

Applies QC filters in sequence: scrublet → MT% → gene count thresholds,
tracking cell removal at each step for a complete audit trail.

Usage:
    python aifi_qc_pipeline.py \
        --input /path/to/labeled_data.h5ad \
        --output /path/to/qc_filtered.h5ad \
        --max_mt_pct 10 \
        --min_genes 200 \
        --max_genes 5000
"""

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce


def run_scrublet_per_sample(adata, sample_key='sample_id',
                             expected_doublet_rate=0.06, random_state=42):
    """Run scrublet independently on each sample."""
    print("Running scrublet per sample...")

    all_scores = pd.Series(dtype=float, name='doublet_score')
    all_calls = pd.Series(dtype=bool, name='predicted_doublet')

    samples = adata.obs[sample_key].unique()
    for i, sample in enumerate(samples):
        mask = adata.obs[sample_key] == sample
        adata_sample = adata[mask].copy()

        try:
            sce.pp.scrublet(
                adata_sample,
                expected_doublet_rate=expected_doublet_rate,
                random_state=random_state
            )
            all_scores = pd.concat([
                all_scores,
                adata_sample.obs['doublet_score']
            ])
            all_calls = pd.concat([
                all_calls,
                adata_sample.obs['predicted_doublet']
            ])
        except Exception as e:
            print(f"  Warning: scrublet failed for {sample}: {e}")
            # Mark all cells as non-doublet if scrublet fails
            all_scores = pd.concat([
                all_scores,
                pd.Series(0.0, index=adata_sample.obs.index, name='doublet_score')
            ])
            all_calls = pd.concat([
                all_calls,
                pd.Series(False, index=adata_sample.obs.index, name='predicted_doublet')
            ])

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(samples)} samples")

    adata.obs['doublet_score'] = all_scores.loc[adata.obs.index]
    adata.obs['predicted_doublet'] = all_calls.loc[adata.obs.index]

    return adata


def serial_qc_filter(adata, max_mt_pct=10.0, min_genes=200, max_genes=5000,
                      skip_scrublet=False, sample_key='sample_id'):
    """
    AIFI-standard serial QC filtering with complete audit trail.

    Filters are applied sequentially in this order:
    1. Scrublet doublet calls
    2. Mitochondrial UMI percentage
    3. Minimum genes detected
    4. Maximum genes detected (likely doublets)
    """
    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, percent_top=[20]
    )

    audit = []
    initial_cells = adata.n_obs
    audit.append({
        'step': 'Pre-QC filtering',
        'cells_removed': 0,
        'cells_remaining': adata.n_obs,
        'pct_of_total': 100.0
    })

    # Step 1: Scrublet
    if not skip_scrublet:
        if 'predicted_doublet' not in adata.obs.columns:
            adata = run_scrublet_per_sample(adata, sample_key=sample_key)

        n_before = adata.n_obs
        adata = adata[~adata.obs['predicted_doublet']].copy()
        removed = n_before - adata.n_obs
        audit.append({
            'step': 'Scrublet doublet call',
            'cells_removed': removed,
            'cells_remaining': adata.n_obs,
            'pct_of_total': removed / initial_cells * 100
        })
        print(f"  Scrublet: removed {removed:,} ({removed/initial_cells:.2%})")

    # Step 2: MT%
    n_before = adata.n_obs
    adata = adata[adata.obs['pct_counts_mt'] <= max_mt_pct].copy()
    removed = n_before - adata.n_obs
    audit.append({
        'step': f'MT UMIs > {max_mt_pct}%',
        'cells_removed': removed,
        'cells_remaining': adata.n_obs,
        'pct_of_total': removed / initial_cells * 100
    })
    print(f"  MT% > {max_mt_pct}: removed {removed:,} ({removed/initial_cells:.2%})")

    # Step 3: Min genes
    n_before = adata.n_obs
    adata = adata[adata.obs['n_genes_by_counts'] >= min_genes].copy()
    removed = n_before - adata.n_obs
    audit.append({
        'step': f'N Genes < {min_genes}',
        'cells_removed': removed,
        'cells_remaining': adata.n_obs,
        'pct_of_total': removed / initial_cells * 100
    })
    print(f"  Genes < {min_genes}: removed {removed:,} ({removed/initial_cells:.2%})")

    # Step 4: Max genes
    n_before = adata.n_obs
    adata = adata[adata.obs['n_genes_by_counts'] <= max_genes].copy()
    removed = n_before - adata.n_obs
    audit.append({
        'step': f'N Genes > {max_genes}',
        'cells_removed': removed,
        'cells_remaining': adata.n_obs,
        'pct_of_total': removed / initial_cells * 100
    })
    print(f"  Genes > {max_genes}: removed {removed:,} ({removed/initial_cells:.2%})")

    # Summary
    total_removed = initial_cells - adata.n_obs
    print(f"\n  Total removed: {total_removed:,} ({total_removed/initial_cells:.2%})")
    print(f"  Cells retained: {adata.n_obs:,} ({adata.n_obs/initial_cells:.2%})")

    return adata, pd.DataFrame(audit)


def main():
    parser = argparse.ArgumentParser(
        description='AIFI-standard serial QC filtering'
    )
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output', required=True, help='Output h5ad file')
    parser.add_argument('--max_mt_pct', type=float, default=10.0)
    parser.add_argument('--min_genes', type=int, default=200)
    parser.add_argument('--max_genes', type=int, default=5000)
    parser.add_argument('--sample_key', default='sample_id')
    parser.add_argument('--skip_scrublet', action='store_true')
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    adata, audit_df = serial_qc_filter(
        adata,
        max_mt_pct=args.max_mt_pct,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        skip_scrublet=args.skip_scrublet,
        sample_key=args.sample_key,
    )

    # Save filtered data
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"\nSaved to {args.output}")

    # Save audit trail
    audit_path = args.output.replace('.h5ad', '_qc_audit.csv')
    audit_df.to_csv(audit_path, index=False)
    print(f"Audit trail saved to {audit_path}")


if __name__ == '__main__':
    main()
