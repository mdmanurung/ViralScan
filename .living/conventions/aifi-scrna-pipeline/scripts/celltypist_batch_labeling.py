#!/usr/bin/env python3
"""
Batch CellTypist labeling at L1/L2/L3 resolution.

Applies AIFI CellTypist models to all samples in a directory,
producing label and score columns for each resolution level.

Usage:
    python celltypist_batch_labeling.py \
        --input_dir /path/to/h5_files \
        --output_dir /path/to/labeled_output \
        --model_l1 /path/to/AIFI_L1_model.pkl \
        --model_l2 /path/to/AIFI_L2_model.pkl \
        --model_l3 /path/to/AIFI_L3_model.pkl
"""

import argparse
import glob
import os
import sys

import celltypist
import numpy as np
import pandas as pd
import scanpy as sc


def label_single_sample(adata, model_files):
    """Apply all CellTypist models to a single AnnData object."""
    for model_name, model_path in model_files.items():
        label_col = f'{model_name}_prediction'
        score_col = f'{model_name}_score'

        predictions = celltypist.annotate(adata, model=model_path)
        labels = predictions.predicted_labels

        adata.obs[label_col] = labels['predicted_labels'].values

        # Extract prediction score for the assigned label
        prob = predictions.probability_matrix
        scores = []
        for i, idx in enumerate(adata.obs.index):
            predicted_type = adata.obs[label_col].iloc[i]
            if predicted_type in prob.columns:
                scores.append(prob.loc[idx, predicted_type])
            else:
                scores.append(0.0)
        adata.obs[score_col] = scores

    return adata


def process_sample(h5_path, model_files, output_dir, save_h5ad=True):
    """Load, normalize, label, and save a single sample."""
    sample_name = os.path.basename(h5_path).replace('.h5', '').replace('.h5ad', '')
    print(f"  Processing {sample_name}...")

    # Load data
    if h5_path.endswith('.h5ad'):
        adata = sc.read_h5ad(h5_path)
    else:
        adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    # Store raw counts before normalization
    adata.layers['counts'] = adata.X.copy()

    # Normalize for CellTypist
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Label
    adata = label_single_sample(adata, model_files)

    # Save outputs
    os.makedirs(output_dir, exist_ok=True)

    # Save labels as CSV (lightweight)
    label_cols = [c for c in adata.obs.columns
                  if 'prediction' in c or 'score' in c]
    labels_df = adata.obs[label_cols]
    labels_df.to_csv(os.path.join(output_dir, f"{sample_name}_labels.csv"))

    # Optionally save full h5ad
    if save_h5ad:
        adata.write_h5ad(os.path.join(output_dir, f"{sample_name}_labeled.h5ad"))

    n_cells = adata.n_obs
    return sample_name, n_cells


def main():
    parser = argparse.ArgumentParser(
        description='Batch CellTypist labeling at L1/L2/L3 resolution'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing .h5 or .h5ad files')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for labeled data')
    parser.add_argument('--model_l1', required=True,
                        help='Path to AIFI L1 CellTypist model')
    parser.add_argument('--model_l2', required=True,
                        help='Path to AIFI L2 CellTypist model')
    parser.add_argument('--model_l3', default=None,
                        help='Path to AIFI L3 CellTypist model (optional)')
    parser.add_argument('--no_h5ad', action='store_true',
                        help='Skip saving full h5ad files (labels CSV only)')
    args = parser.parse_args()

    # Collect model files
    model_files = {
        'AIFI_L1': args.model_l1,
        'AIFI_L2': args.model_l2,
    }
    if args.model_l3:
        model_files['AIFI_L3'] = args.model_l3

    # Validate model files exist
    for name, path in model_files.items():
        if not os.path.exists(path):
            print(f"ERROR: Model file not found: {path}")
            sys.exit(1)

    # Find input files
    h5_files = sorted(
        glob.glob(os.path.join(args.input_dir, '*.h5')) +
        glob.glob(os.path.join(args.input_dir, '*.h5ad'))
    )

    if not h5_files:
        print(f"No .h5 or .h5ad files found in {args.input_dir}")
        sys.exit(1)

    print(f"Found {len(h5_files)} files to process")
    print(f"Models: {', '.join(model_files.keys())}")

    # Process each sample
    summary = []
    for h5_file in h5_files:
        try:
            name, n_cells = process_sample(
                h5_file, model_files, args.output_dir,
                save_h5ad=not args.no_h5ad
            )
            summary.append({'sample': name, 'n_cells': n_cells, 'status': 'OK'})
        except Exception as e:
            print(f"  ERROR: {e}")
            summary.append({'sample': os.path.basename(h5_file),
                           'n_cells': 0, 'status': str(e)})

    # Print summary
    summary_df = pd.DataFrame(summary)
    print(f"\n=== Labeling Summary ===")
    print(f"Processed: {(summary_df['status'] == 'OK').sum()}/{len(summary_df)}")
    print(f"Total cells: {summary_df['n_cells'].sum():,}")
    summary_df.to_csv(os.path.join(args.output_dir, 'labeling_summary.csv'),
                      index=False)


if __name__ == '__main__':
    main()
