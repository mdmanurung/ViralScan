#!/usr/bin/env python3
"""
Pseudobulk aggregation and DESeq2 differential expression via rpy2.

Usage:
    python pseudobulk_deseq2.py \
        --input /path/to/annotated.h5ad \
        --output_dir /path/to/de_results \
        --label_key AIFI_L3 \
        --sample_key sample_id \
        --design "~ Age + CMV + Sex" \
        --contrasts "Age:OlderAdult:YoungAdult,CMV:Positive:Negative,Sex:Male:Female"
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def create_pseudobulk(adata, sample_key='sample_id', label_key='AIFI_L3',
                       layer='counts', min_cells=10):
    """
    Aggregate single-cell counts into pseudobulk (summed UMI counts).

    Returns dict of {cell_type: (count_matrix, metadata)} where
    count_matrix is genes × samples and metadata is sample-level.
    """
    pseudobulk = {}

    for cell_type in sorted(adata.obs[label_key].unique()):
        ct_mask = adata.obs[label_key] == cell_type
        adata_ct = adata[ct_mask]

        samples = adata_ct.obs[sample_key].unique()
        count_dict = {}
        meta_rows = []

        for sample in samples:
            s_mask = adata_ct.obs[sample_key] == sample
            n_cells = s_mask.sum()

            if n_cells < min_cells:
                continue

            X = adata_ct[s_mask].layers[layer] if layer else adata_ct[s_mask].X
            if sparse.issparse(X):
                counts_sum = np.asarray(X.sum(axis=0)).flatten().astype(int)
            else:
                counts_sum = np.asarray(X.sum(axis=0)).flatten().astype(int)

            count_dict[sample] = counts_sum

            # First cell's metadata as representative
            meta = adata_ct.obs[s_mask].iloc[0].to_dict()
            meta['n_cells_pseudobulk'] = int(n_cells)
            meta_rows.append(meta)

        if len(count_dict) < 3:
            continue

        count_matrix = pd.DataFrame(
            count_dict, index=adata.var_names
        )
        metadata = pd.DataFrame(meta_rows, index=list(count_dict.keys()))

        pseudobulk[cell_type] = (count_matrix, metadata)
        print(f"  {cell_type}: {len(count_dict)} samples "
              f"(median {metadata['n_cells_pseudobulk'].median():.0f} cells)")

    return pseudobulk


def filter_genes(count_matrix, adata_sc, cell_type_mask, min_pct=0.10):
    """Filter genes based on detection in single cells (AIFI: ≥10%)."""
    X_sc = adata_sc[cell_type_mask].X
    if sparse.issparse(X_sc):
        pct_detected = np.asarray((X_sc > 0).mean(axis=0)).flatten()
    else:
        pct_detected = (X_sc > 0).mean(axis=0)

    gene_pcts = pd.Series(pct_detected, index=adata_sc.var_names)
    genes_pass = gene_pcts[gene_pcts >= min_pct].index
    genes_in_matrix = count_matrix.index.intersection(genes_pass)

    return count_matrix.loc[genes_in_matrix]


def run_deseq2(count_matrix, metadata, design_formula, contrasts, alpha=0.05):
    """
    Run DESeq2 via rpy2.

    Parameters
    ----------
    count_matrix : DataFrame, genes × samples
    metadata : DataFrame, samples × factors
    design_formula : str, R formula (e.g., '~ Age + CMV + Sex')
    contrasts : list of (factor, numerator, denominator)
    alpha : float

    Returns
    -------
    dict of {contrast_name: DataFrame}
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, Formula
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter
    except ImportError:
        print("ERROR: rpy2 not installed. Install with: pip install rpy2")
        print("Also requires R with DESeq2 installed.")
        return {}

    pandas2ri.activate()
    deseq2 = importr('DESeq2')
    base = importr('base')

    # Ensure factors are properly typed
    for col in metadata.select_dtypes(include='object').columns:
        metadata[col] = pd.Categorical(metadata[col])

    # Convert to R objects
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_counts = ro.conversion.py2rpy(count_matrix.astype(int))
        r_coldata = ro.conversion.py2rpy(metadata)

    # Create DESeqDataSet
    dds = deseq2.DESeqDataSetFromMatrix(
        countData=r_counts,
        colData=r_coldata,
        design=Formula(design_formula)
    )

    # Run DESeq2
    dds = deseq2.DESeq(dds)

    # Extract results for each contrast
    results = {}
    for factor, numerator, denominator in contrasts:
        contrast_name = f"{factor}_{numerator}_vs_{denominator}"
        try:
            res = deseq2.results(
                dds,
                contrast=ro.StrVector([factor, numerator, denominator]),
                alpha=alpha
            )
            with localconverter(ro.default_converter + pandas2ri.converter):
                res_df = ro.conversion.rpy2py(base.as_data_frame(res))
            res_df.index = count_matrix.index
            results[contrast_name] = res_df
        except Exception as e:
            print(f"  Warning: contrast {contrast_name} failed: {e}")

    return results


def parse_contrasts(contrast_string):
    """Parse contrast string: 'Factor:Num:Denom,Factor2:Num2:Denom2'."""
    contrasts = []
    for c in contrast_string.split(','):
        parts = c.strip().split(':')
        if len(parts) == 3:
            contrasts.append(tuple(parts))
        else:
            print(f"  Warning: invalid contrast '{c}', expected Factor:Num:Denom")
    return contrasts


def main():
    parser = argparse.ArgumentParser(
        description='Pseudobulk DE with DESeq2'
    )
    parser.add_argument('--input', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--label_key', default='AIFI_L3')
    parser.add_argument('--sample_key', default='sample_id')
    parser.add_argument('--layer', default='counts')
    parser.add_argument('--design', default='~ Age + CMV + Sex')
    parser.add_argument('--contrasts', required=True,
                        help='Comma-separated Factor:Num:Denom')
    parser.add_argument('--min_cells', type=int, default=10)
    parser.add_argument('--min_pct_detected', type=float, default=0.10)
    parser.add_argument('--alpha', type=float, default=0.05)
    parser.add_argument('--lfc_threshold', type=float, default=0.1)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    contrasts = parse_contrasts(args.contrasts)

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    print("Creating pseudobulk...")
    pseudobulk = create_pseudobulk(
        adata, sample_key=args.sample_key, label_key=args.label_key,
        layer=args.layer, min_cells=args.min_cells
    )

    print(f"\nRunning DESeq2 with design: {args.design}")
    print(f"Contrasts: {contrasts}")

    summary = []
    for cell_type, (count_matrix, metadata) in pseudobulk.items():
        safe_name = cell_type.replace(' ', '_').replace('+', 'pos').replace('-', 'neg')
        print(f"\n  Processing: {cell_type}")

        # Filter genes
        ct_mask = adata.obs[args.label_key] == cell_type
        filtered_counts = filter_genes(
            count_matrix, adata, ct_mask, min_pct=args.min_pct_detected
        )
        print(f"    Genes after filtering: {filtered_counts.shape[0]}")

        # Run DESeq2
        results = run_deseq2(
            filtered_counts, metadata, args.design, contrasts, alpha=args.alpha
        )

        # Save results
        for contrast_name, res_df in results.items():
            out_file = os.path.join(
                args.output_dir, f"{safe_name}_{contrast_name}.csv"
            )
            res_df.to_csv(out_file)

            sig = res_df[
                (res_df['padj'] < args.alpha) &
                (res_df['log2FoldChange'].abs() > args.lfc_threshold)
            ]
            n_up = (sig['log2FoldChange'] > 0).sum()
            n_down = (sig['log2FoldChange'] < 0).sum()
            print(f"    {contrast_name}: {n_up} up, {n_down} down")

            summary.append({
                'cell_type': cell_type,
                'contrast': contrast_name,
                'n_genes_tested': len(res_df),
                'n_sig': len(sig),
                'n_up': n_up,
                'n_down': n_down,
            })

    # Save summary
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(args.output_dir, 'de_summary.csv'), index=False)
    print(f"\nSummary saved to {args.output_dir}/de_summary.csv")


if __name__ == '__main__':
    main()
