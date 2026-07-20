#!/usr/bin/env python3
"""
Cell type frequency analysis with CLR transformation and ALC normalization.

Usage:
    python clr_frequency_analysis.py \
        --input /path/to/annotated.h5ad \
        --output_dir /path/to/frequency_results \
        --label_keys AIFI_L1 AIFI_L2 AIFI_L3 \
        --sample_key sample_id \
        --alc_file /path/to/clinical_blood_counts.csv
"""

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


def tabulate_frequencies(adata, sample_key, label_key):
    """Compute cell type counts and fractions per sample."""
    counts = pd.crosstab(adata.obs[sample_key], adata.obs[label_key])
    totals = counts.sum(axis=1)
    fracs = counts.div(totals, axis=0)

    return counts, fracs, totals


def clr_transform(values_df, pseudocount=0):
    """
    Centered log-ratio transformation.

    Parameters
    ----------
    values_df : DataFrame
        Samples × cell types (counts or fractions)
    pseudocount : int
        Added before log. 0 = standard CLR, 1 = AIFI pseudocount approach
    """
    if pseudocount > 0:
        values_df = values_df + pseudocount

    # Replace zeros with NaN for log, then compute geometric mean
    log_vals = np.log(values_df.replace(0, np.nan))
    geo_mean = np.exp(log_vals.mean(axis=1))

    # CLR = log(value / geometric mean)
    clr = np.log(values_df.div(geo_mean, axis=0))

    # Replace -inf (from log(0)) with NaN
    clr = clr.replace([np.inf, -np.inf], np.nan)

    return clr


def compute_alc_normalized(counts, totals, alc_data, sample_key,
                            alc_column='lymphocyte_count',
                            l1_counts=None):
    """
    Normalize cell counts to absolute lymphocyte counts (ALC).

    Parameters
    ----------
    counts : DataFrame, samples × cell types
    totals : Series, total cells per sample
    alc_data : DataFrame with sample_key and alc_column
    l1_counts : DataFrame, optional L1-level counts for lymphocyte calculation
    """
    # If we have L1 counts, compute scRNA lymphocyte count
    if l1_counts is not None:
        lymph_cols = [c for c in ['T cell', 'NK cell', 'B cell']
                      if c in l1_counts.columns]
        scrna_lymph = l1_counts[lymph_cols].sum(axis=1)
    else:
        # Assume all cells are lymphocyte-derived (rough approximation)
        scrna_lymph = totals

    # Merge with clinical ALC
    alc_series = alc_data.set_index(sample_key)[alc_column]
    common_samples = counts.index.intersection(alc_series.index)

    if len(common_samples) == 0:
        print("  Warning: no matching samples between scRNA and ALC data")
        return None, None

    # ALC ratio
    alc_ratio = alc_series[common_samples] / scrna_lymph[common_samples]

    # Fractions × ALC = estimated absolute count
    fracs = counts.div(totals, axis=0)
    alc_counts = fracs.loc[common_samples].multiply(
        alc_series[common_samples], axis=0
    )

    return alc_counts, alc_ratio


def compare_groups(clr_df, metadata, group_col, group_a, group_b,
                    sample_key='sample_id'):
    """Wilcoxon rank-sum test on CLR values between two groups."""
    meta = metadata.set_index(sample_key) if sample_key in metadata.columns else metadata

    samples_a = meta[meta[group_col] == group_a].index
    samples_b = meta[meta[group_col] == group_b].index

    samples_a = clr_df.index.intersection(samples_a)
    samples_b = clr_df.index.intersection(samples_b)

    results = []
    for cell_type in clr_df.columns:
        vals_a = clr_df.loc[samples_a, cell_type].dropna()
        vals_b = clr_df.loc[samples_b, cell_type].dropna()

        if len(vals_a) >= 3 and len(vals_b) >= 3:
            stat, pval = mannwhitneyu(vals_a, vals_b, alternative='two-sided')
            results.append({
                'cell_type': cell_type,
                f'mean_{group_a}': vals_a.mean(),
                f'mean_{group_b}': vals_b.mean(),
                'diff': vals_b.mean() - vals_a.mean(),
                'pvalue': pval,
                f'n_{group_a}': len(vals_a),
                f'n_{group_b}': len(vals_b),
            })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df['padj'] = multipletests(results_df['pvalue'], method='fdr_bh')[1]
        results_df = results_df.sort_values('padj')

    return results_df


def compare_paired(clr_df, metadata, timepoint_col, tp_a, tp_b,
                    subject_key='subject_id', sample_key='sample_id'):
    """Paired signed-rank Wilcoxon test for longitudinal comparisons."""
    meta = metadata.set_index(sample_key) if sample_key in metadata.columns else metadata

    # Find paired subjects
    samples_a = meta[meta[timepoint_col] == tp_a]
    samples_b = meta[meta[timepoint_col] == tp_b]

    paired_subjects = set(samples_a[subject_key]) & set(samples_b[subject_key])

    results = []
    for cell_type in clr_df.columns:
        vals_a, vals_b = [], []
        for subj in paired_subjects:
            sa = samples_a[samples_a[subject_key] == subj].index
            sb = samples_b[samples_b[subject_key] == subj].index
            sa = clr_df.index.intersection(sa)
            sb = clr_df.index.intersection(sb)
            if len(sa) == 1 and len(sb) == 1:
                va = clr_df.loc[sa[0], cell_type]
                vb = clr_df.loc[sb[0], cell_type]
                if not (np.isnan(va) or np.isnan(vb)):
                    vals_a.append(va)
                    vals_b.append(vb)

        if len(vals_a) >= 5:
            try:
                stat, pval = wilcoxon(vals_b, vals_a)
            except ValueError:
                pval = 1.0
            results.append({
                'cell_type': cell_type,
                f'mean_{tp_a}': np.mean(vals_a),
                f'mean_{tp_b}': np.mean(vals_b),
                'diff': np.mean(vals_b) - np.mean(vals_a),
                'pvalue': pval,
                'n_pairs': len(vals_a),
            })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df['padj'] = multipletests(results_df['pvalue'], method='fdr_bh')[1]
        results_df = results_df.sort_values('padj')

    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='Cell type frequency analysis with CLR and ALC'
    )
    parser.add_argument('--input', required=True, help='Annotated h5ad')
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--label_keys', nargs='+',
                        default=['AIFI_L1', 'AIFI_L2', 'AIFI_L3'])
    parser.add_argument('--sample_key', default='sample_id')
    parser.add_argument('--alc_file', default=None,
                        help='CSV with ALC data per sample')
    parser.add_argument('--alc_column', default='lymphocyte_count')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Compute L1 counts for ALC normalization
    l1_counts = None
    if 'AIFI_L1' in adata.obs.columns:
        l1_counts, _, _ = tabulate_frequencies(adata, args.sample_key, 'AIFI_L1')

    for label_key in args.label_keys:
        if label_key not in adata.obs.columns:
            print(f"  Skipping {label_key} (not in adata.obs)")
            continue

        print(f"\n=== {label_key} ===")

        # Tabulate
        counts, fracs, totals = tabulate_frequencies(
            adata, args.sample_key, label_key
        )
        print(f"  {counts.shape[1]} cell types, {counts.shape[0]} samples")

        # CLR (standard and with pseudocount)
        clr_standard = clr_transform(fracs)
        counts_pseudo = counts + 1
        fracs_pseudo = counts_pseudo.div(counts_pseudo.sum(axis=1), axis=0)
        clr_pseudo = clr_transform(fracs_pseudo)

        # Save all tables
        prefix = os.path.join(args.output_dir, f"{label_key}")
        counts.to_csv(f"{prefix}_counts.csv")
        fracs.to_csv(f"{prefix}_fractions.csv")
        clr_standard.to_csv(f"{prefix}_clr.csv")
        clr_pseudo.to_csv(f"{prefix}_clr_pseudo.csv")

        # ALC normalization
        if args.alc_file:
            alc_data = pd.read_csv(args.alc_file)
            alc_counts, alc_ratio = compute_alc_normalized(
                counts, totals, alc_data, args.sample_key,
                args.alc_column, l1_counts
            )
            if alc_counts is not None:
                alc_counts.to_csv(f"{prefix}_alc_counts.csv")
                print(f"  ALC normalization: {len(alc_counts)} samples")

        # Summary stats
        print(f"  Saved: counts, fractions, CLR, CLR_pseudo")

    # Also save total cells per sample
    totals_df = pd.DataFrame({
        'sample_id': totals.index,
        'total_cells': totals.values,
    })
    totals_df.to_csv(os.path.join(args.output_dir, 'total_cells_per_sample.csv'),
                     index=False)

    print(f"\nAll results saved to {args.output_dir}")


if __name__ == '__main__':
    main()
