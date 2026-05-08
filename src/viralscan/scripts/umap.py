"""
The umap script creates a umap based on the single cell data provided where
each data point represents a single cell. Based on the colour it shows whether
the viral load has been detected in that single cell.
"""

# Importing packages
import os
import logging
import warnings
import numpy as np
import scanpy as sc
import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

from viralscan.constants import VIRUS_NAME_MAP
from viralscan.utils import load_config, setup_script_logging

log = setup_script_logging()

warnings.filterwarnings("ignore")

# Reading Snakefile params and config file
configfile = snakemake.params.configfile
config = load_config(configfile)


def calculate_k_neighbors(n_cells, min_k=10, max_k=200):
    """
    This function dynamically calculates the k-neighbors to determine
    whether or not there is a correlation between virus positive single
    cells.
    ---------------------------------------------------------------------
    Params:
        n_cells (int): total number of cells in data
        min_k (int): minimal k value
        max_k (int): maximal k value
    ---------------------------------------------------------------------
    Returns:
        k (int): k value to use for nearest-neighbor enrichment
    """
    # using square root scaling
    k = int(np.sqrt(n_cells))
    k = max(min_k, min(k, max_k))
    return k


def viral_neighbor_enrichment(coords, labels, k, n_permutations=1000, random_state=0):
    rng = np.random.default_rng(random_state)

    viral_cells = np.where(labels == 1)[0]
    if len(viral_cells) == 0:
        # No viral-positive cells — skip enrichment
        return 0.0, 0.0, 1.0

    nbrs = NearestNeighbors(n_neighbors=k).fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    counts = []
    for i in viral_cells:
        neighbor_idx = indices[i][1:]
        counts.append(np.mean(labels[neighbor_idx]))
    observed = np.mean(counts)

    permuted = []
    for _ in range(n_permutations):
        shuffled = rng.permutation(labels)
        counts_perm = []
        for i in viral_cells:
            neighbor_idx = indices[i][1:]
            counts_perm.append(np.mean(shuffled[neighbor_idx]))
        permuted.append(np.mean(counts_perm))

    expected = np.mean(permuted)
    p_value = (np.sum(np.array(permuted) >= observed) + 1) / (n_permutations + 1)
    return observed, expected, p_value


def umap(adata, found_genes, min_reads_per_cell=2, min_genes_per_cell=1):
    """
    Creates a UMAP of all single cells.
    If viral load is discovered, it highlights and tests for clustering.
    If no viral load is found, it still produces a general UMAP.
    """
    # Violin plot and computing QC metrics
    adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()

    # create violin plot
    p1 = sns.displot(adata.obs["n_counts"], bins=100, kde=False)
    plt.title("Total counts per cell")
    p1.savefig(f"{config['output']}/plots/qc_hist_total_counts.png")
    plt.close()

    # Filtering based on QC threshold (config-driven via PR 11 A4)
    min_counts_threshold = config.get("min_counts", 1000)
    min_genes_threshold = config.get("min_genes", 200)

    adata = adata[
        (adata.obs["n_counts"] >= min_counts_threshold)
        & (adata.obs["n_genes"] >= min_genes_threshold)
    ].copy()

    if adata.n_obs == 0:
        print(
            f"No cells left after QC filtering! "
            f"min_counts={min_counts_threshold}, min_genes={min_genes_threshold}"
        )
        return

    # Resolve counts source once — used for both has_viral check and viral count extraction
    if getattr(adata, "raw", None) and getattr(adata.raw, "X", None) is not None:
        X_counts = adata.raw.X
        var_names = adata.raw.var_names
    elif "counts" in adata.layers:
        X_counts = adata.layers["counts"]
        var_names = adata.var_names
    else:
        X_counts = adata.X
        var_names = adata.var_names

    viral_ids = [g for g in found_genes if g in var_names]
    has_viral = len(viral_ids) > 0

    if not has_viral:
        print("No viral genes detected. Creating standard UMAP (no virus labeling).")

        # Downsample large datasets to avoid memory issues
        if adata.n_obs > 20000:
            print(f"Dataset has {adata.n_obs} cells — subsampling to 10,000 for efficiency.")
            sc.pp.subsample(adata, n_obs=10000, random_state=0)

        # Normalization
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.tl.pca(adata, use_highly_variable=True, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)

        df = pd.DataFrame(
            {
                "UMAP1": adata.obsm["X_umap"][:, 0],
                "UMAP2": adata.obsm["X_umap"][:, 1],
                "barcode": adata.obs_names,
            }
        )

        fig = px.scatter(
            df,
            x="UMAP1",
            y="UMAP2",
            color_discrete_sequence=["#1f77b4"],
            hover_data=["barcode"],
            title="UMAP of All Cells (No Viral Load Detected)",
        )
        fig.update_layout(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white",
        )

        outdir = f"{config['output']}/plots"
        os.makedirs(outdir, exist_ok=True)
        fig.write_html(f"{outdir}/umap_no_virus.html")
        return

    print("Viral genes detected. Computing viral load and annotated UMAP.")
    k_neighbors = calculate_k_neighbors(adata.n_obs)

    viral_counts = np.zeros(adata.n_obs)
    viral_genes_expressed = np.zeros(adata.n_obs)
    viral_presence = {}
    var_name_to_idx = {g: i for i, g in enumerate(var_names)}
    for g in viral_ids:
        idx = var_name_to_idx[g]
        col = X_counts[:, idx]
        arr = col.toarray().flatten() if hasattr(col, "toarray") else np.asarray(col).flatten()
        viral_counts += arr
        viral_genes_expressed += (arr >= 1).astype(int)
        viral_presence[g] = (arr >= 1).astype(int)

    virus_labels = []
    gene_to_virus = {}
    for g in viral_presence:
        for key, virus_name in VIRUS_NAME_MAP.items():
            if g.startswith(key + "_") or g == key:
                gene_to_virus[g] = virus_name
                break
        else:
            gene_to_virus[g] = g  # if the ID is not found, get 'raw'

    for i in range(adata.n_obs):
        detected = list(
            set([gene_to_virus[g] for g in viral_presence if viral_presence[g][i] == 1])
        )
        if len(detected) == 0:
            virus_labels.append("No Virus")
        elif len(detected) == 1:
            virus_labels.append(detected[0])
        else:
            virus_labels.append("Multiple Viruses (" + ", ".join(detected) + ")")

    adata.obs["virus_detected"] = virus_labels

    adata.obs["viral_counts"] = viral_counts
    adata.obs["virus_positive"] = adata.obs["virus_detected"] != "No Virus"

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if "X_umap" not in adata.obsm:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        # Force-include viral genes so they survive HVG filtering
        adata.var["highly_variable"] |= adata.var_names.isin(list(found_genes))
        sc.pp.pca(adata, use_highly_variable=True, random_state=0)
        sc.pp.neighbors(adata, random_state=0)
        sc.tl.umap(adata, random_state=0)

    df = pd.DataFrame(
        {
            "UMAP1": adata.obsm["X_umap"][:, 0],
            "UMAP2": adata.obsm["X_umap"][:, 1],
            "viral_counts": adata.obs["viral_counts"],
            "virus_positive": adata.obs["virus_positive"],
            "barcode": adata.obs_names,
        }
    )

    coords = adata.obsm["X_umap"]
    labels = adata.obs["virus_positive"].astype(int).values
    observed, expected, pval = viral_neighbor_enrichment(coords, labels, k=k_neighbors)

    if pval < 0.05:
        conclusion_text = (
            f"Metric used: Nearest-Neighbor Enrichment (k={k_neighbors}). "
            f"Viral-positive cells show significant clustering (p-value < 0.05) "
            f"(observed={observed:.2f}, expected={expected:.2f}, p={pval:.3f})"
        )
    else:
        conclusion_text = (
            f"Metric used: Nearest-Neighbor Enrichment (k={k_neighbors}). "
            f"No significant clustering detected "
            f"(observed={observed:.2f}, expected={expected:.2f}, p={pval:.3f})"
        )

    df["virus_detected"] = adata.obs["virus_detected"].values

    # Binary UMAP
    df["virus_plot_label"] = np.where(df["virus_positive"], "Viral", "Non-viral")
    df["point_size"] = np.where(df["virus_positive"], 7, 4)

    fig_binary = px.scatter(
        df,
        x="UMAP1",
        y="UMAP2",
        color="virus_detected",
        size="point_size",
        hover_data=["barcode", "viral_counts"],
        size_max=6,
        title="UMAP of Virus-Detected Cells (Binary per Virus)",
    )

    fig_binary.update_traces(marker=dict(opacity=1))
    fig_binary.update_traces(marker=dict(line=dict(width=0)))
    fig_binary.for_each_trace(
        lambda t: t.update(marker_color="lightgray") if t.name == "No Virus" else None
    )

    fig_binary.update_layout(
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        margin=dict(b=100),
    )

    # Add annotation to show the conclusion text
    fig_binary.add_annotation(
        text=conclusion_text,
        xref="paper",
        yref="paper",
        x=0.5,
        y=-0.15,
        showarrow=False,
        font=dict(size=12),
        align="center",
    )

    # Continuous UMAP
    fig_continuous = px.scatter(
        df,
        x="UMAP1",
        y="UMAP2",
        color="viral_counts",
        color_continuous_scale="Reds",
        hover_data=["barcode", "virus_positive"],
        title="UMAP of Viral Load per Cell",
    )
    fig_continuous.update_layout(
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        margin=dict(b=100),
    )

    # Add annotation to show the conclusion text
    fig_continuous.add_annotation(
        text=conclusion_text,
        xref="paper",
        yref="paper",
        x=0.5,
        y=-0.15,
        showarrow=False,
        font=dict(size=12),
        align="center",
    )

    # Save plots to the users output directory
    outdir = f"{config['output']}/plots"
    os.makedirs(outdir, exist_ok=True)
    fig_binary.write_html(f"{outdir}/umap_binary.html")
    fig_continuous.write_html(f"{outdir}/umap_continuous.html")


def main():
    if config["multimapping"]:
        adata = sc.read_h5ad(f"{config['output']}/kb-python/counts_unfiltered/adata_multimap.h5ad")

        if "counts_corrected" in adata.layers and "counts_original" in adata.layers:
            adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
    else:
        adata = sc.read_h5ad(f"{config['output']}/kb-python/counts_unfiltered/adata.h5ad")

    # Load found genes
    found_genes = {}
    with open(f"{config['output']}/log/found_genes.txt") as f:
        for line in f:
            parts = line.strip().split(";")
            if len(parts) == 2:
                gene_id, count = parts
                found_genes[gene_id] = float(count)

    # Check if user wants UMAP
    if config["umap"]:
        print(
            f"You have decided to create a umap. This can take a while before finishing the code. Please wait..."
        )
        umap(adata, found_genes)


main()

# write to output file for Snakemake
with open(snakemake.output[0], "w") as f:
    f.write("done\n")
    if config["umap"]:
        log.info("Umap is done!")
    print(f"All (important) results of ViralScan can be found in {config['output']}summary.txt")
