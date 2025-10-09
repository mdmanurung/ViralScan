"""
The umap script creates a umap based on the single cell data provided where
each data point represents a single cell. Based on the colour it shows whether
the viral load has been detected in that single cell.
"""

# Importing packages
import os
import yaml
import warnings
import numpy as np
import scanpy as sc
import pandas as pd
import plotly.express as px
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import SparseEfficiencyWarning

# Ignore this warning
# warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
warnings.filterwarnings("ignore")

# Reading Snakefile params and config file
configfile = snakemake.params.configfile
with open(configfile, 'r') as f:
    config = yaml.safe_load(f)


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


def umap(adata, found_genes, min_reads_per_cell=2, min_genes_per_cell=1):
    """
    This function creates the umap for the user (if they wanted a umap).
    It creates 2 umaps: a continuous and binary umap.
    ---------------------------------------------------------------------
    Params:
        adata (anndata.AnnData): h5ad file of kb-python used for further
            analysis
        found_genes (dict): dictionary containing information of the gene 
            IDs found and the gene counts
        min_reads_per_cell (int): the minimum amount of reads per cell to
            show in the umap
        min_genes_per_cell (int): the minimum amount of genes per cell to
            show in the umap
    """
    # Filter for empty droplets
    cell_totals = np.array(adata.X.sum(axis=1)).flatten()
    adata = adata[cell_totals >= 1].copy()

    # Calculate dynamic k
    k_neighbors = calculate_k_neighbors(adata.n_obs)

    # Use raw counts if available
    if getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
        X_counts = adata.raw.X
        var_names = adata.raw.var_names
    elif "counts" in adata.layers:
        X_counts = adata.layers["counts"]
        var_names = adata.var_names
    else:
        X_counts = adata.X
        var_names = adata.var_names

    # Keep only viral genes present in adata
    viral_ids = [g for g in found_genes if g in var_names]
    if len(viral_ids) == 0:
        return

    # Compute viral counts and number of viral genes expressed
    viral_counts = np.zeros(adata.n_obs)
    viral_genes_expressed = np.zeros(adata.n_obs)
    for g in viral_ids:
        idx = list(var_names).index(g)
        col = X_counts[:, idx]
        arr = col.toarray().flatten() if hasattr(col, "toarray") else np.asarray(col).flatten()
        viral_counts += arr
        viral_genes_expressed += (arr >= 1).astype(int)

    adata.obs["viral_counts"] = viral_counts
    adata.obs["virus_positive"] = (viral_counts >= min_reads_per_cell) & (viral_genes_expressed >= min_genes_per_cell)

    # Normalize/log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Compute UMAP
    if "X_umap" not in adata.obsm:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata) #  n_neighbors=k_neighbors
        sc.tl.umap(adata)

    # Build dataframe
    df = pd.DataFrame({
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
        "virus_positive": adata.obs["virus_positive"],
        "viral_counts": adata.obs["viral_counts"],
        "barcode": adata.obs_names
    })

    # Compute the Nearest Neighbor Enrichment
    coords = adata.obsm["X_umap"]
    labels = adata.obs["virus_positive"].astype(int).values
    observed, expected, pval = viral_neighbor_enrichment(coords, labels, k=k_neighbors)

    if pval < 0.05:
        conclusion_text = (f"Metric used: Nearest-Neighbor Enrichment (k={k_neighbors}). "
                           f"Viral-positive cells show significant clustering (p-value < 0.05) "
                           f"(observed={observed:.2f}, expected={expected:.2f}, p={pval:.3f})")
    else:
        conclusion_text = (f"Metric used: Nearest-Neighbor Enrichment (k={k_neighbors}). "
                           f"No significant clustering detected "
                           f"(observed={observed:.2f}, expected={expected:.2f}, p={pval:.3f})")



    # Binary UMAP
    fig_binary = px.scatter(
        df, x="UMAP1", y="UMAP2",
        color="virus_positive",
        color_discrete_map={True: "blue", False: "lightgrey"},
        hover_data=["barcode", "viral_counts"],
        title="UMAP of Virus-Positive Cells"
    )
    fig_binary.update_layout(
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        margin=dict(b=100)
    )

    # Add annotation to show the conclusion text
    fig_binary.add_annotation(
        text=conclusion_text,
        xref="paper", yref="paper",
        x=0.5, y=-0.15,
        showarrow=False,
        font=dict(size=12),
        align="center"
    )

    # Continuous UMAP
    fig_continuous = px.scatter(
        df, x="UMAP1", y="UMAP2",
        color="viral_counts",
        color_continuous_scale="Reds",
        hover_data=["barcode","virus_positive"],
        title="UMAP of Viral Load per Cell"
    )
    fig_continuous.update_layout(
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        margin=dict(b=100)
    )

    # Add annotation to show the conclusion text
    fig_continuous.add_annotation(
        text=conclusion_text,
        xref="paper", yref="paper",
        x=0.5, y=-0.15,
        showarrow=False,
        font=dict(size=12),
        align="center"
    )

    # Save plots to the users output directory
    outdir = f"{config['output']}/plots"
    os.makedirs(outdir, exist_ok=True)
    fig_binary.write_html(f"{outdir}/umap_binary.html")
    fig_continuous.write_html(f"{outdir}/umap_continuous.html")


def viral_neighbor_enrichment(coords, labels, k, n_permutations=1000):
    """
    Compute nearest-neighbor enrichment for viral-positive cells to
    determine whether there is a connection between the place in the 
    plot and the viral reads.
    ---------------------------------------------------------------------
    Params:
        coords (numpy.ndarray): the UMAP embedding
        labels (numpy.ndarray): int value showing whether cell is virus 
            positive (1) or not (0) 
        k (int): the amount of neighbors
        n_permutations (int): total amount of permutations for the 
            enrichment
    ---------------------------------------------------------------------
    Returns:
        observed (float): average fraction of neighbors which are virus-
            positive using real labels
        expected (float): average fraction of virus-positive neighbors 
            with random labelling
        p_value (float): probability of observing the results
    """
    # Compute nearest neighbor
    nbrs = NearestNeighbors(n_neighbors=k).fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    # Check for other viral cells -> compute observed fraction
    viral_cells = np.where(labels == 1)[0] # only get viral cells
    counts = []
    for i in viral_cells:
        neighbor_idx = indices[i][1:]
        counts.append(np.mean(labels[neighbor_idx]))
    observed = np.mean(counts)

    # Null distribution -> compute expected fraction
    permuted = []
    for _ in range(n_permutations):
        shuffled = np.random.permutation(labels)
        counts_perm = []
        for i in viral_cells:
            neighbor_idx = indices[i][1:]
            counts_perm.append(np.mean(shuffled[neighbor_idx]))
        permuted.append(np.mean(counts_perm))
    expected = np.mean(permuted)

    p_value = (np.sum(np.array(permuted) >= observed) + 1) / (n_permutations + 1)
    return observed, expected, p_value



def main(): 
    adata = sc.read_h5ad(f"{config['output']}/kb-python/counts_unfiltered/adata_multimap.h5ad") 

    if "counts_corrected" in adata.layers and "counts_original" in adata.layers:
        adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]

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
        print(f"You have decided to create a umap. This can take a while before finishing the code. Please wait...") 
        umap(adata, found_genes)

main()

# write to output file for Snakemake
with open(snakemake.output[0], "w") as f: 
    f.write("done\n") 
    if config["umap"] == "True":
        print("\033[32mUmap is done!\033[0m") 
    print(f"All (important) results of ViralScan can be found in {config["output"]}summary.txt")