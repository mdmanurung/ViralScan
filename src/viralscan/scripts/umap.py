import yaml
import os
import scanpy as sc
import numpy as np
import pandas as pd
import plotly.express as px

# This rule from the snakemake pipeline creates all visuals.

configfile = snakemake.params.configfile
with open(configfile, 'r') as f:
    config = yaml.safe_load(f)

def umap(adata, viral_accessions):
    try:
        # Filtering for empty droplets (now 1, TODO: Check with Mikhael for threshold)
        cell_totals = adata.X.sum(axis=1).A1 if hasattr(adata.X, "toarray") else adata.X.sum(axis=1)
        cells_to_keep = cell_totals >= 1
        adata = adata[cells_to_keep].copy()

        # Create UMAP
        threshold = 1
        all_gene_ids = adata.var_names
        found_genes = {}

        for gene_id in viral_accessions:
            if gene_id in all_gene_ids:
                gene_counts = adata[:, gene_id].X
                if hasattr(gene_counts, 'toarray'):
                    gene_counts = gene_counts.toarray()
                total_count = gene_counts.sum()
                if total_count > threshold:
                    found_genes[gene_id] = total_count

        print(f"Found gene IDs including the gene counts:")
        for g in found_genes:
            print(g, found_genes[g])

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata) 

        # Create a per-cell viral gene count matrix
        viral_counts = np.zeros(adata.n_obs)  # n_obs is the number of cells

        for gene_id in found_genes.keys():
            gene_counts = adata[:, gene_id].X
            if hasattr(gene_counts, 'toarray'):
                gene_counts = gene_counts.toarray().flatten()
            viral_counts += gene_counts

        # Add as observation/cell-level metadata
        adata.obs['viral_counts'] = viral_counts
        adata.obs['virus_positive'] = adata.obs['viral_counts'] > 0  # Boolean label

        # UMAP calculation if not already done
        if 'X_umap' not in adata.obsm:
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

        # Build dataframe for Plotly
        umap = adata.obsm['X_umap']

        # Get metadata
        umap = adata.obsm['X_umap']
        df = pd.DataFrame({
            'UMAP1': umap[:, 0],
            'UMAP2': umap[:, 1],
            'virus_positive': adata.obs['virus_positive'],
            'viral_counts': adata.obs['viral_counts']
        })

        # Create interactive plot
        fig = px.scatter(
            df,
            x='UMAP1',
            y='UMAP2',
            color='virus_positive',
            color_discrete_map={True: 'blue', False: 'lightgrey'},
            hover_data=['viral_counts'],
            title="Interactive UMAP of Virus-Positive Cells"
        )

        fig.update_layout(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white'
        )
        if not os.path.exists(f"{config['output']}/plots/"):
                os.mkdir(f"{config['output']}/plots")
        fig.write_html(f"{config["output"]}/plots/umap.html")
        print("Umap finished!")
    except ValueError:
        print("After filtering, no cells were left. Hence, it is not possible to create a umap with a threshold of 1. Aborting the creation of umap...")


def main():
    adata = sc.read_h5ad(f"{config['output']}/kb-python/counts_unfiltered/adata.h5ad")
    viral_accesssions = []
    file = open(f"{config['output']}log/analysis.txt")
    for line in file:
        viral_accesssions.append(line.strip())

    if config["umap"]:   
        print(f"You have decided to create a umap. This can take a while before finishing the code. Please wait...")
        umap(adata, viral_accesssions)

main()

with open(snakemake.output[0], "w") as f:
    f.write("done\n")
print("\033[32mUmap (if chosen) is done!\033[0m")
print(f"All (important) results of ViralScan can be found in {config["output"]}summary.txt")