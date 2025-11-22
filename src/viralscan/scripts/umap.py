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
import seaborn as sns
import matplotlib.pyplot as plt
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

    # Violin: total_counts
    sc.settings.figdir = f"{config['output']}/plots/"
    sc.pl.violin(adata, "n_counts", save="_n_counts.png")

    sc.pl.violin(adata, "n_genes", save="_n_genes.png")

    # Filtering based on QC threshold
    before = adata.n_obs
    min_counts_threshold = 1000
    min_genes_threshold = 200

    adata = adata[
        (adata.obs["n_counts"] >= min_counts_threshold) &
        (adata.obs["n_genes"] >= min_genes_threshold)
    ].copy()
    after = adata.n_obs
    # print(f"Filtered cells: {before} → {after}")

    # Filtering (for now, different filtering QC have been chosen)
    # cell_totals = np.array(adata.X.sum(axis=1)).flatten()
    # adata = adata[cell_totals >= 1].copy()

    # Identify viral genes and check
    var_names = adata.raw.var_names if getattr(adata, "raw", None) else adata.var_names
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
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)

        df = pd.DataFrame({
            "UMAP1": adata.obsm["X_umap"][:, 0],
            "UMAP2": adata.obsm["X_umap"][:, 1],
            "barcode": adata.obs_names
        })

        fig = px.scatter(
            df, x="UMAP1", y="UMAP2",
            color_discrete_sequence=["#1f77b4"],
            hover_data=["barcode"],
            title="UMAP of All Cells (No Viral Load Detected)"
        )
        fig.update_layout(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white"
        )

        outdir = f"{config['output']}/plots"
        os.makedirs(outdir, exist_ok=True)
        fig.write_html(f"{outdir}/umap_no_virus.html")
        return

    print("Viral genes detected. Computing viral load and annotated UMAP.")
    k_neighbors = calculate_k_neighbors(adata.n_obs)

    # Compute viral counts
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
    viral_counts = np.zeros(adata.n_obs)
    viral_genes_expressed = np.zeros(adata.n_obs)
    viral_presence = {}
    for g in viral_ids:
        idx = list(var_names).index(g)
        col = X_counts[:, idx]
        arr = col.toarray().flatten() if hasattr(col, "toarray") else np.asarray(col).flatten()
        viral_counts += arr
        viral_genes_expressed += (arr >= 1).astype(int)
        viral_presence[g] = (arr >= 1).astype(int)

    map_virus = {"AICHI": "Aichi virus", "AUSBATLYSSA": "Australian Bat Lyssavirus", "BANNA": "Banna virus", "BARMAH": "Barmah forest virus", "BKPOLY": "BK polyomavirus",
        "BUNYAMW":"Bunyamwera virus", "BUNYA": "Bunyavirus La Crosse" , "CERC_HERP": "Cercopithecine herpesvirus", "CHIKUNG": "Chikungunya virus", "COSA_A": "Cosavirus A",
        "COWPOX": "Cowpox virus", "COXSACKIE": "Coxsackievirus", "CRIMEAN": "Crimean Congo hemorrhagic fever virus", "EQUINE_ENCE": "Eastern_equine_encephalitis_virus",
        "EBOLA": "Ebolavirus", "ECHO": "Echovirus", "ENCEPHAL": "Encephalomyocarditis virus", "EPSTEIN": "Epstein-Barr virus", "EURBATLYSSA": "European bat lyssavirus",
        "GB": "GB virus C_Hepatitis G virus", "HANTAAN": "Hantaan virus", "HENDRA": "Hendra virus","HEP_A": "Hepatitis A virus", "HEP_B": "Hepatitis B virus",
        "HEP_C": "Hepatitis C virus", "HEP_DELTA": "Hepatitis delta virus", "HEP_E": "Hepatitis E virus", "HUM_ADENO": "Human adenovirus", "HUM_ASTRO": "Human astrovirus",
        "HUM_COR": "Human coronavirus", "HUM_CYTO": "Human cytomegalovirus", "HUM_ENTERO": "Human enterovirus 68, 70", "HUM_HERP1": "Human herpesvirus 1",
        "HUM_HERP2": "Human herpesvirus 2", "HUM_HERP6B": "Human herpesvirus 6b", "HUM_HERP6": "Human herpesvirus 6", "HUM_HERP7": "Human herpesvirus 7", "HUM_HERP8": "Human herpesvirus 8",
        "HUM_PAP_1618": "Human papillomavirus 16,18", "HUM_PAP_1": "Human papillomavirus 1", "HUM_PAP_2": "Human papillomavirus 2", "HUM_PARA": "Human parainfluenza",
        "HUM_PARVO": "Human parvovirus B19", "HUM_RESP": "Human respiratory syncytial virus", "HUM_RHINO": "Human rhinovirus", "HUM_SARS": "Human SARS coronavirus",
        "INFL_A": "Influenza A virus", "INFL_B": "Influenza B virus", "INFL_C": "Influenza C virus", "JAP_ENCE": "Japanese encephalitis virus", "JC_POLY": "JC polyomavirus",
        "KI_POLY": "KI Polyomavirus", "LAKE_VIC": "Lake Victoria marburgvirus", "LANGAT": "Langat virus", "LASSA": "Lassa virus", "LOUPING": "Louping ill virus",
        'LYMPH': "Lymphocytic choriomeningitis virus", "MAYARO": "Mayaro virus", "MEASLES": "Measles virus", "MERKEL": "Merkel cell polyomavirus", "MERS":"MERS coronavirus",
        "MOLLU": "Molluscum contagiosum virus", "MONKEYPOX": "Monkeypox virus", "MUMPS": "Mumps virus", "MUR_VAL": "Murray valley encephalitis virus", "NIPAH": "Nipah virus",
        "NORWALK": "Norwalk virus", "ORF": "Orf virus", "OROPOU": "Oropouche virus", "ONYONG":"O�nyong-nyong virus", "POLIO": "Poliovirus", "RABIES": "Rabies virus",
        "ROSA_A": "Rosavirus A", "ROSS_RIVER": "Ross river virus", "ROTA_A": "Rotavirus A", "ROTA_B": "Rotavirus B", "ROTA_C": "Rotavirus C", "RUBELLA": "Rubella virus",
        "SALI_A": "Salivirus A", "SAPPORO": "Sapporo virus", "SEMLIKI": "Semliki forest virus", "SEOUL": "Seoul virus", "SINDBIS": "Sindbis virus",
        "ST_LOUIS": "St. louis encephalitis virus", "TICK": "Tick-borne powassan virus", "TTV": "Torque teno virus", "TOSCANA": "Toscana virus", "VACCINIA": "Vaccinia virus",
        "VARICELLA": "Varicella-zoster virus", "VARIOLA": "Variola virus", "VEN_EQU": "Venezuelan equine encephalitis virus", "VES_STOM": "Vesicular stomatitis virus",
        "WES_EQU": "Western equine encephalitis virus", "WES_NILE": "West Nile virus", "WU_POLY": "WU polyomavirus", "YABA": "Yaba-like disease virus", "YELLOW": "Yellow fever virus",
        "ZIKA": "Zika virus",
    }

    virus_labels = []
    gene_to_virus = {}
    for g in viral_presence:
        for key, virus_name in map_virus.items():
            if g.startswith(key + "_") or g == key:
                gene_to_virus[g] = virus_name
                break
        else:
            gene_to_virus[g] = g # if the ID is not found, get 'raw'

    for i in range(adata.n_obs):
        detected = list(set([gene_to_virus[g] for g in viral_presence if viral_presence[g][i] == 1]))
        if len(detected) == 0:
            virus_labels.append("No Virus")
        elif len(detected) == 1:
            virus_labels.append(detected[0])
        else:
            virus_labels.append("Multiple Viruses (" + ", ".join(detected) + ")")

    adata.obs["virus_detected"] = virus_labels

    adata.obs["viral_counts"] = viral_counts
    adata.obs["virus_positive"] = adata.obs["virus_detected"] != "No Virus"

    # HVG selection only when viral genes are present
    # sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=1500)
    # viral_present = [g for g in viral_ids if g in adata.var_names]
    # keep_genes = list(dict.fromkeys(adata.var.index[adata.var.highly_variable].tolist() + viral_present))
    # adata = adata[:, keep_genes].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # sc.pp.scale(adata, max_value=10)
    # sc.tl.pca(adata, svd_solver="arpack")
    # sc.pp.neighbors(adata)
    # sc.tl.umap(adata)

    if "X_umap" not in adata.obsm:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata) #  n_neighbors=k_neighbors
        sc.tl.umap(adata)

    df = pd.DataFrame({
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
        "viral_counts": adata.obs["viral_counts"],
        "virus_positive": adata.obs["virus_positive"],
        "barcode": adata.obs_names
    })

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


    df["virus_detected"] = adata.obs["virus_detected"].values

    # Binary UMAP
    fig_binary = px.scatter(
        df, x="UMAP1", y="UMAP2",
        color="virus_detected",
        hover_data=["barcode", "viral_counts"],
        title="UMAP of Virus-Detected Cells (Binary per Virus)",
        color_discrete_sequence=px.colors.qualitative.Safe
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

# def umap(adata, found_genes, min_reads_per_cell=2, min_genes_per_cell=1):
#     """
#     This function creates the umap for the user (if they wanted a umap).
#     It creates 2 umaps: a continuous and binary umap.
#     ---------------------------------------------------------------------
#     Params:
#         adata (anndata.AnnData): h5ad file of kb-python used for further
#             analysis
#         found_genes (dict): dictionary containing information of the gene 
#             IDs found and the gene counts
#         min_reads_per_cell (int): the minimum amount of reads per cell to
#             show in the umap
#         min_genes_per_cell (int): the minimum amount of genes per cell to
#             show in the umap
#     """
#     # Filter empty droplets
#     cell_totals = np.array(adata.X.sum(axis=1)).flatten()
#     adata = adata[cell_totals >= 1].copy()

#     # Determine data source for counts
#     if getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
#         X_counts = adata.raw.X
#         var_names = adata.raw.var_names
#     elif "counts" in adata.layers:
#         X_counts = adata.layers["counts"]
#         var_names = adata.var_names
#     else:
#         X_counts = adata.X
#         var_names = adata.var_names

#     # Detect viral genes
#     viral_ids = [g for g in found_genes if g in var_names]
#     has_viral = len(viral_ids) > 0

#     # Default annotations
#     adata.obs["virus_detected"] = "No Virus"
#     adata.obs["viral_counts"] = 0.0
#     adata.obs["virus_positive"] = False

#     if has_viral:
#         print(f"Detected {len(viral_ids)} viral genes; computing viral load metrics...")
#         viral_counts = np.zeros(adata.n_obs)
#         viral_presence = {}
#         for g in viral_ids:
#             idx = list(var_names).index(g)
#             col = X_counts[:, idx]
#             arr = col.toarray().flatten() if hasattr(col, "toarray") else np.asarray(col).flatten()
#             viral_counts += arr
#             viral_presence[g] = (arr >= 1).astype(int)

#         virus_labels = []
#         for i in range(adata.n_obs):
#             detected = [g for g in viral_presence if viral_presence[g][i] == 1]
#             virus_labels.append("Virus Detected" if len(detected) > 0 else "No Virus")

#         adata.obs["virus_detected"] = virus_labels
#         adata.obs["viral_counts"] = viral_counts
#         adata.obs["virus_positive"] = adata.obs["virus_detected"] != "No Virus"
#     else:
#         print("No viral genes detected. Proceeding with standard UMAP.")

#     # --- Memory-safe data sanitization ---
#     # Remove all-zero genes
#     adata = adata[:, np.array((adata.X > 0).sum(axis=0)).flatten() > 0].copy()

#     # Clean invalid values if sparse
#     if hasattr(adata.X, "data"):
#         adata.X.data[np.isnan(adata.X.data)] = 0
#         adata.X.data[np.isinf(adata.X.data)] = 0

#     # --- HVG selection ---
#     if has_viral:
#         try:
#             sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
#             keep_genes = adata.var.index[adata.var['highly_variable']].tolist()
#             adata = adata[:, keep_genes].copy()
#         except Exception as e:
#             print(f"Warning: HVG selection failed ({e}). Proceeding without HVG filtering.")
#     else:
#         print("Skipping HVG selection for virus-free dataset.")

#     # --- Normalization and embedding ---
#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)
#     sc.pp.scale(adata, max_value=10)

#     sc.pp.pca(adata)
#     sc.pp.neighbors(adata)
#     sc.tl.umap(adata)

#     # --- DataFrame for plotting ---
#     df = pd.DataFrame({
#         "UMAP1": adata.obsm["X_umap"][:, 0],
#         "UMAP2": adata.obsm["X_umap"][:, 1],
#         "virus_detected": adata.obs["virus_detected"],
#         "viral_counts": adata.obs["viral_counts"],
#         "barcode": adata.obs_names
#     })

#     # --- Viral enrichment test (if applicable) ---
#     k_neighbors = calculate_k_neighbors(adata.n_obs)
#     if has_viral and adata.obs["virus_positive"].sum() > 0:
#         observed, expected, pval = viral_neighbor_enrichment(
#             adata.obsm["X_umap"],
#             adata.obs["virus_positive"].astype(int).values,
#             k=k_neighbors
#         )
#         if pval < 0.05:
#             conclusion_text = (
#                 f"Nearest-Neighbor Enrichment (k={k_neighbors}): "
#                 f"Significant clustering (p={pval:.3f})."
#             )
#         else:
#             conclusion_text = (
#                 f"Nearest-Neighbor Enrichment (k={k_neighbors}): "
#                 f"No significant clustering (p={pval:.3f})."
#             )
#     else:
#         conclusion_text = "No viral load detected; enrichment analysis skipped."

#     # --- Binary coloring (virus detected vs not) ---
#     fig_binary = px.scatter(
#         df, x="UMAP1", y="UMAP2",
#         color="virus_detected",
#         title="UMAP of Single Cells (Virus Presence)",
#         hover_data=["barcode", "viral_counts"],
#         color_discrete_sequence=px.colors.qualitative.Safe
#     )
#     fig_binary.update_layout(
#         xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#         yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#         plot_bgcolor="white",
#         margin=dict(b=100)
#     )
#     fig_binary.add_annotation(
#         text=conclusion_text,
#         xref="paper", yref="paper",
#         x=0.5, y=-0.15, showarrow=False, font=dict(size=12), align="center"
#     )

#     # --- Continuous viral load plot ---
#     fig_continuous = px.scatter(
#         df, x="UMAP1", y="UMAP2",
#         color="viral_counts",
#         color_continuous_scale="Reds",
#         title="UMAP of Viral Load per Cell",
#         hover_data=["barcode"]
#     )
#     fig_continuous.update_layout(
#         xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#         yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#         plot_bgcolor="white",
#         margin=dict(b=100)
#     )
#     fig_continuous.add_annotation(
#         text=conclusion_text,
#         xref="paper", yref="paper",
#         x=0.5, y=-0.15, showarrow=False, font=dict(size=12), align="center"
#     )

#     # --- Save plots ---
#     outdir = f"{config['output']}/plots"
#     os.makedirs(outdir, exist_ok=True)
#     fig_binary.write_html(f"{outdir}/umap_binary.html")
#     fig_continuous.write_html(f"{outdir}/umap_continuous.html")

    # # Filter for empty droplets
    # cell_totals = np.array(adata.X.sum(axis=1)).flatten()
    # adata = adata[cell_totals >= 1].copy()

    # # Calculate dynamic k
    # k_neighbors = calculate_k_neighbors(adata.n_obs)

    # # Use raw counts if available
    # if getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
    #     X_counts = adata.raw.X
    #     var_names = adata.raw.var_names
    # elif "counts" in adata.layers:
    #     X_counts = adata.layers["counts"]
    #     var_names = adata.var_names
    # else:
    #     X_counts = adata.X
    #     var_names = adata.var_names

    # # Keep only viral genes present in adata
    # viral_ids = [g for g in found_genes if g in var_names]

    # if len(viral_ids) == 0: 
    #     print("No viral genes detected. Proceeding with UMAP without virus labels.")
    #     adata.obs["virus_detected"] = ["No Virus"] * adata.n_obs
    #     adata.obs["viral_counts"] = np.zeros(adata.n_obs)
    #     adata.obs["virus_positive"] = np.zeros(adata.n_obs, dtype=bool)
    # else:
    #     # Compute viral counts and number of viral genes expressed
    #     viral_counts = np.zeros(adata.n_obs)
    #     viral_genes_expressed = np.zeros(adata.n_obs)
    #     viral_presence = {}
    #     for g in viral_ids:
    #         idx = list(var_names).index(g)
    #         col = X_counts[:, idx]
    #         arr = col.toarray().flatten() if hasattr(col, "toarray") else np.asarray(col).flatten()
    #         viral_counts += arr
    #         viral_genes_expressed += (arr >= 1).astype(int)
    #         viral_presence[g] = (arr >= 1).astype(int)
    
    #     map_virus = {"AICHI": "Aichi virus", "AUSBATLYSSA": "Australian Bat Lyssavirus", "BANNA": "Banna virus", "BARMAH": "Barmah forest virus", "BKPOLY": "BK polyomavirus",
    #          "BUNYAMW":"Bunyamwera virus", "BUNYA": "Bunyavirus La Crosse" , "CERC_HERP": "Cercopithecine herpesvirus", "CHIKUNG": "Chikungunya virus", "COSA_A": "Cosavirus A",
    #          "COWPOX": "Cowpox virus", "COXSACKIE": "Coxsackievirus", "CRIMEAN": "Crimean Congo hemorrhagic fever virus", "EQUINE_ENCE": "Eastern_equine_encephalitis_virus",
    #          "EBOLA": "Ebolavirus", "ECHO": "Echovirus", "ENCEPHAL": "Encephalomyocarditis virus", "EPSTEIN": "Epstein-Barr virus", "EURBATLYSSA": "European bat lyssavirus",
    #          "GB": "GB virus C_Hepatitis G virus", "HANTAAN": "Hantaan virus", "HENDRA": "Hendra virus","HEP_A": "Hepatitis A virus", "HEP_B": "Hepatitis B virus",
    #          "HEP_C": "Hepatitis C virus", "HEP_DELTA": "Hepatitis delta virus", "HEP_E": "Hepatitis E virus", "HUM_ADENO": "Human adenovirus", "HUM_ASTRO": "Human astrovirus",
    #          "HUM_COR": "Human coronavirus", "HUM_CYTO": "Human cytomegalovirus", "HUM_ENTERO": "Human enterovirus 68, 70", "HUM_HERP1": "Human herpesvirus 1",
    #          "HUM_HERP2": "Human herpesvirus 2", "HUM_HERP6B": "Human herpesvirus 6b", "HUM_HERP6": "Human herpesvirus 6", "HUM_HERP7": "Human herpesvirus 7", "HUM_HERP8": "Human herpesvirus 8",
    #          "HUM_PAP_1618": "Human papillomavirus 16,18", "HUM_PAP_1": "Human papillomavirus 1", "HUM_PAP_2": "Human papillomavirus 2", "HUM_PARA": "Human parainfluenza",
    #          "HUM_PARVO": "Human parvovirus B19", "HUM_RESP": "Human respiratory syncytial virus", "HUM_RHINO": "Human rhinovirus", "HUM_SARS": "Human SARS coronavirus",
    #          "INFL_A": "Influenza A virus", "INFL_B": "Influenza B virus", "INFL_C": "Influenza C virus", "JAP_ENCE": "Japanese encephalitis virus", "JC_POLY": "JC polyomavirus",
    #          "KI_POLY": "KI Polyomavirus", "LAKE_VIC": "Lake Victoria marburgvirus", "LANGAT": "Langat virus", "LASSA": "Lassa virus", "LOUPING": "Louping ill virus",
    #          'LYMPH': "Lymphocytic choriomeningitis virus", "MAYARO": "Mayaro virus", "MEASLES": "Measles virus", "MERKEL": "Merkel cell polyomavirus", "MERS":"MERS coronavirus",
    #          "MOLLU": "Molluscum contagiosum virus", "MONKEYPOX": "Monkeypox virus", "MUMPS": "Mumps virus", "MUR_VAL": "Murray valley encephalitis virus", "NIPAH": "Nipah virus",
    #          "NORWALK": "Norwalk virus", "ORF": "Orf virus", "OROPOU": "Oropouche virus", "ONYONG":"O�nyong-nyong virus", "POLIO": "Poliovirus", "RABIES": "Rabies virus",
    #          "ROSA_A": "Rosavirus A", "ROSS_RIVER": "Ross river virus", "ROTA_A": "Rotavirus A", "ROTA_B": "Rotavirus B", "ROTA_C": "Rotavirus C", "RUBELLA": "Rubella virus",
    #          "SALI_A": "Salivirus A", "SAPPORO": "Sapporo virus", "SEMLIKI": "Semliki forest virus", "SEOUL": "Seoul virus", "SINDBIS": "Sindbis virus",
    #          "ST_LOUIS": "St. louis encephalitis virus", "TICK": "Tick-borne powassan virus", "TTV": "Torque teno virus", "TOSCANA": "Toscana virus", "VACCINIA": "Vaccinia virus",
    #          "VARICELLA": "Varicella-zoster virus", "VARIOLA": "Variola virus", "VEN_EQU": "Venezuelan equine encephalitis virus", "VES_STOM": "Vesicular stomatitis virus",
    #          "WES_EQU": "Western equine encephalitis virus", "WES_NILE": "West Nile virus", "WU_POLY": "WU polyomavirus", "YABA": "Yaba-like disease virus", "YELLOW": "Yellow fever virus",
    #          "ZIKA": "Zika virus",
    #         }

    #     virus_labels = []
    #     gene_to_virus = {}
    #     for g in viral_presence:
    #         for key, virus_name in map_virus.items():
    #             if g.startswith(key + "_") or g == key:
    #                 gene_to_virus[g] = virus_name
    #                 break
    #         else:
    #             gene_to_virus[g] = g # if the ID is not found, get 'raw'


    #     for i in range(adata.n_obs):
    #         detected = list(set([gene_to_virus[g] for g in viral_presence if viral_presence[g][i] == 1]))
    #         if len(detected) == 0:
    #             virus_labels.append("No Virus")
    #         elif len(detected) == 1:
    #             virus_labels.append(detected[0])
    #         else:
    #             virus_labels.append("Multiple Viruses (" + ", ".join(detected) + ")")

    #     adata.obs["virus_detected"] = virus_labels

    #     adata.obs["viral_counts"] = viral_counts
    #     adata.obs["virus_positive"] = adata.obs["virus_detected"] != "No Virus"
    #     # adata.obs["virus_positive"] = (viral_counts >= min_reads_per_cell) & (viral_genes_expressed >= min_genes_per_cell)

    # # added these, might delete later (to avoid circular UMAPs)

    # # Remove genes with all-zero expression
    # adata = adata[:, np.array((adata.X > 0).sum(axis=0)).flatten() > 0].copy()

    # # Convert sparse matrix to dense safely
    # if hasattr(adata.X, "toarray"):
    #     X_tmp = adata.X.toarray()
    # else:
    #     X_tmp = np.asarray(adata.X)

    # # Replace any inf/nan with 0
    # X_tmp[np.isnan(X_tmp)] = 0
    # X_tmp[np.isinf(X_tmp)] = 0
    # adata.X = X_tmp
    # # assume viral_ids is a list of gene IDs present in var_names
    # sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

    # # union HVGs with viral_ids that are actually present in adata.var_names
    # viral_present = [g for g in viral_ids if g in adata.var_names]
    # hvg = adata.var.index[adata.var['highly_variable']].tolist()
    # keep_genes = list(dict.fromkeys(hvg + viral_present))  # preserve order, unique

    # adata = adata[:, keep_genes].copy()
    # # adata = adata[:, adata.var["highly_variable"]].copy()

    # # Normalize/log
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.scale(adata, max_value=10)

    # # Compute UMAP
    # if "X_umap" not in adata.obsm:
    #     sc.pp.pca(adata)
    #     sc.pp.neighbors(adata) #  n_neighbors=k_neighbors
    #     sc.tl.umap(adata)

    # # Build dataframe
    # df = pd.DataFrame({
    #     "UMAP1": adata.obsm["X_umap"][:, 0],
    #     "UMAP2": adata.obsm["X_umap"][:, 1],
    #     "virus_positive": adata.obs["virus_positive"],
    #     "viral_counts": adata.obs["viral_counts"],
    #     "barcode": adata.obs_names
    # })

    # # Compute the Nearest Neighbor Enrichment
    # coords = adata.obsm["X_umap"]
    # labels = adata.obs["virus_positive"].astype(int).values
    # observed, expected, pval = viral_neighbor_enrichment(coords, labels, k=k_neighbors)

    # if pval < 0.05:
    #     conclusion_text = (f"Metric used: Nearest-Neighbor Enrichment (k={k_neighbors}). "
    #                        f"Viral-positive cells show significant clustering (p-value < 0.05) "
    #                        f"(observed={observed:.2f}, expected={expected:.2f}, p={pval:.3f})")
    # else:
    #     conclusion_text = (f"Metric used: Nearest-Neighbor Enrichment (k={k_neighbors}). "
    #                        f"No significant clustering detected "
    #                        f"(observed={observed:.2f}, expected={expected:.2f}, p={pval:.3f})")


    # df["virus_detected"] = adata.obs["virus_detected"].values

    # # Binary UMAP
    # fig_binary = px.scatter(
    #     df, x="UMAP1", y="UMAP2",
    #     color="virus_detected",
    #     hover_data=["barcode", "viral_counts"],
    #     title="UMAP of Virus-Detected Cells (Binary per Virus)",
    #     color_discrete_sequence=px.colors.qualitative.Safe
    # )
    # fig_binary.update_layout(
    #     xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #     yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #     plot_bgcolor="white",
    #     margin=dict(b=100)
    # )

    # # Add annotation to show the conclusion text
    # fig_binary.add_annotation(
    #     text=conclusion_text,
    #     xref="paper", yref="paper",
    #     x=0.5, y=-0.15,
    #     showarrow=False,
    #     font=dict(size=12),
    #     align="center"
    # )

    # # Continuous UMAP
    # fig_continuous = px.scatter(
    #     df, x="UMAP1", y="UMAP2",
    #     color="viral_counts",
    #     color_continuous_scale="Reds",
    #     hover_data=["barcode","virus_positive"],
    #     title="UMAP of Viral Load per Cell"
    # )
    # fig_continuous.update_layout(
    #     xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #     yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #     plot_bgcolor="white",
    #     margin=dict(b=100)
    # )

    # # Add annotation to show the conclusion text
    # fig_continuous.add_annotation(
    #     text=conclusion_text,
    #     xref="paper", yref="paper",
    #     x=0.5, y=-0.15,
    #     showarrow=False,
    #     font=dict(size=12),
    #     align="center"
    # )

    # # Save plots to the users output directory
    # outdir = f"{config['output']}/plots"
    # os.makedirs(outdir, exist_ok=True)
    # fig_binary.write_html(f"{outdir}/umap_binary.html")
    # fig_continuous.write_html(f"{outdir}/umap_continuous.html")


# def viral_neighbor_enrichment(coords, labels, k, n_permutations=1000):
#     """
#     Compute nearest-neighbor enrichment for viral-positive cells to
#     determine whether there is a connection between the place in the 
#     plot and the viral reads.
#     ---------------------------------------------------------------------
#     Params:
#         coords (numpy.ndarray): the UMAP embedding
#         labels (numpy.ndarray): int value showing whether cell is virus 
#             positive (1) or not (0) 
#         k (int): the amount of neighbors
#         n_permutations (int): total amount of permutations for the 
#             enrichment
#     ---------------------------------------------------------------------
#     Returns:
#         observed (float): average fraction of neighbors which are virus-
#             positive using real labels
#         expected (float): average fraction of virus-positive neighbors 
#             with random labelling
#         p_value (float): probability of observing the results
#     """
#     # Compute nearest neighbor
#     nbrs = NearestNeighbors(n_neighbors=k).fit(coords)
#     distances, indices = nbrs.kneighbors(coords)

#     # Check for other viral cells -> compute observed fraction
#     viral_cells = np.where(labels == 1)[0] 
#     counts = []
#     for i in viral_cells:
#         neighbor_idx = indices[i][1:]
#         counts.append(np.mean(labels[neighbor_idx]))
#     observed = np.mean(counts)

#     # Null distribution -> compute expected fraction
#     permuted = []
#     for _ in range(n_permutations):
#         shuffled = np.random.permutation(labels)
#         counts_perm = []
#         for i in viral_cells:
#             neighbor_idx = indices[i][1:]
#             counts_perm.append(np.mean(shuffled[neighbor_idx]))
#         permuted.append(np.mean(counts_perm))
#     expected = np.mean(permuted)

#     p_value = (np.sum(np.array(permuted) >= observed) + 1) / (n_permutations + 1)
#     return observed, expected, p_value

def viral_neighbor_enrichment(coords, labels, k, n_permutations=1000):
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
        print(f"You have decided to create a umap. This can take a while before finishing the code. Please wait...") 
        umap(adata, found_genes)

main()

# write to output file for Snakemake
with open(snakemake.output[0], "w") as f: 
    f.write("done\n") 
    if config["umap"] == "True":
        print("\033[32mUmap is done!\033[0m") 
    print(f"All (important) results of ViralScan can be found in {config["output"]}summary.txt")