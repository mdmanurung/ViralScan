"""
The detection scripts uses the gene IDs found in the analysis script and checks
for presence of the gene IDs in the sample (according to kb count). It also
Creates visualizations (except for the umap), e.g. a barplot and a plot showing
super expressors.
"""

# Importing packages
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse as sparse
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
import os

# Reading Snakefile parameters
file = snakemake.input.file_viral_accessions
snakefile_dir = snakemake.params.snakefile_dir
configfile = snakemake.params.configfile

# Reading config
with open(configfile, 'r') as f:
    config = yaml.safe_load(f)

output = config['output']

def preprocessing():
    """
    This function checks which viral gene IDs have been found in the 
    sample according to the accession list from analysis.py.
    ---------------------------------------------------------------------
    Returns:
        adata (anndata.AnnData): h5ad file of kb-python used for further
            analysis
        found_genes (dict): dictionary containing information of the gene 
            IDs found and the gene counts
        output (str): the path to the output directory defined by the user
    """
    viral_accessions = list()
    viral_file = open(file)
    for f in viral_file:
        viral_accessions.append(f.strip())

    # if os.path.exists(f"{snakefile_dir}/{output}"):
    #     outputpath = f"{snakefile_dir}/{output}"
    # else:
    #     outputpath = output

    adata = sc.read_h5ad(f"{output}/kb-python/counts_unfiltered/adata.h5ad")

    threshold = 1 # TODO: change threshold. (<5 = likely noise), do this after multimapping!
    all_gene_ids = adata.var_names
    found_genes = {}  
    
    # Detect viral gene IDs in samples
    for gene_id in viral_accessions:
        if gene_id in all_gene_ids:
            gene_counts = adata[:, gene_id].X
            if hasattr(gene_counts, 'toarray'):
                gene_counts = gene_counts.toarray()
            total_count = gene_counts.sum()
            if total_count > threshold:
                found_genes[gene_id] = total_count
    return adata, found_genes, output


def histogram(adata, found_genes, map_virus, outputpath):
    """
    This function creates a histogram showing the gene IDs found sorted 
    on the UMI counts.
    ---------------------------------------------------------------------
    Returns:
        adata (anndata.AnnData): h5ad file of kb-python used for further
            analysis
        found_genes (dict): dictionary containing information of the gene 
            IDs found and the gene counts
        output (str): the path to the output directory defined by the user

    """
    if sparse.issparse(adata.X):
        gene_counts = np.array(adata.X.sum(axis=0)).flatten()
    else:
        gene_counts = adata.X.sum(axis=0)

    # Create dataframe with gene IDs and UMI counts
    df = pd.DataFrame({
        'gene_id': adata.var_names,
        'UMI_count': gene_counts
    })

    # Group versions of found viruses in dictionary
    group_by_virus = {}
    detected_viral_genes = set()

    for g in found_genes:
        added = False
        for key in map_virus:
            if key in g:
                detected_viral_genes.add(map_virus[key])
                if map_virus[key] in group_by_virus:
                    group_by_virus[map_virus[key]].append(g) 
                else:
                    group_by_virus[map_virus[key]] = [g] 
                added = True
                break

        # if it is not in the mapping, still add the abbreviation for the entirity of the results.
        if not added:
            if g in group_by_virus:
                group_by_virus[g].append(g)
            else:
                group_by_virus[g] = [g]

    # Check if user wants visualizations
    if config["visual"]:
        for virus in group_by_virus:
            virus_list = group_by_virus[virus]
            df_virus = df[df['gene_id'].isin(virus_list)]
            df_virus_sorted = df_virus.sort_values(by='UMI_count', ascending=False).reset_index(drop=True)
            
            # Plot the UMI counts
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(
                data=df_virus_sorted,
                x='gene_id',
                y='UMI_count'
            )

            # Check for amount of bars before annotating (max of 10 for annotation)
            if len(ax.patches) <= 10:
                for p in ax.patches:
                    ax.annotate(
                        f'{p.get_height():.0f}',  
                        (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='bottom',
                        fontsize=11, color='black',
                        xytext=(0, 3), 
                        textcoords='offset points'
                    )

            # Check if path exists, otherwise create it.
            os.makedirs(f"{outputpath}/plots/", exist_ok=True)

            # Pass variables for bar plot
            plt.xticks(rotation=45, ha='right')
            plt.ylabel('UMI Count')
            plt.xlabel('Gene ID')
            plt.title(f'{virus} Gene UMI Counts (Bar Plot)')
            plt.tight_layout()
            plt.savefig(f"{outputpath}/plots/{virus}_histogram.png")
            plt.close()
    return group_by_virus, detected_viral_genes


def super_expressor(adata, virus, viral_gene_ids, outputpath):
    """
    This function creates a super expressor plot, showing how many data points (single
    cells) have a UMI count > 10. 
    ---------------------------------------------------------------------
    Parameters:
        adata (anndata.AnnData): h5ad file of kb-python used for further
            analysis
        virus (str): virus where plot needs to be made for
        viral_gene_ids (list): gene IDs belonging to the virus param
        outputpath (str): path to the output directory defined by the user
    ---------------------------------------------------------------------
    Raises:
        ValueError: None of the provided viral gene IDs were found in the dataset
    """
    adata.var_names_make_unique()

    # Compute total UMI per cell
    cell_umi_counts = adata.X.sum(axis=1).A1
    adata.obs['nCount_RNA'] = cell_umi_counts

    # Match viral gene IDs to adata and raise ValueError
    viral_mask = adata.var_names.isin(viral_gene_ids)
    matched_genes = adata.var_names[viral_mask]
    if matched_genes.empty:
        raise ValueError("None of the provided viral gene IDs were found in the dataset. No super expressor is therefore found.")

    # Compute viral UMI counts per cell
    adata.obs[virus] = adata[:, viral_mask].X.sum(axis=1).A1

    # Create null model (grey dots)
    total_viral = adata.obs[virus].sum()
    total_rna = adata.obs['nCount_RNA'].sum()
    null_model = total_viral * (adata.obs['nCount_RNA'] / total_rna)

    df_plot = pd.DataFrame({
        'rank': np.arange(adata.n_obs),
        'observed': np.sort(adata.obs[virus].values)[::-1],
        'null': np.sort(null_model.values)[::-1]
    })

    # Count super-expressors (UMI > 10)
    n_SE = (adata.obs[virus] >= 10).sum()
    title = (
        f"Virus {virus}, n_super={n_SE}\n"
        f"n={adata.n_obs}; {virus} max={int(adata.obs[virus].max())}"
    )
    df_long = pd.melt(df_plot, id_vars='rank', value_vars=['observed', 'null'])

    # Start creating plot
    plt.figure(figsize=(5, 4))
    sns.scatterplot(
        data=df_long, x='rank', y='value', hue='variable',
        palette={'observed': 'firebrick', 'null': 'grey'}, s=10
    )
    plt.axhline(10, linestyle='--', color='darkblue', linewidth=1)

    # Pass params for visualization
    plt.yscale('log')
    plt.title(title, fontsize=10, wrap=True)
    plt.xlabel('Cell Rank', fontsize=10)
    plt.ylabel('Viral Counts', fontsize=0)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.legend([], [], frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save plot to output directory
    plt.savefig(f"{outputpath}/plots/nofilter_SuperExpressor_{virus}.png", dpi=500)
    plt.close()


def detect_cells(adata, found_genes, sum):  
    """
    This function detects in which cells (barcodes) the viral genes
    have been found and writes this to the summary in the output
    directory
    ---------------------------------------------------------------------
    Parameters:
        adata (anndata.AnnData): h5ad file of kb-python used for further
            analysis
        found_genes (dict): dictionary containing information of the gene 
            IDs found and the gene counts
        sum (IO[str]): open text file to write the summary to
    """
    # Detect cells and find barcodes for gene IDs
    cells_per_gene = {}
    for viral_gene_name in found_genes:
        gene_counts = adata[:, viral_gene_name].X
        if hasattr(gene_counts, "toarray"):
            gene_counts = gene_counts.toarray()
        expressed_mask = gene_counts.flatten() > 0
        cells_with_gene = adata.obs_names[expressed_mask].tolist()
        cells_per_gene[viral_gene_name] = cells_with_gene

    for gene, barcodes in cells_per_gene.items():
        sum.write(f"{gene} detected in {len(barcodes)} cells. ")
        sum.write(f"Barcodes: {barcodes}\n")
    sum.close()


def main():
    # Dictionary containing all viral load present in the index
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
    adata, found_genes, outputpath = preprocessing()

    # check if user wants visuals in output directory
    group_by_virus, detected_viral_genes = histogram(adata, found_genes, map_virus, outputpath)
    if config["visual"]:
        for virus in group_by_virus:
            super_expressor(adata, virus, group_by_virus[virus], outputpath)

    # Writing results to the summary file
    total_viral_genes = 0
    summary = open(f"{config["output"]}/summary.txt", "w")
    summary.write("Gene ID; Gene Count\n")
    found_genes_file = open(f"{config['output']}log/found_genes.txt", 'w')
    if len(found_genes) > 0:
        for g in found_genes:
            write_to_file = f"{g};{found_genes[g]}\n"
            total_viral_genes += found_genes[g]
            summary.write(write_to_file)
            found_genes_file.write(write_to_file)
        found_genes_file.close()
        summary.write(f"\n\nTotal amount of viral load found: {total_viral_genes}")
        summary.write(f"\n\nOfficial name of viral load detected: {','.join(str(s) for s in detected_viral_genes)}")
    else:
        summary.write("No viral gene IDs found in this sample for the 99 viruses in the index.")
    summary.write(f"\nIf you want to see the cell gene matrix, go to the kb-python/counts_unfiltered/ folder and look for the cells_x_genes.mtx file.\n")
    detect_cells(adata, found_genes, summary)
    

main()

with open(snakemake.output[0], "w") as f:
    f.write("done\n")
print("\033[32mThe detection and visualizations are done!\033[0m")