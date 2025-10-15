# Importing packages
import os
import yaml
import subprocess
import pandas as pd
import anndata as ad
from scipy import sparse

# Get Snakefile params
snakefile_dir = snakemake.params.snakefile_dir
configfile = snakemake.params.configfile

# Read config and define path output
with open(configfile, 'r') as f:
    config = yaml.safe_load(f)
output = config['output']    

def define_paths():
    """
    Define the paths to read for the rest of the code
    ---------------------------------------------------------------------
    Returns:
        str: all paths are in the form of a string
    """
    adata_file = f"{output}/kb-python/counts_unfiltered/adata.h5ad"
    bus_file = f"{output}/kb-python/output.bus"
    ec_file = f"{output}/kb-python/matrix.ec"
    transcript_file = f"{output}/kb-python/transcripts.txt"
    barcodes_file = f"{output}/kb-python/counts_unfiltered/cells_x_genes.barcodes.txt"
    txt_file = f"{output}/kb-python/output.bus.txt"
    genes_file = f"{output}/kb-python/counts_unfiltered/cells_x_genes.genes.txt"
    gene_names_file = f"{output}/kb-python/counts_unfiltered/cells_x_genes.genes.names.txt"
    t2g_file = f"{config['transcripts']}"
    return adata_file, bus_file, ec_file, transcript_file, barcodes_file, txt_file, \
        genes_file, gene_names_file, t2g_file


def load_barcodes(barcodes_file):
    """
    Load the file which contain the barcodes
    ---------------------------------------------------------------------
    Returns:
        barcode_to_idx (dict): dictionary containing information about 
            barcodes
        n_cells (int): the amount of barcodes in the file
    """
    with open(barcodes_file) as f:
        barcodes = [line.strip() for line in f]
    barcodes = [bc.replace("-1", "") for bc in barcodes]
    barcode_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    n_cells = len(barcodes)
    return barcode_to_idx, n_cells


def load_adata(adata_file):
    """
    Load the adata file and get the amount of genes.
    ---------------------------------------------------------------------
    Params:
        adata_file (str): path to the adata file
    ---------------------------------------------------------------------
    Returns:
        adata_orig (anndata.AnnData): the 'original' h5ad (unaltered)
        genes_from_matrix (list): list of genes from h5ad file
        n_genes (int): the total amount of genes
    """
    adata_orig = ad.read_h5ad(adata_file)
    genes_from_matrix = list(adata_orig.var_names)
    n_genes = len(genes_from_matrix)
    return adata_orig, genes_from_matrix, n_genes


def load_genes(genes_file, gene_names_file, n_genes):
    """
    Getting gene IDs and names from the genes file.
    ---------------------------------------------------------------------
    Params:
        genes_file (str): path to the genes file
        gene_names_file (str): path to the gene names file from kb-python
        n_genes (int): the total amount of genes
    ---------------------------------------------------------------------
    Returns:
        gene_ids (list): list containing gene IDs
        gene_names (list): list containing gene names
    """
    with open(genes_file) as f:
        gene_ids = [line.strip() for line in f]
    with open(gene_names_file) as f:
        gene_names = [line.strip() for line in f]
    assert len(gene_ids) == n_genes
    assert len(gene_names) == n_genes
    return gene_ids, gene_names


def load_transcripts(transcript_file, t2g_file):
    """
    Load the transcripts file and the transcripts to genes files.
    ---------------------------------------------------------------------
    Params:
        transcripts_file (str): path to the transcripts file (from kb count)
        t2g_file (str): path to transcripts to genes file (from kb ref)
    ---------------------------------------------------------------------
    Returns:
        transcripts (list): list containing transcripts
        t2g_map (dict): dictionary containing transcripts as key and gene as value
    """
    with open(transcript_file) as f:
        transcripts = [line.strip() for line in f]

    # Load t2g mapping
    t2g = pd.read_csv(t2g_file, sep=r"\s+", header=None, usecols=[0,1], names=["transcript", "gene"])
    t2g_map = dict(zip(t2g["transcript"], t2g["gene"]))
    return transcripts, t2g_map


def read_ec(ec_file, transcripts, t2g_map, gene_ids):
    """
    Read the EC file (from kb count) and create a mapping.
    ---------------------------------------------------------------------
    Params:
        ec_file (str): read the matrix file (from kb count)
        transcripts (list): list containing the transcripts (from kb count)
        t2g_map (dict): dictionary containing transcripts as key and gene as value
        gene_ids (list): list containing gene IDs
    ---------------------------------------------------------------------
    Returns:
        ec_map (dict): dictionary containing EC IDs as key and gene indices as key
    """
    ec_map = {}
    with open(ec_file) as f:
        for i, line in enumerate(f):
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            ec_id = int(parts[0])
            transcript_indices = [int(x) for x in parts[1].split(",") if x.isdigit()]
            transcript_ids = [transcripts[j] for j in transcript_indices if j < len(transcripts)]
            gene_ids_ec = [t2g_map.get(tr) for tr in transcript_ids if tr in t2g_map]
            gene_indices = [gene_ids.index(gid) for gid in gene_ids_ec if gid in gene_ids]
            if gene_indices:
                ec_map[ec_id] = gene_indices
    return ec_map


def normalize_barcodes(bus_df, gene_ids):
    """
    Normalize the barcodes to not have empty values and get the
    viral IDs.
    ---------------------------------------------------------------------
    Params:
        bus_df (pd.DataFrame): DataFrame from output.bus.txt from kb count
        gene_ids (list): list containing gene IDs
    ---------------------------------------------------------------------
    Returns:
        bus_df (pd.DataFrame): DataFrame from output.bus.txt from kb count
        viral_gene_indices (dict): dictionary of viral genes including ID
    """
    bus_df["barcode"] = bus_df["barcode"].str.replace("-1", "", regex=False)
    bus_df["ec"] = pd.to_numeric(bus_df["ec"], errors="coerce").astype("Int64")

    viral_ids_file = os.path.join(output, "log", "analysis.txt")
    viral_gene_indices = set()

    if viral_ids_file and os.path.exists(viral_ids_file):
        with open(viral_ids_file) as f:
            viral_gene_ids = {line.strip() for line in f}
        viral_gene_indices = {i for i, gid in enumerate(gene_ids) if gid in viral_gene_ids}      
    return bus_df, viral_gene_indices


def build_multimap_matrix(bus_df, barcode_to_idx, ec_map, n_cells, n_genes):
    """
    Building the new multimap matrix to (eventually) write to h5ad file.
    ---------------------------------------------------------------------
    Params:
        bus_df (pd.DataFrame): DataFrame from output.bus.txt from kb count
        barcode_to_idx (dict): dictionary containing information about 
            barcodes
        ec_map (dict): dictionary containing EC IDs as key and gene indices as key
        n_cells (int): the total amount of cells
        n_genes (int): the total amount of genes
    ---------------------------------------------------------------------
    Returns:
        corrected_matrix (scipy.sparse.scr_matrix): corrected matrix including 
            data about multimaps to write to h5ad file.
    """
    rows, cols, data = [], [], []
    skipped_no_barcode = 0
    skipped_no_ec = 0

    for idx, row in bus_df.iterrows():
        bc, ec, count = row["barcode"], row["ec"], row["count"]
        if pd.isna(ec):
            skipped_no_ec += 1
            continue
        ec = int(ec)
        if bc not in barcode_to_idx:
            skipped_no_barcode += 1
            continue
        if ec not in ec_map:
            skipped_no_ec += 1
            continue

        cell_idx = barcode_to_idx[bc]
        genes_in_ec = ec_map[ec]
        if not genes_in_ec:
            continue

        share = count / len(genes_in_ec)
        for gid in genes_in_ec:
            rows.append(cell_idx)
            cols.append(gid)
            data.append(share)

    corrected_matrix = sparse.csr_matrix((data, (rows, cols)), shape=(n_cells, n_genes))
    return corrected_matrix


def create_new_h5ad(corrected_matrix, adata_orig, gene_ids, gene_names, genes_from_matrix, viral_gene_indices, n_genes):
    """
    Creating the new h5ad from the corrected matrix.
    """
    adata = ad.AnnData(
        X=corrected_matrix,
        obs=adata_orig.obs.copy(),
        var=pd.DataFrame({"gene_id": gene_ids, "gene_name": gene_names}, index=genes_from_matrix)
    )

    # Mark viral genes
    is_viral = [i in viral_gene_indices for i in range(n_genes)]
    adata.var["is_viral"] = is_viral

    viral_counts = adata[:, adata.var["is_viral"]].X.sum(axis=1).A1
    return adata, viral_counts


def final_results(viral_counts, adata_orig, viral_gene_indices, adata, n_cells):
    """
    Adding the final data to the adata file and write conclusions to the summary file.
    """
    cells_with_virus = (viral_counts > 0).sum()
    total_viral_umis = viral_counts.sum()

    # Extract original viral UMIs per cell
    viral_counts_orig = adata_orig[:, list(viral_gene_indices)].X
    if sparse.issparse(viral_counts_orig):
        viral_counts_orig = viral_counts_orig.toarray()
    
    # Save plots
    adata.layers["counts_corrected"] = adata.X.copy()
    adata.layers["counts_original"] = adata_orig[:, adata.var_names].X.copy()
    adata.layers["counts_combined"] = adata.layers["counts_corrected"] + adata.layers["counts_original"]

    output_file = f"{output}/kb-python/counts_unfiltered/adata_multimap.h5ad"
    adata.write(output_file)

    summary = open(f"{config["output"]}/summary.txt", "w")
    summary.write(f"Viral UMIs in original (not corrected) adata: {adata_orig[:, list(viral_gene_indices)].X.sum()}\n")
    summary.write(f"Total viral UMIs (corrected): {total_viral_umis}\n")
    summary.write(f"Cells with viral reads: {cells_with_virus}/{n_cells}\n\n\n")
    summary.close()


def main():
    if config["multimapping"]:
        adata_file, bus_file, ec_file, transcript_file, barcodes_file, txt_file, \
            genes_file, gene_names_file, t2g_file = define_paths()
        
        # Convert BUS file to text
        subprocess.run(["bustools", "text", "-o", txt_file, bus_file], check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        # Load all data
        adata_orig, genes_from_matrix, n_genes = load_adata(adata_file)
        barcode_to_idx, n_cells = load_barcodes(barcodes_file)
        gene_ids, gene_names = load_genes(genes_file, gene_names_file, n_genes)
        transcripts, t2g_map = load_transcripts(transcript_file, t2g_file)

        # Continue with workflow
        ec_map = read_ec(ec_file, transcripts, t2g_map, gene_ids)
        bus_df = pd.read_csv(txt_file, sep="\t", header=None, names=["barcode", "umi", "ec", "count"])
        bus_df, viral_gene_indices = normalize_barcodes(bus_df, gene_ids)
        corrected_matrix = build_multimap_matrix(bus_df, barcode_to_idx, ec_map, n_cells, n_genes)
        adata, viral_counts = create_new_h5ad(corrected_matrix, adata_orig, gene_ids, gene_names, genes_from_matrix, viral_gene_indices, n_genes)
        final_results(viral_counts, adata_orig, viral_gene_indices, adata, n_cells)

main()

with open(snakemake.output[0], "w") as f:
    f.write("done\n")