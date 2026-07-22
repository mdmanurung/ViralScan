# Importing packages
import gzip
import os
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import pandas as pd

from viralscan.multimapping import build_multimap_layers
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig

# Run-level state, populated by run() from the Run Context. Declared here so the
# helper functions can reference them as module globals; the module imports
# cleanly without Snakemake because nothing reads these at import time.
config: RunConfig = RunConfig()
output: str = ""
kb = None


def strip_10x_suffix(barcode: str) -> str:
    """Remove only the trailing '-1' lane suffix added by 10x Cell Ranger.

    A global ``str.replace("-1", "")`` would corrupt any barcode containing
    '-1' at a non-trailing position; ``removesuffix`` strips just the terminal one.
    """
    return barcode.removesuffix("-1")


def define_paths():
    """
    Define the paths to read for the rest of the code
    ---------------------------------------------------------------------
    Returns:
        str: all paths are in the form of a string

    The kb-python layout is owned by :class:`viralscan.kb_outputs.KbCountOutputs`;
    this just unpacks it (as strings) for the existing call sites. ``t2g_file``
    is the config-supplied t2g, not a kb-python product.
    """
    return (
        str(kb.adata),
        str(kb.bus),
        str(kb.ec),
        str(kb.transcripts_txt),
        str(kb.barcodes),
        str(kb.resolved_bus_txt),
        str(kb.genes),
        str(kb.gene_names),
        config.transcripts,
    )


def prepare_resolved_bus(
    raw_bus: str | Path,
    resolved_bus: str | Path,
    resolved_text: str | Path,
    *,
    whitelist: str | None,
    threads: int,
    corrected_bus: str | Path | None = None,
) -> None:
    """Create the retained corrected/sorted BUS boundary required by v3.

    Kallisto's ``output.bus`` is raw and unsorted. If the run supplied an
    on-list, correction is repeated explicitly so the exact corrected BUS used
    for molecule resolution is retained. Sorting is mandatory in both cases.
    Every tool writes to a staging path before the public artifact is replaced.
    """
    raw_path = Path(raw_bus)
    sorted_path = Path(resolved_bus)
    text_path = Path(resolved_text)
    corrected_path = Path(corrected_bus) if corrected_bus else raw_path.with_name(
        "output.corrected.bus"
    )
    sorted_path.parent.mkdir(parents=True, exist_ok=True)
    sort_input = raw_path
    plain_whitelist: Path | None = None

    try:
        if whitelist:
            whitelist_path = Path(whitelist)
            tool_whitelist = whitelist_path
            if whitelist_path.suffix == ".gz":
                plain_whitelist = sorted_path.with_name("v3_whitelist.txt")
                with gzip.open(whitelist_path, "rb") as source, plain_whitelist.open(
                    "wb"
                ) as target:
                    shutil.copyfileobj(source, target)
                tool_whitelist = plain_whitelist
            corrected_stage = corrected_path.with_suffix(corrected_path.suffix + ".tmp")
            subprocess.run(
                [
                    "bustools",
                    "correct",
                    "-w",
                    str(tool_whitelist),
                    "-o",
                    str(corrected_stage),
                    str(raw_path),
                ],
                check=True,
            )
            corrected_stage.replace(corrected_path)
            sort_input = corrected_path

        sorted_stage = sorted_path.with_suffix(sorted_path.suffix + ".tmp")
        subprocess.run(
            [
                "bustools",
                "sort",
                "-t",
                str(max(1, int(threads))),
                "-o",
                str(sorted_stage),
                str(sort_input),
            ],
            check=True,
        )
        sorted_stage.replace(sorted_path)

        text_stage = text_path.with_suffix(text_path.suffix + ".tmp")
        subprocess.run(
            ["bustools", "text", "-o", str(text_stage), str(sorted_path)],
            check=True,
        )
        text_stage.replace(text_path)
    finally:
        if plain_whitelist is not None:
            plain_whitelist.unlink(missing_ok=True)


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
    barcodes = [strip_10x_suffix(bc) for bc in barcodes]
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
    t2g = pd.read_csv(
        t2g_file, sep=r"\s+", header=None, usecols=[0, 1], names=["transcript", "gene"]
    )
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
    gene_id_to_idx: dict = {gid: i for i, gid in enumerate(gene_ids)}
    with open(ec_file) as f:
        for _i, line in enumerate(f):
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            ec_id = int(parts[0])
            transcript_indices = [int(x) for x in parts[1].split(",") if x.isdigit()]
            transcript_ids = [transcripts[j] for j in transcript_indices if j < len(transcripts)]
            gene_ids_ec = [t2g_map.get(tr) for tr in transcript_ids if tr in t2g_map]
            gene_indices = [gene_id_to_idx[gid] for gid in gene_ids_ec if gid in gene_id_to_idx]
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
    # Strip only the trailing '-1' lane suffix (avoid global replace that
    # would corrupt barcodes with an internal '-1' substring). When barcode is a
    # category, rewrite the (few) distinct labels rather than every row; fall back
    # to a per-value map only if stripping collides two labels into one.
    barcode = bus_df["barcode"]
    if isinstance(barcode.dtype, pd.CategoricalDtype):
        stripped = barcode.cat.categories.map(strip_10x_suffix)
        if stripped.is_unique:
            bus_df["barcode"] = barcode.cat.rename_categories(stripped)
        else:
            bus_df["barcode"] = barcode.astype("string").map(strip_10x_suffix)
    else:
        bus_df["barcode"] = barcode.map(strip_10x_suffix)
    # ec is already int from the typed read; no nullable-Int64 recast needed.

    viral_ids_file = os.path.join(output, "log", "analysis.txt")
    viral_gene_indices = set()

    if viral_ids_file and os.path.exists(viral_ids_file):
        with open(viral_ids_file) as f:
            viral_gene_ids = {line.strip() for line in f}
        viral_gene_indices = {i for i, gid in enumerate(gene_ids) if gid in viral_gene_ids}
    return bus_df, viral_gene_indices


def create_new_h5ad(
    corrected_matrix,
    adata_orig,
    gene_ids,
    gene_names,
    genes_from_matrix,
    viral_gene_indices,
    n_genes,
):
    """
    Creating the new h5ad from the corrected matrix.
    """
    adata = ad.AnnData(
        X=corrected_matrix,
        obs=adata_orig.obs.copy(),
        var=pd.DataFrame({"gene_id": gene_ids, "gene_name": gene_names}, index=genes_from_matrix),
    )

    # Mark viral genes
    is_viral = [i in viral_gene_indices for i in range(n_genes)]
    adata.var["is_viral"] = is_viral

    viral_counts = adata[:, adata.var["is_viral"]].X.sum(axis=1).A1
    return adata, viral_counts


def final_results(viral_counts, adata_orig, viral_gene_indices, adata, n_cells, layers):
    """
    Adding the final data to the adata file and write conclusions to the summary file.
    """
    cells_with_virus = (viral_counts > 0).sum()
    total_viral_molecules = viral_counts.sum()

    # V3 count contract: X is the complete selected-method molecule matrix and
    # these two non-overlapping layers sum to it. Legacy layer names remain only
    # as explicit aliases during the development cycle.
    adata.layers["counts_unique"] = layers.unique
    adata.layers["counts_ambiguous_allocated"] = layers.corrected
    adata.X = adata.layers["counts_unique"] + adata.layers["counts_ambiguous_allocated"]

    # The `layers` object is not used after this function,
    # so assign its sparse matrices directly instead of duplicating each one with
    # .copy() (8 extra full-size sparse copies = a major peak-RSS spike on deep
    # samples). adata.var_names is, by construction, adata_orig.var_names in the
    # same order, so the v3 molecule layers need no reindex/copy.
    adata.layers["counts_corrected"] = layers.corrected
    adata.layers["counts_original"] = layers.unique
    adata.layers["counts_combined"] = adata.X
    adata.layers["counts_multimap_equal"] = layers.equal
    adata.layers["counts_multimap_host_conservative"] = layers.host_conservative
    adata.layers["counts_multimap_unique_weighted"] = layers.unique_weighted
    adata.layers["counts_unique_viral"] = layers.unique_viral
    adata.layers["counts_host_viral_ambiguous"] = layers.host_viral_ambiguous
    adata.layers["counts_host_viral_selected"] = layers.host_viral_selected
    adata.layers["counts_viral_ambiguous_upper"] = layers.viral_ambiguous_upper
    adata.uns["multimap_method"] = config.multimap_method
    adata.uns["multimap_pseudocount"] = config.multimap_pseudocount
    adata.uns["count_schema_version"] = "3.0.0"
    adata.uns["quantification_unit"] = "bustools-resolved-cb-umi-molecule"
    adata.uns["molecule_audit"] = {
        "input_molecules": layers.audit.input_molecules,
        "resolved_molecules": layers.audit.resolved_molecules,
        "unique_molecules": layers.audit.unique_molecules,
        "ambiguous_molecules": layers.audit.ambiguous_molecules,
        "unresolved_molecules": layers.audit.unresolved_molecules,
        "ignored_read_multiplicity": layers.audit.ignored_read_multiplicity,
        "allocated_ambiguous_mass": float(layers.corrected.sum()),
    }
    adata.uns["multimap_diagnostics"] = layers.method_diagnostics

    output_file = str(kb.adata_multimap)
    adata.write(output_file)

    audit = adata.uns["molecule_audit"]
    pd.DataFrame([audit]).to_csv(f"{config.output}/count_audit.tsv", sep="\t", index=False)

    with open(f"{config.output}/summary.txt", "w") as summary:
        summary.write(
            "Viral molecules in unique-count matrix: "
            f"{layers.unique[:, list(viral_gene_indices)].sum()}\n"
        )
        summary.write(f"Total viral molecules (selected method): {total_viral_molecules}\n")
        summary.write(f"Cells with viral reads: {cells_with_virus}/{n_cells}\n\n\n")


def run(ctx, done_file):
    """Entry point: build multimapper layers for one Run, then touch done_file."""
    global config, output, kb
    config = ctx.config
    output = config.output
    kb = ctx.outputs

    if config.multimapping:
        (
            adata_file,
            bus_file,
            ec_file,
            transcript_file,
            barcodes_file,
            txt_file,
            genes_file,
            gene_names_file,
            t2g_file,
        ) = define_paths()

        # Materialize the exact corrected/sorted BUS boundary. Raw output.bus
        # is neither corrected nor sorted and is never a valid v3 count input.
        prepare_resolved_bus(
            bus_file,
            kb.resolved_bus,
            txt_file,
            whitelist=config.whitelist,
            threads=config.cores,
            corrected_bus=kb.corrected_bus,
        )

        # Load all data
        adata_orig, genes_from_matrix, n_genes = load_adata(adata_file)
        barcode_to_idx, n_cells = load_barcodes(barcodes_file)
        gene_ids, gene_names = load_genes(genes_file, gene_names_file, n_genes)
        transcripts, t2g_map = load_transcripts(transcript_file, t2g_file)

        # Continue with workflow
        ec_map = read_ec(ec_file, transcripts, t2g_map, gene_ids)
        # Stream corrected, sorted BUS text. Retaining the UMI is the v3 count
        # boundary; the read-multiplicity column is audit-only.
        viral_ids_file = os.path.join(output, "log", "analysis.txt")
        viral_gene_indices: set[int] = set()
        if os.path.exists(viral_ids_file):
            with open(viral_ids_file) as handle:
                viral_gene_ids = {line.strip() for line in handle}
            viral_gene_indices = {i for i, gid in enumerate(gene_ids) if gid in viral_gene_ids}
        layers = build_multimap_layers(
            bus_df=Path(txt_file),
            barcode_to_idx=barcode_to_idx,
            ec_map=ec_map,
            n_cells=n_cells,
            n_genes=n_genes,
            viral_gene_indices=viral_gene_indices,
            original_counts=adata_orig.X,
            method=config.multimap_method,
            pseudocount=config.multimap_pseudocount,
            em_max_iter=config.multimap_em_max_iter,
            em_tol=config.multimap_em_tol,
        )
        corrected_matrix = layers.unique + layers.corrected
        adata, viral_counts = create_new_h5ad(
            corrected_matrix,
            adata_orig,
            gene_ids,
            gene_names,
            genes_from_matrix,
            viral_gene_indices,
            n_genes,
        )
        final_results(viral_counts, adata_orig, viral_gene_indices, adata, n_cells, layers)

    with open(done_file, "w") as f:
        f.write("done\n")


if "snakemake" in globals():
    run(
        RunContext.from_yaml(snakemake.params.configfile),  # noqa: F821 (snakemake magic global)
        snakemake.output[0],  # noqa: F821
    )
