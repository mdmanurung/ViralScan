"""
The umap script creates a umap based on the single cell data provided where
each data point represents a single cell. Based on the colour it shows whether
the viral load has been detected in that single cell.
"""

# Importing packages
import csv
import logging
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

try:
    import plotly.express as px
except ModuleNotFoundError:  # plotly is only needed to render the UMAP HTML plots
    px = None
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

from viralscan.kb_outputs import KbCountOutputs
from viralscan.multimapping import select_detection_matrix
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scautoqc import add_scautoqc
from viralscan.utils import infected_masks, setup_script_logging, virus_molecules_per_cell
from viralscan.virus_catalog import merged_name_map
from viralscan.virus_grouping import detected_summary_rows, load_run_identity, virus_name_for_gene

log = setup_script_logging()

warnings.filterwarnings("ignore")

# Run-level state, populated by run() from the Run Context. Declared here so the
# helper functions can reference them as module globals; the module imports
# cleanly without Snakemake because nothing reads these at import time.
config: RunConfig = RunConfig()
kb: KbCountOutputs = KbCountOutputs.from_config_output(config.output)
# The Run's Virus Identity table; None (a run directory without one) selects the
# legacy prefix naming.
identity = None


def _gene_to_virus(gene_ids, table=None):
    """gene_id -> virus display name, from the identity table (legacy prefix if none)."""
    if table is not None:
        by_gene = table.by_gene()
        return {g: by_gene[g].virus_name if g in by_gene else g for g in gene_ids}
    name_map = merged_name_map()
    return {g: virus_name_for_gene(g, name_map) for g in gene_ids}


def detected_genes(found_genes, output, table=None):
    """Restrict candidate genes to the summary's detected viruses (AUD-21)."""
    summary = Path(output) / "results" / "viral_summary.tsv"
    if not summary.is_file():
        if table is not None:
            raise FileNotFoundError(f"Viral summary not found: {summary}. Rerun detection first.")
        log.warning("Legacy run has no viral summary; using candidate genes for UMAP labels.")
        return found_genes
    with summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "virus_name" not in reader.fieldnames:
            raise ValueError(f"Viral summary {summary} has no virus_name column.")
        names = {row["virus_name"] for row in detected_summary_rows(list(reader))}
    by_gene = _gene_to_virus(found_genes, table)
    return {gene: count for gene, count in found_genes.items() if by_gene[gene] in names}


def label_infected_cells(adata, matrix, viral_ids, table, min_molecules):
    """Per-cell virus labels and total viral molecules, by the one infected rule (AUD-06).

    *matrix* is the selected-method matrix detection used (``adata.X``), so the UMAP
    cannot disagree with ``viral_summary.tsv``. A cell is labelled with every virus whose
    summed molecules across its genes reach *min_molecules*.
    """
    gene_to_virus = _gene_to_virus(viral_ids, table)
    genes_by_virus: dict = {}
    for g in viral_ids:
        genes_by_virus.setdefault(gene_to_virus[g], []).append(g)
    positive = {
        virus: infected_masks(virus_molecules_per_cell(adata, matrix, genes), min_molecules)[0]
        for virus, genes in genes_by_virus.items()
    }
    labels = []
    for i in range(adata.n_obs):
        hit = sorted(virus for virus, mask in positive.items() if mask[i])
        if not hit:
            labels.append("No Virus")
        elif len(hit) == 1:
            labels.append(hit[0])
        else:
            labels.append("Multiple Viruses (" + ", ".join(hit) + ")")
    return labels, virus_molecules_per_cell(adata, matrix, list(viral_ids))


def _write_scautoqc_umaps(adata):
    """Write UMAPs colored by the in-place scAutoQC cell and cluster flags."""
    if px is None:
        raise RuntimeError("plotly is required to render scAutoQC UMAP HTML plots.")
    if "scautoqc_X_umap_qc" not in adata.obsm:
        log.info("Skipping scAutoQC QC-space UMAP: fewer than six cells were available.")
        return
    coords = adata.obsm["scautoqc_X_umap_qc"]
    df = pd.DataFrame({"UMAP1": coords[:, 0], "UMAP2": coords[:, 1], "barcode": adata.obs_names})
    for key in ("scautoqc_cell_passed_qc", "scautoqc_cluster_passed_qc", "scautoqc_qc_cluster"):
        df[key] = adata.obs[key].astype(str).to_numpy()
        fig = px.scatter(
            df,
            x="UMAP1",
            y="UMAP2",
            color=key,
            hover_data=["barcode"],
            title=f"scAutoQC UMAP: {key.removeprefix('scautoqc_').replace('_', ' ')}",
        )
        fig.update_layout(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white",
        )
        fig.write_html(f"{config.output}/plots/umap_scautoqc_{key.removeprefix('scautoqc_')}.html")


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
    k = min(n_cells, max(min_k, min(k, max_k)))
    return k


def viral_neighbor_enrichment(coords, labels, k, n_permutations=1000, random_state=0):
    rng = np.random.default_rng(random_state)

    viral_cells = np.where(labels == 1)[0]
    if len(viral_cells) == 0 or len(coords) < 2:
        # No viral-positive cells or no other cell to compare — skip enrichment.
        return 0.0, 0.0, 1.0

    k = min(len(coords), max(2, k))
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
    # plots/ is otherwise only created by detection.histogram (visual + a detected
    # virus), so a negative sample or --no-visual run would crash on the QC histogram.
    os.makedirs(f"{config.output}/plots", exist_ok=True)

    # Violin plot and computing QC metrics
    adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()

    if config.scautoqc:
        if "scautoqc_cell_passed_qc" not in adata.obs:
            add_scautoqc(adata, identity)
        _write_scautoqc_umaps(adata)

    # create violin plot
    fig, ax = plt.subplots()
    ax.hist(adata.obs["n_counts"], bins=100)
    ax.set(title="Total counts per cell", xlabel="n_counts", ylabel="Count")
    fig.savefig(f"{config.output}/plots/qc_hist_total_counts.png")
    plt.close(fig)

    # Filtering based on QC threshold (config-driven via PR 11 A4)
    min_counts_threshold = config.min_counts
    min_genes_threshold = config.min_genes
    hvg_min_mean = config.hvg_min_mean
    hvg_max_mean = config.hvg_max_mean
    hvg_min_disp = config.hvg_min_disp
    umap_n_neighbors = config.umap_n_neighbors

    if not config.scautoqc:
        adata = adata[
            (adata.obs["n_counts"] >= min_counts_threshold)
            & (adata.obs["n_genes"] >= min_genes_threshold)
        ].copy()

    if adata.n_obs == 0:
        msg = (
            f"No cells left after QC filtering! "
            f"min_counts={min_counts_threshold}, min_genes={min_genes_threshold}. "
            "Lower --min-counts/--min-genes, or check that the count matrix holds "
            "host genes (a viral-only matrix cannot pass min_counts)."
        )
        log.error(msg)
        raise ValueError(msg)  # fail the Snakemake rule loudly rather than exit 0 with no plots

    # The matrix detection called on (B-23): not raw/counts, which can differ from it.
    # Taken before normalize_total below rewrites adata.X.
    X_counts = select_detection_matrix(adata, config)
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
        sc.pp.highly_variable_genes(
            adata,
            min_mean=hvg_min_mean,
            max_mean=hvg_max_mean,
            min_disp=hvg_min_disp,
        )
        sc.tl.pca(adata, use_highly_variable=True, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors)
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

        outdir = f"{config.output}/plots"
        os.makedirs(outdir, exist_ok=True)
        fig.write_html(f"{outdir}/umap_no_virus.html")
        return

    print("Viral genes detected. Computing viral load and annotated UMAP.")
    k_neighbors = calculate_k_neighbors(adata.n_obs)

    virus_labels, viral_counts = label_infected_cells(
        adata, X_counts, viral_ids, identity, config.infected_min_molecules
    )

    adata.obs["virus_detected"] = virus_labels

    adata.obs["viral_counts"] = viral_counts
    adata.obs["virus_positive"] = adata.obs["virus_detected"] != "No Virus"

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if "X_umap" not in adata.obsm:
        sc.pp.highly_variable_genes(
            adata,
            min_mean=hvg_min_mean,
            max_mean=hvg_max_mean,
            min_disp=hvg_min_disp,
        )
        # Force-include viral genes so they survive HVG filtering
        adata.var["highly_variable"] |= adata.var_names.isin(list(found_genes))
        sc.pp.pca(adata, use_highly_variable=True, random_state=0)
        sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors, random_state=0)
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
    outdir = f"{config.output}/plots"
    os.makedirs(outdir, exist_ok=True)
    fig_binary.write_html(f"{outdir}/umap_binary.html")
    fig_continuous.write_html(f"{outdir}/umap_continuous.html")


def main():
    if config.umap and px is None:
        raise RuntimeError(
            "plotly is required for UMAP HTML plots but is not installed. "
            "Install it with: pip install plotly"
        )

    adata = sc.read_h5ad(str(kb.current_adata(multimapping=config.multimapping)))
    if config.multimapping:
        if adata.uns.get("count_schema_version") != "3.0.0":
            raise ValueError("Multimapping output is not a ViralScan v3 count schema; rebuild it.")

    # Load found genes
    found_genes = {}
    with open(f"{config.output}/log/found_genes.txt") as f:
        for line in f:
            parts = line.strip().split(";")
            if len(parts) == 2:
                gene_id, count = parts
                found_genes[gene_id] = float(count)

    # Check if user wants UMAP
    if config.umap:
        print(
            f"You have decided to create a umap. This can take a while before finishing the code. Please wait..."
        )
        umap(adata, detected_genes(found_genes, config.output, identity))


def run(ctx, done_file):
    """Entry point: optional UMAP for one Run, then touch done_file."""
    global config, kb, identity
    config = ctx.config
    kb = ctx.outputs
    identity = load_run_identity(config.output)

    main()

    with open(done_file, "w") as f:
        f.write("done\n")
        if config.umap:
            log.info("Umap is done!")
    print(f"All (important) results of ViralScan can be found in {config.output}summary.txt")


if "snakemake" in globals():
    run(
        RunContext.from_yaml(snakemake.params.configfile),  # noqa: F821 (snakemake magic global)
        snakemake.output[0],  # noqa: F821
    )
