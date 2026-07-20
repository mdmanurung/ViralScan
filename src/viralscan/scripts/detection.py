"""
The detection scripts uses the gene IDs found in the analysis script and checks
for presence of the gene IDs in the sample (according to kb count). It also
Creates visualizations (except for the umap), e.g. a barplot and a plot showing
super expressors.
"""

# Importing packages
import base64
import datetime
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

from viralscan.anellovirus import merged_name_map
from viralscan.constants import (
    EVE_RISK_GENERA,
    SIBLING_CROSSMAP_RATIO_THRESHOLD,
    SIBLING_VIRUS_PAIRS,
)
from viralscan.enrichment import cell_type_enrichment, write_cell_type_enrichment
from viralscan.multimapping import (
    select_detection_matrix,
    should_write_multimap_evidence,
    summarize_multimap_evidence,
    write_multimap_evidence,
)
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.utils import matrix_for_genes, resolve_count_matrix, setup_script_logging
from viralscan.virus_grouping import group_genes_by_virus

log = setup_script_logging()


# Run-level state, populated by run() from the Run Context. Declared here so the
# helper functions can reference them as module globals; the module imports
# cleanly without Snakemake because nothing reads these at import time.
config: RunConfig = RunConfig()
output: str = ""
kb = None
file: str = ""  # path to the viral-accessions list (the analysis rule's output)


def _gene_counts_from_matrix(matrix, gene_idx):
    gene_counts = matrix[:, gene_idx]
    if hasattr(gene_counts, "toarray"):
        gene_counts = gene_counts.toarray()
    return gene_counts


def _sum_axis1(matrix):
    values = matrix.sum(axis=1)
    if hasattr(values, "A1"):
        return values.A1
    return np.asarray(values).reshape(-1)


def _sum_axis0(matrix):
    values = matrix.sum(axis=0)
    if hasattr(values, "A1"):
        return values.A1
    return np.asarray(values).reshape(-1)


def _count_value(value, ndigits=6):
    """Preserve fractional multimapper counts while keeping whole counts tidy."""
    value = float(value)
    rounded = round(value)
    if np.isclose(value, rounded):
        return int(rounded)
    return round(value, ndigits)


def detect_genes(var_names, matrix, viral_accessions, threshold=1):
    """Return {gene_id: total_count} for viral genes meeting the threshold.

    Pure and importable. ``var_names`` may be a list or a pandas Index; counts
    are summed across all cells and compared with an inclusive ``>=`` threshold.
    """
    pos = {gene_id: i for i, gene_id in enumerate(var_names)}
    found = {}
    for gene_id in viral_accessions:
        if gene_id in pos:
            gene_counts = _gene_counts_from_matrix(matrix, pos[gene_id])
            total_count = gene_counts.sum()
            if total_count >= threshold:
                found[gene_id] = total_count
    return found


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
        viral_accessions (list[str]): viral accessions read from analysis.py output
        detection_matrix: count matrix used for primary viral calls
    """
    viral_accessions = list()
    with open(file) as viral_file:
        for f in viral_file:
            viral_accessions.append(f.strip())

    adata = sc.read_h5ad(str(kb.current_adata(multimapping=config.multimapping)))
    if config.multimapping:
        if "counts_corrected" in adata.layers and "counts_original" in adata.layers:
            adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]

    detection_matrix = select_detection_matrix(adata, config)
    threshold = config.detection_threshold
    found_genes = detect_genes(adata.var_names, detection_matrix, viral_accessions, threshold)
    return adata, found_genes, output, viral_accessions, detection_matrix


def histogram(adata, found_genes, map_virus, outputpath, viral_count_matrix=None):
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
    count_matrix = resolve_count_matrix(viral_count_matrix, adata)
    gene_counts = _sum_axis0(count_matrix)

    # Create dataframe with gene IDs and UMI counts
    df = pd.DataFrame({"gene_id": adata.var_names, "UMI_count": gene_counts})

    # Group versions of found viruses by virus name (boundary-aware match).
    group_by_virus, detected_viral_genes = group_genes_by_virus(found_genes, map_virus)

    # Check if user wants visualizations
    if config.visual:
        for virus in group_by_virus:
            virus_list = group_by_virus[virus]
            if len(virus_list) > 20:
                virus_list = virus_list[:20]
            df_virus = df[df["gene_id"].isin(virus_list)]
            df_virus_sorted = df_virus.sort_values(by="UMI_count", ascending=False).reset_index(
                drop=True
            )

            # Plot the UMI counts
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(data=df_virus_sorted, x="gene_id", y="UMI_count")

            # Check for amount of bars before annotating (max of 10 for annotation)
            if len(ax.patches) <= 10:
                for p in ax.patches:
                    ax.annotate(
                        f"{p.get_height():.0f}",
                        (p.get_x() + p.get_width() / 2.0, p.get_height()),
                        ha="center",
                        va="bottom",
                        fontsize=11,
                        color="black",
                        xytext=(0, 3),
                        textcoords="offset points",
                    )

            # Check if path exists, otherwise create it.
            os.makedirs(f"{outputpath}/plots/", exist_ok=True)

            # Pass variables for bar plot
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("UMI Count")
            plt.xlabel("Gene ID")
            plt.title(f"{virus} Gene UMI Counts (Bar Plot)")
            plt.tight_layout()
            plt.savefig(f"{outputpath}/plots/{virus}_histogram.png")
            plt.close()
    return group_by_virus, detected_viral_genes


def super_expressor(adata, virus, viral_gene_ids, outputpath, viral_count_matrix=None):
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
    adata.obs["nCount_RNA"] = _sum_axis1(adata.X)

    # Match viral gene IDs to adata and raise ValueError
    viral_mask = adata.var_names.isin(viral_gene_ids)
    matched_genes = adata.var_names[viral_mask]

    if matched_genes.empty:
        raise ValueError(
            "None of the provided viral gene IDs were found in the dataset. No super expressor is therefore found."
        )

    # Compute viral UMI counts per cell from the primary-call matrix. Total RNA
    # above remains the full expression matrix for the null model denominator.
    count_matrix = resolve_count_matrix(viral_count_matrix, adata)
    adata.obs[virus] = _sum_axis1(matrix_for_genes(adata, count_matrix, list(matched_genes)))

    # Null Model (grey line)
    total_viral = adata.obs[virus].sum()
    total_rna = adata.obs["nCount_RNA"].sum()
    null_model = total_viral * (adata.obs["nCount_RNA"] / total_rna)

    obs_sorted = np.sort(adata.obs[virus].values)[::-1]
    null_sorted = np.sort(null_model.values)[::-1]

    df_plot = pd.DataFrame(
        {"rank": np.arange(1, len(obs_sorted) + 1), "observed": obs_sorted, "null": null_sorted}
    )

    # Smooth Null model
    window = max(20, int(len(df_plot) / 300))
    df_plot["null_smooth"] = (
        df_plot["null"].rolling(window=window, center=True, min_periods=1).mean()
    )

    # Clip zeros for log plot stability
    df_plot["observed"] = df_plot["observed"].clip(lower=1e-1)

    # Count super-expressors using configurable threshold
    se_threshold = config.se_threshold
    n_SE = (adata.obs[virus] >= se_threshold).sum()

    title = (
        f"Virus {virus}, n_super={n_SE}\nn={adata.n_obs}; {virus} max={int(adata.obs[virus].max())}"
    )

    max_rank_display = 15000  # adjust to visually match paper; 500–20000 is typical
    df_show = df_plot.iloc[:max_rank_display]

    # ---- Plot ----
    plt.figure(figsize=(7, 6))

    plt.plot(
        df_show["rank"],
        df_show["null_smooth"],
        color="grey",
        linewidth=1.5,
        alpha=0.9,
        label="Null model (smoothed)",
    )

    # Mask super-expressors
    super_mask = df_show["observed"] >= se_threshold

    plt.scatter(
        df_show.loc[~super_mask, "rank"],
        df_show.loc[~super_mask, "observed"],
        color="firebrick",
        s=10,
        alpha=0.5,
        label="Observed",
    )

    plt.scatter(
        df_show.loc[super_mask, "rank"],
        df_show.loc[super_mask, "observed"],
        color="red",
        s=25,
        alpha=0.9,
        label=f"Super-expressors (\u2265 {se_threshold} UMI)",
    )

    # Threshold line
    plt.axhline(se_threshold, linestyle="--", color="darkblue", linewidth=1)

    # Log scaling
    plt.yscale("log")
    plt.ylim(1e-1, df_plot["observed"].max() * 1.2)

    # Tidy y-ticks
    ymin, ymax = plt.ylim()
    powers = np.arange(np.floor(np.log10(ymin)), np.ceil(np.log10(ymax)) + 1)
    yticks = 10**powers
    plt.yticks(yticks)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())

    plt.xlabel("Cell Rank", fontsize=11)
    plt.ylabel("Viral UMI Counts", fontsize=11)
    plt.title(title, fontsize=12)

    plt.legend(frameon=True, fontsize=9)
    plt.tight_layout()

    plt.savefig(f"{outputpath}/plots/SuperExpressor_{virus}.png", dpi=500)
    plt.close()


def detect_cells(adata, found_genes, summary, viral_count_matrix=None):
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
        summary (IO[str]): open text file to write the summary to
    """
    # Detect cells and find barcodes for gene IDs
    count_matrix = resolve_count_matrix(viral_count_matrix, adata)
    cells_per_gene = {}
    for viral_gene_name in found_genes:
        gene_counts = matrix_for_genes(adata, count_matrix, [viral_gene_name])
        if hasattr(gene_counts, "toarray"):
            gene_counts = gene_counts.toarray()
        expressed_mask = gene_counts.flatten() > 0
        cells_with_gene = adata.obs_names[expressed_mask].tolist()
        cells_per_gene[viral_gene_name] = cells_with_gene

    for gene, barcodes in cells_per_gene.items():
        summary.write(f"{gene} detected in {len(barcodes)} cells. ")
        summary.write(f"Barcodes: {barcodes}\n")


def compute_stats(
    adata,
    found_genes,
    group_by_virus,
    detected_viral_genes,
    called_mask=None,
    viral_count_matrix=None,
):
    """
    Compute normalized viral detection statistics.

    Parameters
    ----------
    called_mask : np.ndarray[bool] | None
        Boolean mask over ``adata.obs_names`` marking real (non-empty-droplet)
        cells. When provided, per-virus stats are ALSO reported over this subset
        (``*_called`` keys) — the primary, biologically meaningful denominator —
        alongside the all-barcode numbers. ``None`` = all barcodes are cells
        (legacy behaviour; the ``*_called`` values then equal the all-barcode ones).
    viral_count_matrix : matrix | None
        Matrix to use for viral numerator counts and infected-cell masks. ``None``
        preserves legacy behaviour by using ``adata.X``. Total-UMI denominators
        always come from ``adata.X``.

    Returns
    -------
    virus_stats : dict[str, dict]
        Per-virus stats: total_umi, infected_cells, total_cells, pct_infected,
        umi_per_10k, plus n_called_cells, infected_called, pct_infected_called.
    per_cell_df : pd.DataFrame
        One row per cell that carries any viral UMI (with an ``is_called_cell`` flag).
    """
    total_cells = adata.n_obs
    if called_mask is None:
        called_mask = np.ones(total_cells, dtype=bool)
    else:
        called_mask = np.asarray(called_mask, dtype=bool)
    n_called = int(called_mask.sum())

    # Total UMI per cell (sum across all genes)
    total_umi_per_cell = _sum_axis1(adata.X)
    total_umi_all = total_umi_per_cell.sum()

    virus_stats = {}
    cell_rows = []
    count_matrix = resolve_count_matrix(viral_count_matrix, adata)

    for virus, gene_list in group_by_virus.items():
        valid_genes = [g for g in gene_list if g in adata.var_names]
        if not valid_genes:
            continue

        # Per-cell viral UMI for this virus
        viral_matrix = matrix_for_genes(adata, count_matrix, valid_genes)
        if hasattr(viral_matrix, "toarray"):
            viral_matrix = viral_matrix.toarray()
        viral_umi_per_cell = viral_matrix.sum(axis=1)

        total_umi_raw = float(viral_umi_per_cell.sum())
        infected_mask = viral_umi_per_cell > 0
        infected_cells = int(infected_mask.sum())
        pct_infected = round(infected_cells / total_cells * 100, 4) if total_cells else 0.0
        umi_per_10k = round(total_umi_raw / total_umi_all * 10_000, 4) if total_umi_all else 0.0

        # Same stats restricted to called cells (real, non-empty droplets) — the
        # primary denominator. Empty droplets otherwise dilute pct_infected.
        infected_called = int((infected_mask & called_mask).sum())
        pct_infected_called = round(infected_called / n_called * 100, 4) if n_called else 0.0

        # Accession breadth: fraction of reference gene IDs (accessions) for
        # this virus that have ≥1 UMI in any cell. EVE artifacts concentrate on
        # 1-2 host-integrated loci; genuine infection spreads across ORF1/ORF2/ORF3.
        gene_has_count = (viral_matrix > 0).any(axis=0)
        n_acc_detected = int(np.asarray(gene_has_count).flatten().sum())
        n_acc_total = len(valid_genes)
        accession_breadth = round(n_acc_detected / n_acc_total, 4) if n_acc_total else 0.0

        # Host–viral ambiguity fraction: proportion of viral UMI that mapped
        # ambiguously to both host and viral index (written by multimap.py).
        # High fraction indicates reads originating from host genomic regions
        # (e.g. EVE integrations in expressed host genes).
        host_viral_ambig_fraction = None
        if "counts_host_viral_ambiguous" in adata.layers and total_umi_raw > 0:
            ambig_matrix = adata[:, valid_genes].layers["counts_host_viral_ambiguous"]
            if hasattr(ambig_matrix, "toarray"):
                ambig_matrix = ambig_matrix.toarray()
            ambiguity_denominator = float(matrix_for_genes(adata, adata.X, valid_genes).sum())
            host_viral_ambig_fraction = (
                round(float(np.asarray(ambig_matrix).sum()) / ambiguity_denominator, 4)
                if ambiguity_denominator > 0
                else None
            )

        virus_stats[virus] = {
            "total_umi": _count_value(total_umi_raw),
            "infected_cells": infected_cells,
            "total_cells": total_cells,
            "pct_infected": pct_infected,
            "umi_per_10k": umi_per_10k,
            "n_called_cells": n_called,
            "infected_called": infected_called,
            "pct_infected_called": pct_infected_called,
            "accession_breadth": accession_breadth,
            "host_viral_ambig_fraction": host_viral_ambig_fraction,
        }

        # Per-cell rows (only infected cells)
        barcodes = adata.obs_names[infected_mask]
        infected_indices = np.where(infected_mask)[0]
        for i, bc in enumerate(barcodes):
            idx = infected_indices[i]
            cell_total = float(total_umi_per_cell[idx])
            v_umi = float(viral_umi_per_cell[idx])
            cell_rows.append(
                {
                    "barcode": bc,
                    "virus_name": virus,
                    "viral_umi": _count_value(v_umi),
                    "total_umi": _count_value(cell_total),
                    "viral_fraction": round(v_umi / cell_total, 6) if cell_total else 0.0,
                    "is_called_cell": bool(called_mask[idx]),
                }
            )

    per_cell_df = pd.DataFrame(
        cell_rows,
        columns=[
            "barcode",
            "virus_name",
            "viral_umi",
            "total_umi",
            "viral_fraction",
            "is_called_cell",
        ],
    )
    return virus_stats, per_cell_df


def check_sibling_crossmapping(virus_stats):
    """Return {virus_name: note_str} for viruses flagged as likely EM bleed.

    When two viruses that share >80% sequence identity (HHV-6A/6B, HSV-1/2)
    are both detected and the UMI ratio exceeds SIBLING_CROSSMAP_RATIO_THRESHOLD,
    the weaker signal is flagged. The global EM allocates a small fraction of the
    dominant sibling's shared-region multimappers to the other, producing a
    residual that is EM noise rather than genuine co-infection. A log warning is
    also emitted for each flagged pair.
    """
    notes = {}
    checked = set()
    for virus, stats in virus_stats.items():
        sibling = SIBLING_VIRUS_PAIRS.get(virus)
        if sibling is None or virus in checked or sibling not in virus_stats:
            continue
        checked.add(virus)
        checked.add(sibling)
        umi_a = float(stats["total_umi"])
        umi_b = float(virus_stats[sibling]["total_umi"])
        if umi_a <= 0 or umi_b <= 0:
            continue
        ratio = max(umi_a, umi_b) / min(umi_a, umi_b)
        if ratio < SIBLING_CROSSMAP_RATIO_THRESHOLD:
            continue
        weaker, dominant = (virus, sibling) if umi_a < umi_b else (sibling, virus)
        dom_umi = max(umi_a, umi_b)
        wk_umi = min(umi_a, umi_b)
        notes[weaker] = (
            f"possible_em_bleed: {ratio:.0f}:1 ratio vs {dominant} "
            f"({dom_umi} vs {wk_umi} UMI); closely related siblings share "
            f"high k-mer identity — the global EM allocates a small fraction "
            f"of shared-region multimappers to the weaker sibling"
        )
        log.warning(
            "Sibling cross-mapping: %s (%s UMI) vs %s (%s UMI), ratio %.0f:1 "
            "(threshold %.0f). Weaker signal may be EM bleed; see "
            "sibling_crossmap_note in viral_summary.tsv.",
            dominant,
            dom_umi,
            weaker,
            wk_umi,
            ratio,
            SIBLING_CROSSMAP_RATIO_THRESHOLD,
        )
    return notes


def write_tsv_outputs(virus_stats, per_cell_df, outputpath, crossmap_notes=None):
    """Write viral_summary.tsv and per_cell_viral.tsv to results/ sub-folder."""
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    crossmap_notes = crossmap_notes or {}

    # Per-virus summary
    summary_rows = []
    for virus, s in virus_stats.items():
        summary_rows.append(
            {
                "virus_name": virus,
                "total_umi": s["total_umi"],
                # Primary (called-cell) denominator — real, non-empty droplets.
                "infected_called": s.get("infected_called", s["infected_cells"]),
                "n_called_cells": s.get("n_called_cells", s["total_cells"]),
                "pct_infected_called": s.get("pct_infected_called", s["pct_infected"]),
                # Secondary (all-barcode) denominator — kept so the choice is explicit.
                "infected_cells": s["infected_cells"],
                "total_cells": s["total_cells"],
                "pct_infected": s["pct_infected"],
                "umi_per_10k": s["umi_per_10k"],
                "sibling_crossmap_note": crossmap_notes.get(virus, ""),
                # EVE artifact flags
                "accession_breadth": s.get("accession_breadth", 0.0),
                "host_viral_ambig_fraction": s.get("host_viral_ambig_fraction"),
                "eve_risk": any(g in virus for g in EVE_RISK_GENERA),
            }
        )
    virus_df = pd.DataFrame(
        summary_rows,
        columns=[
            "virus_name",
            "total_umi",
            "infected_called",
            "n_called_cells",
            "pct_infected_called",
            "infected_cells",
            "total_cells",
            "pct_infected",
            "umi_per_10k",
            "sibling_crossmap_note",
            "accession_breadth",
            "host_viral_ambig_fraction",
            "eve_risk",
        ],
    )
    virus_df.to_csv(os.path.join(results_dir, "viral_summary.tsv"), sep="\t", index=False)

    # Per-cell viral annotation
    per_cell_df.to_csv(os.path.join(results_dir, "per_cell_viral.tsv"), sep="\t", index=False)
    log.info("Wrote results/viral_summary.tsv and results/per_cell_viral.tsv")


def infected_cell_count(per_cell_df):
    """Return unique virus-positive barcodes represented in per-cell outputs."""
    if per_cell_df is None or per_cell_df.empty or "barcode" not in per_cell_df.columns:
        return 0
    return int(per_cell_df["barcode"].nunique())


def _encode_image(path: str) -> str:
    """Return a base64-encoded PNG string for embedding in HTML."""
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


def generate_html_report(
    virus_stats,
    per_cell_df,
    cell_type_enrichment_df,
    multimap_evidence_df,
    group_by_virus,
    detected_viral_genes,
    outputpath,
    run_date=None,
):
    """Render the Jinja2 HTML report and write it to <outputpath>/report.html."""
    try:
        from jinja2 import Environment, FileSystemLoader, TemplateNotFound
    except ImportError:
        log.warning("jinja2 not installed — skipping HTML report. pip install jinja2 to enable.")
        return

    template_path = os.path.join(os.path.dirname(__file__), "..", "templates")
    env = Environment(loader=FileSystemLoader(template_path), autoescape=True)

    try:
        template = env.get_template("report.html.j2")
    except TemplateNotFound as exc:
        log.warning("Could not load HTML report template: %s", exc)
        return

    # Collect embedded plot images
    plots_dir = os.path.join(outputpath, "plots")
    embedded_plots = {}
    if os.path.isdir(plots_dir):
        for fname in os.listdir(plots_dir):
            if fname.endswith(".png"):
                fpath = os.path.join(plots_dir, fname)
                try:
                    with open(fpath, "rb") as f:
                        embedded_plots[fname] = base64.b64encode(f.read()).decode("utf-8")
                except OSError:
                    pass

    ctx = {
        "run_date": run_date or datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "output_dir": outputpath,
        "virus_stats": virus_stats,
        "detected_viruses": sorted(detected_viral_genes),
        "total_viruses": len(virus_stats),
        "se_threshold": config.se_threshold,
        "detection_threshold": config.detection_threshold,
        "embedded_plots": embedded_plots,
        "any_infected": any(s["infected_cells"] > 0 for s in virus_stats.values()),
        "per_cell_count": infected_cell_count(per_cell_df),
        "cell_type_enrichment": cell_type_enrichment_df.to_dict("records")
        if cell_type_enrichment_df is not None and not cell_type_enrichment_df.empty
        else [],
        "multimap_evidence": multimap_evidence_df.to_dict("records")
        if multimap_evidence_df is not None and not multimap_evidence_df.empty
        else [],
        "multimap_method": config.multimap_method,
        "multimap_primary_call": config.multimap_primary_call,
    }

    html = template.render(**ctx)
    report_path = os.path.join(outputpath, "report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    log.info("HTML report written to %s", report_path)


def main():
    adata, found_genes, outputpath, viral_accessions, detection_matrix = preprocessing()

    # check if user wants visuals in output directory
    group_by_virus, detected_viral_genes = histogram(
        adata, found_genes, merged_name_map(), outputpath, viral_count_matrix=detection_matrix
    )
    if config.visual:
        for virus in group_by_virus:
            super_expressor(
                adata,
                virus,
                group_by_virus[virus],
                outputpath,
                viral_count_matrix=detection_matrix,
            )

    # Cell-calling: label real (non-empty-droplet) barcodes so viral rates are
    # reported over called cells, not over all barcodes (which are mostly empty).
    # Prefers an external list (CellRanger/STARsolo cells); else emptyDrops/knee.
    called_mask = None
    try:
        from viralscan.scripts.cellcalling import call_cells

        counts_dir = os.path.join(config.output, "kb-python", "counts_unfiltered")
        called_mask = call_cells(adata, config, matrix_dir=counts_dir)
    except Exception as exc:  # never let cell-calling break the legacy summary
        log.warning("cell-calling failed (%s); reporting over all barcodes only", exc)

    # Compute normalized statistics (PR 11 A1/A3) over both denominators
    virus_stats, per_cell_df = compute_stats(
        adata,
        found_genes,
        group_by_virus,
        detected_viral_genes,
        called_mask=called_mask,
        viral_count_matrix=detection_matrix,
    )

    # Optional enrichment by cell type labels (PR 11 A5) — restricted to detected viruses.
    detected_groups = {v: genes for v, genes in group_by_virus.items() if v in virus_stats}
    cell_type_df = cell_type_enrichment(
        adata, detected_groups, config, viral_count_matrix=detection_matrix
    )

    # Ambiguity-aware multimapper evidence is additive and does not alter
    # legacy viral_summary.tsv/per_cell_viral.tsv schemas.
    if should_write_multimap_evidence(config):
        evidence_gene_ids = [g for g in viral_accessions if g in adata.var_names]
        evidence_groups, _ = group_genes_by_virus(evidence_gene_ids, merged_name_map())
        multimap_evidence_df = summarize_multimap_evidence(adata, evidence_groups, config)
    else:
        multimap_evidence_df = summarize_multimap_evidence(None, {}, config)

    # Flag sibling pairs with high UMI asymmetry (HHV-6A/6B, HSV-1/2)
    crossmap_notes = check_sibling_crossmapping(virus_stats)

    # Write structured TSV outputs (PR 11 A1)
    write_tsv_outputs(virus_stats, per_cell_df, outputpath, crossmap_notes=crossmap_notes)
    write_cell_type_enrichment(cell_type_df, outputpath)
    if should_write_multimap_evidence(config):
        write_multimap_evidence(multimap_evidence_df, outputpath)

    # Writing results to the legacy summary file (kept for backward-compat)
    found_genes_sorted = dict(sorted(found_genes.items()))
    total_viral_genes = 0
    counts_per_virus = {}
    with (
        open(f"{config.output}/summary.txt", "w") as summary,
        open(f"{config.output}log/found_genes.txt", "w") as found_genes_file,
    ):
        if len(found_genes_sorted) > 0:
            summary.write("Found viral Gene IDs including the count:\n")
            summary.write("Gene ID; Gene Count\n")
            for g in found_genes_sorted:
                key = next((k for k, v in group_by_virus.items() if g in v), None)
                if key not in counts_per_virus:
                    counts_per_virus[key] = found_genes_sorted[g]
                else:
                    counts_per_virus[key] += found_genes_sorted[g]

                write_to_file = f"{g};{found_genes_sorted[g]}\n"
                total_viral_genes += found_genes_sorted[g]
                summary.write(write_to_file)
                found_genes_file.write(write_to_file)
        if len(found_genes_sorted) > 0:
            for virus_name, stats in virus_stats.items():
                summary.write(
                    f"\n{virus_name}: {stats['total_umi']} UMI total, "
                    f"{stats['infected_cells']}/{stats['total_cells']} cells infected "
                    f"({stats['pct_infected']:.2f}%), "
                    f"{stats['umi_per_10k']:.2f} UMI/10k."
                )
            summary.write(f"\n\nTotal amount of viral load found: {total_viral_genes}")
            summary.write(
                f"\n\nOfficial name of viral load detected: {','.join(str(s) for s in detected_viral_genes)}"
            )
        else:
            summary.write("No viral gene IDs found in this sample for the viruses in the index.")
        summary.write(
            f"\nIf you want to see the cell gene matrix, go to the kb-python/counts_unfiltered/ folder and look for the cells_x_genes.mtx file.\n"
        )
        detect_cells(adata, found_genes, summary, viral_count_matrix=detection_matrix)

    # Generate HTML report (PR 11 A2)
    generate_html_report(
        virus_stats,
        per_cell_df,
        cell_type_df,
        multimap_evidence_df,
        group_by_virus,
        detected_viral_genes,
        outputpath,
    )


def run(ctx, viral_accessions_file, done_file):
    """Entry point: detect viruses, write outputs/report for one Run."""
    global config, output, kb, file
    config = ctx.config
    output = config.output
    kb = ctx.outputs
    file = viral_accessions_file

    main()

    with open(done_file, "w") as f:
        f.write("done\n")
    log.info("The detection and visualizations are done!")


if "snakemake" in globals():
    run(
        RunContext.from_yaml(snakemake.params.configfile),  # noqa: F821 (snakemake magic global)
        snakemake.input.file_viral_accessions,  # noqa: F821
        snakemake.output[0],  # noqa: F821
    )
