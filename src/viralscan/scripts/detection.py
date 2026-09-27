"""
The detection scripts uses the gene IDs found in the analysis script and checks
for presence of the gene IDs in the sample (according to kb count). It also
Creates visualizations (except for the umap), e.g. a barplot and a plot showing
super expressors.
"""

# Importing packages
import base64
import datetime
import json
import logging
import os
from pathlib import Path

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
from viralscan.sensitivity import (
    SENSITIVITY_COLUMNS,
    negative_result_statement,
    sensitivity_record,
)
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


def _sum_axis(matrix, axis):
    values = matrix.sum(axis=axis)
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
    # Duplicate accessions in the reference make gene-name -> column lookup
    # ambiguous. Renaming (var_names_make_unique) would silently drop the renamed
    # copy's counts, since the viral accession list still carries the original name.
    # Fail loud with an actionable message instead.
    if not adata.var_names.is_unique:
        dups = adata.var_names[adata.var_names.duplicated()].unique().tolist()
        raise ValueError(
            f"Reference has {len(dups)} duplicate gene IDs (e.g. {dups[:5]}); counts for "
            "duplicates cannot be attributed unambiguously. Rebuild the reference with "
            "unique gene_ids (kb ref collapses one var_name per gene_id)."
        )
    if config.multimapping:
        if adata.uns.get("count_schema_version") != "3.0.0":
            raise ValueError("Multimapping output is not a ViralScan v3 count schema; rebuild it.")

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
    gene_counts = _sum_axis(count_matrix, 0)

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
    adata.obs["nCount_RNA"] = _sum_axis(adata.X, 1)

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
    adata.obs[virus] = _sum_axis(matrix_for_genes(adata, count_matrix, list(matched_genes)), 1)

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


#: Filename patterns detection owns inside ``plots/``. The directory is shared
#: with umap.py, so only these are cleared — never the whole directory.
_OWNED_PLOT_PATTERNS = ("*_histogram.png", "SuperExpressor_*.png")


def clear_stale_virus_plots(outputpath) -> list[str]:
    """Drop last run's per-virus plots before regenerating this run's.

    One plot is written per virus that clears ``detection_threshold``. Which
    viruses clear it depends on the multimap method, so ``rerun-multimap`` can
    legitimately shrink that set — and because the directory was only ever
    created with ``exist_ok=True``, the demoted virus's plot survived.
    ``generate_html_report`` globs this directory, so the regenerated report
    embedded a figure for a virus its own summary table no longer listed.

    Returns the removed paths. umap.py also writes here; its files do not match
    the owned patterns and are left alone.
    """
    plots_dir = Path(outputpath) / "plots"
    if not plots_dir.is_dir():
        return []
    removed = []
    for pattern in _OWNED_PLOT_PATTERNS:
        for path in sorted(plots_dir.glob(pattern)):
            path.unlink()
            removed.append(str(path))
    if removed:
        log.info("cleared %d stale per-virus plot(s) before regenerating", len(removed))
    return removed


def _headline_totals(adata, detected_viral_genes, viral_count_matrix=None) -> dict:
    """Totals that head summary.txt, derived from the matrix currently on disk.

    Deriving these here rather than carrying them forward from multimap keeps
    them consistent with whichever allocation layer the H5AD holds, including
    after ``rerun-multimap`` swaps that layer in place.
    """
    genes = [gene for gene in detected_viral_genes if gene in adata.var_names]
    n_cells = int(adata.n_obs)
    if not genes:
        return {"unique": 0, "selected": 0, "cells_with_virus": 0, "n_cells": n_cells}

    selected_matrix = resolve_count_matrix(viral_count_matrix, adata)
    selected = matrix_for_genes(adata, selected_matrix, genes)
    unique_layer = adata.layers.get("counts_unique")
    unique = matrix_for_genes(adata, unique_layer, genes) if unique_layer is not None else None

    per_cell = np.asarray(_sum_axis(selected, 1)).ravel()
    return {
        "unique": _count_value(np.asarray(_sum_axis(unique, 1)).sum()) if unique is not None else 0,
        "selected": _count_value(per_cell.sum()),
        "cells_with_virus": int((per_cell > 0).sum()),
        "n_cells": n_cells,
    }


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
        viral_count_matrix: count layer to read barcode occupancy from. When
            None, resolve_count_matrix falls back to adata.X; pass the
            multimap-corrected viral matrix to count cells on the same layer
            the reported viral loads come from.
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
        Per-virus molecule-estimate statistics and explicit cell denominators.
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
    total_umi_per_cell = _sum_axis(adata.X, 1)
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

        # A fixed, strategy-independent denominator. See comparable_called_cells.
        comparable = _comparable_called_cells(adata, called_mask)
        infected_comparable = int((infected_mask & comparable).sum())
        pct_infected_comparable = (
            round(infected_comparable / int(comparable.sum()) * 100, 4) if comparable.any() else 0.0
        )

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
                # Clamp to [0,1]: the ambiguous layer is not a strict subset of the
                # full-matrix viral counts, so the raw ratio can slightly exceed 1.
                round(min(1.0, float(np.asarray(ambig_matrix).sum()) / ambiguity_denominator), 4)
                if ambiguity_denominator > 0
                else None
            )

        virus_stats[virus] = {
            "viral_molecules_total_est": _count_value(total_umi_raw),
            "infected_cells": infected_cells,
            "total_cells": total_cells,
            "pct_infected": pct_infected,
            "viral_molecules_per_10k_est": umi_per_10k,
            "n_called_cells": n_called,
            "infected_called": infected_called,
            "pct_infected_called": pct_infected_called,
            "n_comparable_cells": int(comparable.sum()),
            "infected_comparable": infected_comparable,
            "pct_infected_comparable": pct_infected_comparable,
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
                    "viral_molecules_total_est": _count_value(v_umi),
                    "molecules_total_est": _count_value(cell_total),
                    "viral_fraction": round(v_umi / cell_total, 6) if cell_total else 0.0,
                    "is_called_cell": bool(called_mask[idx]),
                }
            )

    per_cell_df = pd.DataFrame(
        cell_rows,
        columns=[
            "barcode",
            "virus_name",
            "viral_molecules_total_est",
            "molecules_total_est",
            "viral_fraction",
            "is_called_cell",
        ],
    )
    return virus_stats, per_cell_df


def measure_positive_control(adata, config, count_matrix=None, depth=None):
    """Recover the k-mer capture term from a spike-in at known abundance.

    A positive control is the only thing in a run that can measure the term
    depth cannot: the fraction of true viral molecules that survive exact
    k-mer matching against this reference. Capture is what decides whether a
    negative is informative — it falls to 0.32 at 15 % divergence and 0.06 at
    20 % — and it is invisible in a count matrix.

    Returns ``(capture, detail_dict)``. ``capture`` is ``None`` when no control
    was configured, in which case the caller must treat every negative as
    uncertifiable.
    """
    gene = getattr(config, "positive_control_gene", None)
    expected = getattr(config, "positive_control_expected_molecules", None)
    if not gene or expected is None:
        return None, {"status": "not-configured"}
    if gene not in adata.var_names:
        # A control gene absent from the index cannot have been measured, and a
        # zero here is the harshest possible signal — but it is a configuration
        # error, not a biological result, so it is reported as such.
        return None, {
            "status": "gene-not-in-reference",
            "gene": gene,
            "detail": (
                f"positive_control_gene {gene!r} is not a column of the count "
                "matrix, so no capture term could be measured. Rebuild the "
                "reference so the control is included."
            ),
        }
    if depth is None:
        depth = float(_sum_axis(adata.X, 1).sum())
    matrix = resolve_count_matrix(count_matrix, adata)
    observed = float(matrix_for_genes(adata, matrix, [gene]).sum())
    capture = observed / float(expected)
    detail = {
        "status": "measured",
        "gene": gene,
        "expected_molecules": float(expected),
        "observed_molecules": observed,
        "capture": capture,
        "implied_divergence": _implied_divergence(capture),
    }
    if capture <= 0:
        detail["status"] = "failed"
        detail["detail"] = (
            f"control gene {gene!r} was planted at {expected:g} molecules but "
            f"{observed:g} were recovered. The reference cannot see this "
            "sequence at this abundance, so every negative in this run is "
            "uninterpretable."
        )
        log.error("Positive control FAILED: %s", detail["detail"])
        return None, detail
    if capture > 1.0:
        # More than was planted. Either the count is not spike-in-specific (the
        # gene is endogenous, or reads leaked from another cell) or the planted
        # estimate was wrong. Either way the control does not cleanly bound
        # capture loss, so it is surfaced rather than quietly clamped.
        detail["status"] = "over-recovered"
        detail["detail"] = (
            f"control gene {gene!r} recovered {observed:g} molecules against "
            f"{expected:g} planted (ratio {capture:.3g}). The control is not "
            "spike-in-specific, or the planted estimate is wrong, so it cannot "
            "bound capture loss. Treated as no measurable loss, which is "
            "optimistic."
        )
        log.warning("Positive control OVER-RECOVERED: %s", detail["detail"])
        return 1.0, detail
    log.info(
        "Positive control %s: %g/%g molecules recovered -> capture=%.4f (implies "
        "~%.1f%% divergence at 90 bp, k=31)",
        gene,
        observed,
        expected,
        capture,
        (detail["implied_divergence"] or 0.0) * 100,
    )
    return capture, detail


#: Divergence beyond which the implied-diversity inversion is abandoned.
#: ``fragment_capture`` underflows to exactly 0.0 in float64 somewhere past ~0.8,
#: and is already below 1e-7 by 0.45, so the curve is numerically flat over the
#: region where a bisection would have to search. A control that recovers less
#: than ``fragment_capture(0.5)`` therefore has no identifiable implied
#: divergence and is reported as such rather than given a confident number.
MAX_IDENTIFIABLE_DIVERGENCE = 0.5


def _implied_divergence(capture: float, read_length: int = 90) -> float | None:
    """Per-base divergence whose capture matches ``capture``, or None if unidentifiable.

    Inverts :func:`viralscan.sensitivity.fragment_capture` by bisection, so the
    result is only meaningful for capture below 1.0. A capture of 1.0 means "no
    loss measurable", which is reported as None rather than 0 % so a reader does
    not over-read it as proof of zero divergence. Likewise a capture at or below
    the value reachable at :data:`MAX_IDENTIFIABLE_DIVERGENCE` is not
    identifiable and also returns None.
    """
    from viralscan.sensitivity import DEFAULT_K, fragment_capture

    if capture >= 0.999 or capture <= 0:
        return None
    lo, hi = 0.0, MAX_IDENTIFIABLE_DIVERGENCE
    if capture <= fragment_capture(hi, read_length=read_length, k=DEFAULT_K):
        return None  # below anything the model can attribute to divergence
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if fragment_capture(mid, read_length=read_length, k=DEFAULT_K) > capture:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


#: Total-UMI floor for the strategy-independent denominator, in host molecules
#: per barcode. Chosen to sit above the empty-droplet mode and below the knee of
#: a real 10x barcode-rank curve, so it selects cells under *any* host-filter
#: strategy. Deliberately not ``defaults.min_counts`` (1000): that is a UMAP QC
#: knob, and after host subtraction most barcodes fall below it.
COMPARABLE_CELL_MIN_UMI = 200.0


def _comparable_called_cells(adata, called_mask):
    """Barcode mask that means the same thing under every host-filter strategy.

    Why this exists
    ---------------
    ``pct_infected_called`` divides by the *called cells of this run*, and
    cell calling is applied after host subtraction, so the denominator changes
    with the strategy. Measured on the same covid PBMC sample:

    ==================  =========  ==================  ===============
    strategy             called      Alphatorquevirus    pct_called
    ==================  =========  ==================  ===============
    no host filter         143,243          1,167,103        56.64 %
    STAR host filter        28,921             57,715        62.89 %
    ==================  =========  ==================  ===============

    Viral molecules fell **20.2x** and the reported prevalence went **up**. The
    denominator collapsed 5.0x faster than the numerator, so the rate inverted.
    Anyone comparing those two runs would conclude host filtering barely changed
    prevalence, which is the opposite of what happened.

    The fix is a denominator that does not move: barcodes whose *host* molecule
    count clears an absolute floor, measured on whatever matrix is in hand. A
    barcode that clears it is a real cell under any strategy, because host
    filtering only ever removes host signal. ``pct_infected_called`` stays as the
    within-run primary; ``pct_infected_comparable`` is the cross-strategy number.
    """
    total = _sum_axis(adata.X, 1)
    return (np.asarray(total) >= COMPARABLE_CELL_MIN_UMI) & np.asarray(called_mask, dtype=bool)


def build_sensitivity_table(adata, virus_stats, config, depth=None, capture=None):
    """Per-virus detection sensitivity for this run.

    ``depth`` is total quantified molecules (sum of ``adata.X``), which is the
    only depth term a run can actually measure: only quantified molecules can
    be detected, so raw read counts would overstate sensitivity.

    Every virus that cleared the detection threshold gets a row, and so does
    every virus present in the reference that did not. The latter is the point:
    a virus that is absent from the table is a virus nobody asked about, and a
    virus with a zero in ``observed_molecules`` is a negative whose meaning
    depends on the LOD columns beside it.
    """
    if depth is None:
        depth = float(_sum_axis(adata.X, 1).sum())
    if capture is None:
        capture, _ = measure_positive_control(adata, config, depth=depth)
    capture_measured = capture is not None
    records = []
    for virus, stats in virus_stats.items():
        observed = float(stats.get("viral_molecules_total_est", 0) or 0)
        records.append(
            sensitivity_record(
                virus,
                observed_molecules=observed,
                depth=depth,
                detection_threshold=int(config.detection_threshold),
                capture=capture if capture_measured else 1.0,
                capture_measured=capture_measured,
            )
        )
    return pd.DataFrame(
        [r.as_row() for r in records],
        columns=list(SENSITIVITY_COLUMNS),
    )


def write_sensitivity_table(sensitivity_df, outputpath):
    """Write results/sensitivity.tsv. Returns its path."""
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, "sensitivity.tsv")
    sensitivity_df.to_csv(path, sep="\t", index=False)
    log.info("Wrote results/sensitivity.tsv (%d virus row(s))", len(sensitivity_df))
    return path


def write_control_report(control_detail, measured_capture, outputpath):
    """Write results/positive_control.json. Returns its path."""
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    payload = dict(control_detail)
    payload["capture_used_for_sensitivity"] = measured_capture
    payload["certifies_negatives"] = measured_capture is not None
    path = os.path.join(results_dir, "positive_control.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    log.info("Wrote results/positive_control.json (status=%s)", payload.get("status"))
    return path


def run_sensitivity_statement(adata, config, depth=None):
    """The caveat text for a run that detected nothing, with its LOD."""
    if depth is None:
        depth = float(_sum_axis(adata.X, 1).sum())
    return negative_result_statement(depth, detection_threshold=int(config.detection_threshold))


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
        umi_a = float(stats["viral_molecules_total_est"])
        umi_b = float(virus_stats[sibling]["viral_molecules_total_est"])
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


def reference_provenance(config, viral_accessions, detected_viruses):
    """Provenance of the viral reference used for a run.

    Viral annotation choices materially change per-virus results (the paper's
    EBV LMP-1/EBNA attribution divergence is annotation-driven), so record the
    exact reference and its viral accessions alongside the results.
    """
    from viralscan import __version__

    return {
        "viralscan_version": __version__,
        "index": config.index or None,
        "transcripts_t2g": config.transcripts or None,
        "gtf": config.gtf,
        "fasta": config.fasta,
        "technology": config.technology,
        "multimapping": bool(config.multimapping),
        "multimap_method": config.multimap_method,
        "multimap_primary_call": config.multimap_primary_call,
        "n_viral_accessions_in_reference": len(viral_accessions),
        "viral_accessions": sorted(viral_accessions),
        "n_viruses_detected": len(detected_viruses),
        "viruses_detected": sorted(detected_viruses),
    }


def write_reference_provenance(config, viral_accessions, detected_viruses, outputpath):
    """Write results/reference_provenance.json; returns its path."""
    prov = reference_provenance(config, viral_accessions, detected_viruses)
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, "reference_provenance.json")
    with open(path, "w") as fh:
        json.dump(prov, fh, indent=2, sort_keys=True)
    return path


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
                "viral_molecules_total_est": s["viral_molecules_total_est"],
                # Primary (called-cell) denominator — real, non-empty droplets.
                "infected_called": s.get("infected_called", s["infected_cells"]),
                "n_called_cells": s.get("n_called_cells", s["total_cells"]),
                "pct_infected_called": s.get("pct_infected_called", s["pct_infected"]),
                # Cross-strategy denominator: an absolute host-UMI floor, so the
                # rate does not move when host filtering changes how many
                # barcodes clear cell calling. See _comparable_called_cells.
                "infected_comparable": s.get("infected_comparable", ""),
                "n_comparable_cells": s.get("n_comparable_cells", ""),
                "pct_infected_comparable": s.get("pct_infected_comparable", ""),
                # Secondary (all-barcode) denominator — kept so the choice is explicit.
                "infected_cells": s["infected_cells"],
                "total_cells": s["total_cells"],
                "pct_infected": s["pct_infected"],
                "viral_molecules_per_10k_est": s["viral_molecules_per_10k_est"],
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
            "viral_molecules_total_est",
            "infected_called",
            "n_called_cells",
            "pct_infected_called",
            "infected_comparable",
            "n_comparable_cells",
            "pct_infected_comparable",
            "infected_cells",
            "total_cells",
            "pct_infected",
            "viral_molecules_per_10k_est",
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


def generate_html_report(
    virus_stats,
    per_cell_df,
    cell_type_enrichment_df,
    multimap_evidence_df,
    group_by_virus,
    detected_viral_genes,
    outputpath,
    run_date=None,
    sensitivity_df=None,
    sensitivity_statement=None,
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
        # Detection sensitivity. `sensitivity_rows` drives the per-virus LOD
        # table; `sensitivity_statement` is the run-level caveat, shown
        # whenever nothing was detected so a zero is never read as an absence.
        "sensitivity_rows": sensitivity_df.to_dict("records")
        if sensitivity_df is not None and not sensitivity_df.empty
        else [],
        "sensitivity_statement": sensitivity_statement or "",
    }

    html = template.render(**ctx)
    report_path = os.path.join(outputpath, "report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    log.info("HTML report written to %s", report_path)


def main():
    adata, found_genes, outputpath, viral_accessions, detection_matrix = preprocessing()

    # Clear last run's per-virus plots first: which viruses clear the detection
    # threshold is method-dependent, so a rerun can demote one, and a surviving
    # stale plot would be re-embedded by generate_html_report (SW-04).
    clear_stale_virus_plots(outputpath)

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
    # Fail closed. Cell calling sets the denominator for every reported viral
    # rate, so a failure that fell back to all barcodes would not lose a number,
    # it would silently change what the number means. Reporting over all
    # barcodes is available, but only by asking for it: --cell-calling none.
    from viralscan.scripts.cellcalling import CellCallingError, call_cells

    counts_dir = os.path.join(config.output, "kb-python", "counts_unfiltered")
    try:
        called_mask = call_cells(adata, config, matrix_dir=counts_dir)
    except CellCallingError:
        raise
    except Exception as exc:
        raise CellCallingError(
            f"cell calling failed: {exc}. Rerun with --cell-calling none to report "
            "over all barcodes deliberately, or fix the caller inputs."
        ) from exc

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

    # Per-virus detection sensitivity (LOD95). Written for every run, not only
    # negatives: a reader who sees Betatorquevirus 1,142 UMI needs the same
    # yardstick as one who sees nothing, and the LOD columns state which of the
    # two limits (depth or reference capture) is actually binding.
    quantified_depth = float(_sum_axis(adata.X, 1).sum())
    measured_capture, control_detail = measure_positive_control(
        adata, config, count_matrix=detection_matrix, depth=quantified_depth
    )
    sensitivity_df = build_sensitivity_table(
        adata, virus_stats, config, depth=quantified_depth, capture=measured_capture
    )

    # Fail closed on an uncertifiable negative when the run was told to require a
    # control. RunConfig already rejects `require_positive_control` with no
    # control, so reaching here with a missing capture means the control was
    # configured and then failed to yield a number (absent gene, or zero
    # recovery). Either way the run cannot certify a negative, and saying so is
    # the whole point of requiring the control.
    nothing_detected = not found_genes
    if nothing_detected and getattr(config, "require_positive_control", False):
        if measured_capture is None:
            raise CellCallingError(
                "require_positive_control is set, this run detected no viral "
                "signal, and no capture term could be measured from the control "
                f"({control_detail.get('status')}). A negative with no "
                "demonstrated ability to see the target is not evidence of "
                f"absence. Control detail: {control_detail.get('detail', control_detail)}"
            )
    if nothing_detected and measured_capture is None:
        log.warning(
            "No viral signal and no positive control: every negative in this run "
            "is a sampling statement, not an absence. The depth limit is "
            "reported in summary.txt and results/sensitivity.tsv."
        )

    # Write structured TSV outputs (PR 11 A1)
    write_tsv_outputs(virus_stats, per_cell_df, outputpath, crossmap_notes=crossmap_notes)
    write_sensitivity_table(sensitivity_df, outputpath)
    write_control_report(control_detail, measured_capture, outputpath)
    write_cell_type_enrichment(cell_type_df, outputpath)
    write_reference_provenance(config, viral_accessions, list(virus_stats.keys()), outputpath)
    if should_write_multimap_evidence(config):
        write_multimap_evidence(multimap_evidence_df, outputpath)

    # Writing results to the legacy summary file (kept for backward-compat).
    #
    # detection is the sole writer of summary.txt. multimap.py used to open the
    # same path with mode "w" and write the three totals below, but it runs
    # earlier in the DAG, so this truncation destroyed them on every run and they
    # were never published. They are recomputed here instead of being passed
    # forward, which also keeps them correct after `rerun-multimap` swaps the
    # selected layer without re-running multimap (SW-04).
    headline = _headline_totals(adata, detected_viral_genes, detection_matrix)
    found_genes_sorted = dict(sorted(found_genes.items()))
    total_viral_genes = 0
    counts_per_virus = {}
    with (
        open(f"{config.output}/summary.txt", "w") as summary,
        open(f"{config.output}log/found_genes.txt", "w") as found_genes_file,
    ):
        summary.write(
            f"Viral molecules in unique-count matrix: {headline['unique']}\n"
            f"Total viral molecules (selected method): {headline['selected']}\n"
            f"Cells with viral reads: {headline['cells_with_virus']}/{headline['n_cells']}\n\n\n"
        )
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
                    f"\n{virus_name}: {stats['viral_molecules_total_est']} viral molecules estimated, "
                    f"{stats['infected_cells']}/{stats['total_cells']} cells infected "
                    f"({stats['pct_infected']:.2f}%), "
                    f"{stats['viral_molecules_per_10k_est']:.2f} viral molecules/10k estimated."
                )
            summary.write(f"\n\nTotal amount of viral load found: {total_viral_genes}")
            summary.write(
                f"\n\nOfficial name of viral load detected: {','.join(str(s) for s in detected_viral_genes)}"
            )
        else:
            summary.write("No viral gene IDs found in this sample for the viruses in the index.\n")
            # A zero is only interpretable next to the limit that produced it.
            # Without this line "no virus found" reads as an absence, which is
            # exactly the inference the depth and capture terms do not support.
            summary.write(
                negative_result_statement(
                    quantified_depth,
                    detection_threshold=int(config.detection_threshold),
                    capture=measured_capture if measured_capture is not None else 1.0,
                    capture_measured=measured_capture is not None,
                )
                + "\n"
            )
            if measured_capture is not None:
                summary.write(
                    "Positive control "
                    f"{control_detail.get('gene')}: "
                    f"{control_detail.get('observed_molecules')} of "
                    f"{control_detail.get('expected_molecules')} planted molecules "
                    f"recovered (capture={measured_capture:.4f}). This negative is "
                    "certifiable at that capture.\n"
                )
            else:
                summary.write(
                    "No positive control was supplied, so no k-mer capture term "
                    "could be measured and this negative CANNOT be read as "
                    "absence. See results/positive_control.json.\n"
                )
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
        sensitivity_df=sensitivity_df,
        sensitivity_statement=(
            negative_result_statement(
                quantified_depth, detection_threshold=int(config.detection_threshold)
            )
        ),
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
