"""
The detection scripts uses the gene IDs found in the analysis script and checks
for presence of the gene IDs in the sample (according to kb count). It also
Creates visualizations (except for the umap), e.g. a barplot and a plot showing
super expressors.
"""

# Importing packages
import base64
import csv
import datetime
import hashlib
import json
import logging
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any, NamedTuple, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.ticker import ScalarFormatter

from viralscan import anello_align
from viralscan.constants import SIBLING_CROSSMAP_RATIO_THRESHOLD
from viralscan.defaults import DEFAULTS
from viralscan.enrichment import cell_type_enrichment, write_cell_type_enrichment
from viralscan.multimapping import (
    select_detection_matrix,
    should_write_multimap_evidence,
    summarize_multimap_evidence,
    write_multimap_evidence,
)
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scripts.cellcalling import COMPARABLE_CELL_MIN_UMI, _strip_suffix
from viralscan.sensitivity import (
    SENSITIVITY_COLUMNS,
    CaptureScope,
    PositiveControl,
    SensitivityRecord,
    scoped_capture,
    sensitivity_record,
    sensitivity_statement,
    verify_exact_control,
)
from viralscan.utils import (
    infected_masks,
    matrix_for_genes,
    resolve_count_matrix,
    setup_script_logging,
    virus_molecules_per_cell,
)
from viralscan.virus_catalog import merged_name_map
from viralscan.virus_grouping import (
    HOST_HOMOLOGY_LOCUS_COLUMNS,
    HOST_HOMOLOGY_NOT_MEASURED,
    group_genes_by_identity,
    group_genes_by_virus,
    host_homology_by_virus,
    host_homology_loci,
    legacy_artifact_risk,
    legacy_claim_scope,
    legacy_eve_risk,
    legacy_sibling_groups,
    load_run_identity,
    virus_facts,
)

log = setup_script_logging()


# Run-level state, populated by run() from the Run Context. Declared here so the
# helper functions can reference them as module globals; the module imports
# cleanly without Snakemake because nothing reads these at import time.
config: RunConfig = RunConfig()
output: str = ""
kb = None
file: str = ""  # path to the viral-accessions list (the analysis rule's output)
# The Run's Virus Identity table (results/virus_identity.tsv); None for a run
# directory that has none, which selects the legacy prefix rules.
identity = None


def group_genes(gene_ids, table=None):
    """Group viral gene IDs by virus: by ``virus_key`` from the table, else legacy prefix."""
    if table is not None:
        return group_genes_by_identity(gene_ids, table)
    return group_genes_by_virus(gene_ids, merged_name_map())


def _in_barcode_set(obs_names, barcodes) -> np.ndarray:
    """Membership of ``obs_names`` in a host barcode set, ignoring a ``-1`` gem-group suffix.

    The host sets are suffix-stripped (cellcalling), so the raw obs_names must be too (B-10).
    """
    wanted = {_strip_suffix(str(b)) for b in barcodes}
    return np.array([_strip_suffix(str(b)) in wanted for b in obs_names], dtype=bool)


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


def detect_genes(var_names, matrix, viral_accessions, threshold=1, groups=None):
    """Return {gene_id: total_count} for viral genes meeting the threshold.

    Pure and importable. ``var_names`` may be a list or a pandas Index; counts
    are summed across all cells and compared with an inclusive ``>=`` threshold.

    With ``groups`` ({virus: [gene_id]}), the threshold applies to each virus's
    summed count (MECH-B) and every nonzero gene of a passing virus is reported,
    so a virus spread thinly over many genes is neither missed nor undercounted.
    Without it, each gene is its own group (the pre-v3 per-gene rule).
    """
    names = list(var_names)
    if len(set(names)) != len(names):
        dupes = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"detect_genes: duplicate var_names {dupes[:5]!r}; make them unique first")
    pos = {gene_id: i for i, gene_id in enumerate(names)}
    totals = {
        g: _gene_counts_from_matrix(matrix, pos[g]).sum() for g in viral_accessions if g in pos
    }
    bad = sorted(g for g, t in totals.items() if not np.isfinite(t))
    if bad:
        raise ValueError(f"detect_genes: non-finite counts (NaN/inf) in viral gene(s) {bad[:5]!r}")
    if groups is None:
        groups = {g: [g] for g in totals}
    found = {}
    for genes in groups.values():
        member = {g: totals[g] for g in genes if g in totals}
        if sum(member.values()) >= threshold:
            found.update({g: c for g, c in member.items() if c > 0})
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
        index_groups (dict): {virus: [gene_id]} over every indexed viral gene
    """
    global identity
    identity = load_run_identity(output)
    if identity is not None:
        # The table is the viral/host partition (and what multimap.py used).
        # Sorted so group order matches the analysis.txt order of older runs.
        viral_accessions = sorted(identity.viral_gene_ids())
    else:
        log.warning(
            "No results/virus_identity.tsv in %s: using log/analysis.txt and the "
            "legacy prefix naming rules.",
            output,
        )
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
    index_groups = group_genes([g for g in viral_accessions if g in adata.var_names], identity)[0]
    found_genes = detect_genes(
        adata.var_names,
        detection_matrix,
        viral_accessions,
        config.detection_threshold,
        groups=index_groups,
    )
    return adata, found_genes, output, viral_accessions, detection_matrix, index_groups


def histogram(adata, found_genes, identity_table, outputpath, viral_count_matrix=None):
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

    # Group found genes by virus (virus_key from the identity table).
    group_by_virus, detected_viral_genes = group_genes(found_genes, identity_table)

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
            ax = plt.gca()
            ax.bar(df_virus_sorted["gene_id"].astype(str), df_virus_sorted["UMI_count"])

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
    cells) have a viral UMI count >= ``config.se_threshold`` (inclusive).
    The caller's ``adata`` is not modified.
    ---------------------------------------------------------------------
    Parameters:
        adata (anndata.AnnData): h5ad file of kb-python used for further
            analysis
        virus (str): virus where plot needs to be made for
        viral_gene_ids (list): gene IDs belonging to the virus param
        outputpath (str): path to the output directory defined by the user
    ---------------------------------------------------------------------
    Raises:
        ValueError: None of the provided viral gene IDs were found in the dataset,
            or the dataset has no cells
    """
    if adata.n_obs == 0:
        raise ValueError("super_expressor: the dataset has no cells, nothing to plot.")
    # Derived columns stay local: the caller's AnnData is not touched (B-09).
    n_count_rna = pd.Series(_sum_axis(adata.X, 1), index=adata.obs_names)

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
    viral_umi = pd.Series(
        _sum_axis(matrix_for_genes(adata, count_matrix, list(matched_genes)), 1),
        index=adata.obs_names,
    )

    # Null Model (grey line)
    total_viral = viral_umi.sum()
    total_rna = n_count_rna.sum()
    null_model = total_viral * (n_count_rna / total_rna)

    obs_sorted = np.sort(viral_umi.values)[::-1]
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
    n_SE = (viral_umi >= se_threshold).sum()

    title = f"Virus {virus}, n_super={n_SE}\nn={adata.n_obs}; {virus} max={viral_umi.max():g}"

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


def _headline_totals(
    adata,
    detected_groups,
    viral_count_matrix=None,
    infected_min_molecules=1.0,
    called_mask=None,
) -> dict:
    """Totals that head summary.txt, for the viruses that passed the gate only.

    *detected_groups* is ``{virus name: [gene IDs]}`` of the viruses whose ``detected``
    flag is true (AUD-07). A candidate-only virus is not headlined: its molecules and
    cells stay in its own ``viral_summary.tsv`` row and per-virus line. "Cells with
    virus" follows the one infected rule (AUD-06, :func:`viralscan.utils.infected_masks`,
    ``--infected-min-molecules``) per virus and takes the union, so it agrees with the
    ``infected_cells`` / ``infected_called`` columns instead of counting any weight > 0.
    ``called_mask=None`` means every barcode is a called cell.

    Deriving these here rather than carrying them forward from multimap keeps
    them consistent with whichever allocation layer the H5AD holds, including
    after ``rerun-multimap`` swaps that layer in place.

    The gene lists are gene IDs (``var_names``), not virus names. The caller
    used to pass the detected virus *names*, none of which is a ``var_name``, so
    every summary.txt read "Viral molecules ... 0" beside a non-zero
    viral_summary.tsv (SW-16). IDs that match nothing now raise instead.
    """
    requested = [gene for genes in detected_groups.values() for gene in genes]
    genes = [gene for gene in dict.fromkeys(requested) if gene in adata.var_names]
    n_cells = int(adata.n_obs)
    empty = {
        "unique": 0,
        "selected": 0,
        "cells_with_virus": 0,
        "cells_with_virus_called": 0,
        "n_cells": n_cells,
    }
    if requested and not genes:
        raise ValueError(
            "_headline_totals expects gene IDs from adata.var_names; none of "
            f"{requested[:3]!r} is one (were virus names passed?)"
        )
    if not genes:
        return empty

    selected_matrix = resolve_count_matrix(viral_count_matrix, adata)
    selected = matrix_for_genes(adata, selected_matrix, genes)
    unique_layer = adata.layers.get("counts_unique")
    unique = matrix_for_genes(adata, unique_layer, genes) if unique_layer is not None else None

    infected_any = np.zeros(n_cells, dtype=bool)
    for virus_genes in detected_groups.values():
        in_matrix = [g for g in virus_genes if g in adata.var_names]
        if in_matrix:
            per_cell = virus_molecules_per_cell(adata, selected_matrix, in_matrix)
            infected_any |= infected_masks(per_cell, infected_min_molecules)[0]
    called = np.ones(n_cells, dtype=bool) if called_mask is None else np.asarray(called_mask, bool)
    return {
        "unique": _count_value(np.asarray(_sum_axis(unique, 1)).sum()) if unique is not None else 0,
        "selected": _count_value(np.asarray(_sum_axis(selected, 1)).sum()),
        "cells_with_virus": int(infected_any.sum()),
        "cells_with_virus_called": int((infected_any & called).sum()),
        "n_cells": n_cells,
    }


def _virus_summary_line(virus_name, stats) -> str:
    """One per-virus line of summary.txt, led by the called-cell numbers (B-05).

    ``pct_infected_called`` is the primary denominator; the all-barcode figure is
    diluted by empty droplets (millions on a raw matrix), so it is labelled as such
    and comes second. No semicolon: compare_legacy_v2_v3 parses any ``x;number``
    line as a gene total.
    """
    line = (
        f"{virus_name}: {stats['viral_molecules_total_est']} viral molecules estimated, "
        f"{stats['infected_called']}/{stats['n_called_cells']} called cells with viral reads "
        f"({stats['pct_infected_called']:.2f}%), "
        f"{stats['infected_cells']}/{stats['total_cells']} all barcodes "
        f"({stats['pct_infected']:.2f}%), "
        f"{stats['viral_molecules_per_10k_est']:.2f} viral molecules/10k estimated."
    )
    # AUD-06/AUD-07: fractional evidence and empty-droplet molecules stay visible.
    if stats.get("infected_subthreshold_cells"):
        line += (
            f" {stats['infected_subthreshold_cells']} cells hold less than the infected-cell "
            "minimum (sub-threshold)."
        )
    if stats.get("background_umis"):
        line += f" {stats['background_umis']} molecules sit in non-called barcodes (background)."
    if stats.get("ambient_assessment") in ("soup_sufficient", "soup_rejected_rho1"):
        line += (
            f" Ambient screen (rho = 1 bound, descriptive): {stats['ambient_assessment']}"
            f", upper tail p = {stats['soup_tail_p_cells']:g}."
        )
    if not stats.get("detected", True):
        line += " CANDIDATE ONLY, the called-cell gate is not met" + (
            " (no infected called cell)." if not stats["infected_called"] else "."
        )
    return line


def _gate_note(config, method, candidate_only) -> str:
    """The summary.txt paragraph that says what "detected" means and what is only a candidate."""
    note = (
        "Detection gate: a virus is detected when at least "
        f"{config.detection_min_cells} called cell(s) each hold at least "
        f"{config.infected_min_molecules:g} of its molecules and the called cells hold at least "
        f"{config.detection_threshold} in total. The gene table above is the all-barcode "
        "candidate screen. Molecules in non-called barcodes are background. A passed gate is "
        "candidate evidence, it does not rule out ambient RNA or barcode swapping."
    )
    if method == "none":
        note += (
            " Cell calling is off, so every barcode counts as a called cell and empty droplets "
            "cannot be told apart from cells."
        )
    if candidate_only:
        note += (
            f"\nCandidate only, not detected (the called-cell gate is not met): "
            f"{', '.join(candidate_only)}"
        )
    return note


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


#: Confidence level of the upper bound on the empty-droplet viral fraction, and the one-sided
#: tail level of the soup screen. Conventional values, not tuned: the screen is descriptive
#: and gates nothing (M-03, docs/plans/ambient-and-barcode-swapping-controls.md section 6).
AMBIENT_UB_LEVEL = 0.95
AMBIENT_TAIL_ALPHA = 0.05

#: The per-virus columns of the ambient screen, in ``viral_summary.tsv`` order.
AMBIENT_COLUMNS = (
    "ambient_viral_fraction",
    "ambient_viral_fraction_ub95",
    "n_ambient_barcodes_with_virus",
    "soup_expected_cells_max",
    "soup_expected_cells_ge2_max",
    "soup_tail_p_cells",
    "viral_depth_spearman",
    "ambient_assessment",
    "ambient_basis",
)

_SWAP_NOT_ASSESSED = (
    "Swapping between libraries is not assessed (scripts/swap_check.py needs a pool manifest)."
)


def _ambient_not_assessed(reason: str) -> dict[str, Any]:
    """The ambient columns of a virus (or run) that cannot be assessed: blank, never 0."""
    row: dict[str, Any] = {c: None for c in AMBIENT_COLUMNS}
    row["ambient_assessment"] = "not_assessed"
    row["ambient_basis"] = f"{reason}. {_SWAP_NOT_ASSESSED}"
    return row


def _poisson_ub(count: float, level: float = AMBIENT_UB_LEVEL) -> float:
    """One-sided exact upper bound on a Poisson mean given *count* events (gamma form).

    ``count = 0`` gives ``-ln(1 - level) = 2.996``, the constant ``sensitivity.py`` uses for LOD95.
    The gamma form also takes a fractional *count* (EM or multimap weight).
    """
    from scipy.stats import gamma

    return float(gamma.ppf(level, count + 1.0))


def soup_screen(
    viral_ambient: float,
    total_ambient: float,
    called_totals,
    infected_called: int,
    infected_min_molecules: float,
    ambient_max_molecules: float,
) -> dict[str, Any]:
    """The calibration-free ``rho = 1`` soup bound for one virus (M-03 design section 3.1).

    ``viral_ambient`` and ``total_ambient`` are the viral and total molecules in the ambient
    barcodes (not called, total at or below ``ambient_max_molecules``). With ``f_U`` the 95 %
    Poisson upper bound on ``viral_ambient / total_ambient``, a called cell of ``N`` molecules
    receives at most ``Poisson(f_U * N)`` soup molecules at any contamination fraction
    ``rho <= 1``. ``soup_expected_cells_max`` (``E1``) sums the chance each called cell reaches
    the infected rule from soup alone and ``soup_expected_cells_ge2_max`` (``E2``) the chance
    of two or more. ``soup_tail_p_cells`` is the upper Poisson tail of the observed
    ``infected_called`` given ``E1``. Descriptive only: ``detected`` never reads it.
    """
    from scipy.stats import poisson

    f_hat = viral_ambient / total_ambient
    f_upper = _poisson_ub(viral_ambient) / total_ambient
    totals = np.asarray(called_totals, dtype=float)
    # Molecules are whole: a cell reaches a fractional infected rule at ceil(rule) molecules.
    k = max(1, int(np.ceil(infected_min_molecules - 1e-9)))
    lam = f_upper * totals
    e1 = float(poisson.sf(k - 1, lam).sum())
    e2 = float(poisson.sf(max(2, k) - 1, lam).sum())
    tail = float(poisson.sf(infected_called - 1, e1)) if infected_called >= 1 else 1.0
    return {
        "ambient_viral_fraction": float(f"{f_hat:.6g}"),
        "ambient_viral_fraction_ub95": float(f"{f_upper:.6g}"),
        "soup_expected_cells_max": float(f"{e1:.6g}"),
        "soup_expected_cells_ge2_max": float(f"{e2:.6g}"),
        "soup_tail_p_cells": float(f"{tail:.6g}"),
        "ambient_assessment": (
            "soup_rejected_rho1" if tail < AMBIENT_TAIL_ALPHA else "soup_sufficient"
        ),
        "ambient_basis": (
            f"rho = 1 soup bound: {AMBIENT_UB_LEVEL:.0%} Poisson upper bound on the viral fraction "
            f"of non-called barcodes with at most {ambient_max_molecules:g} molecules, upper "
            f"tail level {AMBIENT_TAIL_ALPHA:g}. Descriptive: it does not show infection and "
            f"a pass is not a clean result. {_SWAP_NOT_ASSESSED}"
        ),
    }


def _depth_spearman(viral, totals, n_infected_called: int) -> Optional[float]:
    """Spearman correlation of viral with total molecules over the called cells, else None."""
    if n_infected_called < 3:
        return None
    from scipy.stats import spearmanr

    rho = spearmanr(totals, viral)[0]
    return None if not np.isfinite(rho) else round(float(rho), 4)


class _AmbientContext(NamedTuple):
    """Per-barcode totals behind the ambient screen, built once per run."""

    reason: Optional[str]  # why the run cannot be assessed; None when it can
    ambient: Optional[np.ndarray] = None  # mask over adata barcodes
    total_ambient: float = 0.0  # molecules in all ambient barcodes, incl. ones absent from adata
    called_extra: Optional[np.ndarray] = None  # totals of called barcodes absent from adata
    totals: Optional[np.ndarray] = None  # total molecules per adata barcode


def _host_gene_mask(adata) -> Optional[np.ndarray]:
    """Host-gene column mask, or None when the matrix carries no host depth."""
    if (
        adata.uns.get("index_kind") == "virus_only"
        or adata.uns.get("gene_identity_source") == "unavailable"
    ):
        return None
    if "gene_role" in adata.var:
        host = np.asarray(adata.var["gene_role"] == "host")
    elif "is_viral" in adata.var or "viral" in adata.var:
        column = "is_viral" if "is_viral" in adata.var else "viral"
        host = ~np.asarray(adata.var[column], dtype=bool)
    else:
        return None
    return host if host.any() else None


def _ambient_context(
    adata, called_mask, total_per_cell, host_cells, cell_calling_method, ambient_max_molecules
) -> _AmbientContext:
    """Which barcodes are ambient and how deep every barcode is, or why neither is known.

    A barcode's depth is its host molecules plus its row total in the viral matrix. A combined
    run has both in ``adata.X``. A two-step run takes the host part from STARsolo ``GeneFull/raw``
    (``host_cells["host_totals"]``), so barcodes the viral matrix lacks (every called cell, and
    the ambient droplets, with no non-host read) still count. Virus-only indexes, ``--cell-calling
    none``/``knee`` and runs with no host depth are not assessed.
    """
    if cell_calling_method == "none":
        return _AmbientContext("cell calling is off, so there are no empty droplets to measure")
    if cell_calling_method == "knee":
        return _AmbientContext(
            "knee cell calling is sensitivity-only and calls empty droplets as cells"
        )
    totals = np.asarray(total_per_cell, dtype=float)
    extra_called = np.zeros(0)
    extra_ambient = 0.0
    if host_cells:
        host_totals = host_cells.get("host_totals")
        if host_totals is None:
            return _AmbientContext(
                "no STARsolo GeneFull/raw host totals for this run, host depth is unknown"
            )
        names = [_strip_suffix(str(b)) for b in adata.obs_names]
        present = set(names)
        totals = totals + np.array([host_totals.get(b, 0.0) for b in names])
        extra_called = np.array(
            [host_totals[b] for b in host_cells["called"] if b in host_totals and b not in present]
        )
        extra_ambient = sum(
            t
            for b, t in host_totals.items()
            if b not in present and b not in host_cells["called"] and t <= ambient_max_molecules
        )
    elif _host_gene_mask(adata) is None:
        return _AmbientContext(
            "the index has no host genes, so there is no host depth to size empty droplets"
        )
    called = np.asarray(called_mask, dtype=bool)
    ambient = ~called & (totals <= ambient_max_molecules)
    total_ambient = float(totals[ambient].sum()) + float(extra_ambient)
    if total_ambient <= 0:
        return _AmbientContext(
            f"no non-called barcode holds at most {ambient_max_molecules:g} molecules"
        )
    if not (called.any() or len(extra_called)):
        return _AmbientContext("no called cell")
    return _AmbientContext(None, ambient, total_ambient, extra_called, totals)


def compute_stats(
    adata,
    found_genes,
    group_by_virus,
    detected_viral_genes,
    called_mask=None,
    viral_count_matrix=None,
    index_genes_by_virus=None,
    host_cells=None,
    infected_min_molecules=DEFAULTS["infected_min_molecules"],
    detection_threshold=DEFAULTS["detection_threshold"],
    detection_min_cells=DEFAULTS["detection_min_cells"],
    cell_calling_method=None,
    ambient_max_molecules=DEFAULTS["emptydrops_lower"],
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
    index_genes_by_virus : dict[str, list[str]] | None
        Every index gene of each virus (ANDET-01). ``accession_breadth`` is the
        fraction of these with >= 1 molecule. ``None`` falls back to the
        detected genes, where breadth is 1.0 by construction.

    host_cells : dict[str, set[str] | None] | None
        Two-step runs only: ``{"called", "comparable"}`` barcode sets from the host
        matrix. The viral matrix holds only barcodes with non-host reads, so its own
        rows cannot give the cell denominators; ``n_called_cells`` and
        ``n_comparable_cells`` come from these sets (the 263/263 = 100 % trap).
        A missing ``comparable`` set means host depth is unavailable, so the
        comparable statistics are unavailable rather than viral-depth estimates.

    cell_calling_method, ambient_max_molecules :
        The ambient screen (M-03, descriptive, never read by ``detected``): the resolved cell-calling
        method (``none`` and ``knee`` are not assessed) and the total-molecule ceiling of an
        ambient barcode (``--emptydrops-lower``). See :func:`soup_screen`.
        A two-step run needs ``host_cells["host_totals"]`` (per-barcode STARsolo host molecules)
        for it, else the ambient columns read ``not_assessed``.

    infected_min_molecules : float
        The one infected rule (AUD-06, :func:`viralscan.utils.infected_masks`): a cell is
        infected for a virus when its summed allocated molecules reach this value. Cells
        with ``0 < weight < infected_min_molecules`` are counted in
        ``infected_subthreshold_cells`` and do not enter ``per_cell_df``.
    detection_threshold, detection_min_cells : int
        The called-cell gate (AUD-07). ``detected`` is True when at least
        ``detection_min_cells`` *called* cells are infected AND the called cells hold at
        least ``detection_threshold`` viral molecules. Molecules in barcodes that are not
        called cells (empty droplets) are ``background_umis``; they never count toward
        the gate. With ``called_mask=None`` every barcode is a called cell.

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
    n_called = len(host_cells["called"]) if host_cells else int(called_mask.sum())

    # Total UMI per cell (sum across all genes)
    total_umi_per_cell = _sum_axis(adata.X, 1)
    total_umi_all = total_umi_per_cell.sum()

    virus_stats = {}
    cell_rows = []
    count_matrix = resolve_count_matrix(viral_count_matrix, adata)
    ambient_ctx = _ambient_context(
        adata,
        called_mask,
        total_umi_per_cell,
        host_cells,
        cell_calling_method,
        ambient_max_molecules,
    )

    for virus, gene_list in group_by_virus.items():
        valid_genes = [g for g in gene_list if g in adata.var_names]
        if not valid_genes:
            continue

        # Per-cell viral UMI for this virus
        viral_umi_per_cell = virus_molecules_per_cell(adata, count_matrix, valid_genes)

        total_umi_raw = float(viral_umi_per_cell.sum())
        infected_mask, subthreshold_mask = infected_masks(
            viral_umi_per_cell, infected_min_molecules
        )
        infected_cells = int(infected_mask.sum())
        pct_infected = round(infected_cells / total_cells * 100, 4) if total_cells else 0.0
        umi_per_10k = round(total_umi_raw / total_umi_all * 10_000, 4) if total_umi_all else 0.0

        # Same stats restricted to called cells (real, non-empty droplets) — the
        # primary denominator. Empty droplets otherwise dilute pct_infected.
        infected_called = int((infected_mask & called_mask).sum())
        pct_infected_called = round(infected_called / n_called * 100, 4) if n_called else 0.0
        # Called-cell gate (AUD-07): empty-droplet molecules are background, not evidence.
        called_umi = float(viral_umi_per_cell[called_mask].sum())
        background_umi = max(total_umi_raw - called_umi, 0.0)
        detected = bool(
            infected_called >= detection_min_cells and called_umi >= detection_threshold
        )

        # A fixed, strategy-independent denominator. See comparable_called_cells.
        if host_cells and host_cells["comparable"] is not None:
            comparable = _in_barcode_set(adata.obs_names, host_cells["comparable"])
            n_comparable = len(host_cells["comparable"])
        elif host_cells:
            comparable = None
            n_comparable = None
        else:
            comparable = _comparable_called_cells(adata, called_mask)
            n_comparable = int(comparable.sum()) if comparable is not None else None
        infected_comparable = (
            int((infected_mask & comparable).sum()) if comparable is not None else None
        )
        if n_comparable is None or infected_comparable is None:
            pct_infected_comparable = None
        else:
            # no comparable cell: the percentage is undefined, not 0 (a 0 reads as "none infected")
            pct_infected_comparable = (
                round(infected_comparable / n_comparable * 100, 4) if n_comparable else None
            )

        # Accession breadth: fraction of the virus's *index* genes with >= 1
        # molecule in any cell. EVE artefacts concentrate on 1-2 loci; a genuine
        # infection spreads across ORF1/ORF2/ORF3. Over the detected genes alone
        # it was 1.0 by construction (ANDET-01).
        index_genes = (index_genes_by_virus or {}).get(virus) or valid_genes
        in_matrix = [g for g in index_genes if g in adata.var_names]
        n_acc_detected = 0
        if in_matrix:
            index_counts = matrix_for_genes(adata, count_matrix, in_matrix)
            if hasattr(index_counts, "toarray"):
                index_counts = index_counts.toarray()
            n_acc_detected = int(np.asarray((index_counts > 0).any(axis=0)).sum())
        n_acc_total = len(set(index_genes))
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

        if ambient_ctx.reason is not None:
            ambient = _ambient_not_assessed(ambient_ctx.reason)
        else:
            assert ambient_ctx.ambient is not None and ambient_ctx.totals is not None
            assert ambient_ctx.called_extra is not None
            in_a = ambient_ctx.ambient
            called_totals = np.concatenate(
                [ambient_ctx.totals[called_mask], ambient_ctx.called_extra]
            )
            ambient = soup_screen(
                float(viral_umi_per_cell[in_a].sum()),
                ambient_ctx.total_ambient,
                called_totals,
                infected_called,
                infected_min_molecules,
                ambient_max_molecules,
            )
            ambient["n_ambient_barcodes_with_virus"] = int((viral_umi_per_cell[in_a] > 0).sum())
            ambient["viral_depth_spearman"] = _depth_spearman(
                np.concatenate(
                    [viral_umi_per_cell[called_mask], np.zeros(len(ambient_ctx.called_extra))]
                ),
                called_totals,
                infected_called,
            )

        virus_stats[virus] = {
            **ambient,
            "viral_molecules_total_est": _count_value(total_umi_raw),
            "infected_cells": infected_cells,
            "total_cells": total_cells,
            "pct_infected": pct_infected,
            "viral_molecules_per_10k_est": umi_per_10k,
            "n_called_cells": n_called,
            "infected_called": infected_called,
            "pct_infected_called": pct_infected_called,
            "infected_subthreshold_cells": int(subthreshold_mask.sum()),
            "infected_subthreshold_called": int((subthreshold_mask & called_mask).sum()),
            "background_umis": _count_value(background_umi),
            "detected": detected,
            "n_comparable_cells": n_comparable,
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
    negative is informative — under the substitution-only heuristic it falls to
    0.063 at 15 % divergence and 0.013 at 20 % (90 bp, k=31) — and it is
    invisible in a count matrix.

    Returns ``(capture, detail_dict)``. ``capture`` is ``None`` when no control
    was configured, in which case the caller must treat every negative as
    uncertifiable. A measured ratio only speaks for the control's declared scope
    (:func:`control_claim`); use :func:`positive_control_from` and
    :func:`certified_viruses` to learn which rows it covers.
    """
    gene = getattr(config, "positive_control_gene", None)
    expected = getattr(config, "positive_control_expected_molecules", None)
    if not gene or expected is None:
        return None, {"status": "not-configured"}
    if not float(expected) > 0:
        raise ValueError(
            f"positive_control_expected_molecules must be > 0, got {expected!r}: the capture "
            "is observed / expected, so a non-positive plant has no meaning."
        )
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
        "implied_divergence": substitution_model_implied_divergence(capture),
        "implied_divergence_note": IMPLIED_DIVERGENCE_NOTE,
    }
    if getattr(config, "positive_control_scope", None) == "exact_sequence":
        table = load_run_identity(config.output) if config.output else identity
        layer = "X" if matrix is adata.X else "unverified_alternate_matrix"
        detail.update(verify_exact_control(config, table, layer))
        if "is_viral" in adata.var and not bool(adata.var.loc[gene, "is_viral"]):
            detail.update(
                verification_status="failed",
                verification_detail="matrix and VirusIdentityTable disagree on the control gene",
            )
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
    implied = detail["implied_divergence"]
    log.info(
        "Positive control %s: %g/%g molecules recovered -> capture=%.4f "
        "(substitution-only heuristic, not a measurement: %s)",
        gene,
        observed,
        expected,
        capture,
        "no divergence identifiable"
        if implied is None
        else f"i.i.d.-substitution model, 90 bp, k=31, would match ~{implied * 100:.1f}% divergence",
    )
    return capture, detail


#: Divergence beyond which the implied-divergence inversion is abandoned.
#: The exact substitution-only capture is ~1.4e-8 at 0.5 (90 bp, k=31) and keeps
#: falling smoothly, but a recovery ratio that small from a planted control is
#: indistinguishable from zero at any realistic spike-in size, so a control that
#: recovers less than ``fragment_capture_exact(0.5)`` has no identifiable implied
#: divergence and is reported as such rather than given a confident number.
MAX_IDENTIFIABLE_DIVERGENCE = 0.5

#: Wording attached to every reported implied divergence.
IMPLIED_DIVERGENCE_NOTE = (
    "substitution-only heuristic: the per-base divergence at which an i.i.d. "
    "substitution model (90 bp, k=31) would give this recovery ratio. This is not a measured "
    "genome divergence; indels, sequencing error and index competition also lower recovery."
)


def substitution_model_implied_divergence(capture: float, read_length: int = 90) -> Optional[float]:
    """Divergence at which the substitution-only model gives ``capture``, or None.

    A **substitution-only heuristic** (SENS-CORR-02), not a measurement of genome
    divergence: it inverts :func:`viralscan.sensitivity.fragment_capture_exact` by
    bisection and attributes *all* observed loss to i.i.d. substitutions. It never
    sets ``capture_measured``. A capture of 1.0 means "no loss measurable", which is
    reported as None rather than 0 % so a reader does not over-read it as proof of
    zero divergence. A capture at or below the value reachable at
    :data:`MAX_IDENTIFIABLE_DIVERGENCE` is not identifiable and also returns None.
    """
    from viralscan.sensitivity import DEFAULT_K, fragment_capture_exact

    if capture >= 0.999 or capture <= 0:
        return None
    lo, hi = 0.0, MAX_IDENTIFIABLE_DIVERGENCE
    if capture <= fragment_capture_exact(hi, read_length=read_length, k=DEFAULT_K):
        return None  # below anything the model can attribute to divergence
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if fragment_capture_exact(mid, read_length=read_length, k=DEFAULT_K) > capture:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


# ----------------------------------------------------------- capture scope
#: One-time guard for the legacy-control warning (module state, per process).
_LEGACY_SCOPE_WARNED = False


def control_claim(config):
    """``(scope, target)`` a configured positive control is entitled to speak for.

    A control configured without ``positive_control_scope`` (every run before
    SENS-CORR-01) maps to ``panel_mechanics``: it shows the pipeline recovers a
    planted molecule, not that it can see any other virus, so it certifies no
    row. That is logged once per process. No control means no claim and no
    warning.
    """
    global _LEGACY_SCOPE_WARNED
    scope = getattr(config, "positive_control_scope", None)
    if scope:
        return CaptureScope(scope), getattr(config, "positive_control_virus_key", None)
    if getattr(config, "positive_control_gene", None) and not _LEGACY_SCOPE_WARNED:
        _LEGACY_SCOPE_WARNED = True
        log.warning(
            "positive_control_gene is set without positive_control_scope: treating it "
            "as panel_mechanics, which certifies NO virus. Pass --positive-control-scope "
            "exact_sequence with --positive-control-virus-key to let it certify one row."
        )
    return CaptureScope.PANEL_MECHANICS, None


def positive_control_from(config, control_detail):
    """The :class:`PositiveControl` a run measured, or None if there is no usable one.

    Built from the *raw* recovery ratio and only when ``status == "measured"``:
    ``measure_positive_control`` clamps an over-recovered control to 1.0, and that
    clamp must never be mistaken for a measured full recovery.
    """
    if not control_detail or control_detail.get("status") != "measured":
        return None
    scope, target = control_claim(config)
    if scope is CaptureScope.EXACT_SEQUENCE:
        if (
            control_detail.get("verification_status") != "verified"
            or control_detail.get("certification_status") != "single_accession_group"
        ):
            return None
        target = control_detail["resolved_target"]
    return PositiveControl(capture=control_detail.get("capture"), scope=scope, target=target)


def row_capture(control, virus):
    """Capture that applies to row ``virus`` (None outside the control's scope)."""
    if control is None:
        return None
    return scoped_capture(control, virus, control.scope)


def certified_viruses(control, viruses):
    """Sorted rows (of ``viruses``) whose capture was measured in scope."""
    return sorted({v for v in viruses if row_capture(control, v) is not None})


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
    host = _host_gene_mask(adata)
    if host is None:
        return None
    total = _sum_axis(adata.X[:, host], 1)
    return (np.asarray(total) >= COMPARABLE_CELL_MIN_UMI) & np.asarray(called_mask, dtype=bool)


#: ``positive_control.json`` status -> the ``capture_status`` of a row that got no
#: capture from it. A ``measured`` control whose scope does not cover the row is
#: "control-out-of-scope" (unrelated virus or panel_mechanics), never a borrowed value.
_CAPTURE_STATUS_BY_CONTROL = {
    "failed": "control-failed",
    "over-recovered": "control-over-recovered",
    "gene-not-in-reference": "control-gene-missing",
    "measured": "control-out-of-scope",
}


def build_sensitivity_records(
    adata,
    virus_stats,
    config,
    depth=None,
    control_detail=None,
    index_viruses=(),
    observed_below_gate=None,
) -> list[SensitivityRecord]:
    """The canonical per-virus :class:`SensitivityRecord` list for this run.

    ``sensitivity.tsv``, ``summary.txt``, ``report.html`` and the
    ``positive_control.json`` target lists all read this one list (SENS-CORR-03).
    A virus in ``virus_stats`` cleared the detection gate (that membership *is* the
    detection decision); every other indexed virus is a negative, carrying
    ``observed_below_gate[virus]`` molecules when the caller measured a nonzero
    count below the gate.

    ``depth`` is total quantified molecules (sum of ``adata.X``), which is the
    only depth term a run can actually measure: only quantified molecules can
    be detected, so raw read counts would overstate sensitivity.

    ``control_detail`` is the dict from :func:`measure_positive_control` (measured
    here when omitted). Its capture is applied **only to rows inside the control's
    declared scope** (:func:`row_capture`); every other row gets no capture term,
    so ``capture_measured`` is False and the LOD is the labelled depth-only floor
    (effective capture 1.0, never a measured or borrowed value).
    """
    if depth is None:
        depth = float(_sum_axis(adata.X, 1).sum())
    if control_detail is None:
        _, control_detail = measure_positive_control(adata, config, depth=depth)
    control = positive_control_from(config, control_detail)
    below = observed_below_gate or {}
    out_of_scope = _CAPTURE_STATUS_BY_CONTROL.get(control_detail.get("status"))
    undetected = sorted(set(index_viruses) - set(virus_stats))
    rows = [(v, s, True) for v, s in virus_stats.items()] + [(v, {}, False) for v in undetected]
    records = []
    for virus, stats, detected in rows:
        observed = float(stats.get("viral_molecules_total_est", 0) or 0)
        if not detected:
            observed = float(below.get(virus, 0) or 0)
        capture = row_capture(control, virus)
        notes: tuple[str, ...] = ()
        if capture is not None and control is not None:
            notes = (
                f"capture measured on the {control.scope.value} positive control for "
                f"{control.target!r} only, accession {control_detail['accession_version']} "
                f"with pinned sequence SHA-256 {control_detail['sequence_sha256']}; "
                "it does not cover any other virus or untested genome in this group",
            )
        records.append(
            sensitivity_record(
                virus,
                observed_molecules=observed,
                depth=depth,
                detection_threshold=int(config.detection_threshold),
                capture=capture if capture is not None else 1.0,
                capture_measured=capture is not None,
                notes=notes,
                detected=detected,
                capture_status=None if capture is not None else out_of_scope,
            )
        )
    return records


def build_sensitivity_table(
    adata,
    virus_stats,
    config,
    depth=None,
    control_detail=None,
    index_viruses=(),
    observed_below_gate=None,
):
    """The ``sensitivity.tsv`` frame: :func:`build_sensitivity_records`, flattened.

    Every virus that cleared the detection threshold gets a row, and so does every
    virus in the index that did not (``index_viruses``): a virus absent from the
    table is a virus nobody asked about, and a zero in ``observed_molecules`` is a
    negative whose meaning depends on the LOD columns beside it.
    """
    records = build_sensitivity_records(
        adata,
        virus_stats,
        config,
        depth=depth,
        control_detail=control_detail,
        index_viruses=index_viruses,
        observed_below_gate=observed_below_gate,
    )
    return pd.DataFrame([r.as_row() for r in records], columns=list(SENSITIVITY_COLUMNS))


class SensitivityOutputs(NamedTuple):
    """Everything every sensitivity renderer reads, computed once from the records.

    ``table`` is the ``sensitivity.tsv`` frame, ``interpretation`` and ``control_line``
    are the text for ``summary.txt`` and ``report.html``, and the two target tuples
    feed ``positive_control.json``. Nothing downstream re-derives any of them.
    """

    records: list[SensitivityRecord]
    table: pd.DataFrame
    interpretation: str
    control_line: str
    control: Optional[PositiveControl]
    certified_targets: tuple[str, ...]
    informative_targets: tuple[str, ...]


def build_sensitivity_outputs(
    adata, virus_stats, config, depth, control_detail, index_viruses=(), count_matrix=None
) -> SensitivityOutputs:
    """Build the records and every rendering input derived from them (SENS-CORR-03)."""
    records = build_sensitivity_records(
        adata,
        virus_stats,
        config,
        depth=depth,
        control_detail=control_detail,
        index_viruses=index_viruses,
        observed_below_gate=below_gate_counts(adata, count_matrix, index_viruses, virus_stats),
    )
    control = positive_control_from(config, control_detail)
    certified = tuple(sorted(r.virus_name for r in records if r.capture_measured))
    return SensitivityOutputs(
        records=records,
        table=pd.DataFrame([r.as_row() for r in records], columns=list(SENSITIVITY_COLUMNS)),
        interpretation=sensitivity_statement(
            records, depth, detection_threshold=int(config.detection_threshold)
        ),
        control_line=control_summary_line(control_detail, control, certified),
        control=control,
        certified_targets=certified,
        informative_targets=tuple(sorted(r.virus_name for r in records if r.informative)),
    )


def control_summary_line(control_detail, control, certified_targets):
    """One factual line about the positive control, shared by summary.txt and report.html.

    States what was recovered and which rows the capture was applied to. It makes no
    certifiability claim: whether a row is an informative negative is in the records
    (and in the interpretation text built from them), not here.
    """
    status = control_detail.get("status")
    if status == "not-configured":
        return (
            "No positive control was supplied, so no k-mer capture term could be "
            "measured. See results/positive_control.json."
        )
    scope = control.scope.value if control is not None else "no usable scope"
    line = (
        f"Positive control {control_detail.get('gene')}: "
        f"{control_detail.get('observed_molecules')} of "
        f"{control_detail.get('expected_molecules')} planted molecules recovered "
        f"(status: {status}, scope: {scope})."
    )
    if certified_targets:
        line += (
            f" Its capture ({control.capture:.4f}) was applied to {', '.join(certified_targets)} "
            "only; no other virus is covered."
        )
    else:
        line += " It supplied a capture term to no row of this run's table."
    return line + " See results/positive_control.json."


def below_gate_counts(adata, count_matrix, index_groups, detected):
    """Nonzero molecules per *undetected* indexed virus, in the detection count layer.

    Read-only over the existing count matrix (no allocation logic): it lets an
    undetected row carry the honest "N molecules, below the gate" instead of an
    inferred zero.

    ``index_groups`` is ``{virus: [gene_id]}``. A bare name list (accepted by
    :func:`build_sensitivity_records`) names no genes, so it yields ``{}`` (B-06).
    """
    if not isinstance(index_groups, Mapping):
        return {}
    totals = pd.Series(
        _sum_axis(resolve_count_matrix(count_matrix, adata), 0), index=adata.var_names
    )
    counts = {
        v: float(totals.reindex(genes).sum())
        for v, genes in index_groups.items()
        if v not in detected
    }
    return {v: n for v, n in counts.items() if n > 0}


def write_sensitivity_table(sensitivity_df, outputpath):
    """Write results/sensitivity.tsv. Returns its path."""
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, "sensitivity.tsv")
    sensitivity_df.to_csv(path, sep="\t", index=False)
    log.info("Wrote results/sensitivity.tsv (%d virus row(s))", len(sensitivity_df))
    return path


def write_control_report(
    control_detail,
    outputpath,
    control=None,
    certified_targets=(),
    informative_negative_targets=None,
):
    """Write results/positive_control.json. Returns its path.

    ``certified_targets`` lists the rows whose capture was measured in the control's
    scope (``capture_measured`` in the sensitivity records); it is never panel-wide.
    ``informative_negative_targets`` lists the subset that is also an informative
    *negative* (undetected, with adequate adjusted depth), straight from the
    records. ``certifies_negatives`` is true only when that list is non-empty; a
    caller with no records omits it and gets the older "any certified row" meaning.
    ``capture_used_for_sensitivity`` is the control's raw ratio, and only when a
    row actually used it.
    """
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    payload = dict(control_detail)
    payload["schema_version"] = "3.0.0"
    certified = sorted(certified_targets)
    informative = sorted(
        certified if informative_negative_targets is None else informative_negative_targets
    )
    payload["scope"] = control.scope.value if control else None
    payload["target"] = control.target if control else None
    payload["certified_targets"] = certified
    payload["informative_negative_targets"] = informative
    payload["certifies_negatives"] = bool(informative)
    payload["capture_used_for_sensitivity"] = control.capture if certified else None
    if control is not None and control.scope is CaptureScope.PANEL_MECHANICS:
        payload["scope_note"] = (
            "legacy or panel_mechanics control: shows the pipeline recovers a planted "
            "molecule; certifies no virus"
        )
    elif control is not None and control.scope is CaptureScope.EXACT_SEQUENCE:
        payload["scope_note"] = (
            "exact_sequence control: checksum-bound accession in a single-accession "
            "indexed group; no other genome or virus is certified"
        )
    path = os.path.join(results_dir, "positive_control.json")
    from viralscan.validation import require_schema_valid

    require_schema_valid(payload, "positive_control_report.schema.json", path)
    archived_bytes = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    digest = hashlib.sha256(archived_bytes).hexdigest()
    archive = Path(results_dir) / "control_receipts" / f"{digest}.json"
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive_error = None
    try:
        with archive.open("xb") as handle:
            handle.write(archived_bytes)
    except FileExistsError:
        if archive.read_bytes() != archived_bytes:
            archive_error = f"immutable control measurement differs at {archive}"
    except OSError as exc:
        archive_error = f"cannot preserve control measurement at {archive}: {exc}"
    payload["measurement_archive"] = {
        "path": str(archive.relative_to(outputpath)),
        "sha256": digest,
    }
    if archive_error:
        payload["archive_error"] = archive_error
    require_schema_valid(payload, "positive_control_report.schema.json", path)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    if archive_error:
        raise RuntimeError(archive_error)
    log.info(
        "Wrote results/positive_control.json (status=%s, certified targets=%d)",
        payload.get("status"),
        len(certified),
    )
    return path


def check_sibling_crossmapping(virus_stats, sibling_groups=None):
    """Return {virus_name: note_str} for viruses flagged as likely EM bleed.

    ``sibling_groups`` maps virus name -> sibling group (the identity table's
    ``sibling_group``: HHV-6A/6B, HSV-1/2, EBV-1/2). ``None`` uses the retired
    name pairs, for runs with no identity table. Detected viruses of one group
    are compared with the group's dominant member.

    When two viruses that share >80% sequence identity are both detected and
    the UMI ratio exceeds SIBLING_CROSSMAP_RATIO_THRESHOLD, the weaker signal
    is flagged. The global EM allocates a small fraction of the
    dominant sibling's shared-region multimappers to the other, producing a
    residual that is EM noise rather than genuine co-infection. A log warning is
    also emitted for each flagged pair.
    """
    if sibling_groups is None:
        sibling_groups = legacy_sibling_groups()
    members: dict[str, list[str]] = {}
    for virus in virus_stats:
        group = sibling_groups.get(virus)
        if group:
            members.setdefault(group, []).append(virus)
    notes = {}
    for viruses in members.values():
        if len(viruses) < 2:
            continue
        umis = {v: float(virus_stats[v]["viral_molecules_total_est"]) for v in viruses}
        dominant = max(viruses, key=lambda v: umis[v])
        dom_umi = umis[dominant]
        if dom_umi <= 0:
            continue
        for weaker in viruses:
            wk_umi = umis[weaker]
            if weaker == dominant or wk_umi <= 0:
                continue
            ratio = dom_umi / wk_umi
            if ratio < SIBLING_CROSSMAP_RATIO_THRESHOLD:
                continue
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
    from viralscan.run_safety import software_identity

    return {
        "viralscan_version": __version__,
        "software_identity": software_identity(),
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


#: Kallisto count fields an alignment-only row carries as 0 (ANDET-09).
_KALLISTO_COUNT_FIELDS = (
    "viral_molecules_total_est",
    "infected_called",
    "pct_infected_called",
    "infected_comparable",
    "pct_infected_comparable",
    "infected_cells",
    "pct_infected",
    "infected_subthreshold_cells",
    "infected_subthreshold_called",
    "background_umis",
    "viral_molecules_per_10k_est",
)


def _alignment_only_template(summary_rows):
    """Run-level fields for an alignment-only row: denominators kept, counts 0."""
    if not summary_rows:
        return {}
    template: dict[str, Any] = {k: "" for k in summary_rows[0]}
    for k in ("n_called_cells", "n_comparable_cells", "total_cells"):
        template[k] = summary_rows[0][k]
    template.update({k: 0 for k in _KALLISTO_COUNT_FIELDS})
    # Same rule as the kallisto rows: the comparable-cell figures are unavailable when
    # host depth is (blank denominator) and the percentage is undefined at zero cells.
    n_comparable = template["n_comparable_cells"]
    if n_comparable in ("", None):
        template["infected_comparable"] = template["pct_infected_comparable"] = ""
    elif not float(n_comparable):
        template["pct_infected_comparable"] = ""
    template.update(_ambient_not_assessed("alignment-only row: no kallisto counts to screen"))
    return template


def anello_evidence(config, outputpath, identity_table):
    """``(evidence_by_virus, anellovirus_names, status)`` for write_tsv_outputs.

    ``None`` when the identity table is missing (legacy runs keep the old schema).
    """
    if identity_table is None:
        return None
    names = {g.virus_name for g in identity_table.genes if g.viral and g.family == "Anelloviridae"}
    status = anello_align.status_for(
        getattr(config, "anello_align", False),
        getattr(config, "host_index", None),
        getattr(config, "anello_index", None),
    )
    evidence = {}
    path = os.path.join(outputpath, "results", "anello_alignment_by_virus.tsv")
    if status == anello_align.STATUS_OK:
        # The Snakefile makes detection wait for this file, so a missing one is
        # not an empty result: it is a run whose config says the branch was on
        # while the rule never wrote anything -- a resumed pre-ANDET-09 run
        # directory, or a deleted artifact. Say which, rather than reporting
        # zero alignment evidence for every anellovirus.
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"anello_align is enabled but {path} is missing. Rerun the "
                "sample so the anello_align rule writes it, or rerun with "
                "--no-anello-align to report kallisto evidence only."
            )
        with open(path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                evidence[row.pop("virus_name")] = row
    log.info("Anellovirus alignment branch: %s (%d viruses with evidence)", status, len(evidence))
    return evidence, names, status


def write_tsv_outputs(
    virus_stats,
    per_cell_df,
    outputpath,
    crossmap_notes=None,
    facts=None,
    anello=None,
    homology=None,
    homology_loci=None,
):
    """Write viral_summary.tsv and per_cell_viral.tsv to results/ sub-folder.

    *anello* is ``(evidence_by_virus, anellovirus_names, status)`` from
    :func:`anello_evidence`; it adds the ANDET-09 alignment columns and rows.
    *homology* is :func:`host_homology_by_virus` (ANDET-02); a virus it omits is ``not_measured``.
    """
    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    crossmap_notes = crossmap_notes or {}
    facts = facts or {}
    homology = homology or {}

    # Per-virus summary
    summary_rows = []
    for virus, s in virus_stats.items():
        summary_rows.append(
            {
                "virus_name": virus,
                "role": facts[virus].role if virus in facts else "",
                "reference_risk_flags": (
                    "|".join(facts[virus].reference_risk_flags) if virus in facts else ""
                ),
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
                # AUD-06/AUD-07: fractional evidence stays visible, and the called-cell
                # gate is a column, so no candidate row is dropped to apply it.
                "infected_subthreshold_cells": s.get("infected_subthreshold_cells", 0),
                "infected_subthreshold_called": s.get("infected_subthreshold_called", 0),
                "background_umis": s.get("background_umis", 0),
                "detected": s.get("detected", True),
                "viral_molecules_per_10k_est": s["viral_molecules_per_10k_est"],
                "sibling_crossmap_note": crossmap_notes.get(virus, ""),
                # EVE artifact flags
                "accession_breadth": s.get("accession_breadth", 0.0),
                "host_viral_ambig_fraction": s.get("host_viral_ambig_fraction"),
                # From the identity table's risk_class; a virus the table does not
                # describe falls back to the retired genus test (fails closed).
                "eve_risk": facts[virus].eve_risk if virus in facts else legacy_eve_risk(virus),
                # Known low-complexity read-artefact risk (F-019/F-021); a label, not a filter.
                "artifact_risk": (
                    facts[virus].artifact_risk if virus in facts else legacy_artifact_risk(virus)
                ),
                # REF-10: what the call may claim; screening_only for Anelloviridae.
                "claim_scope": (
                    facts[virus].claim_scope if virus in facts else legacy_claim_scope(virus)
                ),
                # ANDET-02: genome-level host homology, measured at build time; no threshold.
                **homology.get(virus, HOST_HOMOLOGY_NOT_MEASURED),
                **{c: s.get(c) for c in AMBIENT_COLUMNS},
            }
        )
    columns = [
        "virus_name",
        "role",
        "reference_risk_flags",
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
        "infected_subthreshold_cells",
        "infected_subthreshold_called",
        "background_umis",
        "detected",
        "viral_molecules_per_10k_est",
        "sibling_crossmap_note",
        "accession_breadth",
        "host_viral_ambig_fraction",
        "eve_risk",
        "artifact_risk",
        "claim_scope",
        "host_homology_status",
        "host_homology_max_identity",
        "host_homology_max_query_coverage",
        "host_homology_max_aligned_bases",
        # M-03: descriptive ambient screen; never read by `detected` (AMBIENT_COLUMNS).
        "ambient_viral_fraction",
        "ambient_viral_fraction_ub95",
        "n_ambient_barcodes_with_virus",
        "soup_expected_cells_max",
        "soup_expected_cells_ge2_max",
        "soup_tail_p_cells",
        "viral_depth_spearman",
        "ambient_assessment",
        "ambient_basis",
    ]
    if anello is not None:
        evidence, anello_names, status = anello
        summary_rows = anello_align.merge_summary_rows(
            summary_rows, evidence, anello_names, status, _alignment_only_template(summary_rows)
        )
        for row in summary_rows:
            if row["detection_source"] == "alignment_only":
                v = str(row["virus_name"])
                row["role"] = facts[v].role if v in facts else ""
                row["reference_risk_flags"] = (
                    "|".join(facts[v].reference_risk_flags) if v in facts else ""
                )
                row["eve_risk"] = facts[v].eve_risk if v in facts else legacy_eve_risk(v)
                row["artifact_risk"] = (
                    facts[v].artifact_risk if v in facts else legacy_artifact_risk(v)
                )
                row["claim_scope"] = facts[v].claim_scope if v in facts else legacy_claim_scope(v)
                row.update(homology.get(v, HOST_HOMOLOGY_NOT_MEASURED))
        columns += list(anello_align.SUMMARY_COLUMNS)
    virus_df = pd.DataFrame(summary_rows, columns=columns)
    virus_df.to_csv(os.path.join(results_dir, "viral_summary.tsv"), sep="\t", index=False)
    locus_rows = []
    for locus in homology_loci or []:
        virus = locus["virus_name"]
        locus_rows.append(
            {
                **locus,
                "eve_risk": facts[virus].eve_risk if virus in facts else legacy_eve_risk(virus),
            }
        )
    pd.DataFrame(locus_rows, columns=HOST_HOMOLOGY_LOCUS_COLUMNS).to_csv(
        os.path.join(results_dir, "host_homology_loci.tsv"), sep="\t", index=False
    )

    # Per-cell viral annotation
    per_cell_df.to_csv(os.path.join(results_dir, "per_cell_viral.tsv"), sep="\t", index=False)
    log.info("Wrote results/viral_summary.tsv and results/per_cell_viral.tsv")


def infected_cell_count(per_cell_df, detected_names=None):
    """Return unique virus-positive barcodes represented in per-cell outputs.

    *detected_names* limits the count to viruses that passed the called-cell gate
    (a candidate-only virus's cells are not "infected cells" in the report headline).
    """
    if per_cell_df is None or per_cell_df.empty or "barcode" not in per_cell_df.columns:
        return 0
    if detected_names is not None:
        per_cell_df = per_cell_df[per_cell_df["virus_name"].isin(list(detected_names))]
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
    facts=None,
    control_line=None,
    detected_names=None,
):
    """Render the Jinja2 HTML report and write it to <outputpath>/report.html.

    ``virus_stats`` lists every candidate; ``detected_names`` (default: the rows whose
    ``detected`` is not False) are those passing the called-cell gate.
    """
    try:
        from jinja2 import Environment, FileSystemLoader, TemplateNotFound
    except ImportError:
        log.warning("jinja2 not installed — skipping HTML report. pip install jinja2 to enable.")
        return

    template_path = os.path.join(os.path.dirname(__file__), "..", "templates")
    env = Environment(loader=FileSystemLoader(template_path), autoescape=True)
    if detected_names is None:
        detected_names = [v for v, s in virus_stats.items() if s.get("detected", True)]

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
        "virus_facts": facts or {},
        "detected_viruses": sorted(detected_names),
        "total_viruses": len(detected_names),
        "se_threshold": config.se_threshold,
        "detection_threshold": config.detection_threshold,
        "embedded_plots": embedded_plots,
        "any_infected": any(
            s["infected_cells"] > 0 for v, s in virus_stats.items() if v in detected_names
        ),
        "per_cell_count": infected_cell_count(per_cell_df, detected_names),
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
        "control_line": control_line or "",
    }

    html = template.render(**ctx)
    report_path = os.path.join(outputpath, "report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    log.info("HTML report written to %s", report_path)


def main():
    adata, found_genes, outputpath, viral_accessions, detection_matrix, index_groups = (
        preprocessing()
    )

    # Clear last run's per-virus plots first: which viruses clear the detection
    # threshold is method-dependent, so a rerun can demote one, and a surviving
    # stale plot would be re-embedded by generate_html_report (SW-04).
    clear_stale_virus_plots(outputpath)

    # check if user wants visuals in output directory
    group_by_virus, detected_viral_genes = histogram(
        adata, found_genes, identity, outputpath, viral_count_matrix=detection_matrix
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
    from viralscan.scripts.cellcalling import (
        CellCallingError,
        HostCellSets,
        call_cells,
        external_host_cells,
        host_called_cells,
        resolve_method,
        solo_raw_dir,
        write_called_cells,
    )

    counts_dir = os.path.join(config.output, "kb-python", "counts_unfiltered")
    # Two-step runs: the kb matrix has no host UMIs, so cells are called on the
    # STARsolo host matrix and the denominators come from that set.
    host_cells: Optional[HostCellSets] = None
    try:
        solo_dir = solo_raw_dir(config)
        if solo_dir is not None and resolve_method(config) == "emptydrops":
            host_cells = host_called_cells(config, solo_dir, COMPARABLE_CELL_MIN_UMI)
            called_mask = _in_barcode_set(adata.obs_names, host_cells["called"])
            if host_cells["called"] and not called_mask.any():
                log.warning(
                    "No barcode of the viral matrix is among the host called cells; every "
                    "candidate will be reported as not detected. Check the barcode suffixes "
                    "of the two matrices."
                )
        else:
            called_mask = call_cells(adata, config, matrix_dir=counts_dir)
            if config.host_index and resolve_method(config) == "external":
                host_cells = external_host_cells(config, COMPARABLE_CELL_MIN_UMI)
        if host_cells:
            os.makedirs(os.path.join(outputpath, "results"), exist_ok=True)
            with open(os.path.join(outputpath, "results", "host_called_cells.tsv"), "w") as fh:
                fh.write("barcode\n")
                fh.writelines(f"{b}\n" for b in sorted(host_cells["called"]))
    except CellCallingError:
        raise
    except Exception as exc:
        raise CellCallingError(
            f"cell calling failed: {exc}. Rerun with --cell-calling none to report "
            "over all barcodes deliberately, or fix the caller inputs."
        ) from exc
    # Layer 2 scores this same set, so it is written out rather than re-derived
    # (PLAN PROG-17).
    write_called_cells(adata.obs_names, called_mask, outputpath)

    # Compute normalized statistics (PR 11 A1/A3) over both denominators
    virus_stats, per_cell_df = compute_stats(
        adata,
        found_genes,
        group_by_virus,
        detected_viral_genes,
        called_mask=called_mask,
        viral_count_matrix=detection_matrix,
        index_genes_by_virus=index_groups,
        host_cells=host_cells,
        infected_min_molecules=config.infected_min_molecules,
        detection_threshold=config.detection_threshold,
        detection_min_cells=config.detection_min_cells,
        cell_calling_method=resolve_method(config),
        ambient_max_molecules=float(
            getattr(config, "emptydrops_lower", DEFAULTS["emptydrops_lower"])
        ),
    )
    # virus_stats holds every candidate of the all-barcode screen (nothing disappears
    # silently); only those passing the called-cell gate are "detected" (AUD-07).
    detected_stats = {v: s for v, s in virus_stats.items() if s["detected"]}
    detected_names = sorted(detected_stats)
    candidate_only = sorted(set(virus_stats) - set(detected_stats))
    if candidate_only:
        log.warning(
            "%d candidate virus(es) fail the called-cell gate (>= %d infected called cell(s) "
            "and >= %d molecules in called cells): %s",
            len(candidate_only),
            config.detection_min_cells,
            config.detection_threshold,
            ", ".join(candidate_only),
        )

    # Optional enrichment by cell type labels (PR 11 A5) — restricted to detected viruses.
    detected_groups = {v: genes for v, genes in group_by_virus.items() if v in detected_stats}
    cell_type_df = cell_type_enrichment(
        adata, detected_groups, config, viral_count_matrix=detection_matrix
    )

    # Ambiguity-aware multimapper evidence is additive and does not alter
    # legacy viral_summary.tsv/per_cell_viral.tsv schemas.
    if should_write_multimap_evidence(config):
        evidence_gene_ids = [g for g in viral_accessions if g in adata.var_names]
        evidence_groups, _ = group_genes(evidence_gene_ids, identity)
        multimap_evidence_df = summarize_multimap_evidence(adata, evidence_groups, config)
    else:
        multimap_evidence_df = summarize_multimap_evidence(None, {}, config)

    # Flag sibling pairs with high UMI asymmetry (HHV-6A/6B, HSV-1/2)
    facts = virus_facts(identity) if identity is not None else {}
    crossmap_notes = check_sibling_crossmapping(
        virus_stats,
        {name: f.sibling_group for name, f in facts.items()} if identity is not None else None,
    )

    # Per-virus detection sensitivity (LOD95). Written for every run, not only
    # negatives: a reader who sees Betatorquevirus 1,142 UMI needs the same
    # yardstick as one who sees nothing, and the LOD columns state which of the
    # two limits (depth or reference capture) is actually binding.
    quantified_depth = float(_sum_axis(adata.X, 1).sum())
    measured_capture, control_detail = measure_positive_control(
        adata, config, count_matrix=detection_matrix, depth=quantified_depth
    )
    # One record list drives the TSV, positive_control.json, summary.txt, the HTML
    # report and the failure message, so no renderer can certify what the others
    # do not (SENS-CORR-03).
    sens = build_sensitivity_outputs(
        adata,
        detected_stats,
        config,
        quantified_depth,
        control_detail,
        index_viruses=index_groups,
        count_matrix=detection_matrix,
    )
    sensitivity_df = sens.table
    control = sens.control
    certified_targets = sens.certified_targets
    # Written before any fail-closed raise below, so a withheld run still leaves
    # its control diagnostics behind.
    write_control_report(
        control_detail,
        outputpath,
        control=control,
        certified_targets=certified_targets,
        informative_negative_targets=sens.informative_targets,
    )

    # Fail closed when the run was told to require a control and none could be
    # measured. This is a mechanics/recovery gate only: passing it does not certify
    # any virus outside the control's scope (see certified_targets). RunConfig already rejects `require_positive_control` with no
    # control, so reaching here with a missing capture means the control was
    # configured and then failed to yield a number (absent gene, or zero
    # recovery). Either way the run cannot certify a negative, and saying so is
    # the whole point of requiring the control.
    nothing_detected = not detected_stats
    if nothing_detected and getattr(config, "require_positive_control", False):
        exact_unverified = (
            getattr(config, "positive_control_scope", None) == "exact_sequence"
            and control_detail.get("verification_status") != "verified"
        )
        if measured_capture is None or exact_unverified:
            raise CellCallingError(
                "require_positive_control is set, this run detected no viral "
                "signal, and no capture term could be measured from the control "
                f"({control_detail.get('status')}; {control_detail.get('verification_detail', '')}). A negative with no "
                "demonstrated ability to see the target is not evidence of "
                f"absence. Control detail: {control_detail.get('detail', control_detail)}"
            )
    if nothing_detected and not certified_targets:
        log.warning(
            "No viral signal and no positive control covering any virus in this run: "
            "every negative is a sampling statement, not an absence. The depth limit "
            "is reported in summary.txt and results/sensitivity.tsv."
        )

    # Write structured TSV outputs (PR 11 A1)
    write_tsv_outputs(
        virus_stats,
        per_cell_df,
        outputpath,
        crossmap_notes=crossmap_notes,
        facts=facts,
        anello=anello_evidence(config, outputpath, identity),
        homology=host_homology_by_virus(identity, config.index) if identity is not None else {},
        homology_loci=host_homology_loci(identity, config.index) if identity is not None else [],
    )
    write_sensitivity_table(sensitivity_df, outputpath)
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
    headline = _headline_totals(
        adata,
        {v: group_by_virus[v] for v in detected_names if v in group_by_virus},
        detection_matrix,
        infected_min_molecules=config.infected_min_molecules,
        called_mask=called_mask,
    )
    n_called_cells = len(host_cells["called"]) if host_cells else int(called_mask.sum())
    found_genes_sorted = dict(sorted(found_genes.items()))
    total_viral_genes = 0
    counts_per_virus = {}
    with (
        open(f"{config.output}/summary.txt", "w") as summary,
        open(Path(config.output) / "log" / "found_genes.txt", "w") as found_genes_file,
    ):
        summary.write(
            "Headline totals cover detected viruses only (called-cell gate met); candidate-only "
            "viruses are listed below.\n"
            f"Viral molecules in unique-count matrix: {headline['unique']}\n"
            f"Total viral molecules (selected method): {headline['selected']}\n"
            f"Cells with viral reads: {headline['cells_with_virus_called']}/{n_called_cells} "
            f"called cells ({headline['cells_with_virus']}/{headline['n_cells']} all barcodes, "
            "empty droplets included), cells holding at least "
            f"{config.infected_min_molecules:g} molecules of a detected virus\n\n\n"
        )
        if len(found_genes_sorted) > 0:
            summary.write(
                "Candidate viral Gene IDs including the count "
                "(all barcodes, empty droplets included):\n"
            )
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
                summary.write("\n" + _virus_summary_line(virus_name, stats))
            summary.write(f"\n\nTotal amount of viral load found: {total_viral_genes}")
            summary.write(f"\n\nOfficial name of viral load detected: {','.join(detected_names)}")
            summary.write("\n\n" + _gate_note(config, resolve_method(config), candidate_only))
            if not detected_names:
                # Candidates only: the negative is a sampling statement, as when nothing is found.
                summary.write("\n" + sens.interpretation + "\n" + sens.control_line)
        else:
            summary.write("No viral gene IDs found in this sample for the viruses in the index.\n")
            # A zero is only interpretable next to the limit that produced it.
            # Without this line "no virus found" reads as an absence, which is
            # exactly the inference the depth and capture terms do not support.
            # Same text as report.html: both come from the records.
            summary.write(sens.interpretation + "\n")
            summary.write(sens.control_line + "\n")
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
        facts=facts,
        sensitivity_statement=sens.interpretation,
        control_line=sens.control_line,
        detected_names=detected_names,
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
