"""Ambiguity-aware multimapper redistribution and evidence summaries."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from viralscan.defaults import (
    DEFAULT_MULTIMAP_METHOD,
    MULTIMAP_METHODS,
)
from viralscan.runconfig import RunConfig

MULTIMAP_EVIDENCE_COLUMNS = [
    "virus_name",
    "gene_id",
    "unique_viral_umi",
    "ambiguous_viral_umi",
    "host_viral_ambiguous_umi",
    "corrected_viral_umi",
    "upper_bound_viral_umi",
    "n_unique_viral_cells",
    "n_ambiguous_viral_cells",
    "multimap_method",
    "call_confidence",
]


@dataclass
class MultimapLayers:
    corrected: sparse.csr_matrix
    equal: sparse.csr_matrix
    host_conservative: sparse.csr_matrix
    unique_weighted: sparse.csr_matrix
    unique_viral: sparse.csr_matrix
    host_viral_ambiguous: sparse.csr_matrix
    host_viral_selected: sparse.csr_matrix
    viral_ambiguous_upper: sparse.csr_matrix


def _empty_matrix(n_cells: int, n_genes: int) -> sparse.csr_matrix:
    return sparse.csr_matrix((n_cells, n_genes), dtype=float)


def _matrix_value(matrix: Any, row: int, col: int) -> float:
    value = matrix[row, col]
    if sparse.issparse(value):
        return float(value.toarray()[0, 0])
    if hasattr(value, "item"):
        return float(value.item())
    return float(value)


def _sum_gene(matrix: Any, gene_idx: int) -> float:
    values = matrix[:, gene_idx]
    return float(values.sum())


def _gene_vector(matrix: Any, gene_idx: int) -> Any:
    values = matrix[:, gene_idx]
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).reshape(-1)


def _csr_from_entries(
    rows: list[int], cols: list[int], data: list[float], n_cells: int, n_genes: int
) -> sparse.csr_matrix:
    if not data:
        return _empty_matrix(n_cells, n_genes)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n_cells, n_genes), dtype=float)


def em_gene_abundances(
    ec_counts: dict[tuple[int, ...], float],
    unique_per_gene: np.ndarray,
    pseudocount: float,
    max_iter: int,
    tol: float,
) -> np.ndarray:
    """Estimate global gene abundances by EM over multi-gene equivalence classes.

    This is the standard RSEM/kallisto-style multimapper model: unique-mapping
    reads are fixed assignments, and each multi-gene EC's mass is fractionally
    allocated to its genes in proportion to the current abundance estimate. One
    pass = one E-step (allocate by current theta) + M-step (theta = unique +
    allocated). ``unique-weighted`` is exactly the *first* E-step of this loop;
    EM iterates it to a fixed point.

    **Global-pool design:** ``ec_counts`` aggregates multi-gene EC masses across
    *all cells* before EM runs. The returned ``theta`` is therefore a
    transcriptome-wide gene-abundance estimate, not a per-cell one. Caller
    (``build_multimap_layers``) uses this single ``theta`` to allocate
    multi-mapping UMIs in every cell. This differs from per-cell EM
    (alevin-fry, STARsolo) and is significantly faster for large datasets at the
    cost of ignoring cell-to-cell abundance variation when resolving ambiguity.

    Parameters
    ----------
    ec_counts:
        Map from a distinct gene-index tuple (the EC's genes) to the total
        multimapper UMI count pooled across all cells.
    unique_per_gene:
        Global unique-mapping UMI per gene (length n_genes) — the fixed mass.
    pseudocount:
        Added to the initial abundances so genes with zero unique support still
        receive non-zero weight on the first E-step.
    max_iter, tol:
        Stop after ``max_iter`` sweeps or when the L1 change in theta (relative
        to its total) drops below ``tol``.

    Returns
    -------
    theta: np.ndarray
        Converged non-negative abundance estimate per gene.
    """
    unique_per_gene = np.asarray(unique_per_gene, dtype=float).reshape(-1)
    n_genes = unique_per_gene.shape[0]
    theta = unique_per_gene + float(pseudocount)
    if not ec_counts:
        return theta

    # Vectorised E/M step: represent the pooled multi-gene ECs as a sparse
    # (n_ec x n_gene) incidence matrix M and a per-EC count vector. Each E-step
    # is then two sparse mat-vecs instead of a Python loop over ECs. This is
    # numerically equivalent to the per-EC allocation `new[genes] += count*w/s`
    # (theta-proportional for well-supported ECs, equal-split for degenerate
    # ones), up to floating-point summation order.
    ec_rows: list[np.ndarray] = []
    ec_cols: list[np.ndarray] = []
    counts = np.empty(len(ec_counts), dtype=float)
    genes_per_ec = np.empty(len(ec_counts), dtype=float)
    for e, (genes, count) in enumerate(ec_counts.items()):
        g = np.asarray(genes, dtype=int)
        ec_rows.append(np.full(g.shape[0], e, dtype=int))
        ec_cols.append(g)
        counts[e] = count
        genes_per_ec[e] = g.shape[0]
    row_idx = np.concatenate(ec_rows)
    col_idx = np.concatenate(ec_cols)
    incidence = sparse.csr_matrix(
        (np.ones(row_idx.shape[0], dtype=float), (row_idx, col_idx)),
        shape=(len(ec_counts), n_genes),
    )
    incidence_t = incidence.T.tocsr()

    for _ in range(int(max_iter)):
        s = incidence @ theta  # per-EC denominator = sum of theta over the EC's genes
        good = s > 1e-12
        # theta-weighted allocation for well-supported ECs; equal split otherwise.
        weighted = theta * (incidence_t @ np.where(good, counts / np.where(good, s, 1.0), 0.0))
        equal = incidence_t @ np.where(good, 0.0, counts / genes_per_ec)
        new = unique_per_gene + weighted + equal
        denom = float(theta.sum()) or 1.0
        if float(np.abs(new - theta).sum()) / denom < tol:
            theta = new
            break
        theta = new
    return theta


def build_multimap_layers(
    bus_df: pd.DataFrame,
    barcode_to_idx: dict[str, int],
    ec_map: dict[int, list[int]],
    n_cells: int,
    n_genes: int,
    viral_gene_indices: set[int],
    original_counts: Any,
    method: str = DEFAULT_MULTIMAP_METHOD,
    pseudocount: float = 1.0,
    em_max_iter: int = 100,
    em_tol: float = 1e-6,
) -> MultimapLayers:
    """Build selected and diagnostic multimapper correction layers.

    Unique ECs are never added to ``corrected`` because they are already present
    in the original kb count matrix. Viral unique ECs are tracked separately for
    confidence reporting.

    ``method="em"`` resolves multimappers by an iterated EM over equivalence
    classes (see :func:`em_gene_abundances`): unique reads are fixed and each
    multi-gene EC's mass is allocated to its genes in proportion to the converged
    global abundance estimate. The three deterministic layers are still built so
    the evidence table and diagnostics remain available under any method.
    """
    if method not in MULTIMAP_METHODS:
        raise ValueError(f"Unknown multimap method: {method}")
    if pseudocount <= 0:
        raise ValueError(f"multimap_pseudocount must be > 0, got {pseudocount}.")

    # EM-only accumulators (collected during the single pass, resolved afterwards)
    use_em = method == "em"
    em_records: list[tuple[int, tuple[int, ...], float]] = []
    em_ec_counts: dict[tuple[int, ...], float] = {}

    equal_rows: list[int] = []
    equal_cols: list[int] = []
    equal_data: list[float] = []
    conservative_rows: list[int] = []
    conservative_cols: list[int] = []
    conservative_data: list[float] = []
    weighted_rows: list[int] = []
    weighted_cols: list[int] = []
    weighted_data: list[float] = []
    unique_rows: list[int] = []
    unique_cols: list[int] = []
    unique_data: list[float] = []
    host_viral_rows: list[int] = []
    host_viral_cols: list[int] = []
    host_viral_data: list[float] = []
    selected_host_viral_rows: list[int] = []
    selected_host_viral_cols: list[int] = []
    selected_host_viral_data: list[float] = []
    upper_rows: list[int] = []
    upper_cols: list[int] = []
    upper_data: list[float] = []

    for row in bus_df.itertuples(index=False):
        bc, ec, count = row.barcode, row.ec, float(row.count)
        if pd.isna(ec):
            continue
        ec = int(ec)
        if bc not in barcode_to_idx or ec not in ec_map:
            continue

        cell_idx = barcode_to_idx[bc]
        genes_in_ec = list(ec_map[ec])
        if not genes_in_ec:
            continue

        distinct_genes = list(dict.fromkeys(genes_in_ec))
        viral_genes = [gid for gid in distinct_genes if gid in viral_gene_indices]
        host_genes = [gid for gid in distinct_genes if gid not in viral_gene_indices]

        if len(genes_in_ec) == 1:
            gid = genes_in_ec[0]
            if gid in viral_gene_indices:
                unique_rows.append(cell_idx)
                unique_cols.append(gid)
                unique_data.append(count)
            continue

        if use_em:
            key = tuple(distinct_genes)
            em_records.append((cell_idx, key, count))
            em_ec_counts[key] = em_ec_counts.get(key, 0.0) + count

        equal_share = count / len(genes_in_ec)
        weights = np.array(
            [_matrix_value(original_counts, cell_idx, gid) + pseudocount for gid in genes_in_ec],
            dtype=float,
        )
        weight_sum = float(weights.sum())

        for i, gid in enumerate(genes_in_ec):
            weighted_share = count * float(weights[i]) / weight_sum
            equal_rows.append(cell_idx)
            equal_cols.append(gid)
            equal_data.append(equal_share)

            if not (viral_genes and host_genes and gid in viral_gene_indices):
                conservative_rows.append(cell_idx)
                conservative_cols.append(gid)
                conservative_data.append(equal_share)

            weighted_rows.append(cell_idx)
            weighted_cols.append(gid)
            weighted_data.append(weighted_share)

            if viral_genes and host_genes and gid in viral_gene_indices:
                selected_host_viral_rows.append(cell_idx)
                selected_host_viral_cols.append(gid)
                if method == "unique-weighted":
                    selected_host_viral_data.append(weighted_share)
                elif method == "host-conservative":
                    selected_host_viral_data.append(0.0)
                else:
                    selected_host_viral_data.append(equal_share)

        if viral_genes and host_genes:
            host_viral_share = count / len(viral_genes)
            for gid in viral_genes:
                host_viral_rows.append(cell_idx)
                host_viral_cols.append(gid)
                host_viral_data.append(host_viral_share)

        if viral_genes:
            upper_share = count / len(viral_genes)
            for gid in viral_genes:
                upper_rows.append(cell_idx)
                upper_cols.append(gid)
                upper_data.append(upper_share)

    equal = _csr_from_entries(equal_rows, equal_cols, equal_data, n_cells, n_genes)
    host_conservative = _csr_from_entries(
        conservative_rows, conservative_cols, conservative_data, n_cells, n_genes
    )
    unique_weighted = _csr_from_entries(
        weighted_rows, weighted_cols, weighted_data, n_cells, n_genes
    )
    host_viral_selected = _csr_from_entries(
        selected_host_viral_rows,
        selected_host_viral_cols,
        selected_host_viral_data,
        n_cells,
        n_genes,
    )

    if use_em:
        # Resolve multimappers by iterated EM over the pooled ECs, then allocate
        # each cell's records by the converged abundances. Also recompute the
        # host-virus-selected diagnostic (viral mass credited from host+viral ECs)
        # from the same allocation so confidence tiers stay consistent.
        unique_per_gene = np.asarray(
            original_counts.sum(axis=0)
            if sparse.issparse(original_counts)
            else original_counts.sum(axis=0)
        ).reshape(-1)
        theta = em_gene_abundances(em_ec_counts, unique_per_gene, pseudocount, em_max_iter, em_tol)
        corr_rows: list[int] = []
        corr_cols: list[int] = []
        corr_data: list[float] = []
        sel_rows: list[int] = []
        sel_cols: list[int] = []
        sel_data: list[float] = []
        for cell_idx, genes, count in em_records:
            gidx = np.asarray(genes, dtype=int)
            w = theta[gidx]
            s = float(w.sum())
            shares = (count * w / s) if s > 0.0 else np.full(len(gidx), count / len(gidx))
            has_viral = any(int(g) in viral_gene_indices for g in gidx)
            has_host = any(int(g) not in viral_gene_indices for g in gidx)
            for g, share in zip(gidx, shares):
                gi = int(g)
                corr_rows.append(cell_idx)
                corr_cols.append(gi)
                corr_data.append(float(share))
                if has_viral and has_host and gi in viral_gene_indices:
                    sel_rows.append(cell_idx)
                    sel_cols.append(gi)
                    sel_data.append(float(share))
        selected = _csr_from_entries(corr_rows, corr_cols, corr_data, n_cells, n_genes)
        host_viral_selected = _csr_from_entries(sel_rows, sel_cols, sel_data, n_cells, n_genes)
    else:
        selected = {
            "equal": equal,
            "host-conservative": host_conservative,
            "unique-weighted": unique_weighted,
        }[method]

    return MultimapLayers(
        corrected=selected,
        equal=equal,
        host_conservative=host_conservative,
        unique_weighted=unique_weighted,
        unique_viral=_csr_from_entries(unique_rows, unique_cols, unique_data, n_cells, n_genes),
        host_viral_ambiguous=_csr_from_entries(
            host_viral_rows, host_viral_cols, host_viral_data, n_cells, n_genes
        ),
        host_viral_selected=host_viral_selected,
        viral_ambiguous_upper=_csr_from_entries(
            upper_rows, upper_cols, upper_data, n_cells, n_genes
        ),
    )


def _confidence(
    unique_umi: float,
    corrected_umi: float,
    host_viral_ambiguous_umi: float,
    selected_host_viral_umi: float,
    threshold: float,
) -> str:
    if unique_umi >= threshold:
        return "strong"
    non_host_ambiguous_umi = max(corrected_umi - selected_host_viral_umi, 0.0)
    if unique_umi + non_host_ambiguous_umi >= threshold:
        return "ambiguous"
    if unique_umi + corrected_umi >= threshold or host_viral_ambiguous_umi >= threshold:
        return "low_confidence"
    return "not_detected"


def summarize_multimap_evidence(
    adata: Any,
    group_by_virus: dict[str, list[str]],
    config: RunConfig,
) -> pd.DataFrame:
    """Summarize unique and ambiguous viral evidence for each grouped gene."""
    if adata is None or not group_by_virus:
        return pd.DataFrame(columns=MULTIMAP_EVIDENCE_COLUMNS)

    n_cells, n_genes = adata.n_obs, adata.n_vars
    zero = _empty_matrix(n_cells, n_genes)
    unique = adata.layers.get("counts_unique_viral", zero)
    corrected = adata.layers.get("counts_corrected", zero)
    host_viral = adata.layers.get("counts_host_viral_ambiguous", zero)
    host_viral_selected = adata.layers.get("counts_host_viral_selected", zero)
    upper = adata.layers.get("counts_viral_ambiguous_upper", corrected)
    method = str(config.multimap_method)
    threshold = float(config.detection_threshold)

    rows = []
    for virus, gene_ids in group_by_virus.items():
        for gene_id in gene_ids:
            if gene_id not in adata.var_names:
                continue
            gene_idx = int(adata.var_names.get_loc(gene_id))
            unique_umi = _sum_gene(unique, gene_idx)
            ambiguous_umi = _sum_gene(corrected, gene_idx)
            host_viral_umi = _sum_gene(host_viral, gene_idx)
            selected_host_viral_umi = _sum_gene(host_viral_selected, gene_idx)
            corrected_umi = unique_umi + ambiguous_umi
            upper_bound_umi = unique_umi + _sum_gene(upper, gene_idx)
            unique_vec = _gene_vector(unique, gene_idx)
            ambiguous_vec = _gene_vector(upper, gene_idx)
            rows.append(
                {
                    "virus_name": virus,
                    "gene_id": gene_id,
                    "unique_viral_umi": round(unique_umi, 6),
                    "ambiguous_viral_umi": round(ambiguous_umi, 6),
                    "host_viral_ambiguous_umi": round(host_viral_umi, 6),
                    "corrected_viral_umi": round(corrected_umi, 6),
                    "upper_bound_viral_umi": round(upper_bound_umi, 6),
                    "n_unique_viral_cells": int((unique_vec > 0).sum()),
                    "n_ambiguous_viral_cells": int((ambiguous_vec > 0).sum()),
                    "multimap_method": method,
                    "call_confidence": _confidence(
                        unique_umi,
                        ambiguous_umi,
                        host_viral_umi,
                        selected_host_viral_umi,
                        threshold,
                    ),
                }
            )

    return pd.DataFrame(rows, columns=MULTIMAP_EVIDENCE_COLUMNS)


def select_detection_matrix(adata: Any, config: RunConfig) -> Any:
    """Return the count matrix used for primary viral calls."""
    if (
        config.multimapping
        and config.multimap_primary_call == "unique-only"
        and "counts_unique_viral" in adata.layers
    ):
        return adata.layers["counts_unique_viral"]
    return adata.X


def should_write_multimap_evidence(config: RunConfig) -> bool:
    """Return whether multimapper evidence outputs should be produced."""
    return bool(config.multimapping)


def write_multimap_evidence(evidence_df: pd.DataFrame, outputpath: str) -> str:
    """Write results/multimap_evidence.tsv and return its path."""
    import os

    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, "multimap_evidence.tsv")
    evidence_df.to_csv(out_path, sep="\t", index=False)
    return out_path
