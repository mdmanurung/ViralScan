"""Ambiguity-aware multimapper redistribution and evidence summaries."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
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
    "viral_molecules_unique",
    "viral_molecules_ambiguous_allocated",
    "host_virus_ambiguous_molecules",
    "viral_molecules_total_est",
    "viral_molecules_upper_bound",
    "n_unique_viral_cells",
    "n_ambiguous_viral_cells",
    "multimap_method",
    "evidence_tier",
]


@dataclass
class MultimapLayers:
    corrected: sparse.csr_matrix
    unique: sparse.csr_matrix
    equal: sparse.csr_matrix
    host_conservative: sparse.csr_matrix
    unique_weighted: sparse.csr_matrix
    unique_viral: sparse.csr_matrix
    host_viral_ambiguous: sparse.csr_matrix
    host_viral_selected: sparse.csr_matrix
    viral_ambiguous_upper: sparse.csr_matrix
    audit: MoleculeAudit
    method_diagnostics: dict[str, Any]


@dataclass(frozen=True)
class MoleculeAudit:
    """Conservation accounting for corrected CB-UMI molecules."""

    input_molecules: int
    resolved_molecules: int
    unique_molecules: int
    ambiguous_molecules: int
    unresolved_molecules: int
    ignored_read_multiplicity: int

    def validate(self, allocated_mass: float) -> None:
        if self.resolved_molecules != self.unique_molecules + self.ambiguous_molecules:
            raise ValueError("Resolved molecule audit does not equal unique plus ambiguous.")
        if (
            self.unique_molecules + self.ambiguous_molecules + self.unresolved_molecules
            != self.input_molecules
        ):
            raise ValueError("Molecule audit does not conserve input molecules.")
        if not np.isfinite(allocated_mass) or allocated_mass < 0:
            raise ValueError("Allocated molecule mass must be finite and non-negative.")
        if not np.isclose(allocated_mass, self.ambiguous_molecules, rtol=0.0, atol=1e-9):
            raise ValueError(
                "Allocated ambiguous mass does not equal the number of ambiguous molecules."
            )


@dataclass
class _AuditCounter:
    input_molecules: int = 0
    resolved_molecules: int = 0
    unique_molecules: int = 0
    ambiguous_molecules: int = 0
    unresolved_molecules: int = 0
    ignored_read_multiplicity: int = 0

    def freeze(self) -> MoleculeAudit:
        return MoleculeAudit(**vars(self))


def _resolved_gene_tuple(ec_ids: Iterable[int], ec_map: dict[int, list[int]]) -> tuple[int, ...]:
    gene_sets = [set(ec_map[ec]) for ec in dict.fromkeys(ec_ids) if ec in ec_map and ec_map[ec]]
    if not gene_sets:
        return ()
    return tuple(sorted(set.intersection(*gene_sets)))


def _iter_bus_text_records(path: Path, buffering: int = -1) -> Iterator[tuple[str, str, int, int]]:
    with path.open(encoding="utf-8", buffering=buffering) as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Malformed BUS text row {line_number}: expected four columns.")
            barcode, umi, ec_raw, count_raw = fields
            yield barcode.removesuffix("-1"), umi, int(ec_raw), int(count_raw)


def iter_cb_umi_molecules(
    records: Iterable[tuple[str, str, int, int]],
    barcode_to_idx: dict[str, int],
    ec_map: dict[int, list[int]],
    audit: _AuditCounter,
) -> Iterator[tuple[int, tuple[int, ...]]]:
    """Stream a corrected, CB-UMI-sorted BUS record iterable."""
    current_key: tuple[str, str] | None = None
    current_cell: int | None = None
    ec_ids: list[int] = []

    def resolve_current() -> tuple[int, tuple[int, ...]] | None:
        if current_key is None or current_cell is None:
            return None
        audit.input_molecules += 1
        genes = _resolved_gene_tuple(ec_ids, ec_map)
        if not genes:
            audit.unresolved_molecules += 1
            return None
        audit.resolved_molecules += 1
        if len(genes) == 1:
            audit.unique_molecules += 1
        else:
            audit.ambiguous_molecules += 1
        return current_cell, genes

    previous_key: tuple[str, str] | None = None
    for barcode, umi, ec, count in records:
        if count < 1:
            raise ValueError("BUS record multiplicity must be a positive integer.")
        cell = barcode_to_idx.get(barcode)
        if cell is None:
            continue
        key = (barcode, umi)
        if previous_key is not None and key < previous_key:
            raise ValueError("BUS text must be corrected and sorted by barcode and UMI.")
        previous_key = key
        audit.ignored_read_multiplicity += count - 1
        if key != current_key:
            resolved = resolve_current()
            if resolved is not None:
                yield resolved
            current_key = key
            current_cell = cell
            ec_ids = [ec]
        else:
            ec_ids.append(ec)
    resolved = resolve_current()
    if resolved is not None:
        yield resolved


def resolve_cb_umi_molecules(
    bus_df: pd.DataFrame,
    barcode_to_idx: dict[str, int],
    ec_map: dict[int, list[int]],
) -> tuple[list[tuple[int, tuple[int, ...]]], MoleculeAudit]:
    """Resolve corrected BUS records into at most one molecule per CB-UMI.

    The input must be the corrected, sorted BUS stream. Record ``count`` is read
    multiplicity and is deliberately excluded from molecule mass. When a
    corrected CB-UMI has multiple EC observations, compatible genes are
    intersected. Disjoint observations are retained in the audit as unresolved
    UMI collisions and do not receive a forced assignment.
    """
    required = {"barcode", "umi", "ec", "count"}
    missing = required.difference(bus_df.columns)
    if missing:
        raise ValueError(f"BUS molecule resolution requires columns: {sorted(missing)}")

    work = bus_df.loc[:, ["barcode", "umi", "ec", "count"]].copy()
    work["barcode"] = work["barcode"].astype(str).str.removesuffix("-1")
    work["cell"] = work["barcode"].map(barcode_to_idx)
    work = work.dropna(subset=["cell", "umi", "ec"])
    if (work["count"] < 1).any():
        raise ValueError("BUS record multiplicity must be a positive integer.")
    work = work.astype({"cell": "int64", "ec": "int64"})

    work = work.sort_values(["barcode", "umi", "ec"], kind="stable")
    counter = _AuditCounter()
    records = (
        (str(row.barcode), str(row.umi), int(row.ec), int(row.count))
        for row in work.itertuples(index=False)
    )
    resolved = list(iter_cb_umi_molecules(records, barcode_to_idx, ec_map, counter))
    return resolved, counter.freeze()


def _empty_matrix(n_cells: int, n_genes: int) -> sparse.csr_matrix:
    return sparse.csr_matrix((n_cells, n_genes), dtype=float)


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
    diagnostics: dict[str, Any] | None = None,
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
        if diagnostics is not None:
            diagnostics.update(
                iterations=0, converged=True, fallback_reason="no_ambiguous_molecules"
            )
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

    converged = False
    iterations = 0
    for iteration in range(1, int(max_iter) + 1):
        iterations = iteration
        s = incidence @ theta  # per-EC denominator = sum of theta over the EC's genes
        good = s > 1e-12
        # theta-weighted allocation for well-supported ECs; equal split otherwise.
        weighted = theta * (incidence_t @ np.where(good, counts / np.where(good, s, 1.0), 0.0))
        equal = incidence_t @ np.where(good, 0.0, counts / genes_per_ec)
        new = unique_per_gene + weighted + equal
        denom = float(theta.sum()) or 1.0
        if float(np.abs(new - theta).sum()) / denom < tol:
            theta = new
            converged = True
            break
        theta = new
    if diagnostics is not None:
        diagnostics.update(
            iterations=iterations,
            converged=converged,
            fallback_reason=None if converged else "max_iterations_reached",
        )
    return theta


def em_cell_abundances(
    ec_counts: dict[tuple[int, ...], float],
    unique_per_gene: np.ndarray,
    global_theta: np.ndarray,
    prior_strength: float,
    max_iter: int,
    tol: float,
    diagnostics: dict[str, Any] | None = None,
) -> np.ndarray:
    """Fit one cell-local abundance model shrunk toward the sample model."""
    unique = np.asarray(unique_per_gene, dtype=float).reshape(-1)
    global_theta = np.asarray(global_theta, dtype=float).reshape(-1)
    compatible = sorted({gene for genes in ec_counts for gene in genes})
    if compatible and float(unique[np.asarray(compatible, dtype=int)].sum()) <= 0.0:
        if diagnostics is not None:
            diagnostics.update(
                iterations=0,
                converged=True,
                fallback_reason="no_local_compatible_unique_support",
            )
        # A positive sample prior alone is not cell-local evidence. Returning
        # zero weights makes the caller retain these molecules as explicit
        # equal-split ambiguity while preserving one unit of molecule mass.
        return np.zeros_like(unique)
    prior = prior_strength * global_theta / (float(global_theta.sum()) or 1.0)
    theta = unique + prior
    for iteration in range(1, int(max_iter) + 1):
        new = unique + prior
        for genes, count in ec_counts.items():
            gidx = np.asarray(genes, dtype=int)
            weights = theta[gidx]
            denom = float(weights.sum())
            if denom > 0:
                new[gidx] += count * weights / denom
            else:
                new[gidx] += count / len(gidx)
        denom = float(theta.sum()) or 1.0
        if float(np.abs(new - theta).sum()) / denom < tol:
            if diagnostics is not None:
                diagnostics.update(iterations=iteration, converged=True, fallback_reason=None)
            return new
        theta = new
    if diagnostics is not None:
        diagnostics.update(
            iterations=int(max_iter), converged=False, fallback_reason="max_iterations_reached"
        )
    return theta


def build_multimap_layers(
    bus_df: pd.DataFrame | Path,
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
    bus_buffer_size: int = -1,
) -> MultimapLayers:
    """Build selected and diagnostic multimapper correction layers.

    Unique ECs are never added to ``corrected`` because they are already present
    in the original kb count matrix. Viral unique ECs are tracked separately for
    confidence reporting.

    ``method="em-global"`` resolves multimappers by an iterated EM over equivalence
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
    use_em = method in {"em-global", "em-cell"}
    ambiguous_records: list[tuple[int, tuple[int, ...], float]] = []
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
    unique_viral_rows: list[int] = []
    unique_viral_cols: list[int] = []
    unique_viral_data: list[float] = []
    host_viral_rows: list[int] = []
    host_viral_cols: list[int] = []
    host_viral_data: list[float] = []
    selected_host_viral_rows: list[int] = []
    selected_host_viral_cols: list[int] = []
    selected_host_viral_data: list[float] = []
    upper_rows: list[int] = []
    upper_cols: list[int] = []
    upper_data: list[float] = []

    # Per-EC invariants are hoisted out of the per-molecule loop: there are far fewer
    # distinct ECs than BUS records (388k vs ~100M on deep samples), and the gene
    # classification / conservative & selected masks depend only on the EC, not the
    # cell or molecule.
    ec_info: dict[int, Any] = {}
    gene_tuple_info: dict[tuple[int, ...], Any] = {}
    for ec_id, genes in ec_map.items():
        genes_in_ec = sorted(set(genes))
        if not genes_in_ec:
            ec_info[ec_id] = None
            continue
        distinct_genes = genes_in_ec
        viral_genes = [gid for gid in distinct_genes if gid in viral_gene_indices]
        has_both = bool(viral_genes) and len(viral_genes) != len(distinct_genes)
        is_viral_pos = [gid in viral_gene_indices for gid in genes_in_ec]
        # conservative keeps a position unless it is a viral gene in a host+viral EC.
        cons_eligible = [not (has_both and v) for v in is_viral_pos]
        sel_eligible = [has_both and v for v in is_viral_pos]
        ec_info[ec_id] = (
            genes_in_ec,
            len(genes_in_ec),
            tuple(distinct_genes),
            viral_genes,
            has_both,
            is_viral_pos,
            cons_eligible,
            sel_eligible,
            np.asarray(genes_in_ec, dtype=np.intp),
        )
        gene_tuple_info[tuple(distinct_genes)] = ec_info[ec_id]

    audit_counter = _AuditCounter()
    if isinstance(bus_df, Path):
        molecule_iter = iter_cb_umi_molecules(
            _iter_bus_text_records(bus_df, bus_buffer_size),
            barcode_to_idx,
            ec_map,
            audit_counter,
        )
    else:
        required = {"barcode", "umi", "ec", "count"}
        missing = required.difference(bus_df.columns)
        if missing:
            raise ValueError(f"BUS molecule resolution requires columns: {sorted(missing)}")
        ordered = bus_df.sort_values(["barcode", "umi", "ec"], kind="stable")
        records = (
            (str(row.barcode).removesuffix("-1"), str(row.umi), int(row.ec), int(row.count))
            for row in ordered.itertuples(index=False)
        )
        molecule_iter = iter_cb_umi_molecules(records, barcode_to_idx, ec_map, audit_counter)

    for cell_idx, distinct_key in molecule_iter:
        # Molecules are already projected to distinct compatible genes. Reuse the
        # cached EC classification where possible, otherwise classify the tuple.
        info = gene_tuple_info.get(distinct_key)
        if info is None:
            genes_in_ec = list(distinct_key)
            viral_genes = [gid for gid in genes_in_ec if gid in viral_gene_indices]
            has_both = bool(viral_genes) and len(viral_genes) != len(genes_in_ec)
            is_viral_pos = [gid in viral_gene_indices for gid in genes_in_ec]
            info = (
                genes_in_ec,
                len(genes_in_ec),
                distinct_key,
                viral_genes,
                has_both,
                is_viral_pos,
                [not (has_both and v) for v in is_viral_pos],
                [has_both and v for v in is_viral_pos],
                np.asarray(genes_in_ec, dtype=np.intp),
            )
            gene_tuple_info[distinct_key] = info
        (
            genes_in_ec,
            n_genes_in_ec,
            distinct_key,
            viral_genes,
            has_both,
            is_viral_pos,
            cons_eligible,
            sel_eligible,
            genes_arr,
        ) = info
        count = 1.0

        if n_genes_in_ec == 1:
            unique_rows.append(cell_idx)
            unique_cols.append(genes_in_ec[0])
            unique_data.append(count)
            if is_viral_pos[0]:
                unique_viral_rows.append(cell_idx)
                unique_viral_cols.append(genes_in_ec[0])
                unique_viral_data.append(count)
            continue

        ambiguous_records.append((cell_idx, distinct_key, count))
        if use_em:
            em_ec_counts[distinct_key] = em_ec_counts.get(distinct_key, 0.0) + count

        equal_share = count / n_genes_in_ec

        for i, gid in enumerate(genes_in_ec):
            equal_rows.append(cell_idx)
            equal_cols.append(gid)
            equal_data.append(equal_share)

            if cons_eligible[i]:
                conservative_rows.append(cell_idx)
                conservative_cols.append(gid)
                conservative_data.append(count / sum(cons_eligible))

            if sel_eligible[i]:
                if method != "unique-weighted":
                    selected_host_viral_rows.append(cell_idx)
                    selected_host_viral_cols.append(gid)
                if method == "host-conservative":
                    selected_host_viral_data.append(0.0)
                elif method != "unique-weighted":
                    selected_host_viral_data.append(equal_share)

        if has_both:
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

    unique = _csr_from_entries(unique_rows, unique_cols, unique_data, n_cells, n_genes)
    unique.sort_indices()
    unique_indptr, unique_indices, unique_values = unique.indptr, unique.indices, unique.data
    for cell_idx, distinct_key, count in ambiguous_records:
        genes_arr = np.asarray(distinct_key, dtype=np.intp)
        start = unique_indptr[cell_idx]
        row_cols = unique_indices[start : unique_indptr[cell_idx + 1]]
        if row_cols.shape[0] == 0:
            weights = np.full(len(genes_arr), pseudocount, dtype=float)
        else:
            pos = np.searchsorted(row_cols, genes_arr)
            safe = np.minimum(pos, row_cols.shape[0] - 1)
            hit = row_cols[safe] == genes_arr
            weights = np.where(hit, unique_values[start + safe], 0.0) + pseudocount
        shares = count * weights / float(weights.sum())
        has_viral = any(int(gene) in viral_gene_indices for gene in genes_arr)
        has_host = any(int(gene) not in viral_gene_indices for gene in genes_arr)
        for gene, share in zip(genes_arr, shares):
            gid = int(gene)
            weighted_rows.append(cell_idx)
            weighted_cols.append(gid)
            weighted_data.append(float(share))
            if method == "unique-weighted" and has_viral and has_host and gid in viral_gene_indices:
                selected_host_viral_rows.append(cell_idx)
                selected_host_viral_cols.append(gid)
                selected_host_viral_data.append(float(share))

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

    method_diagnostics: dict[str, Any] = {
        "method": method,
        "model_scope": "cell" if method == "em-cell" else ("sample" if use_em else "deterministic"),
        "pseudocount": pseudocount,
        "max_iter": em_max_iter if use_em else None,
        "tolerance": em_tol if use_em else None,
    }
    if use_em:
        # Resolve multimappers by iterated EM over the pooled ECs, then allocate
        # each cell's records by the converged abundances. Also recompute the
        # host-virus-selected diagnostic (viral mass credited from host+viral ECs)
        # from the same allocation so confidence tiers stay consistent.
        unique_per_gene = np.asarray(unique.sum(axis=0)).reshape(-1)
        global_diagnostics: dict[str, Any] = {}
        theta = em_gene_abundances(
            em_ec_counts,
            unique_per_gene,
            pseudocount,
            em_max_iter,
            em_tol,
            diagnostics=global_diagnostics,
        )
        method_diagnostics["global_model"] = global_diagnostics
        corr_rows: list[int] = []
        corr_cols: list[int] = []
        corr_data: list[float] = []
        sel_rows: list[int] = []
        sel_cols: list[int] = []
        sel_data: list[float] = []
        cell_theta: dict[int, np.ndarray] = {}
        if method == "em-cell":
            per_cell_ec: dict[int, dict[tuple[int, ...], float]] = {}
            for cell_idx, ec_genes, count in ambiguous_records:
                counts = per_cell_ec.setdefault(cell_idx, {})
                counts[ec_genes] = counts.get(ec_genes, 0.0) + count
            cell_diagnostics: list[dict[str, Any]] = []
            for cell_idx, counts in per_cell_ec.items():
                row = unique.getrow(cell_idx).toarray().reshape(-1)
                diagnostic: dict[str, Any] = {}
                cell_theta[cell_idx] = em_cell_abundances(
                    counts,
                    row,
                    theta,
                    pseudocount,
                    em_max_iter,
                    em_tol,
                    diagnostics=diagnostic,
                )
                cell_diagnostics.append(diagnostic)
            method_diagnostics["cell_models"] = {
                "n_cells": len(cell_diagnostics),
                "n_converged": sum(bool(item["converged"]) for item in cell_diagnostics),
                "max_iterations": max(
                    (int(item["iterations"]) for item in cell_diagnostics), default=0
                ),
                "fallback_reasons": sorted(
                    {
                        str(item["fallback_reason"])
                        for item in cell_diagnostics
                        if item["fallback_reason"]
                    }
                ),
            }

        for cell_idx, ec_genes, count in ambiguous_records:
            gidx = np.asarray(ec_genes, dtype=int)
            w = cell_theta[cell_idx][gidx] if method == "em-cell" else theta[gidx]
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

    audit = audit_counter.freeze()
    audit.validate(float(selected.sum()))
    if not np.isclose(float(unique.sum()), audit.unique_molecules, rtol=0.0, atol=1e-9):
        raise ValueError("Unique molecule layer does not equal the audited unique molecules.")
    return MultimapLayers(
        corrected=selected,
        unique=unique,
        equal=equal,
        host_conservative=host_conservative,
        unique_weighted=unique_weighted,
        unique_viral=_csr_from_entries(
            unique_viral_rows, unique_viral_cols, unique_viral_data, n_cells, n_genes
        ),
        host_viral_ambiguous=_csr_from_entries(
            host_viral_rows, host_viral_cols, host_viral_data, n_cells, n_genes
        ),
        host_viral_selected=host_viral_selected,
        viral_ambiguous_upper=_csr_from_entries(
            upper_rows, upper_cols, upper_data, n_cells, n_genes
        ),
        audit=audit,
        method_diagnostics=method_diagnostics,
    )


def _molecule_evidence_tier(
    unique_umi: float,
    corrected_umi: float,
    host_viral_ambiguous_umi: float,
    selected_host_viral_umi: float,
) -> str:
    """Classify molecule evidence without implying read-level validation."""
    if unique_umi > 0:
        return "candidate_unique"
    non_host_ambiguous_umi = max(corrected_umi - selected_host_viral_umi, 0.0)
    if non_host_ambiguous_umi > 0:
        return "candidate_virus_ambiguous"
    if corrected_umi > 0 or host_viral_ambiguous_umi > 0:
        return "candidate_host_virus_ambiguous"
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
                    "viral_molecules_unique": round(unique_umi, 6),
                    "viral_molecules_ambiguous_allocated": round(ambiguous_umi, 6),
                    "host_virus_ambiguous_molecules": round(host_viral_umi, 6),
                    "viral_molecules_total_est": round(corrected_umi, 6),
                    "viral_molecules_upper_bound": round(upper_bound_umi, 6),
                    "n_unique_viral_cells": int((unique_vec > 0).sum()),
                    "n_ambiguous_viral_cells": int((ambiguous_vec > 0).sum()),
                    "multimap_method": method,
                    "evidence_tier": _molecule_evidence_tier(
                        unique_umi,
                        ambiguous_umi,
                        host_viral_umi,
                        selected_host_viral_umi,
                    ),
                }
            )

    return pd.DataFrame(rows, columns=MULTIMAP_EVIDENCE_COLUMNS)


def select_detection_matrix(adata: Any, config: RunConfig) -> Any:
    """Return v3's complete selected-method molecule matrix."""
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
