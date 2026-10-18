"""Frozen pre-optimization oracles, captured original source.

The multimapping oracles are re-based to ``origin/main`` (post-Q123 "LOCK-01 D1/D2", 329b3dc),
the pre-optimization baseline after the wave-1 rebase; the rest stay at 5e551c9. Per source:

  multimapping.py  origin/main  sha256 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
  detection.py     5e551c9      sha256 6b3fc3f46817e9470106e9c85bd8f5841215ee4d32c922eff83178fce0e4cec2
  umap.py          5e551c9      sha256 802d2c9d3ecb39f7cb485b8ca360ad09da6342d88ce312d93bd10faa7e0d6255

Only ``MULTIMAP_EVIDENCE_COLUMNS`` and ``summarize_multimap_evidence`` differ from 5e551c9;
the Q123 columns (``locus_id``, ``n_unique_positive_cells_called``) are carried here so the
identity check is against the same feature set, not against the feature's absence. main's
only other change under these oracles is a docstring in ``viral_neighbor_enrichment``.

Helpers and schema are independent of the candidate. Full original sources and
independent file goldens are retained in the execution packet.
"""

from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from viralscan.runconfig import RunConfig

# multimapping.py @ origin/main SHA-256: 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
# (5e551c9 was 2b06374abcbee607b9a31023c50ccd15d5662e752e8cbcc809effa81aed67b1b)
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
    "locus_id",
    "n_unique_positive_cells_called",
]


# multimapping.py @ origin/main SHA-256: 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
# (5e551c9 was 2b06374abcbee607b9a31023c50ccd15d5662e752e8cbcc809effa81aed67b1b)
def _empty_matrix(n_cells: int, n_genes: int) -> sparse.csr_matrix:
    return sparse.csr_matrix((n_cells, n_genes), dtype=float)


# multimapping.py @ origin/main SHA-256: 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
# (5e551c9 was 2b06374abcbee607b9a31023c50ccd15d5662e752e8cbcc809effa81aed67b1b)
def _sum_gene(matrix: Any, gene_idx: int) -> float:
    values = matrix[:, gene_idx]
    return float(values.sum())


# multimapping.py @ origin/main SHA-256: 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
# (5e551c9 was 2b06374abcbee607b9a31023c50ccd15d5662e752e8cbcc809effa81aed67b1b)
def _gene_vector(matrix: Any, gene_idx: int) -> Any:
    values = matrix[:, gene_idx]
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).reshape(-1)


# multimapping.py @ origin/main SHA-256: 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
# (5e551c9 was 2b06374abcbee607b9a31023c50ccd15d5662e752e8cbcc809effa81aed67b1b)
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


# Original file SHA-256: 6b3fc3f46817e9470106e9c85bd8f5841215ee4d32c922eff83178fce0e4cec2
def _count_value(value, ndigits=6):
    """Preserve fractional multimapper counts while keeping whole counts tidy."""
    value = float(value)
    rounded = round(value)
    if np.isclose(value, rounded):
        return int(rounded)
    return round(value, ndigits)


# Original file SHA-256: 802d2c9d3ecb39f7cb485b8ca360ad09da6342d88ce312d93bd10faa7e0d6255
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


# multimapping.py @ origin/main SHA-256: 076172803de664ceac8177d019d7c215a8b9fbdeae4a5aec239bfa049bace902
# (5e551c9 was 2b06374abcbee607b9a31023c50ccd15d5662e752e8cbcc809effa81aed67b1b)
def summarize_multimap_evidence(
    adata: Any,
    group_by_virus: dict[str, list[str]],
    config: RunConfig,
    called_mask: Any = None,
    locus_by_gene: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Summarize unique and ambiguous viral evidence for each grouped gene.

    ``locus_id`` is the gene's independent locus (``viral_loci``; empty when the run has no annotation).
    ``n_unique_positive_cells_called`` is per virus, repeated on each of its rows: the number of called
    barcodes (``called_mask``; None = every barcode) holding at least one unique-layer molecule at any gene
    of the virus, the union and not the maximum over genes (LOCK-01 D2)."""
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
    called = (
        np.ones(n_cells, dtype=bool) if called_mask is None else np.asarray(called_mask, dtype=bool)
    )
    locus_by_gene = locus_by_gene or {}

    rows = []
    for virus, gene_ids in group_by_virus.items():
        idx = [int(adata.var_names.get_loc(g)) for g in gene_ids if g in adata.var_names]
        positive = (
            np.asarray(unique[:, idx].sum(axis=1)).reshape(-1) > 0
            if idx
            else np.zeros(n_cells, bool)
        )
        n_positive = int((positive & called).sum())
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
                    "locus_id": locus_by_gene.get(gene_id, ""),
                    "n_unique_positive_cells_called": n_positive,
                }
            )

    return pd.DataFrame(rows, columns=MULTIMAP_EVIDENCE_COLUMNS)
