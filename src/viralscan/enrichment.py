"""Cell-type enrichment utilities.

Extracted from ``scripts/detection.py`` so these functions can be imported and
tested independently without triggering the Snakemake magic-globals that are
evaluated at ``detection.py`` module level.
"""

from __future__ import annotations

import logging
import os
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import fisher_exact

log = logging.getLogger("viralscan")


def _bh_adjust(pvals: list[float] | npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Benjamini-Hochberg correction for a 1D list/array of p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return np.array([])

    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = min(prev, ranked[i] * n / rank)
        adjusted[i] = value
        prev = value

    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def cell_type_enrichment(
    adata: Any, group_by_virus: dict[str, list[str]], cfg: dict[str, Any]
) -> pd.DataFrame:
    """Compute per-virus enrichment by cell type using Fisher exact tests."""
    cell_types_path = cfg.get("cell_types")
    if not cell_types_path:
        return pd.DataFrame()

    if not os.path.exists(cell_types_path):
        log.warning("cell_types CSV not found at %s; skipping enrichment.", cell_types_path)
        return pd.DataFrame()

    try:
        labels = pd.read_csv(cell_types_path)
    except Exception as exc:
        log.warning("Failed to read cell_types CSV (%s); skipping enrichment.", exc)
        return pd.DataFrame()

    if labels.empty:
        return pd.DataFrame()

    if "barcode" not in labels.columns or "cell_type" not in labels.columns:
        # Fallback for files with unnamed first two columns.
        if len(labels.columns) >= 2:
            labels = labels.rename(
                columns={labels.columns[0]: "barcode", labels.columns[1]: "cell_type"}
            )
        else:
            log.warning("cell_types CSV must contain barcode and cell_type columns; skipping.")
            return pd.DataFrame()

    labels = labels[["barcode", "cell_type"]].dropna()
    if labels.empty:
        return pd.DataFrame()

    labels["barcode"] = labels["barcode"].astype(str)
    labels["cell_type"] = labels["cell_type"].astype(str)
    labels = labels.drop_duplicates(subset=["barcode"], keep="first")

    obs_labels = pd.DataFrame({"barcode": adata.obs_names.astype(str)})
    merged = obs_labels.merge(labels, on="barcode", how="left")
    merged = merged[merged["cell_type"].notna()].copy()
    if merged.empty:
        log.warning("No barcodes overlap between AnnData and cell_types CSV; skipping.")
        return pd.DataFrame()

    barcode_to_idx = {bc: i for i, bc in enumerate(adata.obs_names.astype(str))}
    merged["idx"] = merged["barcode"].map(barcode_to_idx)
    merged = merged[merged["idx"].notna()].copy()
    merged["idx"] = merged["idx"].astype(int)
    if merged.empty:
        return pd.DataFrame()

    idx = merged["idx"].to_numpy()
    cell_types = merged["cell_type"].to_numpy()
    unique_cell_types = sorted(pd.unique(cell_types))
    type_masks = {ct: (cell_types == ct) for ct in unique_cell_types}

    rows = []
    for virus, gene_list in group_by_virus.items():
        valid_genes = [g for g in gene_list if g in adata.var_names]
        if not valid_genes:
            continue

        viral_matrix = adata[:, valid_genes].X
        if hasattr(viral_matrix, "toarray"):
            viral_matrix = viral_matrix.toarray()

        infected_mask = viral_matrix.sum(axis=1) > 0
        infected_in_labeled = infected_mask[idx]

        for ct in unique_cell_types:
            ct_mask = type_masks[ct]
            n_total = int(ct_mask.sum())
            if n_total == 0:
                continue

            a = int(np.sum(infected_in_labeled[ct_mask]))
            b = int(n_total - a)
            c = int(np.sum(infected_in_labeled[~ct_mask]))
            d = int(np.sum((~infected_in_labeled)[~ct_mask]))

            odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append(
                {
                    "virus": virus,
                    "cell_type": ct,
                    "n_infected": a,
                    "n_total": n_total,
                    "pct": round((a / n_total) * 100, 4),
                    "OR": float(odds_ratio),
                    "pvalue": float(pvalue),
                }
            )

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    result["padj"] = _bh_adjust(result["pvalue"].to_numpy())
    return (
        result[["virus", "cell_type", "n_infected", "n_total", "pct", "OR", "pvalue", "padj"]]
        .sort_values(["virus", "padj", "pvalue", "cell_type"])
        .reset_index(drop=True)
    )


def write_cell_type_enrichment(cell_type_df: pd.DataFrame, outputpath: str) -> str | None:
    """Write cell-type enrichment table when data are available."""
    if cell_type_df is None or cell_type_df.empty:
        return None

    results_dir = os.path.join(outputpath, "results")
    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, "cell_type_enrichment.tsv")
    cell_type_df.to_csv(out_path, sep="\t", index=False)
    log.info("Wrote results/cell_type_enrichment.tsv")
    return out_path
