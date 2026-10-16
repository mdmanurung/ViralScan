"""Sum viral molecules up the NCBI lineage when reads cannot tell related viruses apart (PLAN ``TAXO-01``).

A molecule compatible with EBV and EBV type 2 counts for neither virus alone, but it is certainly
*Lymphocryptovirus humangamma4*. The multimap step stores, per cell, the molecules unique to each virus,
species and genus and those merely compatible with it (``stamp_viral_groups``). This module turns them
into ``results/viral_summary_by_species.tsv`` / ``viral_summary_by_genus.tsv`` and the per-virus
``unique_fraction`` of ``viral_summary.tsv``. Nothing here changes a count.

For one group, ``molecules_unique <= molecules_package_total <= molecules_compatible`` holds under every
multimap method; ``unique_fraction = molecules_unique / molecules_compatible``. A low fraction means the
group's own total is an allocation, and the next level up is the number to report.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pandas as pd

from viralscan.constants import ASSIGNMENT_UNIQUE_FRACTION_MIN
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

LEVELS = ("virus", "species", "genus")
GROUP_COLUMNS = [
    "level",
    "group_key",
    "group_name",
    "members",
    "n_members",
    "molecules_unique",
    "molecules_package_total",
    "molecules_compatible",
    "molecules_host_shared",
    "unique_fraction",
    "molecules_unique_called",
    "molecules_package_total_called",
    "molecules_compatible_called",
    "n_called_cells_unique",
    "n_called_cells_compatible",
    "multimap_method",
]
FILES = {"species": "viral_summary_by_species.tsv", "genus": "viral_summary_by_genus.tsv"}


def _col_sum(m: Any, rows: np.ndarray | None = None) -> np.ndarray:
    m = m if rows is None else m[rows]
    return np.asarray(m.sum(axis=0)).reshape(-1)


def _col_nnz(m: Any, rows: np.ndarray) -> np.ndarray:
    return np.asarray((m[rows] > 0).sum(axis=0)).reshape(-1)


def group_summaries(
    adata: Any, identity: VirusIdentityTable | None, called_mask: Any = None
) -> dict[str, pd.DataFrame]:
    """level -> one row per group with any compatible molecule; ``{}`` when the run has no groups."""
    meta = adata.uns.get("viral_groups") if adata is not None else None
    if identity is None or meta is None or "viral_groups_compatible" not in adata.obsm:
        return {}
    called = (
        np.ones(adata.n_obs, bool) if called_mask is None else np.asarray(called_mask, dtype=bool)
    )
    unique, compat, host = (
        adata.obsm[f"viral_groups_{k}"] for k in ("unique", "compatible", "host_shared")
    )
    sums = {
        "molecules_unique": _col_sum(unique),
        "molecules_compatible": _col_sum(compat),
        "molecules_host_shared": _col_sum(host),
        "molecules_unique_called": _col_sum(unique, called),
        "molecules_compatible_called": _col_sum(compat, called),
        "n_called_cells_unique": _col_nnz(unique, called),
        "n_called_cells_compatible": _col_nnz(compat, called),
    }
    var_index = {g: i for i, g in enumerate(adata.var_names)}
    members: dict[tuple[str, str], set[str]] = {}
    genes: dict[tuple[str, str], list[int]] = {}
    for gene in identity.genes:
        if gene.viral and gene.gene_id in var_index:
            for level in LEVELS:
                key = (level, gene.group_key(level))
                members.setdefault(key, set()).add(gene.virus_name)
                genes.setdefault(key, []).append(var_index[gene.gene_id])
    x_all, x_called = _col_sum(adata.X), _col_sum(adata.X, called)
    method = str(adata.uns.get("multimap_method", ""))
    rows: dict[str, list[dict[str, Any]]] = {level: [] for level in LEVELS}
    for j, (level, key, name) in enumerate(
        zip(meta["level"], meta["key"], meta["name"], strict=True)
    ):
        if sums["molecules_compatible"][j] <= 0:
            continue
        idx = genes.get((level, key), [])
        row: dict[str, Any] = {"level": level, "group_key": key, "group_name": name}
        row["members"] = "|".join(sorted(members.get((level, key), ())))
        row["n_members"] = len(members.get((level, key), ()))
        row.update(
            {k: float(v[j]) if k.startswith("molecules") else int(v[j]) for k, v in sums.items()}
        )
        row["molecules_package_total"] = round(float(x_all[idx].sum()), 6)
        row["molecules_package_total_called"] = round(float(x_called[idx].sum()), 6)
        row["unique_fraction"] = round(row["molecules_unique"] / row["molecules_compatible"], 6)
        row["multimap_method"] = method
        rows[str(level)].append(row)
    return {
        level: pd.DataFrame(r, columns=GROUP_COLUMNS).sort_values(
            "molecules_compatible", ascending=False, kind="stable"
        )
        for level, r in rows.items()
    }


def write_group_summaries(frames: dict[str, pd.DataFrame], outputpath: str) -> None:
    """Write the species and genus tables (header only when the run has no groups)."""
    results = os.path.join(outputpath, "results")
    os.makedirs(results, exist_ok=True)
    for level, name in FILES.items():
        frames.get(level, pd.DataFrame(columns=GROUP_COLUMNS)).to_csv(
            os.path.join(results, name), sep="\t", index=False
        )


def assignment_columns(
    frames: dict[str, pd.DataFrame], identity: VirusIdentityTable | None
) -> dict[str, dict[str, Any]]:
    """virus_name -> the ``viral_summary.tsv`` columns ``molecules_unique_to_virus``,
    ``molecules_compatible``, ``unique_fraction`` and ``aggregate_note``.

    The note names the first level up whose own unique fraction clears
    :data:`~viralscan.constants.ASSIGNMENT_UNIQUE_FRACTION_MIN`, or says none does.
    """
    if identity is None or "virus" not in frames:
        return {}
    by_key = {lvl: f.set_index("group_key") for lvl, f in frames.items()}
    first: dict[str, GeneIdentity] = {}
    for row in identity.genes:
        if row.viral:
            first.setdefault(row.virus_name, row)
    out: dict[str, dict[str, Any]] = {}
    for name, gene in first.items():
        v = by_key["virus"].loc[gene.virus_key] if gene.virus_key in by_key["virus"].index else None
        if v is None:
            continue
        frac = float(v["unique_fraction"])
        note = ""
        if frac < ASSIGNMENT_UNIQUE_FRACTION_MIN:
            note = (
                f"unique_fraction {frac:.3f} < {ASSIGNMENT_UNIQUE_FRACTION_MIN:g}: "
                "total is mostly an allocation"
            )
            for level in ("species", "genus"):
                key = gene.group_key(level)
                g = by_key[level].loc[key] if key in by_key[level].index else None
                if g is None or int(g["n_members"]) < 2:
                    continue
                g_frac = float(g["unique_fraction"])
                note += (
                    f"; {level} {g['group_name']}: unique_fraction {g_frac:.3f} ({FILES[level]})"
                )
                if g_frac >= ASSIGNMENT_UNIQUE_FRACTION_MIN:
                    break
        out[name] = {
            "molecules_unique_to_virus": float(v["molecules_unique"]),
            "molecules_compatible": float(v["molecules_compatible"]),
            "unique_fraction": frac,
            "aggregate_note": note,
        }
    return out
