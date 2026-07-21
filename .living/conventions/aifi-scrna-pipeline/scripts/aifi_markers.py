"""
AIFI marker annotation library — programmatic access.

Loads the marker panels, L2 doublet thresholds, and L3 refinement rules extracted
from aifimmunology/sound-life-scrna-analysis, and provides the cluster-level
filtering primitives AIFI used to apply them.

Data files live in ../assets/ relative to this script.

Usage
-----
    from aifi_markers import (
        broad_markers, class_markers, markers_for, marker_frac_df,
        remove_cl, extract_cl, l2_rules_for, l3_rules_for,
    )

    panel = class_markers('nk')                  # NK class panel
    genes = markers_for('CD4 MAIT')              # genes AIFI checked for this type
    adata, n = remove_cl(adata, 'IL7R', 'above', 0.4, 'leiden_2')
"""

import json
import csv
import os
from functools import lru_cache

_ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets")

LINEAGES = [
    "b_cell",
    "cd4_t_cell",
    "cd8_t_cell",
    "other_t_cell",
    "myeloid",
    "nk",
    "other",
]


# ---------------------------------------------------------------- loading

@lru_cache(maxsize=1)
def _library():
    with open(os.path.join(_ASSETS, "aifi_marker_library.json")) as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def _l2_thresholds():
    with open(os.path.join(_ASSETS, "aifi_l2_doublet_thresholds.json")) as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def _l3_rules():
    with open(os.path.join(_ASSETS, "aifi_l3_refinement_rules.csv")) as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- panels

def broad_markers(lineage):
    """Cross-lineage contamination markers for a lineage. Returns list of gene symbols."""
    _check_lineage(lineage)
    return [g for g, _ in _library()["lineage_panels"][lineage]["broad_markers"]]


def class_markers(lineage):
    """Within-lineage identity markers. Returns list of gene symbols."""
    _check_lineage(lineage)
    return [g for g, _ in _library()["lineage_panels"][lineage]["class_markers"]]


def annotated_panel(lineage, which="class_markers"):
    """Return [(gene, comment), ...] preserving AIFI's inline annotations."""
    _check_lineage(lineage)
    return [tuple(x) for x in _library()["lineage_panels"][lineage][which]]


def markers_for(cell_type):
    """Genes AIFI inspected for a specific AIFI_L3 population.

    Returns [] when the population had no ad-hoc marker check recorded, which is
    not the same as having no markers — the lineage panel still applies.
    """
    return list(_library()["per_celltype_markers"].get(cell_type, []))


def l3_roster(lineage=None):
    """AIFI_L3 labels, optionally restricted to one lineage partition."""
    roster = _library()["l3_roster_by_lineage"]
    if lineage is None:
        return {k: list(v) for k, v in roster.items()}
    _check_lineage(lineage)
    return list(roster[lineage])


def lineage_of(cell_type):
    """Which clustering partition an AIFI_L3 label belongs to, or None."""
    for lin, types in _library()["l3_roster_by_lineage"].items():
        if cell_type in types:
            return lin
    return None


def _check_lineage(lineage):
    if lineage not in LINEAGES:
        raise KeyError(f"unknown lineage {lineage!r}; expected one of {LINEAGES}")


# ---------------------------------------------------------------- rules

def l2_rules_for(population, kind="fraction_based"):
    """Doublet rules for an AIFI_L2 class.

    Returns a list of {'reason', 'gene', 'gt'} in source order. Duplicate reasons
    are preserved deliberately: three populations (CD14 monocyte, cDC1, Effector
    B cell) carry two T cell doublet rules combined by OR, which a dict keyed by
    reason would silently collapse. See references/l2_doublet_filter_rules.md.
    """
    table = _l2_thresholds()[kind]
    val = table.get(population, [] if kind == "fraction_based" else None)
    if kind == "mean_expression_based":
        return [val] if val else []
    return list(val)


def l3_rules_for(cell_type):
    """Refinement rules for an AIFI_L3 label: gene, direction, cutoff, action."""
    return [r for r in _l3_rules() if r["cell_type"] == cell_type]


# ---------------------------------------------------------------- primitives

def marker_frac_df(adata, markers, clusters="leiden_2"):
    """Fraction of cells per cluster with detected expression of each marker.

    This reproduces AIFI's approach of reading the dot_size_df off a scanpy
    dotplot rather than computing the fraction directly.
    """
    import scanpy as sc

    if isinstance(markers, str):
        markers = [markers]
    fig = sc.pl.dotplot(
        adata, groupby=clusters, var_names=markers, return_fig=True
    )
    return fig.dot_size_df


def _select_clusters(adata, gene, direction, cutoff, clusters):
    frac = marker_frac_df(adata, gene, clusters)
    col = frac[gene]
    if direction == "above":
        return frac.index[col > cutoff].tolist()
    elif direction == "below":
        return frac.index[col < cutoff].tolist()
    raise ValueError(f"direction must be 'above' or 'below', got {direction!r}")


def remove_cl(adata, gene, direction, cutoff, clusters="leiden_2", verbose=True):
    """Drop whole clusters passing a marker threshold. Returns (adata, n_removed)."""
    sel = _select_clusters(adata, gene, direction, cutoff, clusters)
    if verbose:
        print(f"removing clusters {sel} ({gene} {direction} {cutoff})")
    remove_idx = adata.obs[clusters].isin(sel)
    n_removed = int(remove_idx.sum())
    return adata[~remove_idx].copy(), n_removed


def extract_cl(adata, gene, direction, cutoff, clusters="leiden_2", verbose=True):
    """Pull out whole clusters passing a threshold. Returns (subset, n_extracted)."""
    sel = _select_clusters(adata, gene, direction, cutoff, clusters)
    if verbose:
        print(f"extracting clusters {sel} ({gene} {direction} {cutoff})")
    sel_idx = adata.obs[clusters].isin(sel)
    return adata[sel_idx].copy(), int(sel_idx.sum())


def apply_l2_filter(adata, population, clusters="leiden_2", verbose=True):
    """Apply every published fraction-based doublet rule for an L2 class in order.

    Returns (adata, log) where log is a list of {'reason','gene','gt','n_removed'}.
    """
    log = []
    for rule in l2_rules_for(population):
        adata, n = remove_cl(
            adata, rule["gene"], "above", rule["gt"], clusters, verbose=verbose
        )
        log.append({**rule, "n_removed": n})
    return adata, log


if __name__ == "__main__":
    print(f"lineages: {len(LINEAGES)}")
    roster = l3_roster()
    print(f"L3 populations: {sum(len(v) for v in roster.values())}")
    print(f"L2 classes with rules: {len(_l2_thresholds()['fraction_based'])}")
    print(f"L3 refinement rules: {len(_l3_rules())}")
    print(f"\nNK class panel: {class_markers('nk')}")
    print(f"CD4 MAIT checks: {markers_for('CD4 MAIT')}")
    print(f"CD14 monocyte L2 rules: {l2_rules_for('CD14 monocyte')}")
