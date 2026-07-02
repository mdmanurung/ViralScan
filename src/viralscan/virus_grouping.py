"""Map viral gene IDs to virus names — the prefix-grouping domain rule.

This logic was previously duplicated (and had drifted) across ``detection.py``
(substring match) and ``umap.py`` (prefix match). The two rules disagreed on 151
of the 2692 bundled gene IDs, and both were buggy: substring over-matched
(``EPSTEIN_HHV4_BORF1`` matched key ``ORF``), prefix under-matched
(``TTV7_gp2`` matched nothing because the key is ``TTV``).

The unified rule here is **boundary-aware**: a key matches a gene ID when the ID
equals the key, or starts with the key and the next character is ``_`` or a
digit. Keys are tried longest-first so the most specific prefix wins (e.g.
``HUM_HERP6B`` over ``HUM_HERP6``). See ``CONTEXT.md`` ("Viral Accession").
"""

from __future__ import annotations

from collections.abc import Iterable

from viralscan.constants import VIRUS_NAME_MAP


def virus_name_for_gene(gene_id: str, name_map: dict[str, str] | None = None) -> str:
    """Return the virus name for ``gene_id``, or the gene ID itself if unmatched.

    Boundary-aware prefix match (see module docstring). ``name_map`` defaults to
    :data:`VIRUS_NAME_MAP`.
    """
    if name_map is None:
        name_map = VIRUS_NAME_MAP
    for key in sorted(name_map, key=len, reverse=True):
        if gene_id == key:
            return name_map[key]
        if gene_id.startswith(key):
            nxt = gene_id[len(key)]
            if nxt == "_" or nxt.isdigit():
                return name_map[key]
    return gene_id


def group_genes_by_virus(
    gene_ids: Iterable[str], name_map: dict[str, str] | None = None
) -> tuple[dict[str, list[str]], set[str]]:
    """Group gene IDs by virus name.

    Returns ``(group_by_virus, detected_viruses)`` where:

    - ``group_by_virus`` maps each virus name (or the raw gene ID, when unmatched)
      to the list of gene IDs assigned to it, preserving input order.
    - ``detected_viruses`` is the set of virus names that matched a map key
      (raw-fallback gene IDs are excluded).
    """
    group_by_virus: dict[str, list[str]] = {}
    detected: set[str] = set()
    for gene_id in gene_ids:
        name = virus_name_for_gene(gene_id, name_map)
        group_by_virus.setdefault(name, []).append(gene_id)
        if name != gene_id:  # matched a map key (unmatched genes map to themselves)
            detected.add(name)
    return group_by_virus, detected
