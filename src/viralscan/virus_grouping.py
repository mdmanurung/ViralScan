"""Map viral gene IDs to virus names.

Primary path (PLAN ``MECH-A`` step 4): a Run's gene -> virus assignment comes
from its Virus Identity table (``results/virus_identity.tsv``, see
:mod:`viralscan.virus_identity`): :func:`load_run_identity`,
:func:`group_genes_by_identity` and :func:`virus_facts`. Nothing is re-derived
from gene-ID strings.

Legacy fallback: the prefix-grouping domain rule below is used only when a Run
has no table (a run directory made before the table existed, or a caller that
passes none). It is kept so those runs still report.

The prefix rule, historically:

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

import json
import logging
import os
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from viralscan.constants import (
    EVE_RISK_GENERA,
    LOW_COMPLEXITY_RISK_GENERA,
    SIBLING_VIRUS_PAIRS,
    VIRUS_GENE_ID_ALIASES,
    VIRUS_NAME_MAP,
)
from viralscan.virus_identity import (
    ANELLOVIRIDAE,
    TABLE_FILENAME,
    UNCATALOGUED,
    VirusIdentityTable,
    deversion,
)

log = logging.getLogger(__name__)

PathLike = Union[str, "os.PathLike[str]"]

#: ``risk_class`` value marking endogenous-viral-element (EVE) risk.
RISK_EVE = "eve"
#: ``risk_class`` value marking a known low-complexity read-artefact risk (F-019, F-021).
RISK_LOW_COMPLEXITY = "low_complexity"
#: ``claim_scope`` for a family with no orthogonally confirmed positive sample
#: (REF-10): a call may be reported as a screen, never as a confirmed infection.
CLAIM_SCOPE_SCREENING = "screening_only"


def virus_name_for_gene(
    gene_id: str,
    name_map: dict[str, str] | None = None,
    aliases: dict[str, str] | None = None,
) -> str:
    """Return the virus name for ``gene_id``, or the gene ID itself if unmatched.

    Two tiers, tried in order:

    1. **Strict boundary-aware prefix match** over ``name_map`` (default
       :data:`VIRUS_NAME_MAP`): a key matches when the ID equals the key, or
       starts with the key and the next character is ``_`` or a digit. Keys are
       tried longest-first so the most specific prefix wins (``HUM_HERP6B`` over
       ``HUM_HERP6``). See :data:`VIRUS_NAME_MAP`.
    2. **Alias prefix match** over ``aliases`` (default
       :data:`VIRUS_GENE_ID_ALIASES`): plain ``startswith``, longest key first.
       This tier exists only for panel schemes that write the virus token and
       the gene token with no separator (``Ydvgp129``, ``TTVgp1``), which tier 1
       rejects on purpose. Because it runs only after tier 1 has matched
       nothing, it can never change the name of a gene ID that already resolves.

    Returns the raw ``gene_id`` when neither tier matches.
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
    if aliases is None:
        aliases = VIRUS_GENE_ID_ALIASES
    for key in sorted(aliases, key=len, reverse=True):
        if gene_id.startswith(key):
            return aliases[key]
    return gene_id


def group_genes_by_virus(
    gene_ids: Iterable[str],
    name_map: dict[str, str] | None = None,
    aliases: dict[str, str] | None = None,
) -> tuple[dict[str, list[str]], set[str]]:
    """Group gene IDs by virus name.

    Returns ``(group_by_virus, detected_viruses)`` where:

    - ``group_by_virus`` maps each virus name (or the raw gene ID, when unmatched)
      to the list of gene IDs assigned to it, preserving input order.
    - ``detected_viruses`` is the set of virus names that matched a map key or an
      alias (raw-fallback gene IDs are excluded).
    """
    group_by_virus: dict[str, list[str]] = {}
    detected: set[str] = set()
    for gene_id in gene_ids:
        name = virus_name_for_gene(gene_id, name_map, aliases)
        group_by_virus.setdefault(name, []).append(gene_id)
        if name != gene_id:  # matched a map key or alias (unmatched genes map to themselves)
            detected.add(name)
    return group_by_virus, detected


# ----------------------------------------------------------- identity-table path


def identity_path(output: PathLike) -> Path:
    """Where a Run's Virus Identity table lives under its output directory."""
    return Path(output) / "results" / TABLE_FILENAME


def load_run_identity(output: PathLike) -> VirusIdentityTable | None:
    """The Run's Virus Identity table, or ``None`` when the Run has none.

    ``None`` selects the legacy prefix rules (an old run directory). A table
    that exists but cannot be parsed raises: silently falling back would change
    which genes count as viral.
    """
    path = identity_path(output)
    if not path.is_file():
        return None
    return VirusIdentityTable.read_tsv(path)


def group_genes_by_identity(
    gene_ids: Iterable[str], table: VirusIdentityTable
) -> tuple[dict[str, list[str]], set[str]]:
    """Group viral gene IDs by virus: one group per ``virus_key``, named ``virus_name``.

    Same contract as :func:`group_genes_by_virus`: ``(group_by_virus,
    detected_viruses)``, groups keyed by display name in first-seen gene order.
    Display names are unique across keys (the table guarantees it), so a name
    identifies exactly one virus_key. A gene that is host or absent from the
    table is not viral and is skipped, never grouped under its own ID.
    """
    by_gene = table.by_gene()
    group_by_virus: dict[str, list[str]] = {}
    for gene_id in gene_ids:
        identity = by_gene.get(gene_id)
        if identity is None or not identity.viral:
            continue
        group_by_virus.setdefault(identity.virus_name, []).append(gene_id)
    return group_by_virus, set(group_by_virus)


@dataclass(frozen=True)
class VirusFacts:
    """Per-virus facts a report needs, resolved from the identity table."""

    virus_key: str
    virus_name: str
    sibling_group: str
    eve_risk: bool
    #: ``"low_complexity"`` or ``""``: a diagnostic label, never a filter (ANELLO-PRIOR.4).
    artifact_risk: str = ""
    #: ``"screening_only"`` or ``""`` (no restriction): what a call may claim (ANDET-03).
    claim_scope: str = ""
    #: Curated reference role; empty means uncurated/unknown, never a target assertion.
    role: str = ""


def legacy_eve_risk(virus_name: str) -> bool:
    """The retired substring test: does the name contain an EVE-risk genus?"""
    return any(genus in virus_name for genus in EVE_RISK_GENERA)


def legacy_artifact_risk(virus_name: str) -> str:
    """Genus-name fallback for ``artifact_risk`` when no catalogue risk_class exists."""
    return RISK_LOW_COMPLEXITY if any(g in virus_name for g in LOW_COMPLEXITY_RISK_GENERA) else ""


def legacy_claim_scope(virus_name: str) -> str:
    """Genus-name fallback for ``claim_scope`` when no catalogue family exists.

    ponytail: reuses the low-complexity genus list, which today is exactly the
    Anelloviridae plus "Torque teno virus"; give it its own list if they diverge.
    """
    return CLAIM_SCOPE_SCREENING if legacy_artifact_risk(virus_name) else ""


def virus_facts(table: VirusIdentityTable) -> dict[str, VirusFacts]:
    """virus display name -> :class:`VirusFacts`.

    ``eve_risk`` fails closed:

    - a virus is at risk when *any* of its genes has ``risk_class == "eve"``, so
      a catalogue that disagrees across one virus's genomes keeps the flag;
    - a virus whose identity is not catalogued (``uncatalogued``, or
      ``legacy_prefix`` still keyed ``name:``) has no ``risk_class`` to read, so the retired genus
      substring test decides, and the Run warns. Never silently ``False``.
    """
    risk: dict[str, bool] = {}
    artifact: dict[str, bool] = {}
    anello: dict[str, bool] = {}
    sibling: dict[str, str] = {}
    key_of: dict[str, str] = {}
    unknown: dict[str, bool] = {}
    roles: dict[str, set[str]] = {}
    for g in table.genes:
        if not g.viral:
            continue
        name = g.virus_name
        roles.setdefault(name, set())
        if g.role:
            roles[name].add(g.role)
        key_of.setdefault(name, g.virus_key)
        risk[name] = risk.get(name, False) or g.risk_class == RISK_EVE
        artifact[name] = artifact.get(name, False) or g.risk_class == RISK_LOW_COMPLEXITY
        anello[name] = anello.get(name, False) or g.family == ANELLOVIRIDAE
        if not sibling.get(name):
            sibling[name] = g.sibling_group
        # A legacy prefix name that adopted a catalogue row has a real risk_class.
        if g.status == UNCATALOGUED or g.virus_key.startswith("name:"):
            unknown[name] = True
    facts: dict[str, VirusFacts] = {}
    for name, key in key_of.items():
        if len(roles[name]) > 1:
            raise ValueError(f"{name} has conflicting catalogue roles: {sorted(roles[name])}")
        eve = risk[name]
        art = RISK_LOW_COMPLEXITY if artifact[name] else ""
        scope = CLAIM_SCOPE_SCREENING if anello[name] else ""
        if unknown.get(name) and not scope:
            scope = legacy_claim_scope(name)
        if unknown.get(name) and not (eve or art):
            eve, art = legacy_eve_risk(name), legacy_artifact_risk(name)
            log.warning(
                "%s has no catalogue risk_class (status uncatalogued/legacy); "
                "eve_risk=%s and artifact_risk=%r come from the genus-name fallback.",
                name,
                eve,
                art,
            )
        facts[name] = VirusFacts(
            key, name, sibling[name], eve, art, scope, next(iter(roles[name]), "")
        )
    return facts


#: Genome-level host homology columns of ``viral_summary.tsv`` (PLAN ``ANDET-02``).
HOST_HOMOLOGY_COLUMNS = (
    "host_homology_status",
    "host_homology_max_identity",
    "host_homology_max_query_coverage",
    "host_homology_max_aligned_bases",
)
HOST_HOMOLOGY_NOT_MEASURED = dict.fromkeys(HOST_HOMOLOGY_COLUMNS, "") | {
    "host_homology_status": "not_measured"
}


def host_homology_by_virus(table: VirusIdentityTable, index: PathLike | None) -> dict[str, dict]:
    """virus display name -> :data:`HOST_HOMOLOGY_COLUMNS`, from the index's reference manifest.

    ``viralscan build-ref --genome-dlist`` measures each viral genome against the host
    genome and records the maximum alignment per genome in ``reference_manifest.json``
    beside the index. A virus reports the maximum over its genomes. ``status`` is
    ``measured`` (every genome), ``partial`` (some) or ``not_measured`` (none, or no
    manifest): absence is never reported as zero homology. These are genome-level
    measurements with no threshold; where the *reads* land is the evidence route's job.
    """
    measured: dict[str, dict] = {}
    path = Path(index).parent / "reference_manifest.json" if index else None
    if path is not None and path.is_file():
        try:
            records = json.loads(path.read_text(encoding="utf-8")).get("sequences", [])
        except (OSError, ValueError) as exc:
            log.warning(
                "Unreadable reference manifest %s (%s): host homology not measured.", path, exc
            )
            records = []
        for rec in records:
            if rec.get("host_homology_status") == "measured":
                measured[deversion(rec["accession_version"].upper())] = rec
    genomes: dict[str, set[str]] = {}
    for g in table.genes:
        if g.viral:
            genomes.setdefault(g.virus_name, set()).add(deversion(g.genome_accession.upper()))
    out: dict[str, dict] = {}
    for name, accessions in genomes.items():
        hits = [measured[a] for a in accessions if a in measured]
        if not hits:
            out[name] = dict(HOST_HOMOLOGY_NOT_MEASURED)
            continue
        out[name] = {
            "host_homology_status": "measured" if len(hits) == len(accessions) else "partial",
            "host_homology_max_identity": max(h["host_homology_max_identity"] for h in hits),
            "host_homology_max_query_coverage": max(
                h["host_homology_max_query_coverage"] for h in hits
            ),
            "host_homology_max_aligned_bases": max(
                h["host_homology_max_aligned_bases"] for h in hits
            ),
        }
    return out


def legacy_sibling_groups() -> dict[str, str]:
    """``SIBLING_VIRUS_PAIRS`` as ``{virus_name: group}`` for runs with no table."""
    groups: dict[str, str] = {}
    for a, b in SIBLING_VIRUS_PAIRS.items():
        group = groups.get(a) or groups.get(b) or f"{a}|{b}"
        groups[a] = groups[b] = group
    return groups
