"""The per-Run Virus Identity table: which indexed genes are viral, and which virus.

PLAN ``MECH-A``. Before this module every consumer re-derived a gene's viral
status and virus name from string formats: ``log/analysis.txt`` (GTF gene IDs)
for the viral/host partition, and prefix matching of gene IDs against name maps
for the virus. The two disagreed, and the prefix rule cannot see a genome
accession, so one genome could appear under two names and 1,912 max-panel
genomes would have been reported as raw gene IDs.

The table is built once per Run from the index's own t2g, whose fifth column is
the genome accession (the GTF seqname) of each transcript, joined to the
packaged catalogue (``virus_catalog.tsv``). One row per indexed gene_id.

Resolution, in order (the first that applies wins):

1. **catalogued** — column 5 is a catalogued accession (versioned or bare).
   Viral, whether or not the gene is in a GTF.
2. **uncatalogued** — the gene is in the GTF gene set (bundled panel, ``--gtf``,
   anellovirus candidates) and its column 5 is not a transcript of the index.
   Viral, keyed and named by its genome accession; the Run warns with a count.
3. **host** — everything else.
4. **legacy_prefix** — a t2g without a fifth column (a pre-v3 prebuilt index)
   carries no genome accession, so only the GTF gene set marks a gene viral and
   the legacy prefix maps name it.

An index built by ``viralscan build-ref`` or ``viralscan --reference`` also
carries a **build manifest** (:func:`write_build_manifest`): the host and viral
gene sets of the index, in ``<index>.build_manifest.json``. When the Run finds
it, the manifest replaces the ``--gtf`` gene set, and a Run whose ``--gtf``
would classify the index's genes differently is refused
(:class:`BuildManifestContradiction`). Gene IDs are compared de-versioned, so
Ensembl release drift in the host set is accepted. The structural guard below
stays as a backstop.

The rule-2 guard exists because a combined index built by ``kb ref`` from a host
cDNA FASTA writes the transcript ID, not a chromosome, into column 5 of every
host row. A combined GTF passed as ``--gtf`` would otherwise turn every host
gene viral.

Virus key (what one "virus" row of the outputs is):

- ``taxid:<n>`` — the NCBI taxid. A segmented virus's segments share their
  isolate's taxid, so they group into one row;
- ``taxid:<n>|<strain>`` — one isolate of a segmented virus whose taxid is
  shared by several isolates (the species-level taxid 11320 "Influenza A
  virus" holds four). A taxid is split this way only when it holds the same
  segment twice. The strain text alone is not a key: one RefSeq set leaves it
  empty on some segments and spells it "A/goose/Guangdong/1/96(H5N1)" or
  "A/goose/Guangdong/1/1996" on others, while the isolate taxid (93838) is
  the same on all eight;
- ``genus:<genus>`` — Anelloviridae, reported per genus (``ANDET-05``). NCBI
  taxids for anelloviruses are catch-all bins: taxid 2055263 "Anelloviridae
  sp." spans four genera, so a taxid cannot carry the genus rollup;
- ``accession:<acc>`` — an uncatalogued genome, or a catalogue row with no taxid;
- ``name:<name>`` — the legacy prefix fallback.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Union

from viralscan import __version__ as _VIRALSCAN_VERSION
from viralscan import virus_catalog

log = logging.getLogger(__name__)

TABLE_FILENAME = "virus_identity.tsv"

CATALOGUED = "catalogued"
UNCATALOGUED = "uncatalogued"
HOST = "host"
LEGACY_PREFIX = "legacy_prefix"

#: Index Kind: ``combined`` holds host genes, ``virus_only`` holds none.
COMBINED = "combined"
VIRUS_ONLY = "virus_only"

ANELLOVIRIDAE = "Anelloviridae"
#: Display name of anelloviruses with no assigned genus (user decision Q8b = A2).
ANELLO_UNASSIGNED = "Anelloviridae (genus unassigned)"

PathLike = Union[str, "os.PathLike[str]"]


@dataclass(frozen=True)
class GeneIdentity:
    """One indexed gene and the virus it belongs to (empty fields for host)."""

    gene_id: str
    genome_accession: str
    status: str
    viral: bool
    virus_key: str = ""
    virus_name: str = ""
    organism: str = ""
    taxid: str = ""
    species: str = ""
    genus: str = ""
    family: str = ""
    segment: str = ""
    strain: str = ""
    sibling_group: str = ""
    risk_class: str = ""
    role: str = ""


COLUMNS: tuple[str, ...] = tuple(f.name for f in fields(GeneIdentity))


# --------------------------------------------------------------------------- t2g


@dataclass(frozen=True)
class T2gGenes:
    """The per-gene facts the table needs from a kb ``t2g``.

    ``accession`` maps gene_id to the column-5 value of its first transcript
    (``""`` in a legacy t2g). ``structural_host`` holds genes whose column 5
    names a transcript of the index rather than a genome sequence.
    """

    accession: dict[str, str]
    structural_host: frozenset[str]
    legacy: bool


def read_t2g(path: PathLike) -> T2gGenes:
    """Parse a kb t2g with a strict tab split.

    Host rows of a combined t2g leave columns 3 and 4 empty
    (``ENST…\\tENSG…\\t\\t\\tENST…``); a whitespace split would collapse them and
    shift column 5. A line with no tab at all (hand-written legacy t2g files) is
    split on whitespace.
    """
    rows: list[list[str]] = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t") if "\t" in line else line.split()
            if len(cols) >= 2:
                rows.append(cols)
    legacy = any(len(cols) < 5 for cols in rows)
    transcripts = {cols[0] for cols in rows}
    accession: dict[str, str] = {}
    structural_host: set[str] = set()
    for cols in rows:
        gene_id = cols[1]
        genome = "" if legacy else cols[4].strip()
        accession.setdefault(gene_id, genome)
        if genome and genome in transcripts:
            structural_host.add(gene_id)
    return T2gGenes(accession, frozenset(structural_host), legacy)


# --------------------------------------------------------------------- catalogue

_SUBTYPE_SUFFIX = re.compile(r"\s*\((?:H\d+N\d+|H\d+|N\d+|mixed)\)\s*$", re.IGNORECASE)


def normalise_strain(text: str) -> str:
    """Strain/isolate text as a grouping key: subtype suffix dropped, casefolded."""
    return " ".join(_SUBTYPE_SUFFIX.sub("", text or "").split()).casefold()


def anellovirus_genus(row: Mapping[str, str], anello_genus: Mapping[str, str]) -> str:
    """The ``ANDET-05`` genus label of an anellovirus catalogue row.

    The anellovirus accession table is authoritative; the catalogue's
    lineage-derived genus is empty on 704 anellovirus rows. Max-panel genomes
    absent from the table fall back to the catalogue genus when it is an ICTV
    genus name, then the family.
    """
    label = anello_genus.get(row.get("accession_version", "")) or anello_genus.get(
        row.get("accession", "")
    )
    if label:
        return label
    genus = (row.get("genus") or "").strip()
    # A lineage token such as "Small anellovirus" is not an ICTV genus; genus
    # names are one word ending in "virus".
    return genus if re.fullmatch(r"[A-Z][a-z]+virus", genus) else ANELLOVIRIDAE


def _segment_label(row: Mapping[str, str]) -> str:
    """Segment name as a comparison key ("L RNA" and "L" are one segment)."""
    label = (row.get("segment") or "").strip().casefold()
    label = re.sub(r"^segment\s+", "", label)
    return re.sub(r"\s+rna$", "", label)


def assign_virus_keys(
    rows: Sequence[Mapping[str, str]], anello_genus: Mapping[str, str]
) -> list[str]:
    """The Virus key of each catalogue row (see the module docstring).

    Keys are assigned over the whole catalogue, not row by row, because whether
    a segmented taxid is one isolate depends on the other rows with that taxid.
    """
    keys: list[str] = []
    for row in rows:
        taxid = (row.get("taxid") or "").strip()
        if (row.get("family") or "").strip() == ANELLOVIRIDAE:
            keys.append(f"genus:{anellovirus_genus(row, anello_genus)}")
        elif taxid:
            keys.append(f"taxid:{taxid}")
        else:
            keys.append(f"accession:{(row.get('accession') or '').strip()}")

    segmented: dict[str, list[int]] = defaultdict(list)
    for i, (row, key) in enumerate(zip(rows, keys)):
        if key.startswith("taxid:") and _segment_label(row):
            segmented[key].append(i)
    for key, members in segmented.items():
        if max(Counter(_segment_label(rows[i]) for i in members).values()) > 1:
            for i in members:
                strain = rows[i].get("strain") or rows[i].get("isolate") or ""
                keys[i] = f"{key}|{normalise_strain(strain)}"
    return keys


def _display_names(members_by_key: Mapping[str, list[Mapping[str, str]]]) -> dict[str, str]:
    """One display name per key, unique across keys.

    The curated ``common_name`` wins (most frequent among the key's rows, so
    existing outputs keep "Epstein-Barr virus"); otherwise the NCBI organism of
    the key's first RefSeq row, with the strain appended when a segmented
    isolate's organism names only its species. Anellovirus keys are named by their genus. A name still shared by
    two keys is suffixed with the key so no two viruses merge in a name-keyed
    output; the packaged catalogue is tested to need no suffix.
    """
    names: dict[str, str] = {}
    for key, members in members_by_key.items():
        if key.startswith("genus:"):
            genus = key.split(":", 1)[1]
            # A bare family label beside genus rows reads as their total (Q8b = A2).
            names[key] = ANELLO_UNASSIGNED if genus == ANELLOVIRIDAE else genus
            continue
        curated = Counter(
            (r.get("common_name") or "").strip()
            for r in members
            if (r.get("common_name") or "").strip()
        )
        if curated:
            names[key] = min(curated.items(), key=lambda kv: (-kv[1], kv[0]))[0]
            continue
        refseq = [r for r in members if (r.get("refseq") or "").lower() == "true"] or list(members)
        first = min(refseq, key=lambda r: r.get("accession_version") or "")
        name = (
            (first.get("organism") or "").strip()
            or (first.get("species") or "").strip()
            or key.split(":", 1)[1]
        )
        # A segmented key is one isolate; name it so, when the organism names
        # only the species (a species-level taxid).
        strain = (first.get("strain") or first.get("isolate") or "").strip()
        if (
            any(_segment_label(r) for r in members)
            and strain
            and normalise_strain(strain) not in normalise_strain(name)
        ):
            name = f"{name} ({strain})"
        names[key] = name
    shared = {name for name, n in Counter(names.values()).items() if n > 1}
    return {key: (f"{name} [{key}]" if name in shared else name) for key, name in names.items()}


class _Catalogue:
    """Catalogue rows indexed by accession, with each row's key and display name."""

    def __init__(
        self,
        rows: Iterable[Mapping[str, str]],
        anello_genus: Mapping[str, str],
    ) -> None:
        self.anello_genus = anello_genus
        self._by_accession: dict[str, Mapping[str, str]] = {}
        self._key: dict[str, str] = {}
        members: dict[str, list[Mapping[str, str]]] = defaultdict(list)
        usable = [
            row
            for row in rows
            if (row.get("accession_version") or "").strip() or (row.get("accession") or "").strip()
        ]
        for row, key in zip(usable, assign_virus_keys(usable, anello_genus)):
            versioned = (row.get("accession_version") or "").strip()
            bare = (row.get("accession") or "").strip() or versioned.split(".")[0]
            members[key].append(row)
            for acc in (versioned, bare):
                if acc:
                    self._by_accession.setdefault(acc, row)
                    self._key.setdefault(acc, key)
        self.names = _display_names(members)

    def lookup(self, accession: str) -> tuple[Mapping[str, str], str] | None:
        for acc in (accession, accession.split(".")[0]):
            row = self._by_accession.get(acc)
            if row is not None:
                return row, self._key[acc]
        return None


def _default_anello_genus() -> dict[str, str]:
    try:
        from viralscan.anellovirus import anello_name_map

        return anello_name_map()
    except Exception as exc:  # noqa: BLE001 - the table is an enrichment
        log.warning("Anellovirus accession table unavailable (%s); using catalogue genus.", exc)
        return {}


# ---------------------------------------------------------------- build manifest

MANIFEST_SUFFIX = ".build_manifest.json"
MANIFEST_SCHEMA_VERSION = "1.0"

_VERSION_SUFFIX = re.compile(r"\.\d+$")
_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


class BuildManifestContradiction(ValueError):
    """The Run's ``--gtf`` or index contradicts the index's build manifest."""


def deversion(gene_id: str) -> str:
    """``ENSG00000123.4`` -> ``ENSG00000123`` (a trailing ``.N`` only)."""
    return _VERSION_SUFFIX.sub("", gene_id)


def manifest_path_for_index(index: PathLike) -> Path:
    """Where the build manifest of the kb index at ``index`` lives (next to it)."""
    index = Path(index)
    return index.with_name(index.name + MANIFEST_SUFFIX)


def gtf_gene_ids(path: PathLike) -> set[str]:
    """Every ``gene_id`` of a GTF (plain or ``.gz``)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    ids: set[str] = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 9:
                match = _GENE_ID_RE.search(cols[8])
                if match:
                    ids.add(match.group(1))
    return ids


@dataclass(frozen=True)
class BuildManifest:
    """The host and viral gene sets of one kb index, de-versioned."""

    host_gene_ids: frozenset[str]
    viral_gene_ids: frozenset[str]
    provenance: Mapping[str, object]
    path: Path | None = None


def write_build_manifest(
    index: PathLike,
    host_gene_ids: Iterable[str],
    viral_gene_ids: Iterable[str],
    provenance: Mapping[str, object] | None = None,
) -> Path:
    """Write ``<index>.build_manifest.json``; returns its path.

    Gene IDs are stored sorted and as given (versioned IDs stay versioned; the
    reader de-versions). A gene in both sets is a caller bug and raises.
    """
    host = sorted(set(host_gene_ids))
    viral = sorted(set(viral_gene_ids))
    both = {deversion(g) for g in host} & {deversion(g) for g in viral}
    if both:
        raise ValueError(
            f"{len(both)} gene ID(s) are listed as both host and viral in the build "
            f"manifest (e.g. {', '.join(sorted(both)[:5])})."
        )
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "viralscan_version": _VIRALSCAN_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "index": Path(index).name,
        "n_host_genes": len(host),
        "n_viral_genes": len(viral),
        "provenance": dict(provenance or {}),
        "host_gene_ids": host,
        "viral_gene_ids": viral,
    }
    out = manifest_path_for_index(index)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def write_build_manifest_from_t2g(
    index: PathLike,
    t2g_path: PathLike,
    viral_gtf_gene_ids: Iterable[str],
    provenance: Mapping[str, object] | None = None,
) -> Path:
    """Manifest of an index built from a user-supplied FASTA/GTF (``--reference``).

    The GTF does not say which genes are host, so the split follows the index's
    own t2g: a GTF gene whose column 5 names a transcript of the index is a host
    cDNA row (the structural guard), every other GTF gene is viral, and every
    indexed gene outside the GTF is host.
    """
    t2g = read_t2g(t2g_path)
    gtf_genes = set(viral_gtf_gene_ids)
    viral = {g for g in t2g.accession if g in gtf_genes and g not in t2g.structural_host}
    host = set(t2g.accession) - viral
    return write_build_manifest(index, host, viral, provenance)


def load_build_manifest(path: PathLike) -> BuildManifest:
    """Read a build manifest; raises :class:`ValueError` when it is malformed."""
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path}: unreadable index build manifest ({exc})") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{path}: not an index build manifest")
    major = str(data.get("schema_version", "")).split(".")[0]
    if major != MANIFEST_SCHEMA_VERSION.split(".")[0]:
        raise ValueError(
            f"{path}: unsupported build manifest schema_version "
            f"{data.get('schema_version')!r} (this ViralScan reads {MANIFEST_SCHEMA_VERSION})."
        )
    for key in ("host_gene_ids", "viral_gene_ids"):
        if not isinstance(data.get(key), list):
            raise ValueError(f"{path}: build manifest has no {key} list")
    return BuildManifest(
        host_gene_ids=frozenset(deversion(str(g)) for g in data["host_gene_ids"]),
        viral_gene_ids=frozenset(deversion(str(g)) for g in data["viral_gene_ids"]),
        provenance=data.get("provenance") or {},
        path=path,
    )


def _examples(ids: Iterable[str], limit: int = 8) -> str:
    ordered = sorted(ids)
    more = f" (+{len(ordered) - limit} more)" if len(ordered) > limit else ""
    return ", ".join(ordered[:limit]) + more


# ------------------------------------------------------------------------- table


@dataclass(frozen=True)
class VirusIdentityTable:
    """The Virus Identity table of one Run, one :class:`GeneIdentity` per gene."""

    genes: tuple[GeneIdentity, ...]

    @property
    def index_kind(self) -> str:
        return COMBINED if any(g.status == HOST for g in self.genes) else VIRUS_ONLY

    def by_gene(self) -> dict[str, GeneIdentity]:
        return {g.gene_id: g for g in self.genes}

    def viral_gene_ids(self) -> list[str]:
        return [g.gene_id for g in self.genes if g.viral]

    def status_counts(self) -> Counter[str]:
        return Counter(g.status for g in self.genes)

    def groups(self) -> dict[str, list[str]]:
        """virus_key → its viral gene_ids, in table order."""
        out: dict[str, list[str]] = {}
        for g in self.genes:
            if g.viral:
                out.setdefault(g.virus_key, []).append(g.gene_id)
        return out

    def virus_names(self) -> dict[str, str]:
        """virus_key → display name."""
        return {g.virus_key: g.virus_name for g in self.genes if g.viral}

    def write_tsv(self, path: PathLike) -> Path:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(COLUMNS)
            for g in self.genes:
                writer.writerow(
                    [str(getattr(g, c)).lower() if c == "viral" else getattr(g, c) for c in COLUMNS]
                )
        return out

    @classmethod
    def read_tsv(cls, path: PathLike) -> VirusIdentityTable:
        with open(path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"{path}: not a Virus Identity table (missing columns {missing})")
            genes = tuple(
                GeneIdentity(
                    **{c: (row[c] == "true") if c == "viral" else (row[c] or "") for c in COLUMNS}
                )
                for row in reader
            )
        return cls(genes)


def _catalogued(
    gene_id: str, accession: str, row: Mapping[str, str], key: str, cat: _Catalogue
) -> GeneIdentity:
    family = (row.get("family") or "").strip()
    genus = (
        anellovirus_genus(row, cat.anello_genus)
        if family == ANELLOVIRIDAE
        else (row.get("genus") or "").strip()
    )
    return GeneIdentity(
        gene_id=gene_id,
        genome_accession=accession,
        status=CATALOGUED,
        viral=True,
        virus_key=key,
        virus_name=cat.names[key],
        organism=(row.get("organism") or "").strip(),
        taxid=(row.get("taxid") or "").strip(),
        species=(row.get("species") or "").strip(),
        genus=genus,
        family=family,
        segment=(row.get("segment") or "").strip(),
        strain=(row.get("strain") or row.get("isolate") or "").strip(),
        sibling_group=(row.get("sibling_group") or "").strip(),
        risk_class=(row.get("risk_class") or "").strip(),
        role=(row.get("role") or "").strip(),
    )


def _resolve_genes(
    t2g: T2gGenes,
    gtf_genes: set[str],
    cat: _Catalogue | None,
    t2g_path: PathLike,
    *,
    warn_legacy: bool = True,
) -> list[GeneIdentity]:
    """Apply the resolution rules of the module docstring to every indexed gene."""
    genes: list[GeneIdentity] = []
    if t2g.legacy:
        from viralscan.virus_grouping import virus_name_for_gene

        name_map = virus_catalog.merged_name_map()
        for gene_id in t2g.accession:
            if gene_id in gtf_genes:
                name = virus_name_for_gene(gene_id, name_map)
                genes.append(GeneIdentity(gene_id, "", LEGACY_PREFIX, True, f"name:{name}", name))
            else:
                genes.append(GeneIdentity(gene_id, "", HOST, False))
        if warn_legacy:
            log.warning(
                "t2g %s has no genome-accession column (pre-v3 index): viral genes come "
                "from the GTF gene set and are named by legacy prefix maps.",
                t2g_path,
            )
    else:
        assert cat is not None
        for gene_id, accession in t2g.accession.items():
            hit = cat.lookup(accession) if accession else None
            if hit is not None:
                genes.append(_catalogued(gene_id, accession, hit[0], hit[1], cat))
            elif gene_id in gtf_genes and gene_id not in t2g.structural_host:
                key = f"accession:{accession or gene_id}"
                genes.append(
                    GeneIdentity(gene_id, accession, UNCATALOGUED, True, key, accession or gene_id)
                )
            else:
                genes.append(GeneIdentity(gene_id, accession, HOST, False))
    return genes


def _manifest_gene_set(
    manifest: BuildManifest,
    t2g: T2gGenes,
    gtf_genes: set[str],
    cat: _Catalogue | None,
    t2g_path: PathLike,
) -> set[str]:
    """The viral gene set the manifest dictates, after checking it against the Run.

    Raises :class:`BuildManifestContradiction` when the viral genes ``--gtf``
    would yield (de-versioned, among the indexed genes) differ from the
    manifest's, or when the index holds genes the manifest does not list.
    """
    where = manifest.path or "the index build manifest"
    indexed = {deversion(g) for g in t2g.accession}
    unlisted = indexed - manifest.viral_gene_ids - manifest.host_gene_ids
    if unlisted:
        raise BuildManifestContradiction(
            f"The index {t2g_path} holds {len(unlisted)} gene(s) that the build manifest "
            f"{where} lists as neither host nor viral (e.g. {_examples(unlisted)}). "
            "The manifest does not describe this index; rebuild the index or remove the manifest."
        )
    # Compare de-versioned: a GTF written for another release names the same genes.
    gtf_bare = {deversion(g) for g in gtf_genes}
    gtf_in_index = {g for g in t2g.accession if deversion(g) in gtf_bare}
    from_gtf = {
        deversion(g.gene_id)
        for g in _resolve_genes(t2g, gtf_in_index, cat, t2g_path, warn_legacy=False)
        if g.viral
    }
    from_manifest = manifest.viral_gene_ids & indexed
    if from_gtf != from_manifest:
        gtf_only = from_gtf - from_manifest
        manifest_only = from_manifest - from_gtf
        parts = []
        if gtf_only:
            parts.append(
                f"{len(gtf_only)} gene(s) that --gtf marks viral are not viral in the "
                f"manifest (e.g. {_examples(gtf_only)})"
            )
        if manifest_only:
            parts.append(
                f"{len(manifest_only)} viral gene(s) of the manifest are not viral under "
                f"--gtf (e.g. {_examples(manifest_only)})"
            )
        raise BuildManifestContradiction(
            f"--gtf contradicts the build manifest {where} of the index {t2g_path}: "
            + "; ".join(parts)
            + ". Pass the GTF the index was built from, or rebuild the index."
        )
    return {g for g in t2g.accession if deversion(g) in manifest.viral_gene_ids}


def build_identity_table(
    t2g_path: PathLike,
    gtf_gene_ids: Iterable[str],
    *,
    catalogue_rows: Iterable[Mapping[str, str]] | None = None,
    anello_genus: Mapping[str, str] | None = None,
    build_manifest: Path | None = None,
) -> VirusIdentityTable:
    """Resolve every gene of the index at ``t2g_path`` (rules in the module docstring).

    ``gtf_gene_ids`` is the GTF gene set ``analysis.obtain_gtf`` collects.
    ``catalogue_rows`` defaults to the packaged catalogue, ``anello_genus`` to
    :func:`viralscan.anellovirus.anello_name_map`.

    ``build_manifest`` is the index's :func:`write_build_manifest` file. When
    given, its viral gene set replaces ``gtf_gene_ids`` in rule 2, after a check
    that the two agree (see :func:`_manifest_gene_set`); ``None`` keeps the
    ``gtf_gene_ids`` behaviour exactly.

    Raises :class:`ValueError` when no gene resolves as viral: every later step
    would then report a clean negative for the wrong reason, and
    :class:`BuildManifestContradiction` (a ``ValueError``) when the manifest and
    the Run disagree.
    """
    t2g = read_t2g(t2g_path)
    gtf_genes = set(gtf_gene_ids)
    cat = None
    if not t2g.legacy:
        rows = virus_catalog.load_catalogue() if catalogue_rows is None else catalogue_rows
        cat = _Catalogue(rows, _default_anello_genus() if anello_genus is None else anello_genus)
    if build_manifest is not None:
        manifest = load_build_manifest(build_manifest)
        gtf_genes = _manifest_gene_set(manifest, t2g, gtf_genes, cat, t2g_path)
        log.info("Using index build manifest %s (its gene sets override --gtf).", build_manifest)
    genes = _resolve_genes(t2g, gtf_genes, cat, t2g_path)

    table = VirusIdentityTable(tuple(genes))
    counts = table.status_counts()
    if not table.viral_gene_ids():
        guarded = len(gtf_genes & t2g.structural_host)
        detail = (
            f" {guarded} GTF gene(s) were treated as host because their t2g column 5 "
            "names a transcript, not a genome (a host cDNA row)."
            if guarded
            else ""
        )
        raise ValueError(
            f"No viral gene in the index {t2g_path}: none of its {len(genes)} genes has a "
            f"catalogued genome accession or appears in the viral GTF gene set.{detail}"
        )
    if counts[UNCATALOGUED]:
        examples = sorted({g.genome_accession for g in genes if g.status == UNCATALOGUED})
        log.warning(
            "%d viral gene(s) on %d genome(s) are not in the virus catalogue and are "
            "named by accession (e.g. %s). Add them with extras/build_virus_catalog.py.",
            counts[UNCATALOGUED],
            len(examples),
            ", ".join(examples[:5]),
        )
    log.info(
        "Virus Identity table: %s index, %s",
        table.index_kind,
        ", ".join(f"{k}={v}" for k, v in sorted(counts.items())),
    )
    return table


def write_identity_table(config: Any, gtf_gene_ids: set[str]) -> Path | None:
    """Build the Run's Virus Identity table and write ``results/virus_identity.tsv``.

    Resolution runs against the index t2g (``config.transcripts``). A Snakemake
    run always has one (``kb count`` needs it), so the skip below only serves a
    direct call without an index. Raises :class:`ValueError` when the index has
    no viral gene, and :class:`~viralscan.virus_identity.BuildManifestContradiction`
    when the index's build manifest (``<index>.build_manifest.json``) disagrees
    with ``--gtf``.

    ``config`` is a :class:`~viralscan.runconfig.RunConfig` (``transcripts``,
    ``index`` and ``output`` are read).
    """
    t2g = (config.transcripts or "").strip()
    if not t2g or not Path(t2g).exists():
        log.warning("No t2g file at %r; %s was not written.", t2g, TABLE_FILENAME)
        return None
    manifest = manifest_path_for_index(config.index) if (config.index or "").strip() else None
    if manifest is not None and manifest.is_file():
        table = build_identity_table(t2g, gtf_gene_ids, build_manifest=manifest)
    else:
        log.warning(
            "No index build manifest at %s: host and viral genes come from the --gtf "
            "gene set (the index was built outside `viralscan build-ref`/`--reference`).",
            manifest or "(no index path)",
        )
        table = build_identity_table(t2g, gtf_gene_ids)
    return table.write_tsv(Path(config.output) / "results" / TABLE_FILENAME)


def backfill_identity_table(config: Any) -> Path | None:
    """Build ``results/virus_identity.tsv`` for a Run made before the table existed.

    Runs made before 35940ec have ``log/analysis.txt`` but no table. Once the
    consumer rules declare the table as an input, snakemake would rerun
    ``analysis`` and everything after it (PLAN MECH-A step 4c). The table is
    built from that Run's own ``analysis.txt`` gene set, not a fresh GTF glob,
    so catalogue or GTF changes since the Run cannot shift its viral/host
    partition. The file takes ``analysis.txt``'s mtime, so snakemake sees the
    ``analysis`` step as complete and nothing downstream turns stale.

    Returns the written path, or ``None`` when a table already exists or the
    Run has no ``analysis.txt``.
    """
    out = Path(config.output)
    gene_list = out / "log" / "analysis.txt"
    if (out / "results" / TABLE_FILENAME).exists() or not gene_list.is_file():
        return None
    gene_ids = {
        line.strip() for line in gene_list.read_text(encoding="utf-8").splitlines() if line.strip()
    }
    written = write_identity_table(config, gene_ids)
    if written is not None:
        stat = gene_list.stat()
        os.utime(written, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        log.info("Backfilled %s from %s", written, gene_list)
    return written
