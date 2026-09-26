"""Named HPV ORFs: catalogue loader, gene classes, and panel-ID resolution.

Why a catalogue rather than a coordinate table
---------------------------------------------
The reference ViralScan ships today represents HPV with four accessions whose
gene IDs are ``HpV16gp1`` ... ``HpV16gp8`` — RefSeq ``locus_tag`` values. Nothing
in the kallisto ``t2g`` says which of those is E6 and which is L1, so an
oncogene-versus-capsid contrast is not expressible against that index. The
bundled RefSeq GTFs do carry the answer in a ``gene`` attribute, but the packager
that produced the index kept the ``locus_tag`` as the ID and dropped ``gene``.

The obvious workaround — map ORFs to names by genomic coordinate against a
published HPV16 coordinate table — is exactly the wrong thing to do here, for
two reasons that are measurable rather than stylistic:

1. **The coordinate table would have to be trusted blindly anyway.** Papillomavirus
   genomes are submitted as linearised circles with the linearisation point
   chosen by the submitter, inside E1. HPV16 E6 therefore sits at 7125-7601 in
   ``NC_001526.4`` and at 105-581 in ``NC_001357.1`` (HPV18) — the same ORF, at
   opposite ends of the record, because the two records start at different points
   in their circles. Any single hard-coded coordinate map silently mis-assigns
   every type but the one it was copied from.
2. **It is unnecessary.** Every HPV complete-genome record examined carries the
   name as a feature qualifier: ``/gene="E6"`` in RefSeq and in most INSDC
   submissions, ``/product="transforming protein E6"`` in the rest. The names are
   in the record; a coordinate table is a second, weaker source for something
   the first source already states.

So this catalogue stores *the record's own* annotation, with the qualifier each
name came from recorded in ``name_source`` so a reader can audit it.

What the catalogue can and cannot support
-----------------------------------------
It supports an **oncogene-versus-capsid** contrast per genotype: E6/E7 are
early-region oncoproteins transcribed from the same locus in every
carcinogen-driven HPV-positive oropharyngeal tumour, whereas L1/L2 are
late-region capsid genes transcribed only in productive infection. E6/E7
transcripts therefore report viral gene expression where L1 reports virion
production.

It does **not** support confident per-genotype attribution of L1 signal. L1 is the
most conserved coding region in the papillomavirus genus — it is the region
pan-HPV PCR primers are designed against — so L1 reads cross-map freely between
the 16 genotypes here, and kallisto's multimapping will distribute them. Genotype
labels on L1 counts are not independent evidence of which type is present. E6/E7
are far more type-divergent and so are better behaved, but the count of a
cross-mapped transcript is not the count of a transcript. See the ``note`` column
and ``PLAN.md`` before quoting a per-type number.

Scope: this is a *transcriptomic* catalogue. A transcriptionally silent
integrated HPV genome — the common state in tonsillar crypt epithelium, and the
state that drives HPV-positive oropharyngeal carcinoma — produces no reads at
all and is invisible here, not "negative".
"""

from __future__ import annotations

import csv
import importlib.resources
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

#: Functional class of a canonical HPV open reading frame.
#:
#: ``oncogene``        E6 and E7 themselves — the viral oncoproteins.
#: ``oncogene_locus``  E6* and E7*, transcribed from the E6/E7 locus but *not*
#:                    encoding the oncoprotein: E6* is a truncated E6 and E7*
#:                    lacks the LXCXE pRb-binding motif, so neither is
#:                    transforming. Kept separate so that E6* expression is never
#:                    silently counted as E6 activity.
#: ``early``           the other early-region ORFs (E1, E2, E4, E5, the E1^E4
#:                    fusion, and the E8/E9/E10/E11 family). E5 is deliberately
#:                    *not* an oncogene here: it is a transforming protein, but
#:                    it is not part of the E6/E7 axis this catalogue exists to
#:                    measure, and lumping it in would overstate what a
#:                    positive call means.
#: ``late_capsid``     L1 and L2, the late-region structural genes.
GENE_CLASSES = frozenset({"oncogene", "oncogene_locus", "early", "late_capsid", "other"})

#: The viral oncoproteins: the E6/E7 locus.
ONCOGENE_GENES = frozenset({"E6", "E6*", "E7", "E7*"})

#: Genes whose transcripts come from the E6/E7 locus without encoding the
#: oncoprotein. Reported separately from :data:`ONCOGENE_GENES` on purpose.
ONCOGENE_LOCUS_GENES = frozenset({"E6*", "E7*"})

#: Late-region structural genes. Present in productive infection; absent from a
#: silently integrated genome, which is why their absence is not evidence of
#: absence of virus.
CAPSID_GENES = frozenset({"L1", "L2"})

#: The 14 oncogenic (high-risk) genotypes: those HPV types that drive
#: HPV-associated oropharyngeal squamous cell carcinoma.
HIGH_RISK_GENOTYPES = (
    "16",
    "18",
    "31",
    "33",
    "35",
    "39",
    "45",
    "51",
    "52",
    "56",
    "58",
    "59",
    "66",
    "69",
)

#: Columns the loader requires to be present and non-empty for every row.
REQUIRED_COLUMNS = (
    "accession",
    "accession_version",
    "genotype",
    "high_risk",
    "genome_length",
    "topology",
    "gene_id",
    "canonical_gene_name",
    "gene_class",
    "genome_start",
    "genome_end",
    "strand",
    "n_exons",
    "exon_blocks",
    "spans_origin",
    "cds_length_nt",
    "name_source",
    "product_as_in_ncbi",
    "source",
)

_TSV_FILENAME = "hpv_genes.tsv"
_PACKAGE_DATA = "viralscan.data"

Row = dict[str, Any]


def _default_tsv_path() -> Path:
    """Return the packaged catalogue via importlib.resources.

    Mirrors ``gene_programs._default_tsv_path`` so every packaged table resolves
    the same way regardless of how the package was installed.
    """
    try:
        ref = importlib.resources.files(_PACKAGE_DATA).joinpath(_TSV_FILENAME)
        return Path(str(ref))
    except AttributeError:  # pragma: no cover - Python 3.8 fallback
        import importlib.resources as pkg_resources

        with pkg_resources.path(_PACKAGE_DATA, _TSV_FILENAME) as path:
            return Path(path)


def load_catalogue(path: str | os.PathLike[str] | None = None) -> list[Row]:
    """Load the HPV ORF catalogue.

    Returns one dict per ORF, with every TSV column as a string. The boolean-ish
    columns (``high_risk``, ``spans_origin``) are kept as the strings
    ``"true"``/``"false"`` rather than coerced, so a malformed value is visible in
    the raw row instead of silently becoming ``False``.
    """
    tsv_path = Path(path) if path is not None else _default_tsv_path()
    rows: list[Row] = []
    with open(tsv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(dict(row))
    return rows


def validate_catalogue(rows: list[Row]) -> list[str]:
    """Return a list of problems with ``rows``; empty means valid.

    Checks the controlled vocabularies, the boolean columns, integer columns,
    coordinate ordering, and that no ``(accession, canonical_gene_name)`` pair
    repeats. Deliberately does not check the sequences: the FASTA those
    coordinates refer to is not packaged, and that check belongs to the build
    script and to the tests that run against the NCBI cache.
    """
    errors: list[str] = []
    if not rows:
        return ["HPV gene catalogue is empty"]
    seen: set[tuple[str, str]] = set()
    for index, row in enumerate(rows, start=2):
        where = f"row {index}"
        blank = [
            column
            for column in REQUIRED_COLUMNS
            if column not in row or not str(row.get(column, "")).strip()
        ]
        if blank:
            errors.append(f"{where}: empty or missing column(s) {blank}")
            continue
        if row.get("gene_class") not in GENE_CLASSES:
            errors.append(
                f"{where}: gene_class {row.get('gene_class')!r} not in {sorted(GENE_CLASSES)}"
            )
        for column in ("high_risk", "spans_origin"):
            if row.get(column) not in {"true", "false"}:
                errors.append(f"{where}: {column}={row.get(column)!r} is not true/false")
        for column in ("genome_length", "genome_start", "genome_end", "n_exons", "cds_length_nt"):
            try:
                int(str(row.get(column, "")))
            except ValueError:
                errors.append(f"{where}: {column}={row.get(column)!r} is not an integer")
        try:
            start = int(str(row.get("genome_start")))
            end = int(str(row.get("genome_end")))
            length = int(str(row.get("genome_length")))
            if start < 1 or end < start:
                errors.append(f"{where}: coordinates {start}..{end} are not a valid span")
            elif end > length:
                errors.append(f"{where}: end {end} exceeds genome_length {length}")
        except ValueError:
            pass
        if row.get("strand") not in {"+", "-"}:
            errors.append(f"{where}: strand {row.get('strand')!r} is not + or -")
        if row.get("topology") not in {"circular", "linear"}:
            errors.append(f"{where}: topology {row.get('topology')!r} is invalid")
        if row.get("name_source") not in {"gene", "product"}:
            errors.append(f"{where}: name_source {row.get('name_source')!r} is not gene or product")
        if row.get("source") not in {"refseq", "insdc"}:
            errors.append(f"{where}: source {row.get('source')!r} is invalid")
        key = (str(row.get("accession_version", "")), str(row.get("canonical_gene_name", "")))
        if key in seen:
            errors.append(f"{where}: duplicate (accession, gene) {key}")
        seen.add(key)
    return errors


@dataclass(frozen=True)
class HpvOrf:
    """One catalogued HPV open reading frame."""

    gene_id: str
    genotype: str
    canonical_gene_name: str
    gene_class: str
    genome_start: int
    genome_end: int
    strand: str
    n_exons: int
    high_risk: bool

    @property
    def is_oncogene(self) -> bool:
        """True for E6/E7 themselves, false for E6*/E7*."""
        return self.gene_class == "oncogene"

    @property
    def is_capsid(self) -> bool:
        return self.gene_class == "late_capsid"

    @property
    def is_from_oncogene_locus(self) -> bool:
        return self.gene_class in {"oncogene", "oncogene_locus"}


def orfs(rows: list[Row]) -> list[HpvOrf]:
    """Coerce catalogue rows to typed :class:`HpvOrf` records."""
    return [
        HpvOrf(
            gene_id=row["gene_id"],
            genotype=row["genotype"],
            canonical_gene_name=row["canonical_gene_name"],
            gene_class=row["gene_class"],
            genome_start=int(row["genome_start"]),
            genome_end=int(row["genome_end"]),
            strand=row["strand"],
            n_exons=int(row["n_exons"]),
            high_risk=row["high_risk"] == "true",
        )
        for row in rows
    ]


def genotypes(rows: list[Row], *, high_risk_only: bool = False) -> list[str]:
    """Genotypes present in ``rows``, sorted, deduplicated."""
    found = {
        row["genotype"]
        for row in rows
        if row["genotype"] and not (high_risk_only and row["high_risk"] != "true")
    }
    return sorted(found, key=_genotype_sort_key)


def _genotype_sort_key(genotype: str) -> tuple[int, str]:
    match = re.match(r"^(\d+)([A-Za-z]*)$", genotype)
    if not match:
        return (10**6, genotype)
    return (int(match.group(1)), match.group(2))


def _catalogued_genotypes(rows: list[Row]) -> list[str]:
    """:func:`genotypes` under a name that survives a parameter of the same name.

    Several functions here take a ``genotypes`` filter argument, which shadows the
    module-level function inside their own body.
    """
    return genotypes(rows)


def gene_ids_for(
    rows: list[Row],
    gene_class: str,
    genotypes: Iterable[str] | None = None,
) -> list[str]:
    """Panel gene IDs whose ORFs fall in ``gene_class``.

    With ``genotypes`` given, only those genotypes are considered; without it,
    every genotype in the catalogue is. A caller comparing an E6/E7 signal against
    an L1 signal must pass the *same* genotype subset to both, otherwise the two
    sides are summed over different genomes and the ratio is meaningless.
    """
    wanted = None if genotypes is None else set(genotypes)
    return [
        row["gene_id"]
        for row in rows
        if row["gene_class"] == gene_class and (wanted is None or row["genotype"] in wanted)
    ]


def oncogene_gene_ids(rows: list[Row], genotypes: Iterable[str] | None = None) -> list[str]:
    """Gene IDs for E6 and E7 — the transcripts that report oncogene expression.

    E6* and E7* are excluded on purpose. They come from the same locus, so their
    reads are indistinguishable from E6/E7 at the sequence level, and reporting
    them as separate oncoprotein evidence would double-count one locus. Use
    :func:`oncogene_locus_gene_ids` if the wider locus is wanted.
    """
    return [
        gene_id
        for gene_id in [
            row["gene_id"]
            for row in rows
            if row["gene_class"] == "oncogene"
            and (genotypes is None or row["genotype"] in set(genotypes))
        ]
    ]


def oncogene_locus_gene_ids(rows: list[Row], genotypes: Iterable[str] | None = None) -> list[str]:
    """Gene IDs for every transcript of the E6/E7 locus, E6* and E7* included."""
    wanted = None if genotypes is None else set(genotypes)
    return [
        row["gene_id"]
        for row in rows
        if row["gene_class"] in {"oncogene", "oncogene_locus"}
        and (wanted is None or row["genotype"] in wanted)
    ]


def capsid_gene_ids(rows: list[Row], genotypes: Iterable[str] | None = None) -> list[str]:
    """Gene IDs for the late-region capsid genes L1 and L2."""
    return gene_ids_for(rows, "late_capsid", genotypes)


def early_gene_ids(rows: list[Row], genotypes: Iterable[str] | None = None) -> list[str]:
    """Gene IDs for the early-region ORFs other than E6/E7."""
    return gene_ids_for(rows, "early", genotypes)


def split_oncogene_vs_capsid(
    counts: dict[str, float],
    rows: list[Row],
    genotypes: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Summarise ``counts`` as oncogene-locus versus capsid signal.

    ``counts`` maps panel gene ID to a count, so this is usable on any count
    matrix whose columns the caller has already resolved. Returns total counts,
    the number of ORFs contributing on each side, and the ratio **only** when both
    sides are non-zero — a one-sided result is reported as ``None`` rather than as
    an infinite or zero ratio, because "capsid but no oncogene transcript" is a
    biologically ordinary state (productive infection) and not a ratio of 0.

    The ratio is an ORF-length-normalised-agnostic bulk count and should be read
    as *which region the transcripts came from*, not as a fold-change in
    oncogene activity. It inherits every cross-mapping caveat in the module
    docstring, and L1 cross-mapping across genotypes is the one that bites.
    """
    wanted = None if genotypes is None else set(genotypes)
    oncogene_total = 0.0
    capsid_total = 0.0
    oncogene_n = 0
    capsid_n = 0
    for row in rows:
        if wanted is not None and row["genotype"] not in wanted:
            continue
        value = float(counts.get(row["gene_id"], 0.0) or 0.0)
        if row["gene_class"] in {"oncogene", "oncogene_locus"}:
            oncogene_total += value
            oncogene_n += 1
        elif row["gene_class"] == "late_capsid":
            capsid_total += value
            capsid_n += 1
    ratio: float | None = None
    if oncogene_total > 0 and capsid_total > 0:
        ratio = oncogene_total / capsid_total
    return {
        "oncogene_locus_total": oncogene_total,
        "capsid_total": capsid_total,
        "oncogene_locus_orfs": oncogene_n,
        "capsid_orfs": capsid_n,
        "oncogene_to_capsid_ratio": ratio,
        "genotypes": sorted(wanted) if wanted is not None else _catalogued_genotypes(rows),
    }
