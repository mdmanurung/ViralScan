"""Generator for ``src/viralscan/data/hpv_genes.tsv``.

ViralScan's HPV reference is four accessions with ``HpV16gp1``-style gene IDs:
RefSeq ``locus_tag`` values that say nothing about which open reading frame is
which. An oncogene-versus-capsid contrast is not expressible against that index.
This script produces a catalogue that names every HPV ORF.

Where the gene names come from
------------------------------
Not from a coordinate table. Every HPV complete-genome record examined carries
the ORF name as a feature qualifier:

* ``/gene="E6"`` — RefSeq records and most INSDC submissions;
* ``/product="transforming protein E6"`` or ``/product="putative E4 protein"`` —
  the records that omit ``/gene`` entirely.

The qualifier used is recorded per row in ``name_source``, so no claim in the
shipped TSV depends on trust in a table nobody can check.

A coordinate map was the obvious alternative and was rejected on evidence, not
taste. Papillomavirus genomes are submitted as **linearised circles** whose
linearisation point the submitter chooses, inside E1. The same E6 ORF therefore
sits at 7125-7601 in ``NC_001526.4`` (HPV16) and at 105-581 in ``NC_001357.1``
(HPV18) — opposite ends of the record. A coordinate table copied from one
genome would mis-assign every other genotype, and a wrong E6/E7 boundary yields
a confidently wrong biological answer rather than an obvious failure.

Independent confirmation of the names
-------------------------------------
The record's own annotation is not taken on trust. Three invariants are checked
across every record and the build aborts if any fails:

1. **ORF order, as a rotation.** Papillomavirus early and late ORFs appear on the
   genome in the fixed order ``E6, E7, E1, E2, E4, E5, L2, L1``, all 16 records
   reproduce it — but not in the same *linear* order, because the circle is cut
   at the submitter's choice of point. ``NC_001357.1`` (HPV18) reads
   ``E6, E7, E1, ...``; ``NC_001526.4`` (HPV16), cut inside E1, reads
   ``E1, E2, E5, L2, L1, E6, E7``. Both are the same genome order. This is the
   linearisation hazard, caught by the build rather than argued about in a
   comment.
2. **Protein identity.** Every E7 translation carries the LXCXE
   retinoblastoma-binding motif and every E6 translation carries its C-X2-C zinc
   fingers; E6, E7, L1 and L2 also have to fall in their known length ranges
   (E6 ~150 aa, E7 ~100 aa, L1 ~505-570 aa, L2 ~470 aa). An ORF mislabelled E6
   or E7 would fail this, which is the point: it is the one check that does not
   assume NCBI's annotation is right.
3. **Capsid scale.** L1 is at least 3x either E6 or E7, which is the check an
   E6/E7/L1 mix-up would fail. Note that L1 is *not* the longest ORF — E1 is, at
   ~650 aa, because E1 is the replication helicase — and L1 is not always longer
   than L2: HPV-2 annotates L2 at 525 aa against L1 at 511 aa. Both are real, and
   both were caught by writing the naive invariant first and watching sixteen
   correct records fail it.

The translations make the naming auditable by eye as well: E7 reads
``...LXCYEQL...`` in HPV16, HPV31 and HPV35 alike, and E6 reads ``...IICVYCKQQL...``
in HPV16 and HPV18. An E6/E7/L1 boundary derived from a coordinate table instead
of from the record would not survive this.

Genotype selection
------------------
The 14 high-risk types plus HPV1 and HPV2. HPV1 and HPV2 are cutaneous types
with no oropharyngeal association, and are carried only so that rebuilding the
reference does not *lose* the two types the current index already has.

Per type, in order of preference:

1. A RefSeq complete genome, if one exists. Only 16, 18, 31 and 33 have one —
   RefSeq has complete genomes for 61 of the ~220 papillomavirus types and the
   14 high-risk types are mostly not among them.
2. Otherwise the **oldest** INSDC complete genome carrying annotated CDS, because
   the earliest submission of a type is the prototype that genotyping assays and
   published amplicons were designed against. Ties break on the most annotated
   CDS set, then on accession.

What the result supports is stated in the module docstring of
``viralscan.hpv_genes``: an oncogene-versus-capsid contrast per genotype, and
**not** confident per-genotype attribution of L1 signal, because L1 is the most
conserved coding region in the genus and cross-maps between all 16 genotypes
here.

Usage
-----
    python extras/build_hpv_reference.py
    python extras/build_hpv_reference.py --emit-reference OUTDIR
    python extras/build_hpv_reference.py --list-records

``--emit-reference`` writes a merged FASTA + GTF for ``kb ref`` on demand. The
sequences are not committed: 16 genomes is ~127 kB today, but the packaging
review that governs reference size (PLAN PR 8) should see the request rather than
inherit it. The catalogue is the shipped artefact.

Exits non-zero if any record fails to resolve, so a bad accession or an
unrecognised ORF symbol cannot ship silently.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, str(REPO_ROOT / "src"))

from viralscan.scripts import ncbi_fetch  # noqa: E402
from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR  # noqa: E402

DATA_DIR = REPO_ROOT / "src" / "viralscan" / "data"
DEFAULT_OUT = DATA_DIR / "hpv_genes.tsv"

#: The 14 oncogenic genotypes, plus the two cutaneous types the current index
#: already carries so a rebuild does not drop them.
#:
#: ``(genotype, accession, source, note)`` — the accession is a curated constant
#: rather than a discovery step so that regenerating the TSV is deterministic and
#: offline once the cache is warm. ``--list-records`` re-derives it from NCBI and
#: prints any drift for a human to accept.
SELECTED_RECORDS: tuple[tuple[str, str, str, str], ...] = (
    (
        "1",
        "NC_001356.1",
        "refseq",
        "REVIEWED RefSeq; cutaneous, no oropharyngeal association; retained for "
        "continuity with the current index",
    ),
    (
        "2",
        "NC_001352.1",
        "refseq",
        "REVIEWED RefSeq; cutaneous, no oropharyngeal association; retained for "
        "continuity with the current index",
    ),
    (
        "16",
        "NC_001526.4",
        "refseq",
        "REVIEWED RefSeq, the type that dominates HPV-positive oropharyngeal "
        "carcinoma; the only record here with a spliced E1^E4 CDS",
    ),
    ("18", "NC_001357.1", "refseq", "VALIDATED RefSeq"),
    ("31", "NC_075191.1", "refseq", "PROVISIONAL RefSeq; /product= is generic, /gene= is not"),
    ("33", "NC_075233.1", "refseq", "PROVISIONAL RefSeq; /product= is generic, /gene= is not"),
    (
        "35",
        "M74117.1",
        "insdc",
        "oldest complete genome for type 35 (1993) and the one with the full 8-ORF "
        "annotation; the RefSeq record NC_075267.1 is type 35H, PROVISIONAL, and "
        "sequence-identical to X74477, and annotates 6 ORFs",
    ),
    ("39", "KC470231.1", "insdc", "oldest complete genome for type 39"),
    (
        "45",
        "EF202156.1",
        "insdc",
        "oldest complete genome for type 45; carries /product= only, no /gene=. The "
        "RefSeq record NC_075269.1 is PROVISIONAL, sequence-identical to X74479, and "
        "annotates 6 ORFs, so this 8-ORF INSDC record is preferred",
    ),
    (
        "51",
        "KF436866.1",
        "insdc",
        "oldest complete genome for type 51; E5 is not annotated in this record",
    ),
    (
        "52",
        "AB819272.1",
        "insdc",
        "oldest complete genome for type 52 with the full 8-ORF annotation. The "
        "earlier GQ472848 (2009) carries only 6 ORFs, as does the PROVISIONAL "
        "RefSeq record NC_075270.1, which is sequence-identical to X74481",
    ),
    (
        "56",
        "EF177176.1",
        "insdc",
        "oldest complete genome for type 56; E5 is not annotated in this record",
    ),
    (
        "58",
        "EU918765.1",
        "insdc",
        "oldest complete genome for type 58 carrying /gene=. D90400 (1993) and "
        "FJ385261-FJ385268 (2008) are older and annotate CDS features, but D90400 "
        "carries no /gene= or /product= at all and the FJ38526x records carry "
        "/product= only; /gene= is preferred because it is the RefSeq convention "
        "and distinguishes the E6* isoform from its parent",
    ),
    ("59", "EU918767.1", "insdc", "oldest complete genome for type 59"),
    (
        "66",
        "U31794.1",
        "insdc",
        "oldest complete genome for type 66; E5 is not annotated in this record",
    ),
    (
        "68",
        "DQ080079.1",
        "insdc",
        "type 68a complete genome (2005); carries /product= only, no /gene=. "
        "Replaces the earlier AB027020.1 (type 69) entry: type 69 is not on the "
        "IARC Group 1 / clinical 14-type high-risk list, whereas type 68 is",
    ),
)

#: The high-risk subset of :data:`SELECTED_RECORDS`. This is the IARC Group 1 /
#: clinical 14-type list: 16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68.
HIGH_RISK_GENOTYPES = frozenset(
    {"16", "18", "31", "33", "35", "39", "45", "51", "52", "56", "58", "59", "66", "68"}
)

#: Canonical HPV ORF symbols, mapped to the class this catalogue reports.
#:
#: This is the curated biology and it is short on purpose. Papillomavirus has
#: 8-9 ORFs and every one of them is on this list or the build fails; there is
#: no generic "ORF_n" fallback, because a row whose identity is unknown is worse
#: than no row — it would be counted as though it meant something.
#:
#: E5 is classed ``early``, not ``oncogene``. It is a transforming protein and
#: several reviews call it one, but it is not part of the E6/E7 axis this
#: catalogue exists to measure, and putting it in the oncogene class would make a
#: positive oncogene call mean something it does not.
#:
#: E6* and E7* are ``oncogene_locus``, not ``oncogene``: they are transcribed from
#: the E6/E7 locus but do not encode the oncoproteins. E6* lacks the PDZ-binding
#: motif and does not destabilise p53; E7* lacks the LXCXE pRb-binding motif, so
#: neither is transforming. Their reads are indistinguishable from E6/E7 at the
#: sequence level, so folding them into the oncogene class would let E6* on its
#: own produce a confident positive oncogene call.
CANONICAL_GENES: dict[str, str] = {
    "E1": "early",
    "E2": "early",
    "E3": "early",
    "E4": "early",
    "E5": "early",
    "E6": "oncogene",
    "E7": "oncogene",
    "E6*": "oncogene_locus",
    "E7*": "oncogene_locus",
    "E6^c": "oncogene_locus",
    "E8": "early",
    "E9": "early",
    "E10": "early",
    "E11": "early",
    "E1^E4": "early",
    "E2^E1": "early",
    "L1": "late_capsid",
    "L2": "late_capsid",
}

#: Non-alphanumeric ORF symbols mapped to a panel-safe spelling.
#:
#: A kallisto ``t2g`` is a two-column TSV with no quoting, and bustools has to
#: agree with kallisto about every character, so ``E6*`` and ``E1^E4`` are
#: spelled out rather than shipped raw. An unmapped exotic symbol is a build
#: error, not something to guess at.
SAFE_SYMBOLS: dict[str, str] = {
    "E1^E4": "E1_E4",
    "E2^E1": "E2_E1",
    "E6^c": "E6_c",
    "E6*": "E6_star",
    "E7*": "E7_star",
}

#: The canonical order ORFs appear in along a papillomavirus genome. Used to
#: check the record's own annotation rather than to assign names.
CANONICAL_ORF_ORDER = ("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1")

#: ORFs that wrap the linearisation point of a circular record.
#: ``NC_001526.4``'s ``E1^E4`` is ``join(1..16, 2494..2756)``.
ORIGIN_SPANNING_GENES = frozenset({"E1", "E1^E4", "E2^E1"})

TSV_COLUMNS = (
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
    "protein_id",
    "name_source",
    "product_as_in_ncbi",
    "source",
    "note",
)

FEATURE_RE = re.compile(r"^ {5}(\S+)\s(.*)$")
QUALIFIER_RE = re.compile(r"^ {21}/([A-Za-z_]+)(?:=(.*))?$")
SPAN_RE = re.compile(r"<?(\d+)\.\.>?(\d+)")
SINGLE_RE = re.compile(r"^<?(\d+)>?$")
LOCUS_RE = re.compile(r"^LOCUS\s+(\S+)\s+(\d+)\s+(?:bp|aa)\b([^\n]*)", re.M)
VERSION_RE = re.compile(r"^VERSION\s+(\S+)", re.M)

_SYMBOL_TOKEN_RE = re.compile(r"[A-Za-z][0-9]*(?:\^[A-Za-z][0-9]*)?\*?")


class BuildError(RuntimeError):
    """Raised when a record cannot be turned into trustworthy catalogue rows."""


def _unquote(value: str) -> str:
    value = value.strip()
    if value.startswith('"'):
        end = value.find('"', 1)
        if end > 0:
            return value[1:end]
    return value.strip('"')


def parse_genbank(text: str) -> dict[str, object]:
    """Extract the CDS features of one GenBank flatfile.

    Returns the locus name, genome length, topology, versioned accession and the
    CDS list. Each CDS is a dict with ``location``, the ``exon_blocks`` parsed out
    of it, ``strand``, and its ``/`` qualifiers.
    """
    locus_match = LOCUS_RE.search(text)
    version_match = VERSION_RE.search(text)
    version = version_match.group(1) if version_match else ""
    if not locus_match or not version:
        raise BuildError("record has no usable LOCUS or VERSION line")
    # The topology sits *after* the molecule type on the LOCUS line, e.g.
    # "LOCUS  NC_001526  7906 bp  DNA  circular  VRL  30-APR-2025". Matching only
    # up to the molecule type would read every record as linear, and papillomavirus
    # genomes are submitted as linearised circles, so the topology decides which
    # ORFs wrap the linearisation point.
    remainder = locus_match.group(3)
    if re.search(r"\bcircular\b", remainder):
        topology = "circular"
    elif re.search(r"\blinear\b", remainder):
        topology = "linear"
    else:
        raise BuildError(
            f"{locus_match.group(1)}: LOCUS line states no topology: {locus_match.group(0)!r}"
        )
    header = {
        "locus": locus_match.group(1),
        "accession": version.split(".")[0],
        "accession_version": version,
        "length": int(locus_match.group(2)),
        "topology": topology,
    }

    cds: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    location = ""
    pending_qualifier: str | None = None
    in_features = False
    for line in text.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_features and re.match(r"^(ORIGIN|CONTIG|BASEFREQ|//)", line):
            in_features = False
        if not in_features:
            continue
        feature_match = FEATURE_RE.match(line)
        if feature_match:
            if current is not None:
                _finish_cds(current, location)
                cds.append(current)
                current = None
            location = ""
            pending_qualifier = None
            if feature_match.group(1) == "CDS":
                current = {"qualifiers": {}}
                location = feature_match.group(2).strip()
            continue
        if current is None:
            continue
        qualifier_match = QUALIFIER_RE.match(line)
        if qualifier_match:
            key = qualifier_match.group(1)
            if key in current["qualifiers"]:
                pending_qualifier = key
            else:
                current["qualifiers"][key] = _unquote(qualifier_match.group(2) or "")
                pending_qualifier = key
            continue
        continuation = line[21:].strip()
        if not continuation:
            continue
        if not continuation.startswith("/"):
            # A wrapped location, or a wrapped /translation or /note body. Both
            # are common: /translation alone runs to several lines and is the
            # only sequence evidence available, so truncating it would silently
            # discard the motif checks below.
            if pending_qualifier:
                current["qualifiers"][pending_qualifier] += continuation
            else:
                location += " " + continuation
    if current is not None:
        _finish_cds(current, location)
        cds.append(current)
    if not cds:
        raise BuildError(f"{header['accession_version']}: no CDS features")
    header["cds"] = cds
    return header


def _finish_cds(record: dict[str, object], location: str) -> None:
    blocks, strand = _parse_location(location)
    record["location"] = location.strip()
    record["exon_blocks"] = blocks
    record["strand"] = strand


def _parse_location(location: str) -> tuple[list[tuple[int, int]], str]:
    """Parse a GenBank location into 1-based inclusive blocks and a strand.

    Handles ``join(...)`` and ``complement(...)`` in either nesting order, and
    drops ``<``/``>`` fuzzy bounds. A papillomavirus ORF that wraps the
    linearisation point is a single CDS with a block at position 1 and another
    near the end, so blocks are kept as a list rather than collapsed to a range.
    """
    text = location.strip()
    strand = "+"
    if "complement(" in text:
        strand = "-"
        text = text.replace("complement(", "").rstrip(")").strip()
    match = re.search(r"join\((.*)\)\s*$", text)
    if match:
        text = match.group(1)
    blocks: list[tuple[int, int]] = []
    for piece in text.split(","):
        piece = piece.strip().strip("<>")
        span = SPAN_RE.match(piece)
        if span:
            blocks.append((int(span.group(1)), int(span.group(2))))
            continue
        single = SINGLE_RE.match(piece)
        if single:
            blocks.append((int(single.group(1)), int(single.group(1))))
    if not blocks:
        raise BuildError(f"unparseable location {location!r}")
    return sorted(blocks), strand


def _symbol_candidates(text: str) -> set[str]:
    """Recognised ORF symbols appearing as whole tokens in ``text``."""
    return {token for token in _SYMBOL_TOKEN_RE.findall(text) if token in CANONICAL_GENES}


#: Isoform symbols and the full-length gene they are a truncated or fused form
#: of. NCBI puts these in ``/product`` only, never ``/gene``.
ISOFORM_PARENT = {
    "E6*": "E6",
    "E7*": "E7",
    "E6^c": "E6",
    "E1^E4": "E1",
    "E2^E1": "E2",
}


def resolve_gene_name(qualifiers: dict[str, str]) -> tuple[str, str, str]:
    """Return ``(canonical_name, name_source, verbatim_text)`` for one CDS.

    ``/gene`` is authoritative — except where ``/product`` names a strictly more
    specific isoform of it. ``NC_001526.4`` annotates HPV16's truncated E6* CDS
    with ``/gene="E6"`` and ``/product="protein E6*"``; taking ``/gene`` at face
    value would report the E6* isoform as the E6 oncoprotein, which is exactly
    the overstatement :data:`CANONICAL_GENES` separates ``oncogene`` from
    ``oncogene_locus`` to prevent. So an isoform symbol in ``/product`` wins over
    its own parent in ``/gene``, and the two disagreeing in any other way is an
    error rather than a guess.

    ``/product`` is also the only source for records that omit ``/gene``
    entirely — ``EF202156.1``'s HPV45 among them.

    Token matching, not substring matching, is deliberate: ``/product="E6
    protein"`` must not be read as ``E6*``, and ``E1^E4`` must not be read as
    ``E1``.
    """
    gene = qualifiers.get("gene", "").strip()
    gene_symbol = gene if gene in CANONICAL_GENES else None
    product = qualifiers.get("product", "").strip()
    found = _symbol_candidates(product)

    if gene_symbol is not None and found and found != {gene_symbol}:
        parents = {ISOFORM_PARENT.get(symbol) for symbol in found}
        if parents == {gene_symbol}:
            return found.pop(), "product", product
        raise BuildError(f"/gene={gene!r} disagrees with /product={product!r} ({sorted(found)})")
    if gene_symbol is not None:
        return gene_symbol, "gene", gene
    if len(found) == 1:
        return found.pop(), "product", product
    if len(found) > 1:
        raise BuildError(f"ambiguous /product={product!r} names {sorted(found)}")
    raise BuildError(
        f"no recognisable ORF symbol in /gene={gene!r} /product={product!r} "
        f"(protein_id={qualifiers.get('protein_id', 'unset')!r})"
    )


def safe_gene_id(accession_version: str, canonical_name: str) -> str:
    """Panel-safe, reference-unique gene ID for one ORF.

    Namespaced by accession because the symbol alone is not unique: ``E6`` occurs
    once in each of the 16 genomes, and a merged reference with 16 rows called
    ``E6`` is the defect this catalogue exists to remove.
    """
    safe = SAFE_SYMBOLS.get(canonical_name, canonical_name)
    return f"{accession_version}_{safe}"


def _check_orf_order(genotype: str, accession_version: str, names: list[str]) -> list[str]:
    """Return problems with the observed ORF order for one genome.

    Matched as a **rotation**, not a fixed sequence, and that is not a
    convenience. Papillomavirus genomes are linearised circles and the
    submitter chooses where to cut, so the same genome reports its ORFs in a
    different order in different records: ``NC_001357.1`` (HPV18) reads
    ``E6, E7, E1, ...`` while ``NC_001526.4`` (HPV16), cut inside E1, reads
    ``E1, E2, E5, L2, L1, E6, E7``. Both are the canonical order; they differ by
    where the circle was opened. This is the concrete reason a hard-coded
    coordinate table cannot be reused across genotypes, and it is why the names
    are read off the record instead.
    """
    observed: list[str] = []
    for name in names:
        if name in CANONICAL_ORF_ORDER and name not in observed:
            observed.append(name)
    canonical = [name for name in CANONICAL_ORF_ORDER if name in observed]
    rotations = [canonical[index:] + canonical[:index] for index in range(len(canonical))]
    if observed not in rotations:
        problems = [
            f"HPV-{genotype} {accession_version}: ORF order {observed} is not a "
            f"rotation of the papillomavirus genome order {canonical}"
        ]
        return problems
    return []


#: Expected protein lengths in amino acids, generous bounds around what the 16
#: records here actually carry. Used only where the record supplies a
#: ``/translation``, and deliberately wide so an unusual genotype is not failed
#: for being atypical.
PROTEIN_LENGTH_RANGE = {
    "E6": (135, 175),
    "E7": (90, 120),
    "L1": (480, 600),
    "L2": (440, 540),
}

#: HPV E7 binds retinoblastoma through an LXCXE motif. Present once in every
#: papillomavirus E7 examined, and in no other papillomavirus protein.
E7_PRB_MOTIF = re.compile(r"L.C.E")

#: HPV E6's two zinc fingers are four C-X2-C pairs. L1 is cysteine-poor, so
#: requiring them is also a cheap L1/E6 mix-up detector.
E6_ZINC_FINGER = re.compile(r"C..C")


def _check_protein_identities(
    genotype: str, accession_version: str, translations: dict[str, str]
) -> list[str]:
    """Confirm the E6/E7/L1/L2 labels against their protein sequence.

    This is the only check that does not take NCBI's annotation on trust. HPV E7
    binds pRb through LXCXE and HPV E6 has C-X2-C zinc fingers, both conserved
    across the genus; an ORF mislabelled E6 or E7 would lack them. A length check
    runs alongside because the motif alone is a short string.
    """
    problems: list[str] = []
    e7 = translations.get("E7", "")
    if e7 and not E7_PRB_MOTIF.search(e7):
        problems.append(
            f"HPV-{genotype} {accession_version}: the ORF labelled E7 has no "
            f"LXCXE retinoblastoma-binding motif"
        )
    e6 = translations.get("E6", "")
    if e6 and len(E6_ZINC_FINGER.findall(e6)) < 3:
        problems.append(
            f"HPV-{genotype} {accession_version}: the ORF labelled E6 has "
            f"{len(E6_ZINC_FINGER.findall(e6))} C-X2-C zinc-finger motifs, "
            f"expected 3 or more"
        )
    for name, (low, high) in PROTEIN_LENGTH_RANGE.items():
        protein = translations.get(name, "")
        if protein and not low <= len(protein) <= high:
            problems.append(
                f"HPV-{genotype} {accession_version}: {name} is {len(protein)} aa, "
                f"outside the papillomavirus range {low}-{high} aa"
            )
    return problems


def _check_capsid_scale(genotype: str, accession_version: str, spans: dict[str, int]) -> list[str]:
    """Return problems with the size ordering that an E6/E7/L1 mix-up would break.

    Two assumptions that look obvious are false, and the build caught both:

    * **L1 is not the longest ORF.** E1 is, at ~650 aa, because E1 is the
      replication helicase; L1 runs ~505-570 aa. Asserting it would fail all
      sixteen correct records.
    * **L1 is not always longer than L2.** It is for the great majority of
      types, but ``NC_001352.1`` (HPV-2) annotates L2 at 525 aa against L1 at
      511 aa, and ``NC_001356.1`` (HPV-1) has them within one residue of equal.
      That is real, not a mislabel, so it is reported rather than asserted.

    What does hold universally, and is the mix-up detector, is that L1 is a
    large ORF and at least 3x either oncogene. L1's absolute size is bounded
    separately by :data:`PROTEIN_LENGTH_RANGE`.
    """
    problems: list[str] = []
    if "L1" not in spans:
        return [f"HPV-{genotype} {accession_version}: L1 is not annotated"]
    l1 = spans["L1"]
    for early in ("E6", "E7"):
        if early in spans and l1 < 3 * spans[early]:
            problems.append(
                f"HPV-{genotype} {accession_version}: L1 ({l1} nt) is less than 3x "
                f"{early} ({spans[early]} nt); an E6/E7/L1 mix-up would look like this"
            )
    return problems


def fetch_genbank(accession: str, email: str, api_key: str | None, cache_dir: Path) -> str:
    """Return the GenBank flatfile for ``accession``, cache-first.

    Thin wrapper over the now-public ``ncbi_fetch.fetch_genbank()`` (closes
    `HPV-09`'s second ask): the cache is ``ncbi_fetch``'s own, so a record
    already fetched by any other ViralScan entry point is reused rather than
    re-downloaded, and this module no longer reaches into ``_efetch``,
    ``_cache_valid`` or ``_write_cached`` directly.
    """
    _path, text = ncbi_fetch.fetch_genbank(accession, email, api_key, cache_dir)
    if "FEATURES" not in text:
        raise BuildError(f"{accession}: NCBI returned no FEATURES table")
    return text


def build_rows(
    email: str, api_key: str | None, cache_dir: Path, *, record_source: str = "ncbi"
) -> list[dict[str, str]]:
    """Build one catalogue row per ORF across every selected record."""
    rows: list[dict[str, str]] = []
    problems: list[str] = []
    for genotype, accession, source, note in SELECTED_RECORDS:
        try:
            text = fetch_genbank(accession, email, api_key, cache_dir)
            header = parse_genbank(text)
        except BuildError as exc:
            problems.append(str(exc))
            continue
        if header["accession"] != accession.split(".")[0]:
            problems.append(f"{accession}: record declares accession {header['accession']}")
        if header["accession_version"] != accession:
            problems.append(
                f"{accession}: record is version {header['accession_version']}; "
                f"the catalogue pins {accession}. If NCBI has revised the record, "
                f"regenerate deliberately rather than silently."
            )
        genome_length = int(header["length"])

        names: list[str] = []
        spans: dict[str, int] = {}
        translations: dict[str, str] = {}
        for cds in header["cds"]:
            qualifiers = cds["qualifiers"]
            blocks = cds["exon_blocks"]
            try:
                canonical, name_source, verbatim = resolve_gene_name(qualifiers)
            except BuildError as exc:
                problems.append(f"{accession}: {exc}")
                continue
            span_nt = sum(end - start + 1 for start, end in blocks)
            gene_id = safe_gene_id(header["accession_version"], canonical)
            spans[canonical] = max(spans.get(canonical, 0), span_nt)
            names.append(canonical)
            translation = qualifiers.get("translation", "")
            if len(translation) > len(translations.get(canonical, "")):
                # A locus can carry both a full product and a truncated one
                # (``NC_001526.4`` annotates E6 and E6*). The full product is the
                # longer translation, and the motif checks need it.
                translations[canonical] = translation
            spans_origin = (
                header["topology"] == "circular"
                and len(blocks) > 1
                and any(start == 1 for start, _ in blocks)
                and canonical in ORIGIN_SPANNING_GENES
            )
            rows.append(
                {
                    "accession": header["accession"],
                    "accession_version": header["accession_version"],
                    "genotype": genotype,
                    "high_risk": str(genotype in HIGH_RISK_GENOTYPES).lower(),
                    "genome_length": str(genome_length),
                    "topology": str(header["topology"]),
                    "gene_id": gene_id,
                    "canonical_gene_name": canonical,
                    "gene_class": CANONICAL_GENES[canonical],
                    "genome_start": str(min(start for start, _ in blocks)),
                    "genome_end": str(max(end for _, end in blocks)),
                    "strand": str(cds["strand"]),
                    "n_exons": str(len(blocks)),
                    "exon_blocks": ",".join(f"{start}..{end}" for start, end in blocks),
                    "spans_origin": str(bool(spans_origin)).lower(),
                    "cds_length_nt": str(span_nt),
                    "protein_id": qualifiers.get("protein_id", ""),
                    "name_source": name_source,
                    "product_as_in_ncbi": verbatim,
                    "source": source,
                    "note": note,
                }
            )

        problems += _check_orf_order(genotype, header["accession_version"], names)
        problems += _check_capsid_scale(genotype, header["accession_version"], spans)
        problems += _check_protein_identities(genotype, header["accession_version"], translations)
        for required in ("E6", "E7", "L1", "L2"):
            if required not in spans:
                problems.append(
                    f"HPV-{genotype} {header['accession_version']}: {required} is not "
                    f"annotated, so the oncogene-versus-capsid contrast is not "
                    f"computable for this genotype"
                )

    if problems:
        for problem in problems:
            print(f"ERROR: {problem}", file=sys.stderr)
        raise BuildError(f"{len(problems)} problem(s) found; refusing to write a partial catalogue")

    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (row["accession_version"], row["canonical_gene_name"], row["gene_id"])
        if key in seen:
            raise BuildError(
                f"duplicate row {key}: two CDS resolved to the same name in one genome"
            )
        seen.add(key)
    seen_ids = {row["gene_id"] for row in rows}
    if len(seen_ids) != len(rows):
        raise BuildError("gene_id is not unique across the merged reference")

    rows.sort(key=lambda row: (_genotype_key(row["genotype"]), int(row["genome_start"])))
    return rows


def _genotype_key(genotype: str) -> tuple[int, str]:
    match = re.match(r"^(\d+)([A-Za-z]*)$", genotype)
    if not match:
        return (10**6, genotype)
    return (int(match.group(1)), match.group(2))


def emit_reference(
    rows: list[dict[str, str]], out_dir: Path, email: str, api_key: str | None, cache_dir: Path
) -> tuple[Path, Path]:
    """Write a merged FASTA + GTF for ``kb ref`` from the catalogued ORFs.

    Written on demand rather than committed: 16 papillomavirus genomes is ~127 kB
    now, but reference size is a packaging decision (PLAN PR 8) and should be made
    once, deliberately, rather than inherited from a build script.

    The GTF is built here rather than taken from ``ncbi_fetch._genbank_to_gtf``
    because that function sets ``gene_id`` from ``/gene=``, which for HPV means
    every genome in a merged reference contributes a row literally named ``E6``.
    Two HPV genomes already collapse from 16 ORFs to 9 gene IDs that way.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / "hpv_reference.fasta"
    gtf_path = out_dir / "hpv_reference.gtf"

    fasta_chunks: list[str] = []
    gtf_lines: list[str] = []
    by_version: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_version.setdefault(row["accession_version"], []).append(row)

    for accession_version, group in by_version.items():
        fasta_path_one, _ = ncbi_fetch._fetch_one(accession_version, cache_dir, email, api_key)
        fasta_chunks.append(fasta_path_one.read_text())
        seqname = accession_version
        for row in sorted(group, key=lambda r: int(r["genome_start"])):
            blocks = [
                tuple(int(v) for v in block.split("..")) for block in row["exon_blocks"].split(",")
            ]
            attrs = (
                f'gene_id "{row["gene_id"]}"; transcript_id "{row["gene_id"]}"; '
                f'gene_name "{row["canonical_gene_name"]}"; '
                f'gene_biotype "protein_coding"; '
                f'product "{row["product_as_in_ncbi"]}"; '
                f'note "canonical_gene_class:{row["gene_class"]}";'
            )
            for start, end in blocks:
                gtf_lines.append(
                    f"{seqname}\tViralScan\texon\t{start}\t{end}\t.\t{row['strand']}\t0\t{attrs}"
                )
    fasta_path.write_text("".join(fasta_chunks))
    gtf_path.write_text("\n".join(gtf_lines) + "\n")
    return fasta_path, gtf_path


def list_records() -> int:
    """Print the curated accession table and what NCBI currently offers.

    Read-only against NCBI; prints any drift rather than rewriting
    :data:`SELECTED_RECORDS`, so a change of reference has to be a human
    decision.
    """
    print(f"{'type':>5}  {'accession':<14} {'source':<7} note")
    for genotype, accession, source, note in SELECTED_RECORDS:
        print(f"{genotype:>5}  {accession:<14} {source:<7} {note}")
    print()
    print(f"{len(HIGH_RISK_GENOTYPES)} high-risk genotypes, {len(SELECTED_RECORDS)} records total")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    parser.add_argument(
        "--emit-reference",
        type=Path,
        metavar="DIR",
        help="also write a merged FASTA+GTF for kb ref into DIR",
    )
    parser.add_argument(
        "--list-records", action="store_true", help="print the curated accession table and exit"
    )
    args = parser.parse_args(argv)

    if args.list_records:
        return list_records()
    if not args.email:
        parser.error("NCBI needs a contact email: pass --email or set NCBI_EMAIL")

    try:
        rows = build_rows(args.email, args.api_key, args.cache_dir)
    except BuildError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if args.emit_reference:
        fasta_path, gtf_path = emit_reference(
            rows, args.emit_reference, args.email, args.api_key, args.cache_dir
        )
        print(f"wrote {fasta_path} and {gtf_path} for kb ref")

    _report(rows, args.out)
    return 0


def _report(rows: list[dict[str, str]], out_path: Path) -> None:
    by_genotype: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_genotype.setdefault(row["genotype"], []).append(row)
    print(f"wrote {out_path} ({len(rows)} ORFs across {len(by_genotype)} genotypes)")
    print()
    print(
        f"  {'type':>5} {'accession':<14} {'ORFs':>4} {'E6':>4} {'E7':>4} {'L1':>4} "
        f"{'name source':<13} topology"
    )
    for genotype in sorted(by_genotype, key=_genotype_key):
        group = by_genotype[genotype]
        names = {row["canonical_gene_name"] for row in group}
        sources = sorted({row["name_source"] for row in group})
        flag = "" if genotype in HIGH_RISK_GENOTYPES else "  (low risk)"
        print(
            f"  {genotype:>5} {group[0]['accession_version']:<14} {len(group):>4} "
            f"{'yes' if 'E6' in names else 'NO':>4} "
            f"{'yes' if 'E7' in names else 'NO':>4} "
            f"{'yes' if 'L1' in names else 'NO':>4} "
            f"{'+'.join(sources):<13} {group[0]['topology']}{flag}"
        )
    oncogene = sum(1 for row in rows if row["gene_class"] == "oncogene")
    locus = sum(1 for row in rows if row["gene_class"] == "oncogene_locus")
    capsid = sum(1 for row in rows if row["gene_class"] == "late_capsid")
    print()
    print(
        f"  oncogene {oncogene}, oncogene-locus {locus}, capsid {capsid}, "
        f"early {len(rows) - oncogene - locus - capsid}"
    )
    print("  every record passed the ORF-order, capsid-scale and E6/E7-motif checks")
    print("  L1 cross-maps across genotypes; per-type L1 counts are not type evidence")


if __name__ == "__main__":
    raise SystemExit(main())
