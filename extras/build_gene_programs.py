"""One-time generator for ``src/viralscan/data/gene_programs.tsv``.

ViralScan's second reporting layer infers a *gene programme* (latent vs
productive) for viruses the first layer detected. That inference needs three
things this script produces:

1. **Curated biology.** Which ORF belongs to which programme. This cannot be
   derived from a file -- it is textbook virology plus genome-specific
   nomenclature -- so it lives in :data:`CURATED_MARKERS` below.
2. **Panel attributes.** ``gene``/``product``/``gene_biotype``/``has_cds`` for
   each marker. These *are* derivable, but only from the bundled panel: the
   merged and STARsolo GTFs do not carry them (the merged panel has no
   ``gene`` or ``product`` attributes at all, and 1-bp ``exon 1 1`` stubs for
   every curated gene).
3. **Overlap groups.** Computed here from real GTF intervals rather than
   hand-assigned, because this is the column the whole design rests on.

Why overlap groups exist
------------------------
Measured on the bundled EBV LCL run (``SRR12682296``), which is latently
infected by construction, the per-gene totals were:

    LATENT (10 genes)  236,342     LYTIC (12 genes)  247,633
      BNLF2a   49,662                 BHLF1   142,954   <- highest EBV gene
      LMP-1    50,237                 BMRF1    45,005
      BaRF1.1  46,956                 BBLF2/3  13,246
      BNLF2b   45,108                 BRLF1     9,560
      EBNA-1.1    920  <-             BZLF1     9,308

A latent:lytic ratio of 1.15 in a cell line defined by latency is not biology.
EBNA-1 is expressed from every latent episome and must be present in every
infected cell, yet it sits ~155x below BHLF1. The cause is pervasive
overlapping-ORF cross-mapping: EBV's latent transcripts are transcribed from a
region densely packed with nested and antisense lytic ORFs, so reads
cross-map in both directions and inflate the lytic side to match the latent
side. Summing 94 EBV genes into one row hides it completely.

The defence is to require *breadth across distinct overlap groups* rather than
a count of genes. Two nested lytic ORFs are one piece of evidence, not two, so
groups are derived by interval intersection within a genome and breadth counts
groups. Hand-assigning them across 94 EBV and 173 CMV ORFs would be both
tedious and unverifiable; computing them makes the catalogue reproducible and
lets a reviewer check the arithmetic.

Gene-ID forms
-------------
Gene IDs are not portable across panels:

============  ==========================  ==========================  ==================
panel         EBV                        CMV                         HHV-6B
============  ==========================  ==========================  ==================
bundled       ``EPSTEIN_HHV4_BZLF1``     ``HUM_CYTO_HHV5wtgp048``    ``HUM_HERP6B_U94``
merged        ``EPSTEIN_HHV4_BZLF1``     ``HUM_CYTO_HHV5wtgp048``    ``HUM_HERP6B_U94``
STARsolo      ``HHV4_BZLF1``             ``HHV5wtgp048``             ``HUM_HERP6B_U94``
============  ==========================  ==========================  ==================

STARsolo strips the ``PREFIX_`` from most curated genes but **retains** it for
HHV-6B, and collapses ``*_unassigned_gene_N`` to a bare ``unassigned_gene_N``
that is shared by seven or more genomes. Those two facts are why
``do_not_normalise`` and ``available_in_starsolo`` exist as columns.

Usage
-----
    python extras/build_gene_programs.py [--out PATH] [--report-groups]

Exits non-zero if any curated marker fails to resolve, so a typo in
:data:`CURATED_MARKERS` cannot ship silently.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "src", "viralscan", "data")
DEFAULT_OUT = os.path.join(DATA_DIR, "gene_programs.tsv")

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
ATTR_RE_CACHE: dict[str, dict[str, re.Pattern[str]]] = {}


# ---------------------------------------------------------------------------
# Curated biology. Hand-written; see module docstring.
#
# `programme` is "latent" or "productive". A marker may be listed once; where a
# gene serves as an anchor for one programme it is not repeated for the other.
# `note` is carried into the shipped TSV so a reader never has to guess why a
# gene is in the set.
# ---------------------------------------------------------------------------
CURATED_MARKERS: tuple[dict[str, str], ...] = (
    # ── Epstein-Barr virus (NC_007605) — panel_completeness: complete ──────
    # Latency programme. EBNA-1 is the anchor: it is expressed from every
    # latent episome, so it is present in every latently infected cell. RefSeq
    # annotates the locus twice -- once as an mRNA record (no CDS,
    # gene_biotype "other") and once as the protein -- and _resolve returns
    # both, so a single curated entry covers the whole locus.
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "EBNA-1",
        "programme": "latent",
        "note": "anchor for every latency type; annotated as both an mRNA record "
        "(gene_biotype other, no CDS) and a protein record",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "EBNA-2",
        "programme": "latent",
        "note": "latent; co-recruits TBP on the EBNA-LP promoter",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "EBNA-LP",
        "programme": "latent",
        "note": "leader protein; keeps the EBNA2B promoter open during latency",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "LMP-1",
        "programme": "latent",
        "note": "latency II/III; oncogenic; NF-kB independent",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "LMP-2A",
        "programme": "latent",
        "note": "latency II/III; blocks BCR signalling",
    },
    # BARF1 is deliberately absent: it is latent only in epithelial cancers
    # (NPC, EBV-gastric; PMID 32708965, 39329759), so in B cells and LCLs it is not
    # a latency marker. The old entry also pulled in BaRF1 (BamHI-a, the lytic
    # ribonucleotide reductase) through a case-insensitive match.
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BNLF2a",
        "programme": "latent",
        "note": "BamHI-W latent transcript; heavily overlapped by lytic ORFs",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BNLF2b",
        "programme": "latent",
        "note": "BamHI-W latent transcript; heavily overlapped by lytic ORFs",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BHRF1",
        "programme": "latent",
        "note": "latent; Becl-2-like; absent from the lytic programme",
    },
    # Productive programme. Deliberately restricted to the lytic DNA
    # polymerase machinery plus the two immediate-early transactivators: these
    # are the least overlapped ORFs in the genome and therefore the ones where
    # unique-placing molecules carry the most information.
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BcLF1",
        "programme": "productive",
        "note": "vlf1, lytic DNA polymerase; master lytic-cycle switch",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BALF5",
        "programme": "productive",
        "note": "lytic polymerase processivity subunit; B* family, heavily overlapped",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BBLF4",
        "programme": "productive",
        "note": "lytic polymerase catalytic subunit",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BBLF2",
        "programme": "productive",
        "note": "lytic polymerase helicase; BBLF2/BBLF3 locus",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BBLF1",
        "programme": "productive",
        "note": "lytic polymerase priming subunit",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BTRF1",
        "programme": "productive",
        "note": "lytic polymerase RNB subunit",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BZLF1",
        "programme": "productive",
        "note": "Zta; lytic immediate-early transactivator; the canonical lytic marker",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BRLF1",
        "programme": "productive",
        "note": "Rta; lytic immediate-early transactivator",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BGLF4",
        "programme": "productive",
        "note": "viral thymidine kinase; lytic",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BMRF1",
        "programme": "productive",
        "note": "lytic polymerase processivity factor",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BSRF1",
        "programme": "productive",
        "note": "lytic transactivator; binds RTA",
    },
    # ── Human cytomegalovirus (NC_006273) — partial ───────────────────────
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL123",
        "programme": "productive",
        "note": "IE1; immediate-early; marks reactivation, not latency (PMID 29535194)",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL122",
        "programme": "productive",
        "note": "IE2; immediate-early; marks reactivation, not latency (PMID 29535194)",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL111A",
        "programme": "latent",
        "note": "cmvIL-10; latency-associated; NOT UL138 despite the shared locus",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL138",
        "programme": "latent",
        "note": "latency; mRNA in the UL138-UL139-UL140 region",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL44",
        "programme": "productive",
        "note": "DNA polymerase processivity factor; lytic",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL84",
        "programme": "productive",
        "note": "UL84; lytic",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL97",
        "programme": "productive",
        "note": "UL97 kinase; lytic",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL83",
        "programme": "productive",
        "note": "pp65 tegument protein; lytic; the standard CMV lytic marker",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL54",
        "programme": "productive",
        "note": "UL54; lytic",
    },
    # ── Human herpesvirus 6A (NC_001664) — complete ───────────────────────
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U90",
        "programme": "latent",
        "note": "IE1; HHV-6A IE1 is U90, NOT U94 (HHV-6B nomenclature differs)",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U86",
        "programme": "latent",
        "note": "IE2; HHV-6A IE2 is U86, not U95",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U91",
        "programme": "latent",
        "note": "spliced antisense IE1; latency-associated",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U94",
        "programme": "latent",
        "note": "Rep protein; latency",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U41",
        "programme": "productive",
        "note": "U41; lytic",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U42",
        "programme": "productive",
        "note": "U42 transactivator; lytic",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U50",
        "programme": "productive",
        "note": "U50 DNA packaging protein; lytic",
    },
    {
        "virus": "Human herpesvirus 6",
        "refseq_gene": "U38",
        "programme": "productive",
        "note": "U38 DNA polymerase; lytic",
    },
    # ── Human herpesvirus 7 (NC_001716) — complete ────────────────────────
    {"virus": "Human herpesvirus 7", "refseq_gene": "U90", "programme": "latent", "note": "IE1"},
    {"virus": "Human herpesvirus 7", "refseq_gene": "U86", "programme": "latent", "note": "IE2"},
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U91",
        "programme": "latent",
        "note": "U91; latency-associated",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U17",
        "programme": "latent",
        "note": "US22 family member; latency",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U42",
        "programme": "productive",
        "note": "multifunctional expression regulator; lytic",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U75",
        "programme": "productive",
        "note": "gH; lytic",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U57",
        "programme": "productive",
        "note": "major capsid; lytic",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U37",
        "programme": "productive",
        "note": "UL37 homolog; lytic",
    },
    # ── HSV-1 / HSV-2 (NC_001806 / NC_001798) — partial ───────────────────
    # Latency is represented by LAT alone. A single transcript is too thin an
    # anchor set: absence of LAT is not evidence of latency, so
    # latency_observable_in_rna=false and the latent state is unreachable.
    {
        "virus": "Human herpesvirus 1",
        "refseq_gene": "LAT",
        "programme": "latent",
        "note": "latency-associated transcript; sole latency transcript, too thin an anchor set",
    },
    {
        "virus": "Human herpesvirus 2",
        "refseq_gene": "LAT",
        "programme": "latent",
        "note": "latency-associated transcript; sole latency transcript, too thin an anchor set",
    },
    {
        "virus": "Human herpesvirus 1",
        "refseq_gene": "UL44",
        "programme": "productive",
        "note": "gC; lytic",
    },
    {
        "virus": "Human herpesvirus 1",
        "refseq_gene": "UL54",
        "programme": "productive",
        "note": "UL54; lytic",
    },
    {
        "virus": "Human herpesvirus 1",
        "refseq_gene": "UL9",
        "programme": "productive",
        "note": "helicase/primase; lytic",
    },
    {
        "virus": "Human herpesvirus 1",
        "refseq_gene": "UL41",
        "programme": "productive",
        "note": "vhs; lytic",
    },
    {
        "virus": "Human herpesvirus 2",
        "refseq_gene": "UL44",
        "programme": "productive",
        "note": "gC; lytic",
    },
    {
        "virus": "Human herpesvirus 2",
        "refseq_gene": "UL54",
        "programme": "productive",
        "note": "UL54; lytic",
    },
    {
        "virus": "Human herpesvirus 2",
        "refseq_gene": "UL9",
        "programme": "productive",
        "note": "helicase/primase; lytic",
    },
    {
        "virus": "Human herpesvirus 2",
        "refseq_gene": "UL41",
        "programme": "productive",
        "note": "vhs; lytic",
    },
    # ── VZV (NC_001348) — partial ─────────────────────────────────────────
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF4",
        "programme": "latent",
        "note": "IE62; VZV latency transcript; not explicitly annotated in the panel",
    },
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF61",
        "programme": "productive",
        "note": "ICP0 homologue; lytic",
    },
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF62",
        "programme": "productive",
        "note": "ICP4 homologue; lytic",
    },
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF33",
        "programme": "productive",
        "note": "maturase; lytic",
    },
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF31",
        "programme": "productive",
        "note": "gB; lytic",
    },
    # ── KSHV / HHV-8 (NC_009333) — partial ───────────────────────────────
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF73",
        "programme": "latent",
        "note": "LANA; principal latency antigen. RefSeq annotates it as description "
        "ORF73 with no gene/product attribute",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "K1",
        "programme": "latent",
        "note": "K1; latency",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF16",
        "programme": "productive",
        "note": "vBcl-2; lytic (PMID 20860481); not vGPCR, which is ORF74",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF50",
        "programme": "productive",
        "note": "RTA; master lytic-cycle transactivator",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF17",
        "programme": "productive",
        "note": "MTA; lytic; RTA target",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF37",
        "programme": "productive",
        "note": "lytic",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF43",
        "programme": "productive",
        "note": "lytic",
    },
    # ── HHV-6B (AF157706) — partial ──────────────────────────────────────
    # The bundled HHV-6B GTF carries no attributes at all; these IDs are
    # resolvable only because they are standard HHV-6B nomenclature.
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U94",
        "programme": "latent",
        "note": "ID-only: the bundled GTF has no product text for HHV-6B",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U95",
        "programme": "latent",
        "note": "ID-only: the bundled GTF has no product text for HHV-6B",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U91",
        "programme": "latent",
        "note": "ID-only: the bundled GTF has no product text for HHV-6B",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U41",
        "programme": "productive",
        "note": "ID-only",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U50",
        "programme": "productive",
        "note": "ID-only",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U39",
        "programme": "productive",
        "note": "gB; ID-only",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U48",
        "programme": "productive",
        "note": "gH; ID-only",
    },
)

#: Per-virus facts that are properties of the *virus*, not of any one gene.
#:
#: ``complete``  both programmes are anchored by genes that are plausibly
#:               independent, so a latent call can be supported.
#: ``partial``   one programme's anchor set is too thin to support an absence
#:               claim (typically a single latency transcript), so the latent
#:               state is unreachable.
#:
#: ``latency_observable_in_rna`` is False wherever the latency anchor set is
#: thinner than two independent genes. It is also conceptually False for
#: DNA-level latency in general: scRNA-seq measures transcripts, so a provirus
#: or a silent integrated genome is invisible rather than latent.
VIRUS_FACTS: dict[str, dict[str, object]] = {
    "Epstein-Barr virus": {"panel_completeness": "complete", "latency_observable_in_rna": True},
    # Single-cell HCMV latency shows no restricted latency programme: it mirrors
    # a late-lytic programme at much lower levels (PMID 29535194), so no marker's
    # presence separates latent from lytic cells.
    "Human cytomegalovirus": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 6": {"panel_completeness": "complete", "latency_observable_in_rna": True},
    "Human herpesvirus 7": {"panel_completeness": "complete", "latency_observable_in_rna": True},
    "Human herpesvirus 1": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 2": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 6b": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Varicella-zoster virus": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 8": {"panel_completeness": "partial", "latency_observable_in_rna": False},
}

#: Genomes whose gene IDs must never have zero-padding normalised.
#:
#: HHV-6A uses two padding widths for genuinely different ORFs --
#: ``HHV6gp041`` is U43 while ``HHV6gp41`` is U42, ``HHV6gp050`` is U51 while
#: ``HHV6gp50`` is U57, and six further pairs differ. Normalising padding would
#: silently merge distinct genes.
DO_NOT_NORMALISE_VIRUSES = frozenset({"Human herpesvirus 6"})

#: Which bundled GTF holds each virus's genes.
VIRUS_GTF: dict[str, str] = {
    "Epstein-Barr virus": "Epstein_Barr_virus_NC_007605.gtf",
    "Human cytomegalovirus": "Human_cytomegalovirus_NC_006273.gtf",
    "Human herpesvirus 6": "Human_herpesvirus_6_NC_001664.gtf",
    "Human herpesvirus 6b": "Human_herpesvirus6B.gtf",
    "Human herpesvirus 7": "Human_herpesvirus_7_NC_001716.gtf",
    "Human herpesvirus 1": "Human_herpesvirus_1_NC_001806.gtf",
    "Human herpesvirus 2": "Human_herpesvirus_2_NC_001798.gtf",
    "Varicella-zoster virus": "Varicella_zoster_virus_NC_001348.gtf",
    "Human herpesvirus 8": "Human_herpesvirus_8_NC_009333.gtf",
}

#: STARsolo keeps the prefix for these viruses; every other curated virus has
#: it stripped.
STARsolo_KEEPS_PREFIX = frozenset({"Human herpesvirus 6b"})

#: Panel prefixes stripped by the STARsolo packager.
PANEL_PREFIXES = (
    "EPSTEIN_",
    "HUM_CYTO_",
    "HUM_HERP6_",
    "HUM_HERP7_",
    "HUM_HERP1_",
    "HUM_HERP2_",
    "HUM_HERP8_",
    "VARICELLA_",
    "CERC_HERP_",
)

TSV_COLUMNS = (
    "virus",
    "programme",
    "panel_completeness",
    "latency_observable_in_rna",
    "gene_id_bundled",
    "gene_id_merged",
    "gene_id_starsolo",
    "refseq_gene",
    "product",
    "gene_biotype",
    "has_cds",
    "overlap_group",
    "non_overlapping",
    "available_in_starsolo",
    "do_not_normalise",
    "note",
)


def _attr_regex(field: str) -> re.Pattern[str]:
    if field not in ATTR_RE_CACHE:
        ATTR_RE_CACHE[field] = re.compile(rf'{field} "([^"]*)"')
    return ATTR_RE_CACHE[field]


def parse_gtf(path: str) -> dict[str, dict[str, object]]:
    """Return ``{gene_id: record}`` for one bundled viral GTF.

    ``record`` carries the gene's **exonic blocks** (used for overlap
    grouping), the attribute text, and whether any ``CDS`` feature exists.

    Blocks are collected from ``exon`` rows, falling back to ``CDS`` and then to
    ``gene``. A bounding box is deliberately *not* used. EBV makes the reason
    concrete: ``EPSTEIN_HHV4_LMP-2A`` has exons at 58-272 and again at
    166103-166458, because LMP-2 is spliced across the origin and the genome
    carries terminal repeats. Its bounding box is therefore 1-171,823 -- the
    whole genome -- which would place it in an overlap group with all 95 other
    EBV genes. Cross-mapping happens where read sequence is actually shared, so
    the interval set is the correct primitive and the bounding box is not.
    """
    records: dict[str, dict[str, object]] = {}
    seqname: str | None = None
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature, start, end, attrs = cols[2], int(cols[3]), int(cols[4]), cols[8]
            if feature == "gene":
                seqname = cols[0]
            m = GENE_ID_RE.search(attrs)
            if not m:
                continue
            gene_id = m.group(1)
            rec = records.setdefault(
                gene_id,
                {
                    "blocks": set(),
                    "span_start": start,
                    "span_end": end,
                    "seqname": seqname,
                    "gene": "",
                    "product": "",
                    "description": "",
                    "gene_biotype": "",
                    "has_cds": False,
                    "has_exon": False,
                },
            )
            rec["span_start"] = min(rec["span_start"], start)
            rec["span_end"] = max(rec["span_end"], end)
            if feature == "exon":
                rec["has_exon"] = True
                rec["blocks"].add((start, end))
            elif feature == "CDS":
                rec["has_cds"] = True
                rec["blocks"].add((start, end))
            elif feature == "gene":
                for field in ("gene", "product", "description", "gene_biotype"):
                    hit = _attr_regex(field).search(attrs)
                    if hit and hit.group(1):
                        rec[field] = hit.group(1)
    # Genes with only a `gene` feature (no exon/CDS) fall back to their span.
    for rec in records.values():
        if not rec["blocks"]:
            rec["blocks"] = {(int(rec["span_start"]), int(rec["span_end"]))}
    return records


def build_overlap_groups(records: dict[str, dict[str, object]]) -> dict[str, str]:
    """Group gene IDs that share at least one exonic interval.

    Exact, via union-find over the atomic intervals induced by every block
    boundary. For each atomic interval, every gene covering it is unioned with
    the rest, so genes that overlap transitively (A overlaps B, B overlaps C, A
    does not overlap C) end up in one group -- which is the correct semantics,
    because a read anywhere in that span is ambiguous across all of them.

    Groups are numbered in order of each group's first member's genomic start,
    so identifiers are stable across runs and can be diffed in review.
    """
    boundaries: set[int] = set()
    for rec in records.values():
        for start, end in rec["blocks"]:
            boundaries.add(int(start))
            boundaries.add(int(end) + 1)
    ordered = sorted(boundaries)

    parent: dict[str, str] = {g: g for g in records}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for lo, hi in zip(ordered, ordered[1:]):
        if hi <= lo:
            continue
        covering = [
            g
            for g, rec in records.items()
            if any(int(s) <= lo and hi <= int(e) + 1 for s, e in rec["blocks"])
        ]
        for other in covering[1:]:
            union(covering[0], other)

    members: dict[str, list[str]] = defaultdict(list)
    for gene_id in records:
        members[find(gene_id)].append(gene_id)
    order = sorted(
        members,
        key=lambda root: min(int(records[g]["span_start"]) for g in members[root]),
    )
    return {g: f"g{index + 1}" for index, root in enumerate(order) for g in members[root]}


def starsolo_form(gene_id: str, virus: str) -> str:
    """Gene ID as the STARsolo packager writes it."""
    if virus in STARsolo_KEEPS_PREFIX:
        return gene_id
    for prefix in PANEL_PREFIXES:
        if gene_id.startswith(prefix):
            return gene_id[len(prefix) :]
    return gene_id


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument(
        "--report-groups",
        action="store_true",
        help="print the overlap group size distribution per virus",
    )
    args = parser.parse_args(argv)

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    for virus, gtf_name in sorted(VIRUS_GTF.items()):
        gtf_path = os.path.join(DATA_DIR, gtf_name)
        if not os.path.exists(gtf_path):
            errors.append(f"{virus}: bundled GTF missing: {gtf_path}")
            continue
        records = parse_gtf(gtf_path)
        groups = build_overlap_groups(records)
        facts = VIRUS_FACTS[virus]
        for entry in CURATED_MARKERS:
            if entry["virus"] != virus:
                continue
            refseq = entry["refseq_gene"]
            matches = _resolve(records, refseq)
            if not matches:
                errors.append(f"{virus}: no gene in {gtf_name} matches {refseq!r}")
                continue
            for gene_id in sorted(matches):
                rec = records[gene_id]
                group = groups[gene_id]
                note = entry["note"]
                # A curated entry can match several records for one locus (EBV
                # annotates EBNA-1 twice, as an mRNA and as a protein). Only
                # claim "no CDS" on the record that actually lacks one, or the
                # shipped TSV asserts something false about its sibling.
                if "no CDS" in note and rec["has_cds"]:
                    note = note.split("; anchor")[0].split("; mRNA")[0].strip(" ,;")
                    if not note:
                        note = f"second record of the {refseq} locus"
                # A marker is an independent breadth unit only if no other
                # curated marker of the same virus shares its overlap group.
                rows.append(
                    {
                        "virus": virus,
                        "programme": entry["programme"],
                        "panel_completeness": facts["panel_completeness"],
                        "latency_observable_in_rna": str(
                            facts["latency_observable_in_rna"]
                        ).lower(),
                        "gene_id_bundled": gene_id,
                        "gene_id_merged": gene_id,
                        "gene_id_starsolo": starsolo_form(gene_id, virus),
                        "refseq_gene": refseq,
                        "product": rec["product"] or refseq,
                        "gene_biotype": rec["gene_biotype"] or "",
                        "has_cds": str(bool(rec["has_cds"])).lower(),
                        "overlap_group": group,
                        "non_overlapping": "",  # filled in below
                        "available_in_starsolo": "true",
                        "do_not_normalise": str(virus in DO_NOT_NORMALISE_VIRUSES).lower(),
                        "note": note,
                    }
                )

    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        print(
            f"\n{len(errors)} curated marker(s) failed to resolve; "
            "refusing to write a partial catalogue.",
            file=sys.stderr,
        )
        return 1

    _fill_non_overlapping(rows)
    _flag_starsolo_unavailable(rows)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"wrote {args.out} ({len(rows)} rows, {len({r['virus'] for r in rows})} viruses)")
    _report(rows, verbose=args.report_groups)
    return 0


def _resolve(records: dict[str, dict[str, object]], refseq: str) -> list[str]:
    """Gene IDs in ``records`` that encode ``refseq``.

    Three conventions coexist in the bundled panel, so all three are tried in
    order of decreasing confidence:

    1. **The gene ID is or ends with the name** -- EBV, the clearest case:
       ``EPSTEIN_HHV4_BZLF1`` *is* BZLF1.
    2. **The ``gene`` attribute equals the name** -- the herpesvirus convention:
       ``HUM_HERP1_HHV1gp00p01`` carries ``gene "UL44"``. This is what makes
       HSV-1/2, VZV, CMV, HHV-6A/7 and KSHV resolvable despite their opaque
       ``gpNN`` identifiers.
    3. **The name is a slash-delimited component of the ``gene`` value** --
       RefSeq fuses genuinely overlapping ORFs into one record, e.g.
       ``EPSTEIN_HHV4_BBLF2/BBLF3``. Both components are real genes.
    4. **The ``description`` attribute equals or contains the name** -- needed
       for KSHV LANA, which RefSeq annotates as ``description "ORF73"`` and
       which carries no ``gene`` or ``product`` attribute at all.

    Matching is case-sensitive and anchored -- EBV ``BARF1`` (BamHI-A) and
    ``BaRF1`` (BamHI-a, the lytic ribonucleotide reductase) differ only in case
    and are different genes. A name must match the whole
    ``gene`` value, a whole slash-delimited component of it, or appear as a
    whitespace/punctuation-delimited token of ``description``. Substring
    matching on free text would make ``UL4`` match ``UL41`` and ``UL44``, which
    are different genes.
    """
    by_id, by_gene, by_fused, by_desc = [], [], [], []
    token = re.compile(rf"(?:^|[^A-Za-z0-9]){re.escape(refseq)}(?:[^A-Za-z0-9]|$)")
    for gene_id, rec in records.items():
        if gene_id == refseq or gene_id.endswith("_" + refseq):
            by_id.append(gene_id)
            continue
        gene_attr = str(rec["gene"]).strip()
        if gene_attr == refseq:
            by_gene.append(gene_id)
            continue
        if refseq in [part.strip() for part in gene_attr.split("/")]:
            by_fused.append(gene_id)
            continue
        if token.search(str(rec["description"])):
            by_desc.append(gene_id)
    return sorted(by_id) or sorted(by_gene) or sorted(by_fused) or sorted(by_desc)


def _fill_non_overlapping(rows: list[dict[str, object]]) -> None:
    """A marker is non-overlapping iff no *other* marker of its virus shares
    its overlap group."""
    by_virus_group: dict[tuple[str, str], list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        by_virus_group[(str(row["virus"]), str(row["overlap_group"]))].append(idx)
    for indices in by_virus_group.values():
        # A lone marker in a group is independent by construction.
        shared = len(indices) > 1
        for idx in indices:
            rows[idx]["non_overlapping"] = "false" if shared else "true"


def _flag_starsolo_unavailable(rows: list[dict[str, object]]) -> None:
    """Mark markers unusable under STARsolo.

    The STARsolo packager rewrites ``*_unassigned_gene_N`` as a bare
    ``unassigned_gene_N``, which seven or more genomes share. Those rows carry
    the most biologically important transcripts in EBV (the EBERs) and cannot be
    attributed to a genome there, so they are marked unavailable rather than
    silently mis-assigned.
    """
    for row in rows:
        gene_id = str(row["gene_id_starsolo"])
        if re.match(r"^unassigned_gene_\d+$", gene_id):
            row["available_in_starsolo"] = "false"
            row["note"] = (
                f"{row['note']}; STARsolo collapses this to a bare {gene_id} shared "
                "by >=7 genomes, so it cannot be attributed to this virus there"
            ).strip("; ")


def _report(rows: list[dict[str, object]], *, verbose: bool) -> None:
    by_virus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_virus[str(row["virus"])].append(row)
    print()
    for virus, vrows in sorted(by_virus.items()):
        for programme in ("latent", "productive"):
            prows = [r for r in vrows if r["programme"] == programme]
            if not prows:
                continue
            groups = {str(r["overlap_group"]) for r in prows if r["non_overlapping"] == "true"}
            complete = vrows[0]["panel_completeness"] == "complete"
            flag = "" if complete else "   [partial]"
            print(
                f"  {virus:<26}{programme:<11} {len(prows):>3} genes, "
                f"{len(groups):>2} independent overlap groups{flag}"
            )
    if verbose:
        print()
        for virus, vrows in sorted(by_virus.items()):
            sizes = defaultdict(int)
            for row in vrows:
                sizes[str(row["overlap_group"])] += 1
            biggest = max(sizes.values())
            print(f"  {virus:<26} {len(sizes):>3} groups, largest {biggest} genes")


if __name__ == "__main__":
    raise SystemExit(main())
