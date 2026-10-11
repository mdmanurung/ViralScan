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
3. **Overlap groups.** Computed here from the 31-mers the index actually holds (shared
   sequence between genes, either strand) rather than hand-assigned or inferred from exon
   coordinates, because this is the column the whole design rests on (AUD-01).

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
EBNA-1 is expressed from the latent episome in latency I-III, so a latency I-III cell
line should show it well above a lytic gene, yet it sits ~155x below BHLF1. The cause is pervasive
overlapping-ORF cross-mapping: EBV's latent transcripts are transcribed from a
region densely packed with nested and antisense lytic ORFs, so reads
cross-map in both directions and inflate the lytic side to match the latent
side. Summing 94 EBV genes into one row hides it completely.

The defence is to require *breadth across distinct overlap groups* rather than
a count of genes. Two nested lytic ORFs are one piece of evidence, not two, so
groups are derived from shared 31-mers (cliques of the link graph, see
:func:`build_overlap_groups`) and breadth counts groups. Hand-assigning them across 94 EBV and 173 CMV ORFs would be both
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
    python extras/build_gene_programs.py [--out PATH] [--gtf-dir DIR] [--cache-dir DIR]
        [--kmer-threshold F] [--sensitivity TSV] [--report-groups]

Needs the bundled GTFs and the cached NCBI FASTA/GenBank files (``ncbi_fetch`` cache). Exits
non-zero if any curated marker fails to resolve, if a genome FASTA is not cached, or if two
curated markers of different programmes share >= 50 % of their 31-mers while ending in
different groups (or in one group), so none of these can ship silently.
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
    # Latency programme. EBNA-1 is the anchor: it is expressed from the latent
    # episome in latency I-III (not necessarily in latency 0). RefSeq
    # annotates the locus twice -- once as an mRNA record (no CDS,
    # gene_biotype "other") and once as the protein -- and _resolve returns
    # both, so a single curated entry covers the whole locus.
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "EBNA-1",
        "programme": "latent",
        "note": "anchor for every latency type (latency I-III; not necessarily detected in a "
        "cell in latency 0); annotated as both an mRNA record "
        "(gene_biotype other, no CDS) and a protein record",
        "kinetic_class": "latent",
        "kinetic_pmid": "31266868",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "EBNA-2",
        "programme": "latent",
        "note": "latent; co-recruits TBP on the EBNA-LP promoter",
        "kinetic_class": "latent",
        "kinetic_pmid": "31266868",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "EBNA-LP",
        "programme": "latent",
        "note": "leader protein; keeps the EBNA2B promoter open during latency",
        "kinetic_class": "latent",
        "kinetic_pmid": "31266868",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "LMP-1",
        "programme": "latent",
        "note": "latency II/III; oncogenic; RefSeq product: latent membrane protein LMP-1",
        "kinetic_class": "latent",
        "kinetic_pmid": "31266868",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "LMP-2A",
        "programme": "latent",
        "note": "latency II/III; blocks BCR signalling",
        "kinetic_class": "latent",
        "kinetic_pmid": "31266868",
    },
    # BNLF2a, BNLF2b and BHRF1 are deliberately absent. All three are early
    # lytic (CAGE kinetic class early; PMID 29864140), so they were never latency
    # markers, but they cannot join the productive set either: BNLF2a/b lie
    # wholly inside LMP-1's interval (overlap group g50, antisense), so a
    # productive count there is indistinguishable from latent LMP-1, and the
    # miR-BHRF1-2/3 stem-loops sit in the 3'UTR of latent EBNA-LP transcripts
    # (PMID 28950226), so BHRF1-locus reads are not guaranteed lytic.
    # BARF1 is deliberately absent: it is latent only in epithelial cancers
    # (NPC, EBV-gastric; PMID 32708965, 39329759), so in B cells and LCLs it is not
    # a latency marker. The old entry also pulled in BaRF1 (BamHI-a, the lytic
    # ribonucleotide reductase) through a case-insensitive match.
    # Productive programme. Mostly the lytic DNA polymerase machinery and the two
    # immediate-early transactivators, plus BcLF1 (major capsid, late), BMRF1,
    # BSRF1 and BGLF4 (early); the aim is the least overlapped ORFs, where
    # unique-placing molecules carry the most information.
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BcLF1",
        "programme": "productive",
        "note": "major capsid protein; late (CAGE class late, PMID 29864140); not a "
        "polymerase subunit and not a lytic switch",
        "kinetic_class": "late",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BALF5",
        "programme": "productive",
        "note": "RefSeq product: DNA polymerase catalytic subunit; lytic; B* family, heavily overlapped",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BBLF4",
        "programme": "productive",
        "note": "RefSeq product: helicase-primase helicase subunit; lytic",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BBLF2",
        "programme": "productive",
        "note": "RefSeq product: helicase-primase subunit; lytic; BBLF2/BBLF3 locus",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BBLF1",
        "programme": "productive",
        "note": "RefSeq product: myristylated tegument protein; lytic",
        "kinetic_class": "leaky_late",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BTRF1",
        "programme": "productive",
        "note": "RefSeq product: tegument protein UL88; lytic",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BZLF1",
        "programme": "productive",
        "note": "Zta; lytic immediate-early transactivator; the canonical lytic marker",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BRLF1",
        "programme": "productive",
        "note": "Rta; lytic immediate-early transactivator",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BGLF4",
        "programme": "productive",
        "note": "RefSeq product: tegument serine/threonine protein kinase; lytic; not the EBV "
        "thymidine kinase (BXLF1)",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BMRF1",
        "programme": "productive",
        "note": "RefSeq product: DNA polymerase processivity subunit; lytic",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    {
        "virus": "Epstein-Barr virus",
        "refseq_gene": "BSRF1",
        "programme": "productive",
        "note": "RefSeq product: tegument protein UL51; lytic",
        "kinetic_class": "early",
        "kinetic_pmid": "29864140",
    },
    # ── Human cytomegalovirus (NC_006273) — partial ───────────────────────
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL123",
        "programme": "productive",
        "note": "IE1; immediate-early; marks reactivation, not latency (PMID 29535194)",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "29535194",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL122",
        "programme": "productive",
        "note": "IE2; immediate-early; marks reactivation, not latency (PMID 29535194)",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "29535194",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL111A",
        "programme": "latent",
        "note": "cmvIL-10; latency-associated; NOT UL138 despite the shared locus",
        "kinetic_class": "latent",
        "kinetic_pmid": "29535194",
    },
    {
        "virus": "Human cytomegalovirus",
        "refseq_gene": "UL138",
        "programme": "latent",
        "note": "latency; mRNA in the UL138-UL139-UL140 region",
        "kinetic_class": "latent",
        "kinetic_pmid": "29535194",
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
    # ── Human herpesvirus 6A (NC_001664) — partial ────────────────────────
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U90",
        "programme": "productive",
        "note": "IE1; HHV-6A IE1 is U90, NOT U94 (HHV-6B nomenclature differs); immediate-"
        "early in productive infection (PMID 12388818). Latency-associated "
        "transcripts share this coding region (PMID 11907257), so reads here cannot "
        "be called latent",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "12367753;12388818",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U86",
        "programme": "productive",
        "note": "IE2; HHV-6A IE2 is U86, not U95; IE2 mRNA is expressed under immediate-early"
        " conditions, 2-4 h (PMID 12706083). Latency-associated transcripts share the"
        " IE locus (PMID 11907257), so reads here cannot be called latent",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "12706083",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U91",
        "programme": "latent",
        "note": "spliced antisense IE1; latency-associated",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U94",
        "programme": "latent",
        "note": "Rep protein; latency",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U41",
        "programme": "productive",
        "note": "U41; lytic",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U42",
        "programme": "productive",
        "note": "U42 transactivator; lytic",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U50",
        "programme": "productive",
        "note": "U50 DNA packaging protein; lytic",
    },
    {
        "virus": "Human herpesvirus 6A",
        "refseq_gene": "U38",
        "programme": "productive",
        "note": "U38 DNA polymerase; lytic",
    },
    # ── Human herpesvirus 7 (NC_001716) — partial ─────────────────────────
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U90",
        "programme": "productive",
        "note": "IE1; U89/90 is an immediate-early (alpha) transcript in productive infection"
        " (PMID 10573164); no HHV-7 latency transcript is known",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "10573164",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U86",
        "programme": "productive",
        "note": "IE2; kinetic class not verified for HHV-7. Not a latency marker: no HHV-7 "
        "transcript was detected in latently infected PBMC (PMID 10573164)",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U91",
        "programme": "latent",
        "note": "U91; UNVERIFIED as latency-associated: no HHV-7 latency transcript was "
        "detected in PBMC (PMID 10573164, PMC4944276). Inert while latency is not "
        "observable",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U17",
        "programme": "latent",
        "note": "US22 family member; UNVERIFIED as latency-associated: no HHV-7 latency "
        "transcript was detected in PBMC (PMID 10573164, PMC4944276). Inert while "
        "latency is not observable",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U42",
        "programme": "productive",
        "note": "multifunctional expression regulator; lytic",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "10573164",
    },
    {
        "virus": "Human herpesvirus 7",
        "refseq_gene": "U75",
        "programme": "productive",
        "note": "RefSeq product: tegument protein UL7; lytic",
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
        "note": "RefSeq product: nuclear egress lamina protein; lytic",
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
        "kinetic_class": "latent",
        "kinetic_pmid": "29875240",
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
        "note": "RefSeq product: DNA replication origin-binding helicase; lytic",
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
        "note": "RefSeq product: DNA replication origin-binding helicase; lytic",
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
        "note": "IE4, the ICP27 homologue (not IE62, which is ORF62); RNA and protein are "
        "detected in latently infected ganglia (PMID 15890936); not explicitly "
        "annotated in the panel",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "15890936",
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
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "16537613",
    },
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF33",
        "programme": "productive",
        "note": "RefSeq product: capsid maturation protease; lytic",
    },
    {
        "virus": "Varicella-zoster virus",
        "refseq_gene": "ORF31",
        "programme": "productive",
        "note": "gB; lytic",
    },
    # ── KSHV / HHV-8 (NC_009333) — complete (PROG-08) ───────────────────────────────
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF73",
        "programme": "latent",
        "note": "LANA; principal latency antigen. RefSeq annotates it as description "
        "ORF73 with no gene/product attribute",
        "kinetic_class": "latent",
        "kinetic_pmid": "24453964;9733875",
    },
    # PROG-08 (2026-10-03): K1 removed from the latent set -- its role is in lytic
    # replication (PMID 27307571) and no source supports it as a latency anchor.
    # ORF72/ORF71 share LANA's promoter and 3' end (PMID 9733875), so they are
    # one breadth unit with ORF73 (CO_TRANSCRIBED below).
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF72",
        "programme": "latent",
        "note": "v-cyclin; latency cluster, co-transcribed with LANA (one breadth unit)",
        "kinetic_class": "latent",
        "kinetic_pmid": "9733875",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF71",
        "programme": "latent",
        "note": "vFLIP (K13); latency cluster, 3' end of the LANA/v-cyclin mRNA (one breadth unit)",
        "kinetic_class": "latent",
        "kinetic_pmid": "9733875",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "K12",
        "programme": "latent",
        "note": "kaposin; the most abundant latent transcript, but also induced in "
        "lytic replication (PMID 17913828), so its presence alone does not exclude "
        "lytic activity",
        "kinetic_class": "latent",
        "kinetic_pmid": "9733875;17913828",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "vIRF-3",
        "programme": "latent",
        "note": "LANA2 (ORF K10.5); B-cell-specific latent protein (PEL, Castleman), "
        "not induced by phorbol ester; absent from KS lesions",
        "kinetic_class": "latent",
        "kinetic_pmid": "11119611",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF16",
        "programme": "productive",
        "note": "vBcl-2; lytic (PMID 20860481, a review: cited content unverified); not vGPCR, "
        "which is ORF74",
        "kinetic_class": "early",
        "kinetic_pmid": "24453964",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF50",
        "programme": "productive",
        "note": "RTA; master lytic-cycle transactivator",
        "kinetic_class": "early",
        "kinetic_pmid": "24453964",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF17",
        "programme": "productive",
        "note": "capsid maturation protease (assemblin); late, 48 h (PMID 24453964; PMID "
        "36733461 unverified for the 48 h timing). Not MTA, which is ORF57",
        "kinetic_class": "late",
        "kinetic_pmid": "24453964;36733461",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF37",
        "programme": "productive",
        "note": "lytic",
        "kinetic_class": "early",
        "kinetic_pmid": "24453964",
    },
    {
        "virus": "Human herpesvirus 8",
        "refseq_gene": "ORF43",
        "programme": "productive",
        "note": "lytic",
        "kinetic_class": "late",
        "kinetic_pmid": "24453964;30735904",
    },
    # ── HHV-6B (AF157706) — partial ──────────────────────────────────────
    # The bundled HHV-6B GTF carries no attributes at all; these IDs are
    # resolvable only because they are standard HHV-6B nomenclature.
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U94",
        "programme": "latent",
        "note": "ID-only: the bundled GTF has no product text for HHV-6B",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "33627386",
    },
    {
        "virus": "Human herpesvirus 6b",
        "refseq_gene": "U95",
        "programme": "productive",
        "note": "ID-only: the bundled GTF has no product text for HHV-6B. HHV-6B immediate-"
        "early gene in productive infection, not a latency marker (PMID 17928352, "
        "33627386)",
        "kinetic_class": "immediate_early",
        "kinetic_pmid": "33627386;17928352",
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
    # U90/U86 are immediate-early, so the remaining latent rows are U91
    # (unverified) and U94, which is also immediate-early in productive HHV-6B
    # infection (PMID 33627386): neither separates latency from productive
    # infection, so latency is not observable (same reasoning as HHV-7).
    "Human herpesvirus 6A": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    # No HHV-7 latency-associated transcript is established: viral transcripts
    # were not detected in latently infected PBMC (PMID 10573164), and
    # PBMC latency is described as "very limited transcriptional activity"
    # (PMC4944276). U90/U86 are immediate-early (lytic), so the remaining
    # "latent" rows (U91, U17) are unverified and latency is not observable.
    "Human herpesvirus 7": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 1": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 2": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Human herpesvirus 6b": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    "Varicella-zoster virus": {"panel_completeness": "partial", "latency_observable_in_rna": False},
    # PROG-08 (2026-10-03): three independent latency units -- the LANA/v-cyclin/
    # vFLIP cluster (one mRNA family, PMID 9733875), kaposin K12, and LANA2/vIRF-3
    # (B cells only, PMID 11119611). K12 is also lytic-induced (PMID 17913828).
    "Human herpesvirus 8": {"panel_completeness": "complete", "latency_observable_in_rna": True},
}

#: Markers transcribed as one mRNA family, so they count as one breadth unit even
#: when their CDS models do not overlap (overlap groups are exonic-only).
CO_TRANSCRIBED: dict[str, list[tuple[str, ...]]] = {
    "Human herpesvirus 8": [("ORF73", "ORF72", "ORF71")],  # PMID 9733875
}

#: Genomes whose gene IDs must never have zero-padding normalised.
#:
#: HHV-6A uses two padding widths for genuinely different ORFs --
#: ``HHV6gp041`` is U43 while ``HHV6gp41`` is U42, ``HHV6gp050`` is U51 while
#: ``HHV6gp50`` is U57, and six further pairs differ. Normalising padding would
#: silently merge distinct genes.
DO_NOT_NORMALISE_VIRUSES = frozenset({"Human herpesvirus 6A"})

#: Which bundled GTF holds each virus's genes.
VIRUS_GTF: dict[str, str] = {
    "Epstein-Barr virus": "Epstein_Barr_virus_NC_007605.gtf",
    "Human cytomegalovirus": "Human_cytomegalovirus_NC_006273.gtf",
    "Human herpesvirus 6A": "Human_herpesvirus_6_NC_001664.gtf",
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

#: Allowed ``kinetic_class`` values: the marker's expression class in
#: productive infection, or ``latent`` where it is a latency-restricted
#: transcript, or ``unclassified`` where no retrieved primary source classifies
#: it (never a guess). Classes follow the cited source's own assay, which
#: differs between sources (cycloheximide/PAA inhibitors, CAGE-seq, ribosome
#: profiling time points, or convention), so ``kinetic_pmid`` is part of the
#: meaning. Nominally: ``immediate_early`` needs no protein synthesis, ``early``
#: needs it but not viral DNA replication, ``late`` needs replication,
#: ``leaky_late`` starts early and rises with replication. ``kinetic_class`` is
#: catalogue metadata only; no output table carries it yet.
KINETIC_CLASSES = (
    "latent",
    "immediate_early",
    "early",
    "leaky_late",
    "late",
    "unclassified",
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
    "kinetic_class",
    "kinetic_pmid",
    "note",
)


def _attr_regex(field: str) -> re.Pattern[str]:
    if field not in ATTR_RE_CACHE:
        ATTR_RE_CACHE[field] = re.compile(rf'{field} "([^"]*)"')
    return ATTR_RE_CACHE[field]


def parse_gtf(path: str, cache_dir: str | None = None) -> dict[str, dict[str, object]]:
    """Return ``{gene_id: record}`` for one bundled viral GTF.

    ``record`` carries the gene's **exonic blocks** (used for overlap
    grouping), the attribute text, and whether any ``CDS`` feature exists.

    The file is read through the REF-13 normaliser (``viralscan.gtf_normalise``), so the
    ``exon`` rows are the targets ``kb ref`` indexes: the gene span, or the spliced
    exons of a single-protein CDS. Blocks come from those ``exon`` rows only (``CDS``
    rows would put introns and UTRs back, and the index has neither), falling back to
    the gene span. A bounding box is deliberately *not* used. EBV makes the reason
    concrete: ``EPSTEIN_HHV4_LMP-2A`` has exons at 58-272 and again at
    166103-166458, because LMP-2 is spliced across the origin and the genome
    carries terminal repeats. Its bounding box is therefore 1-171,823 -- the
    whole genome -- which would place it in an overlap group with all 95 other
    EBV genes. Cross-mapping happens where read sequence is actually shared, so
    the interval set is the correct primitive and the bounding box is not.
    """
    from viralscan.gtf_normalise import normalise_gtf_file

    records: dict[str, dict[str, object]] = {}
    for line in normalise_gtf_file(path, cache_dir):
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        feature, start, end, attrs = cols[2], int(cols[3]), int(cols[4]), cols[8]
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
                "seqname": cols[0],
                "gene": "",
                "product": "",
                "cds_product": "",
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
            # RefSeq puts `product` on the CDS rows, not the gene row.
            hit = _attr_regex("product").search(attrs)
            if hit and hit.group(1) and not rec["cds_product"]:
                rec["cds_product"] = hit.group(1)
        elif feature == "gene":
            for field in ("gene", "product", "description", "gene_biotype"):
                hit = _attr_regex(field).search(attrs)
                if hit and hit.group(1):
                    rec[field] = hit.group(1)
    # Genes with only a `gene` feature (no exon/CDS) fall back to their span.
    for rec in records.values():
        cds_product = rec.pop("cds_product")
        rec["product"] = rec["product"] or cds_product
        if not rec["blocks"]:
            rec["blocks"] = {(int(rec["span_start"]), int(rec["span_end"]))}
    return records


KMER_SIZE = 31

#: Two genes are *linked* (a read from one can pseudoalign to the other) when the shared 31-mers
#: are at least this fraction of the smaller gene's 31-mers. Chosen from
#: ``analysis/audit_2026-10/aud01_sensitivity.tsv``; the generator refuses a catalogue in which
#: two curated markers of different programmes share >= :data:`CROSS_PROGRAMME_FLOOR` of their
#: 31-mers whatever this value is.
DEFAULT_KMER_THRESHOLD = 0.5
CROSS_PROGRAMME_FLOOR = 0.5
SENSITIVITY_THRESHOLDS = (0.25, 0.5, 0.9)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _canonical_kmers(seq: str, k: int = KMER_SIZE) -> set[str]:
    """Canonical (min of forward, reverse complement) k-mers of ``seq``, skipping any with N.

    Canonical because ``kb count`` runs unstranded by default (see
    ``analysis/audit_2026-10/M01_strand_note.md``): a read from the antisense strand of gene X
    that matches gene Y's sense sequence is ambiguous, so strand cannot be used to separate
    genes. This is the conservative choice; an explicit ``--strand forward/reverse`` run is a
    subset of it.
    """
    out: set[str] = set()
    for i in range(len(seq) - k + 1):
        word = seq[i : i + k]
        if set(word) <= set("ACGT"):
            out.add(min(word, word.translate(_COMPLEMENT)[::-1]))
    return out


def gene_kmers(
    gtf_path: str, cache_dir: str | None = None, k: int = KMER_SIZE
) -> dict[str, set[str]]:
    """``{gene_id: canonical k-mers of the gene's indexed transcripts}``.

    Sequence comes from the REF-13 normalised ``exon`` rows (what ``kb ref`` indexes, one
    transcript per ``transcript_id``, exons joined in genomic order) and the genome FASTA in the
    ``ncbi_fetch`` cache (``<cache>/<seqname>/<seqname>.fasta``). A gene with several
    transcripts (repeat copies, ``-span`` isoforms) contributes the union. Raises
    ``FileNotFoundError`` when a sequence is not cached: a catalogue regenerated without
    sequence would silently fall back to exon geometry, which is the defect AUD-01 fixes.
    """
    from pathlib import Path

    from viralscan.gtf_normalise import normalise_gtf_file
    from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR

    root = Path(cache_dir) if cache_dir else Path(DEFAULT_CACHE_DIR)
    exons: dict[tuple[str, str], list[tuple[str, int, int]]] = defaultdict(list)
    for line in normalise_gtf_file(gtf_path, cache_dir):
        cols = line.split("\t")
        if line.startswith("#") or len(cols) < 9 or cols[2] != "exon":
            continue
        attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', cols[8]))
        key = (attrs["gene_id"], attrs.get("transcript_id", attrs["gene_id"]))
        exons[key].append((cols[0], int(cols[3]), int(cols[4])))
    genomes: dict[str, str] = {}
    kmers: dict[str, set[str]] = defaultdict(set)
    for (gene_id, _tx), blocks in sorted(exons.items()):
        for seqname in {b[0] for b in blocks}:
            if seqname not in genomes:
                fasta = root / seqname / f"{seqname}.fasta"
                if not fasta.is_file():
                    raise FileNotFoundError(f"cached FASTA missing for {seqname}: {fasta}")
                genomes[seqname] = "".join(fasta.read_text().split("\n", 1)[1].split()).upper()
        seq = "".join(
            genomes[sq][start - 1 : end] for sq, start, end in sorted(blocks, key=lambda b: b[1])
        )
        kmers[gene_id] |= _canonical_kmers(seq, k)
    return dict(kmers)


def kmer_sharing(kmers: dict[str, set[str]]) -> dict[tuple[str, str], float]:
    """``{(a, b): shared / min(|a|, |b|)}`` for every pair (a < b) sharing >= 1 k-mer."""
    holders: dict[str, list[str]] = defaultdict(list)
    for gene_id in sorted(kmers):
        for kmer in kmers[gene_id]:
            holders[kmer].append(gene_id)
    shared: dict[tuple[str, str], int] = defaultdict(int)
    for genes in holders.values():
        for i, a in enumerate(genes):
            for b in genes[i + 1 :]:
                shared[(a, b)] += 1
    return {(a, b): n / min(len(kmers[a]), len(kmers[b])) for (a, b), n in shared.items()}


def coordinate_sharing(records: dict[str, dict[str, object]]) -> dict[tuple[str, str], float]:
    """Fallback when no sequence is available: pairs whose exon blocks overlap on one sequence.

    Keyed by ``seqname`` (genes on different sequences never share coordinates). Strand is
    ignored, as in :func:`_canonical_kmers`. Every overlapping pair gets sharing 1.0.
    """
    out: dict[tuple[str, str], float] = {}
    ids = sorted(records)
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if records[a]["seqname"] != records[b]["seqname"]:
                continue
            if any(
                int(s1) <= int(e2) and int(s2) <= int(e1)
                for s1, e1 in records[a]["blocks"]
                for s2, e2 in records[b]["blocks"]
            ):
                out[(a, b)] = 1.0
    return out


def build_overlap_groups(
    records: dict[str, dict[str, object]],
    sharing: dict[tuple[str, str], float] | None = None,
    threshold: float = DEFAULT_KMER_THRESHOLD,
) -> dict[str, str]:
    """Partition genes into read-ambiguity groups in which **every pair is linked**.

    Genes X and Y are linked when ``sharing[(X, Y)] >= threshold`` (:func:`kmer_sharing`);
    with ``sharing=None`` the fallback is exon-block overlap within one sequence
    (:func:`coordinate_sharing`). A group is a *clique* of the link graph, never a connected
    component: two genes that both overlap a third but not each other (A-B-C chains, the
    EBV BMRF1 ~ BaRF1.1 ~ BORF1 ~ BOLF1 ~ EBNA-1 chain of audit A-01) are not merged.

    Groups are built greedily and deterministically: genes are ordered by
    ``(seqname, span_start, gene_id)``; each ungrouped gene seeds a group, which then absorbs
    every ungrouped gene linked to all current members, strongest partner of the seed first
    (ties by the same order). The result is a partition into
    *maximal* cliques of what remains, not a *minimum* clique partition (NP-hard). Counting
    distinct groups is therefore an upper bound on the independent evidence: it can over-count
    when linked genes land in different cliques, which :func:`audit_markers` reports.

    Groups are numbered ``g1..`` in seed order, so identifiers are stable and diffable.
    """
    if sharing is None:
        sharing = coordinate_sharing(records)

    def linked(a: str, b: str) -> bool:
        return sharing.get((a, b) if a < b else (b, a), 0.0) >= threshold

    order = sorted(
        records,
        key=lambda g: (str(records[g]["seqname"]), int(records[g]["span_start"]), g),  # type: ignore[call-overload]
    )
    groups: dict[str, str] = {}
    n_groups = 0
    for seed in order:
        if seed in groups:
            continue
        n_groups += 1
        label = f"g{n_groups}"
        members = [seed]
        groups[seed] = label
        # strongest partner of the seed first, so two records of one locus stay together
        for cand in sorted(
            order, key=lambda g: -sharing.get((seed, g) if seed < g else (g, seed), 0.0)
        ):
            if cand not in groups and all(linked(cand, m) for m in members):
                members.append(cand)
                groups[cand] = label
    return groups


def audit_markers(
    programmes: dict[str, str],
    groups: dict[str, str],
    sharing: dict[tuple[str, str], float],
    threshold: float,
) -> dict[str, list[tuple[str, str]]]:
    """Check the curated markers (``programmes``: gene id -> programme) against the grouping.

    * ``mixed_group``: a latent and a productive marker in one group;
    * ``cross_programme_split``: markers of different programmes sharing
      >= :data:`CROSS_PROGRAMME_FLOOR` of their 31-mers (smaller gene) yet in different groups,
      which would let one molecule count as breadth for both programmes;
    * ``same_programme_split``: linked (>= ``threshold``) markers of one programme in different
      groups; each such pair can over-count breadth by one.
    """
    out: dict[str, list[tuple[str, str]]] = {
        "mixed_group": [],
        "cross_programme_split": [],
        "same_programme_split": [],
    }
    ids = sorted(programmes)
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            same_group = groups[a] == groups[b]
            f = sharing.get((a, b), 0.0)
            if programmes[a] != programmes[b]:
                if same_group:
                    out["mixed_group"].append((a, b))
                elif f >= CROSS_PROGRAMME_FLOOR:
                    out["cross_programme_split"].append((a, b))
            elif not same_group and f >= threshold:
                out["same_programme_split"].append((a, b))
    return out


SENSITIVITY_COLUMNS = (
    "virus",
    "threshold",
    "n_genes",
    "n_gene_groups",
    "largest_gene_group",
    "n_markers",
    "n_marker_groups",
    "productive_marker_groups",
    "latent_marker_groups",
    "markers_in_mixed_groups",
    "cross_programme_split_pairs",
    "same_programme_split_pairs",
)


def marker_groups(
    virus: str,
    records: dict[str, dict[str, object]],
    programmes: dict[str, str],
    sharing: dict[tuple[str, str], float],
    threshold: float,
) -> dict[str, str]:
    """Overlap groups of the curated markers of ``virus`` (``programmes``: gene id -> programme).

    Only markers are partitioned: breadth is counted over markers, so a non-marker gene must not
    seed (and thereby split) a marker clique. :data:`CO_TRANSCRIBED` families are then forced
    into one group; that is a transcript-level fact, not a sequence-sharing one, so those groups
    are the one place a group is not a clique of the link graph.
    """
    groups = build_overlap_groups({g: records[g] for g in programmes}, sharing, threshold)
    for names in CO_TRANSCRIBED.get(virus, []):
        ids = sorted(
            {g for n in names for g in _resolve(records, n)}, key=lambda g: int(groups[g][1:])
        )
        for g in ids:
            groups[g] = groups[ids[0]]
    return groups


def sensitivity_rows(
    virus: str,
    records: dict[str, dict[str, object]],
    programmes: dict[str, str],
    sharing: dict[tuple[str, str], float],
) -> list[dict[str, object]]:
    """One row per :data:`SENSITIVITY_THRESHOLDS` value for ``analysis/.../aud01_sensitivity.tsv``."""
    out: list[dict[str, object]] = []
    for threshold in SENSITIVITY_THRESHOLDS:
        gene_groups = build_overlap_groups(records, sharing, threshold)
        sizes: dict[str, int] = defaultdict(int)
        for label in gene_groups.values():
            sizes[label] += 1
        groups = marker_groups(virus, records, programmes, sharing, threshold)
        audit = audit_markers(programmes, groups, sharing, threshold)
        by_group: dict[str, set[str]] = defaultdict(set)
        for gene_id, label in groups.items():
            by_group[label].add(programmes[gene_id])
        mixed = {label for label, progs in by_group.items() if len(progs) > 1}
        out.append(
            {
                "virus": virus,
                "threshold": threshold,
                "n_genes": len(records),
                "n_gene_groups": len(sizes),
                "largest_gene_group": max(sizes.values()),
                "n_markers": len(programmes),
                "n_marker_groups": len(by_group),
                "productive_marker_groups": len(
                    {groups[g] for g, p in programmes.items() if p == "productive"}
                ),
                "latent_marker_groups": len(
                    {groups[g] for g, p in programmes.items() if p == "latent"}
                ),
                "markers_in_mixed_groups": sum(1 for g in programmes if groups[g] in mixed),
                "cross_programme_split_pairs": len(audit["cross_programme_split"]),
                "same_programme_split_pairs": len(audit["same_programme_split"]),
            }
        )
    return out


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
        "--gtf-dir",
        default=DATA_DIR,
        help="directory holding the bundled viral GTFs (default: src/viralscan/data)",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="ncbi_fetch cache holding <accession>/<accession>.fasta and .gb (default: "
        "~/.cache/viralscan/ncbi)",
    )
    parser.add_argument(
        "--kmer-threshold",
        type=float,
        default=DEFAULT_KMER_THRESHOLD,
        help="link two genes when they share this fraction of the smaller gene's 31-mers",
    )
    parser.add_argument(
        "--sensitivity",
        default=None,
        help="also write the per-virus group counts at thresholds 0.25/0.5/0.9 to this TSV",
    )
    parser.add_argument(
        "--report-groups",
        action="store_true",
        help="print the overlap group size distribution per virus",
    )
    args = parser.parse_args(argv)

    sys.path.insert(0, os.path.join(REPO_ROOT, "src"))
    from viralscan.repeat_merges import load_repeat_merges

    repeat_merges = load_repeat_merges()
    rows: list[dict[str, object]] = []
    errors: list[str] = []
    warnings: list[str] = []
    linked_by_virus: dict[str, set[frozenset[str]]] = {}
    sensitivity: list[dict[str, object]] = []
    for virus, gtf_name in sorted(VIRUS_GTF.items()):
        gtf_path = os.path.join(args.gtf_dir, gtf_name)
        if not os.path.exists(gtf_path):
            errors.append(f"{virus}: bundled GTF missing: {gtf_path}")
            continue
        records = parse_gtf(gtf_path, args.cache_dir)
        try:
            kmers = gene_kmers(gtf_path, args.cache_dir)
        except FileNotFoundError as exc:
            errors.append(f"{virus}: {exc}")
            continue
        # PROG-13: a repeat copy is one gene with its survivor in the index's t2g, so it is
        # one record here: its exons and k-mers join the survivor's before groups are computed.
        for copy, survivor in repeat_merges.items():
            if copy in records and survivor in records:
                gone = records.pop(copy)
                records[survivor]["blocks"] |= gone["blocks"]
                records[survivor]["span_start"] = min(
                    records[survivor]["span_start"], gone["span_start"]
                )
                records[survivor]["span_end"] = max(records[survivor]["span_end"], gone["span_end"])
                kmers.setdefault(survivor, set()).update(kmers.pop(copy, set()))
        sharing = kmer_sharing(kmers)
        matched: dict[int, list[str]] = {}
        programmes: dict[str, str] = {}
        for idx, entry in enumerate(CURATED_MARKERS):
            if entry["virus"] != virus:
                continue
            matches = _resolve(records, entry["refseq_gene"])
            if not matches:
                errors.append(f"{virus}: no gene in {gtf_name} matches {entry['refseq_gene']!r}")
                continue
            matched[idx] = sorted(matches)
            for gene_id in matches:
                if programmes.setdefault(gene_id, entry["programme"]) != entry["programme"]:
                    errors.append(f"{virus}: {gene_id} is curated into two programmes")
        if not programmes:
            continue
        groups = marker_groups(virus, records, programmes, sharing, args.kmer_threshold)
        report = audit_markers(programmes, groups, sharing, args.kmer_threshold)
        for pair in report["mixed_group"]:
            errors.append(f"{virus}: latent and productive markers share a group: {pair}")
        for pair in report["cross_programme_split"]:
            errors.append(
                f"{virus}: markers of different programmes share >= "
                f"{CROSS_PROGRAMME_FLOOR:.0%} of 31-mers but sit in different groups: {pair}"
            )
        for pair in report["same_programme_split"]:
            warnings.append(f"{virus}: linked same-programme markers in different groups: {pair}")
        linked_by_virus[virus] = {
            frozenset(pair)
            for pair, f in sharing.items()
            if f >= args.kmer_threshold and pair[0] in programmes and pair[1] in programmes
        }
        sensitivity += sensitivity_rows(virus, records, programmes, sharing)
        facts = VIRUS_FACTS[virus]
        for idx, gene_ids in matched.items():
            entry = CURATED_MARKERS[idx]
            refseq = entry["refseq_gene"]
            for gene_id in gene_ids:
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
                        "kinetic_class": entry.get("kinetic_class", "unclassified"),
                        "kinetic_pmid": entry.get("kinetic_pmid", ""),
                        "note": note,
                    }
                )

    if args.sensitivity:
        with open(args.sensitivity, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=SENSITIVITY_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(sensitivity)
        print(f"wrote {args.sensitivity} ({len(sensitivity)} rows)")

    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        print(
            f"\n{len(errors)} catalogue error(s); refusing to write a partial or inconsistent "
            "catalogue.",
            file=sys.stderr,
        )
        return 1

    bad = [r for r in rows if r["kinetic_class"] not in KINETIC_CLASSES]
    if bad:
        raise SystemExit(f"invalid kinetic_class: {bad[0]['refseq_gene']}")
    _fill_non_overlapping(rows, linked_by_virus)
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


def _fill_non_overlapping(
    rows: list[dict[str, object]], linked_by_virus: dict[str, set[frozenset[str]]]
) -> None:
    """A marker is non-overlapping iff no *other* marker of its virus shares its overlap group
    or is linked to it (>= the k-mer threshold). Groups are cliques, so a marker can be linked to
    a marker of a neighbouring group; that must not read as independent."""
    by_virus_group: dict[tuple[str, str], list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        by_virus_group[(str(row["virus"]), str(row["overlap_group"]))].append(idx)
    for indices in by_virus_group.values():
        for idx in indices:
            row = rows[idx]
            partners = {
                other
                for pair in linked_by_virus.get(str(row["virus"]), set())
                if row["gene_id_bundled"] in pair
                for other in pair
                if other != row["gene_id_bundled"]
            }
            siblings = {
                rows[i]["gene_id_bundled"]
                for i in indices
                if rows[i]["gene_id_bundled"] != row["gene_id_bundled"]
            }
            shared = bool(partners or siblings)
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
