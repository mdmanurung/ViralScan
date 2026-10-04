"""Anellovirus alignment branch that does not depend on kallisto (PLAN ANDET-09).

kallisto needs an exact 31-mer, and the panel's anellovirus set was clustered
at 95 % ANI, so a divergent strain shares few k-mers with its nearest reference
(F-013). The evidence replay (``viralscan evidence``) only re-aligns reads that
already sit in a kallisto equivalence class. This branch aligns *every*
host-unmapped read pair to the panel's own anellovirus genomes with STARsolo,
then reports read-level evidence per accession and per virus.

The evidence is a label, never a filter (ANELLO-PRIOR). Thresholds are plan
section 11 starting values; ANDET-09e calibrates them.

The pure helpers here are shared by ``scripts/build_bundled_panel_ref.py``
(index build), ``scripts/anello_align.py`` (Snakemake rule) and
``scripts/detection.py`` (merge into ``viral_summary.tsv``).
"""

from __future__ import annotations

import csv
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

#: Directory, next to the kb index, that holds the STAR anellovirus index.
INDEX_DIRNAME = "anello_star"

#: Per-accession table written by the Snakemake rule (run-relative path).
ACCESSION_TSV = "results/anello_alignment_by_accession.tsv"

#: Values of the ``alignment_status`` column.
STATUS_OK = "ok"
STATUS_DISABLED = "disabled"
STATUS_NO_HOST_FILTER = "skipped_no_host_filter"
STATUS_NO_INDEX = "skipped_no_anello_index"

#: A read with a homopolymer run this long is the poly-G / poly-A artefact
#: class (F-019 used >=15 nt).
HOMOPOLYMER_RUN = 15
_HOMOPOLYMER_RE = re.compile(r"A{%d,}|C{%d,}|G{%d,}|T{%d,}" % ((HOMOPOLYMER_RUN,) * 4))

#: The read-side artefact measures (F-019 update 2026-10-04). The pileup
#: position, the poly-A fraction and NM are each equally consistent with a
#: genuine 3'-end read, so none of them discriminates. These three do:
#:
#: * ``complex_body_fraction`` — a genuine 10x R2 read of a viral mRNA 3' end is
#:   ``[complex viral sequence][untemplated poly-A]``. A TSO/poly-A chimera has
#:   no complex viral body at all. This is the decisive measure.
#: * ``tso_fraction`` — in 10x 3' chemistry only fragments carrying the
#:   bead-oligo end are sequenced, so the TSO end of the cDNA is physically
#:   discarded. TSO sequence inside R2 cannot come from a genuine molecule.
#: * ``median_query_coverage`` — identity without the length it was measured
#:   over is not evidence: a soft-clipped 34 nt perfect match reports identity
#:   1.0 exactly like a full-length one.
#:
#: All three are reported, never used to drop a read. Anelloviruses are a
#: commensal virome, so a call is biologically expected (``ANELLO-PRIOR``), and
#: any homopolymer-based filter removes precisely the genuine 3'-end reads it
#: would take to prove one. Flagging bounds the artefact; filtering would hide a
#: real infection (user, 2026-10-04).

#: Shortest templated body that can carry usable evidence. The review used
#: mlen >= 25; this allows 20 because the body/tail boundary is not exact — a
#: genuine 3'UTR ending in a few A's merges them into the tail, costing up to
#: 6 nt of body (measured: 75 % lose nothing, worst case 6 nt over 5,000 draws).
#: At 25 that boundary effect alone discarded 22 % of genuine 25 nt bodies.
MIN_BODY_LEN = 20

#: Dinucleotide Shannon entropy, in bits, above which a body counts as complex.
#: Measured, not assumed. Over 5,000 random ACGT draws the minimum is 2.21 bits
#: at 20 nt and 2.49 at 25 nt, while the artefact classes sit at 0.00 (poly-A,
#: poly-G), 1.00 (AC repeat), 1.58 (CAG repeat) and 1.70 (A-rich). 2.0 sits in
#: the gap between 1.70 and 2.21, so no genuine body of any length >= 20 is lost
#: and every measured artefact class is still flagged.
#:
#: The review's 3.5 does not transfer: it was measured on 90 nt reads (median
#: 3.88), and at 25 nt the median genuine body is 3.49 — it would discard half
#: of them, which the "do not miss a real infection" setting rules out.
#:
#: Entropy alone does not catch the TSO (H = 3.56); ``has_tso`` does. That is
#: why both are reported.
MIN_BODY_ENTROPY = 2.0

#: 10x template-switch oligo, screened in both orientations. The 25 nt core is
#: the part that appears in chimeric reads; ``ATGGG`` completes the full oligo.
TSO = "AAGCAGTGGTATCAACGCAGAGTAC"
TSO_MAX_MISMATCH = 2

#: STAR genomeGenerate settings for a ~6 Mb, ~2,000-contig reference:
#: SAindexNbases = min(14, log2(6e6)/2 - 1) ~= 10; ChrBinNbits =
#: log2(genome / n_contigs) ~= 11 (plan section 10).
GENOME_GENERATE_ARGS: tuple[str, ...] = (
    "--genomeSAindexNbases",
    "10",
    "--genomeChrBinNbits",
    "11",
)

#: Alignment settings (plan section 11, single-end: the cDNA read is the only
#: aligned read in 10x, so the mate-gap option is dropped). Both multimap caps
#: are 100 so every secondary record reaches the BAM and NH is the true count.
ALIGN_ARGS: tuple[str, ...] = (
    "--twopassMode",
    "Basic",
    "--alignIntronMax",
    "5000",
    "--outFilterMultimapNmax",
    "100",
    "--outSAMmultNmax",
    "100",
    "--winAnchorMultimapNmax",
    "200",
    "--outFilterMismatchNmax",
    "999",
    "--outFilterMismatchNoverLmax",
    "0.08",
    "--outFilterMatchNminOverLread",
    "0.80",
    "--outFilterScoreMinOverLread",
    "0.80",
    "--outSAMattributes",
    "NH",
    "HI",
    "AS",
    "nM",
    "NM",
    "MD",
    "jM",
    "jI",
    "CB",
    "UB",
    "--outSAMtype",
    "BAM",
    "SortedByCoordinate",
    # Strand: the forward-strand default loses most 5' reads (F-020).
    "--soloStrand",
    "Unstranded",
    "--soloFeatures",
    "GeneFull",
    "--soloMultiMappers",
    "Unique",
    "EM",
)

ACCESSION_COLUMNS: tuple[str, ...] = (
    "accession",
    "length",
    "reads",
    "unique_reads",
    "weighted_reads",
    "median_nh",
    "median_identity",
    "median_query_coverage",
    "complex_body_fraction",
    "tso_fraction",
    "breadth",
    "breadth_unique",
    "start_sites",
    "homopolymer_fraction",
    "splice_reads",
    "sense_fraction",
)

#: Columns merged into viral_summary.tsv. NA on rows that are not anellovirus.
SUMMARY_COLUMNS: tuple[str, ...] = (
    "detection_source",
    "alignment_status",
    "alignment_reads",
    "alignment_unique_reads",
    "alignment_molecules_unique",
    "alignment_cells_unique",
    "alignment_median_identity",
    "alignment_median_query_coverage",
    "alignment_complex_body_fraction",
    "alignment_tso_fraction",
    "alignment_accessions",
    "alignment_start_sites",
    "alignment_homopolymer_fraction",
    "alignment_splice_reads",
)


# ── Index ─────────────────────────────────────────────────────────────────────
def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    """Yield ``(first header token, sequence)`` from a plain FASTA."""
    name: Optional[str] = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def write_star_reference(
    viral_fasta: Path, accessions: set[str], out_dir: Path
) -> tuple[Path, Path, int]:
    """Write the anellovirus FASTA and a one-gene-per-contig GTF into *out_dir*.

    A whole-contig gene is deliberate: the panel's anellovirus gene models are
    CDS-only and stop short of polyA (learning 2026-10-01), so a GTF built from
    them would drop the 3' reads STARsolo should count.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / "anello.fa"
    gtf = out_dir / "anello.gtf"
    bare = {a.split(".")[0] for a in accessions}
    n = 0
    with open(fasta, "w", encoding="utf-8") as fa, open(gtf, "w", encoding="utf-8") as gt:
        for name, seq in iter_fasta(viral_fasta):
            if name not in accessions and name.split(".")[0] not in bare:
                continue
            fa.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i : i + 80] + "\n")
            attrs = f'gene_id "{name}"; transcript_id "{name}"; gene_name "{name}";'
            for feature in ("gene", "transcript", "exon"):
                gt.write(f"{name}\tviralscan\t{feature}\t1\t{len(seq)}\t.\t+\t.\t{attrs}\n")
            n += 1
    if n == 0:
        raise ValueError(f"no anellovirus records found in {viral_fasta}")
    return fasta, gtf, n


def genome_generate_cmd(star: str, out_dir: Path, fasta: Path, gtf: Path, threads: int) -> list[str]:
    """STAR genomeGenerate command for the anellovirus reference (pure)."""
    return [
        star,
        "--runMode",
        "genomeGenerate",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(out_dir),
        "--genomeFastaFiles",
        str(fasta),
        "--sjdbGTFfile",
        str(gtf),
        "--sjdbOverhang",
        "100",
        *GENOME_GENERATE_ARGS,
    ]


def resolve_index(kb_index: str) -> Optional[str]:
    """Return ``<kb index dir>/anello_star`` when a built STAR index is there."""
    candidate = Path(kb_index).resolve().parent / INDEX_DIRNAME
    return str(candidate) if (candidate / "SA").is_file() else None


def align_cmd(
    star: str,
    index: str,
    cdna: str,
    barcode: str,
    barcode_args: Sequence[str],
    prefix: str,
    threads: int,
) -> list[str]:
    """STARsolo alignment command; cDNA read first, barcode read second (pure)."""
    read_cmd = "zcat" if cdna.endswith(".gz") else "-"
    return [
        star,
        "--runThreadN",
        str(threads),
        "--genomeDir",
        index,
        "--readFilesIn",
        cdna,
        barcode,
        "--readFilesCommand",
        read_cmd,
        "--outFileNamePrefix",
        prefix,
        *barcode_args,
        *ALIGN_ARGS,
    ]


def status_for(enabled: bool, host_index: Optional[str], anello_index: Optional[str]) -> str:
    """The run-level ``alignment_status`` (pure)."""
    if not enabled:
        return STATUS_DISABLED
    if not host_index:
        return STATUS_NO_HOST_FILTER
    if not anello_index:
        return STATUS_NO_INDEX
    return STATUS_OK


# ── SAM parsing ───────────────────────────────────────────────────────────────
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class Alignment:
    qname: str
    flag: int
    rname: str
    pos: int  # 1-based leftmost
    cigar: str
    seq: str
    tags: dict[str, str]

    @property
    def secondary(self) -> bool:
        return bool(self.flag & 0x100)

    @property
    def reverse(self) -> bool:
        return bool(self.flag & 0x10)

    @property
    def nh(self) -> int:
        return int(self.tags.get("NH", "1"))

    def ops(self) -> list[tuple[int, str]]:
        return [(int(n), op) for n, op in _CIGAR_RE.findall(self.cigar)]

    def aligned_bases(self) -> int:
        return sum(n for n, op in self.ops() if op in "MI=X")

    def ref_blocks(self) -> list[tuple[int, int]]:
        """0-based half-open reference intervals covered by aligned bases."""
        blocks, ref = [], self.pos - 1
        for n, op in self.ops():
            if op in "M=X":
                blocks.append((ref, ref + n))
                ref += n
            elif op in "DN":
                ref += n
        return blocks

    def identity(self) -> Optional[float]:
        """Edit-distance identity ``1 - NM / aligned`` (NM counts indels too).

        Normalised by *aligned* length, so this says nothing about how much of
        the read took part: a soft-clipped 34 nt perfect match scores 1.0 just
        like a full-length one. Always read it beside ``query_coverage``.
        """
        aligned = self.aligned_bases()
        if "NM" not in self.tags or aligned == 0:
            return None
        return 1.0 - int(self.tags["NM"]) / aligned

    def read_seq(self) -> str:
        """SEQ put back into sequencing orientation.

        SAM stores SEQ reverse-complemented on a reverse-strand record, which
        silently inverts every read-side measure: a genuine ``[body][poly-A]``
        arrives as ``[poly-T][rc body]``, so its body reads as empty and it is
        flagged an artefact, while a ``[poly-A][TSO]`` chimera arrives with the
        TSO in front and reads as complex. With ``--soloStrand Unstranded`` that
        is about half of all records. Anything looking at read *sequence* must
        go through here; NH, NM and the CIGAR are orientation-free.
        """
        return _revcomp(self.seq) if self.reverse else self.seq

    def query_coverage(self) -> Optional[float]:
        """Fraction of the read that took part in the alignment.

        The soft-clipped remainder is where the falsifying evidence hides — in
        the covid case the entire TSO sat in it (F-019, update 2026-10-04).
        """
        if not self.seq or self.seq == "*":
            return None
        return min(1.0, self.aligned_bases() / len(self.seq))

    def spliced(self) -> bool:
        return "N" in self.cigar


def parse_sam_line(line: str) -> Optional[Alignment]:
    """Parse one SAM body line; ``None`` for headers and unmapped records."""
    if line.startswith("@"):
        return None
    f = line.rstrip("\n").split("\t")
    flag = int(f[1])
    if flag & 0x4:
        return None
    tags = {}
    for t in f[11:]:
        key, _, value = t.split(":", 2)
        tags[key] = value
    return Alignment(f[0], flag, f[2], int(f[3]), f[5], f[9], tags)


def has_homopolymer(seq: str) -> bool:
    return bool(_HOMOPOLYMER_RE.search(seq.upper()))


# ── Read-side artefact measures ───────────────────────────────────────────────
def read_body(seq: str) -> str:
    """The templated part of a read: everything 5' of its first long run.

    A genuine 3'-end read is ``[body][untemplated poly-A]``; an artefact is a
    long run with nothing informative in front of it. Splitting at the run is
    what lets the two be told apart — on the body, not on the run.
    """
    match = _HOMOPOLYMER_RE.search(seq.upper())
    return seq[: match.start()] if match else seq


def dinucleotide_entropy(seq: str) -> float:
    """Shannon entropy in bits over the 16 dinucleotides (0 for a homopolymer)."""
    pairs = [seq[i : i + 2] for i in range(len(seq) - 1)]
    if not pairs:
        return 0.0
    total = len(pairs)
    return -sum(
        (n / total) * math.log2(n / total) for n in Counter(pairs).values()
    )


def is_complex_body(seq: str) -> bool:
    """Does this read carry enough complex templated sequence to mean anything?

    Note this is deliberately *not* a virus check: it asks whether the read
    could carry evidence at all. Whether that body is viral is settled by the
    alignment, which is why the two are reported as separate columns.
    """
    body = read_body(seq)
    return len(body) >= MIN_BODY_LEN and dinucleotide_entropy(body) >= MIN_BODY_ENTROPY


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _contains(haystack: str, needle: str, max_mismatch: int) -> bool:
    for i in range(len(haystack) - len(needle) + 1):
        if sum(a != b for a, b in zip(haystack[i:], needle)) <= max_mismatch:
            return True
    return False


def has_tso(seq: str) -> bool:
    """Is the 10x TSO present in this read, in either orientation?

    It cannot be there in a genuine molecule: 10x 3' chemistry only sequences
    fragments carrying the bead-oligo end, so the TSO end of the cDNA is
    discarded. Its presence is positive evidence of a chimera.
    """
    upper = seq.upper()
    return _contains(upper, TSO, TSO_MAX_MISMATCH) or _contains(
        upper, _revcomp(TSO), TSO_MAX_MISMATCH
    )


def _valid(tag: Optional[str]) -> bool:
    return bool(tag) and tag != "-"


# ── Per-accession metrics ─────────────────────────────────────────────────────
@dataclass
class _Acc:
    reads: set = field(default_factory=set)
    unique_reads: int = 0
    weighted: float = 0.0
    nh: list = field(default_factory=list)
    identity: list = field(default_factory=list)
    covered: set = field(default_factory=set)
    covered_unique: set = field(default_factory=set)
    starts: set = field(default_factory=set)
    homopolymer: int = 0
    splice: int = 0
    sense: int = 0
    coverage: list = field(default_factory=list)
    complex_body: int = 0
    tso: int = 0


def accession_metrics(
    alignments: Iterable[Alignment], lengths: Mapping[str, int]
) -> list[dict[str, object]]:
    """Per-accession evidence. Reads are counted once per accession (by QNAME).

    *weighted_reads* is the sum of 1/NH over records, so a read placed on k
    accessions contributes 1/k to each; it describes ambiguity, it is not EM.
    Breadth uses every record, and *breadth_unique* only NH == 1 records, so
    100 secondary placements cannot make every related genome look covered.
    """
    acc: dict[str, _Acc] = defaultdict(_Acc)
    for a in alignments:
        s = acc[a.rname]
        s.weighted += 1.0 / a.nh
        positions = [p for lo, hi in a.ref_blocks() for p in range(lo, hi)]
        s.covered.update(positions)
        if a.nh == 1:
            s.covered_unique.update(positions)
        if a.qname in s.reads:
            continue
        s.reads.add(a.qname)
        s.nh.append(a.nh)
        s.unique_reads += a.nh == 1
        ident = a.identity()
        if ident is not None:
            s.identity.append(ident)
        s.starts.add(a.pos)
        s.homopolymer += has_homopolymer(a.seq)
        s.splice += a.spliced()
        s.sense += not a.reverse
        cov = a.query_coverage()
        if cov is not None:
            s.coverage.append(cov)
        read = a.read_seq()  # never a.seq: see Alignment.read_seq
        s.complex_body += is_complex_body(read)
        s.tso += has_tso(read)

    rows = []
    for name in sorted(acc):
        s = acc[name]
        length = lengths.get(name, 0)
        rows.append(
            {
                "accession": name,
                "length": length,
                "reads": len(s.reads),
                "unique_reads": s.unique_reads,
                "weighted_reads": round(s.weighted, 4),
                "median_nh": statistics.median(s.nh),
                "median_identity": round(statistics.median(s.identity), 4)
                if s.identity
                else "",
                "median_query_coverage": round(statistics.median(s.coverage), 4)
                if s.coverage
                else "",
                "complex_body_fraction": round(s.complex_body / len(s.reads), 4)
                if s.reads
                else "",
                "tso_fraction": round(s.tso / len(s.reads), 4) if s.reads else "",
                "breadth": round(len(s.covered) / length, 4) if length else "",
                "breadth_unique": round(len(s.covered_unique) / length, 4) if length else "",
                "start_sites": len(s.starts),
                "homopolymer_fraction": round(s.homopolymer / len(s.reads), 4),
                "splice_reads": s.splice,
                "sense_fraction": round(s.sense / len(s.reads), 4),
            }
        )
    return rows


def virus_molecules(
    alignments: Iterable[Alignment], acc_to_virus: Mapping[str, str]
) -> dict[str, dict[str, int]]:
    """Molecules and cells per virus, from reads whose every placement is one virus.

    A read placed on two accessions of the same genus is still unique to that
    genus, so this is the genus-level analogue of STARsolo's gene-unique count.
    Molecules are distinct (CB, UB). Reads STARsolo left without a CB or UB
    (``-``: barcode off the list, or a filtered UMI such as a homopolymer) are
    skipped.
    """
    viruses: dict[str, set[str]] = defaultdict(set)
    keys: dict[str, tuple[str, str]] = {}
    for a in alignments:
        viruses[a.qname].add(acc_to_virus.get(a.rname, a.rname))
        cb, ub = a.tags.get("CB"), a.tags.get("UB")
        if _valid(cb) and _valid(ub):
            keys[a.qname] = (cb, ub)
    molecules: dict[str, set] = defaultdict(set)
    for qname, vs in viruses.items():
        if len(vs) == 1 and qname in keys:
            molecules[next(iter(vs))].add(keys[qname])
    return {
        v: {"molecules": len(m), "cells": len({cb for cb, _ in m})}
        for v, m in molecules.items()
    }


# A table written before these columns existed simply lacks them, so both
# helpers treat a missing column as "not measured" rather than failing. That
# keeps a resumed run readable instead of crashing on its own older output.
def _read_weighted(rows: Sequence[Mapping[str, object]], column: str, reads: int) -> object:
    """Read-weighted mean of a per-accession fraction; "" when nothing was measured."""
    present = [r for r in rows if r.get(column, "") != ""]
    if not reads or not present:
        return ""
    weighed = sum(int(r["reads"]) for r in present)
    if not weighed:
        return ""
    return round(
        sum(float(r[column]) * int(r["reads"]) for r in present) / weighed, 4
    )


def _median_of(rows: Sequence[Mapping[str, object]], column: str) -> object:
    values = [float(r[column]) for r in rows if r.get(column, "") != ""]
    return round(statistics.median(values), 4) if values else ""


def virus_summary(
    acc_rows: Sequence[Mapping[str, object]],
    molecules: Mapping[str, Mapping[str, int]],
    acc_to_virus: Mapping[str, str],
) -> dict[str, dict[str, object]]:
    """Aggregate per-accession rows to virus-level ``alignment_*`` columns.

    Fractions are weighted by reads, not averaged over accessions, so one
    accession with three reads cannot outvote one with three thousand.
    """
    groups: dict[str, list] = defaultdict(list)
    for r in acc_rows:
        groups[acc_to_virus.get(str(r["accession"]), str(r["accession"]))].append(r)
    out = {}
    for virus in sorted(set(groups) | set(molecules)):
        rows = groups.get(virus, [])
        reads = sum(int(r["reads"]) for r in rows)
        idents = [float(r["median_identity"]) for r in rows if r["median_identity"] != ""]
        mol = molecules.get(virus, {})
        out[virus] = {
            "alignment_reads": reads,
            "alignment_unique_reads": sum(int(r["unique_reads"]) for r in rows),
            "alignment_molecules_unique": mol.get("molecules", 0),
            "alignment_cells_unique": mol.get("cells", 0),
            "alignment_median_identity": round(statistics.median(idents), 4) if idents else "",
            "alignment_median_query_coverage": _median_of(rows, "median_query_coverage"),
            "alignment_accessions": len(rows),
            "alignment_start_sites": sum(int(r["start_sites"]) for r in rows),
            # Read-weighted across accessions.
            "alignment_homopolymer_fraction": _read_weighted(
                rows, "homopolymer_fraction", reads
            ),
            "alignment_complex_body_fraction": _read_weighted(
                rows, "complex_body_fraction", reads
            ),
            "alignment_tso_fraction": _read_weighted(rows, "tso_fraction", reads),
            "alignment_splice_reads": sum(int(r["splice_reads"]) for r in rows),
        }
    return out


# ── I/O ───────────────────────────────────────────────────────────────────────
def write_accession_tsv(rows: Sequence[Mapping[str, object]], path: Path, meta: Mapping) -> None:
    """Write the per-accession table with a one-line ``# {json}`` provenance header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("# " + json.dumps(dict(meta), sort_keys=True) + "\n")
        writer = csv.DictWriter(fh, fieldnames=ACCESSION_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_accession_tsv(path: Path) -> tuple[dict, list[dict[str, str]]]:
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
        meta = json.loads(first[2:]) if first.startswith("# ") else {}
        if not first.startswith("# "):
            fh.seek(0)
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return meta, rows


# ── Merge into viral_summary.tsv ──────────────────────────────────────────────
def merge_summary_rows(
    rows: list[dict[str, object]],
    evidence: Mapping[str, Mapping[str, object]],
    anello_viruses: set[str],
    status: str,
    row_template: Mapping[str, object],
) -> list[dict[str, object]]:
    """Add ``SUMMARY_COLUMNS`` to every row, and rows for alignment-only viruses.

    * Non-anellovirus rows: ``detection_source = kallisto``, other columns NA.
    * Anellovirus rows: evidence columns filled (0 when the branch ran and saw
      nothing); ``kallisto+alignment`` when it holds >= 1 unique molecule.
    * An anellovirus with >= 1 unique alignment molecule but no kallisto row gets
      a new row, ``detection_source = alignment_only``, with *row_template*
      supplying the run-level fields (denominators) and kallisto counts at 0.
    """
    ran = status == STATUS_OK
    zero = {c: 0 for c in SUMMARY_COLUMNS[2:]}
    seen = set()
    out = []
    for row in rows:
        virus = str(row["virus_name"])
        row = dict(row)
        if virus in anello_viruses:
            seen.add(virus)
            ev = evidence.get(virus, zero) if ran else {}
            row.update({c: ev.get(c, "") for c in SUMMARY_COLUMNS[2:]})
            hit = ran and int(ev.get("alignment_molecules_unique", 0)) >= 1
            row["detection_source"] = "kallisto+alignment" if hit else "kallisto"
            row["alignment_status"] = status
        else:
            row.update({c: "" for c in SUMMARY_COLUMNS})
            row["detection_source"] = "kallisto"
        out.append(row)
    if ran:
        for virus in sorted(evidence):
            ev = evidence[virus]
            if virus in seen or int(ev.get("alignment_molecules_unique", 0)) < 1:
                continue
            new = dict(row_template)
            new["virus_name"] = virus
            new.update({c: ev.get(c, "") for c in SUMMARY_COLUMNS[2:]})
            new["detection_source"] = "alignment_only"
            new["alignment_status"] = status
            out.append(new)
    return out
