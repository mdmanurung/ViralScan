"""Splice junctions in the reads behind a viral call (PLAN ``EVID-02``, ``viralscan evidence --splice``).

``viralscan evidence`` aligns with ``minimap2 -ax sr``, which never splits a read, so its BAMs hold no
junctions. This module realigns the **member reads of counted virus-best molecules** to the viral genome alone
with ``minimap2 -ax splice:sr`` and tallies every intron (an ``N`` in a CIGAR). The unit is fixed in
``docs/evidence_unit.md``: a molecule is a corrected (CB, UMI) of ``molecule_verdicts.tsv.gz``; its member
reads are the ``read_lineage.tsv.gz`` rows with that CB and UB, joined to the competitive BAM on
``read_number``. ``viral_reads.umi_dedup.fasta`` is never an input: its representative read favours the virus.

A molecule supports a junction when any member read does. Read counts are reported beside molecule counts.

The alignment is viral-only, so a junction supports viral RNA processing. It does not establish the origin of
the molecule, infection or reactivation: ambient RNA and integrated genomes (ciHHV-6) also give spliced viral
reads (Peddu 2019, doi:10.1128/JVI.01418-19).

Outputs, in ``<evidence out>/splice/``: ``splice.bam`` (+ ``.bai``), ``junctions.tsv``, ``splice_summary.tsv``,
``target_genome.fa`` and ``splice_manifest.json``.
"""

import csv
import gzip
import json
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any

from viralscan.evidence import _run
from viralscan.molecule_verdict import VERDICTS, _read_number
from viralscan.run_safety import sha256_file
from viralscan.utils import open_maybe_gzip

# ---------------------------------------------------------------------------
# Junction rules. Fixed before any real data were looked at; change them only with a PLAN row.
# ---------------------------------------------------------------------------

#: The ``splice:sr`` preset appeared in minimap2 2.29; older versions have no short-read splice preset and
#: plain ``splice`` maps none of the 64 bp HSV-1 reads, so there is no fallback.
MIN_MINIMAP2_VERSION = (2, 29)
#: Longest intron minimap2 may open (``-G``). EBV latency transcripts have introns of tens of kb.
MAX_INTRON_LENGTH = 100_000
#: Shortest intron kept. Shorter N gaps are deletions that minimap2 chose to call splices.
MIN_INTRON_LENGTH = 20
#: Aligned bases (M/=/X) needed on each side of an N, up to the next N or the read end. An anchor of 4-7 bp
#: matches by chance somewhere in a 100 kb window, so a spliced call on it is not evidence.
MIN_ANCHOR = 8
#: Largest NM / aligned query bases of an alignment that may report a junction. A read with many mismatches
#: is more likely to be a poor alignment placed across a gap than a spliced read.
MAX_MISMATCH_FRACTION = 0.08
#: Largest soft-clipped share of the read for an alignment that may report a junction. A read that is
#: mostly clipped is mostly unexplained, and a short aligned part carries the junction by chance.
MAX_CLIP_FRACTION = 0.25
#: Extra SAM FLAG bits that disqualify a record: secondary 0x100, supplementary 0x800. Unmapped 0x4 has no CIGAR.
SKIP_FLAGS = 0x904

MINIMAP2_SPLICE_OPTIONS: tuple[str, ...] = (
    "-ax",
    "splice:sr",
    "--secondary=no",
    "-G",
    str(MAX_INTRON_LENGTH),
)

#: Intron end dinucleotides, in transcript orientation, that spliceosomes use (U2 GT-AG, GC-AG; U12 AT-AC).
CANONICAL_MOTIFS = frozenset({"GT-AG", "GC-AG", "AT-AC"})

VERDICTS_FILE = "molecule_verdicts.tsv.gz"
LINEAGE_FILE = "read_lineage.tsv.gz"
RAW_BAM = "competitive_reads.raw.bam"
SPLICE_DIR = "splice"

NOTE = (
    "Viral-only alignment of member reads of counted virus-best molecules. A junction supports viral RNA "
    "processing. It does not establish the origin of the molecule, infection or reactivation."
)

#: ``annotation_status`` values.
ANNOT_AVAILABLE = "annotated_introns_available"
ANNOT_NO_INTRONS = "no_annotated_introns"
ANNOT_NOT_EVALUABLE = "intron_matching_not_evaluable"
ANNOT_PARTIAL = "partly_evaluable"

JUNCTION_FIELDS = (
    "contig",
    "start",
    "end",
    "strand",
    "motif",
    "motif_canonical",
    "supporting_molecules",
    "supporting_reads",
    "annotated",
    "annotated_transcript",
    "counted_in_spliced_molecules",
)
CATEGORY_FIELDS = (
    "counted_virus_best",
    "counted_host_best",
    "counted_tie",
    "counted_unaligned",
    "not_counted",
)
SUMMARY_FIELDS = (
    "status",
    "reason",
    *CATEGORY_FIELDS,
    "reads_realigned",
    "reads_aligned",
    "spliced_reads",
    "spliced_molecules",
    "spliced_molecules_annotated",
    "unique_junctions",
    "unique_junctions_excluded_noncanonical",
    "molecules_excluded_noncanonical",
    "junction_calls_rejected",
    "annotation_status",
    "evidence_manifest_sha256",
    "note",
)


class SpliceError(ValueError):
    """The splice check cannot run as specified (old minimap2, inconsistent inputs)."""


# ---------------------------------------------------------------------------
# minimap2 version
# ---------------------------------------------------------------------------


def parse_minimap2_version(text: str) -> tuple[int, int]:
    """``(major, minor)`` from ``minimap2 --version`` output such as ``2.31-r1307``."""
    match = re.match(r"\s*(\d+)\.(\d+)", text)
    if match is None:
        raise SpliceError(f"Cannot read the minimap2 version from {text!r}")
    return int(match.group(1)), int(match.group(2))


def minimap2_version() -> str:
    """Raw ``minimap2 --version`` string; ``SpliceError`` when it is not installed."""
    try:
        proc = subprocess.run(["minimap2", "--version"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SpliceError(f"--splice needs minimap2 on PATH: {exc}") from exc
    return proc.stdout.strip()


def require_splice_minimap2() -> str:
    """Version string of a minimap2 that has ``splice:sr``; refuse an older one, never fall back."""
    version = minimap2_version()
    if parse_minimap2_version(version) < MIN_MINIMAP2_VERSION:
        need = ".".join(map(str, MIN_MINIMAP2_VERSION))
        raise SpliceError(
            f"--splice needs minimap2 >= {need} (the splice:sr preset); found {version}. "
            "Upgrade minimap2. There is no fallback to plain `splice`, which maps no short reads."
        )
    return version


# ---------------------------------------------------------------------------
# Annotated introns
# ---------------------------------------------------------------------------


def annotated_introns(
    gtfs: Iterable[Path], contigs: Iterable[str]
) -> tuple[dict[str, dict[tuple[int, int], list[str]]], set[str]]:
    """``(contig -> {(start, end): transcripts}, contigs with exon rows)`` from GTFs.

    An intron is the gap between two consecutive ``exon`` rows of one ``transcript_id`` (1-based, inclusive).
    CDS rows are never used: the gap between two CDS rows of one gene can be a UTR exon, an unannotated
    exon or a real intron. A contig with no exon row is not evaluable; one with exon rows and no gap
    has an empty intron set.
    """
    wanted = set(contigs)
    exons: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for path in gtfs:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or line.startswith("#") or fields[2] != "exon":
                    continue
                if fields[0] not in wanted:
                    continue
                attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
                if "transcript_id" not in attrs:
                    continue
                exons[fields[0]][attrs["transcript_id"]].append((int(fields[3]), int(fields[4])))
    introns: dict[str, dict[tuple[int, int], list[str]]] = {}
    for contig, by_transcript in exons.items():
        found: dict[tuple[int, int], list[str]] = defaultdict(list)
        for transcript, rows in by_transcript.items():
            rows.sort()
            for (_, left_end), (right_start, _) in zip(rows[:-1], rows[1:], strict=True):
                if right_start - left_end > 1:
                    found[left_end + 1, right_start - 1].append(transcript)
        introns[contig] = {k: sorted(v) for k, v in found.items()}
    return introns, set(exons)


# ---------------------------------------------------------------------------
# Junctions from one alignment
# ---------------------------------------------------------------------------

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def junctions_of(
    flag: int, pos: int, cigar: str, nm: int | None, read_length: int
) -> tuple[list[tuple[int, int]], Counter[str]]:
    """Accepted ``(start, end)`` introns (1-based, inclusive) of one SAM record, and why others were refused.

    Only a primary, non-supplementary, mapped record may report a junction. A record that fails the
    mismatch or clip rule reports none; an N that fails the anchor or length rule is refused alone.
    """
    rejected: Counter[str] = Counter()
    if flag & (SKIP_FLAGS | 0x4) or cigar == "*":
        return [], rejected
    ops = [(int(n), op) for n, op in _CIGAR.findall(cigar)]
    if not any(op == "N" for _, op in ops):
        return [], rejected
    n_introns = sum(op == "N" for _, op in ops)
    aligned = sum(n for n, op in ops if op in "M=X")
    clipped = sum(n for n, op in ops if op in "SH")
    if nm is not None and aligned and nm / aligned > MAX_MISMATCH_FRACTION:
        rejected["mismatch_fraction"] += n_introns
        return [], rejected
    if read_length and clipped / read_length > MAX_CLIP_FRACTION:
        rejected["clip_fraction"] += n_introns
        return [], rejected
    # M/=/X bases between neighbouring Ns (or the read ends) are the anchors
    blocks, current = [], 0
    for n, op in ops:
        if op == "N":
            blocks.append(current)
            current = 0
        elif op in "M=X":
            current += n
    blocks.append(current)
    accepted: list[tuple[int, int]] = []
    ref, k = pos, 0
    for n, op in ops:
        if op == "N":
            start, end = ref, ref + n - 1
            if min(blocks[k], blocks[k + 1]) < MIN_ANCHOR:
                rejected["anchor"] += 1
            elif not MIN_INTRON_LENGTH <= n <= MAX_INTRON_LENGTH:
                rejected["intron_length"] += 1
            else:
                accepted.append((start, end))
            k += 1
        if op in "MDN=X":
            ref += n
    return accepted, rejected


def motif_of(genome: str, start: int, end: int) -> tuple[str, bool, str]:
    """``(motif, canonical, strand)`` of the intron ``start..end`` (1-based) in transcript orientation.

    The forward reading is ``donor-acceptor`` of the genome; the reverse reading is the reverse complement
    of the intron's two ends. A canonical reading wins; with none, the forward reading is reported with
    strand ``.``.
    """
    head, tail = genome[start - 1 : start + 1].upper(), genome[end - 2 : end].upper()
    forward = f"{head}-{tail}"
    reverse = f"{revcomp(tail)}-{revcomp(head)}"
    if forward in CANONICAL_MOTIFS:
        return forward, True, "+"
    if reverse in CANONICAL_MOTIFS:
        return reverse, True, "-"
    return forward, False, "."


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def _molecule_categories(path: Path) -> tuple[Counter[str], set[tuple[str, str]]]:
    """Counts of the five categories of ``docs/evidence_unit.md`` and the counted virus-best molecules."""
    counts: Counter[str] = Counter()
    members: set[tuple[str, str]] = set()
    seen: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["cb"], row["ub"])
            if key in seen:
                raise SpliceError(f"{path} lists molecule {key} twice")
            seen.add(key)
            if row["counted"] != "True":
                counts["not_counted"] += 1
                continue
            if row["verdict"] not in VERDICTS:
                raise SpliceError(f"{path}: unknown verdict {row['verdict']!r}")
            counts[f"counted_{row['verdict']}"] += 1
            if row["verdict"] == "virus_best":
                members.add(key)
    return counts, members


def _member_reads(path: Path, molecules: set[tuple[str, str]]) -> dict[int, tuple[str, str]]:
    """``read_number -> (cb, ub)`` of every lineage row of a selected molecule."""
    reads: dict[int, tuple[str, str]] = {}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["cb"], row["ub"])
            if key in molecules:
                reads[int(row["read_number"])] = key
    return reads


def _stream_bam(bam: Path) -> Iterator[list[str]]:
    """SAM fields of every non-supplementary record, streamed from ``samtools view``."""
    with subprocess.Popen(
        ["samtools", "view", "-F", "0x800", str(bam)], stdout=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line.rstrip("\n").split("\t")
    if proc.returncode:
        raise SpliceError(f"samtools view failed (exit {proc.returncode}) on {bam}")


def _read_sequences(
    bam: Path, reads: dict[int, tuple[str, str]]
) -> tuple[dict[str, str], set[str]]:
    """``(qname -> sequence as sequenced, viral contigs the reads aligned to)`` for the member reads.

    The sequence comes from the primary record (reverse-complemented back when it is stored reversed).
    Secondary records only name the viral contigs. Supplementary records are never read.
    """
    sequences: dict[str, str] = {}
    contigs: set[str] = set()
    for fields in _stream_bam(bam):
        if _read_number(fields[0]) not in reads:
            continue
        flag, rname = int(fields[1]), fields[2]
        if rname.startswith("VIRUS|"):
            contigs.add(rname.removeprefix("VIRUS|"))
        if flag & 0x100:
            continue
        seq = fields[9]
        if seq == "*":
            raise SpliceError(f"Primary record of {fields[0]} has no sequence")
        sequences[fields[0]] = revcomp(seq) if flag & 0x10 else seq
    missing = len(reads) - len(sequences)
    if missing:
        raise SpliceError(
            f"{missing} member reads of {LINEAGE_FILE} have no primary record in {bam}"
        )
    return sequences, contigs


def _write_target_genome(viral_fasta: str, contigs: set[str], out: Path) -> dict[str, str]:
    """Write the named contigs of *viral_fasta* (first header token) to *out*; return their sequences."""
    genome: dict[str, str] = {}
    name: str | None = None
    with open_maybe_gzip(viral_fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in genome:
                    raise SpliceError(f"Duplicate FASTA ID {name!r} in {viral_fasta}")
                if name in contigs:
                    genome[name] = ""
            elif name in genome:
                genome[name] += line.strip()
    absent = contigs - set(genome)
    if absent:
        raise SpliceError(f"{sorted(absent)} are in the competitive BAM but not in {viral_fasta}")
    out.write_text("".join(f">{n}\n{genome[n]}\n" for n in sorted(genome)), encoding="utf-8")
    return genome


# ---------------------------------------------------------------------------
# Alignment and tally
# ---------------------------------------------------------------------------


def _align(reads_fasta: Path, genome_fasta: Path, out_bam: Path, threads: int) -> None:
    """minimap2 splice:sr -> sorted, indexed BAM (list-form subprocesses, errors surfaced)."""
    with tempfile.TemporaryFile() as sam:
        _run(
            [
                "minimap2",
                *MINIMAP2_SPLICE_OPTIONS,
                "-t",
                str(threads),
                str(genome_fasta),
                str(reads_fasta),
            ],
            stdout_file=sam,
        )
        sam.seek(0)
        _run(["samtools", "sort", "-@", str(threads), "-o", str(out_bam), "-"], stdin_file=sam)
    _run(["samtools", "index", str(out_bam)])


def _tally(
    bam: Path, reads: dict[int, tuple[str, str]]
) -> tuple[int, int, dict[tuple[str, int, int], dict[str, set[Any]]], Counter[str]]:
    """``(reads aligned, reads spliced, junction -> {molecules, reads}, rejections)`` from the splice BAM."""
    n_aligned = 0
    junctions: dict[tuple[str, int, int], dict[str, set[Any]]] = defaultdict(
        lambda: {"molecules": set(), "reads": set()}
    )
    spliced: set[str] = set()
    rejected: Counter[str] = Counter()
    for fields in _stream_bam(bam):
        flag = int(fields[1])
        if flag & SKIP_FLAGS or flag & 0x4:
            continue
        n_aligned += 1
        nm = next((int(t[5:]) for t in fields[11:] if t.startswith("NM:i:")), None)
        accepted, why = junctions_of(flag, int(fields[3]), fields[5], nm, len(fields[9]))
        rejected.update(why)
        number = _read_number(fields[0])
        assert number is not None
        for start, end in accepted:
            entry = junctions[fields[2], start, end]
            entry["molecules"].add(reads[number])
            entry["reads"].add(fields[0])
            spliced.add(fields[0])
    return n_aligned, len(spliced), junctions, rejected


# ---------------------------------------------------------------------------
# Public entry
# ---------------------------------------------------------------------------


def _write_tsv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _finish(
    out: Path,
    summary: dict[str, object],
    junction_rows: Sequence[dict[str, object]],
    manifest: dict[str, Any],
) -> dict[str, object]:
    """Write ``junctions.tsv``, ``splice_summary.tsv`` and the receipt (hashes of every other output)."""
    _write_tsv(out / "junctions.tsv", junction_rows, JUNCTION_FIELDS)
    _write_tsv(out / "splice_summary.tsv", [summary], SUMMARY_FIELDS)
    manifest["status"] = summary["status"]
    manifest["reason"] = summary["reason"]
    manifest["outputs"] = {
        p.name: sha256_file(p)
        for p in sorted(out.iterdir())
        if p.is_file() and p.name != "splice_manifest.json"
    }
    (out / "splice_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def run_splice(
    evidence_dir: Path,
    *,
    viral_fasta: str | None,
    annotation_gtfs: Sequence[Path] = (),
    competitor_mode: bool = False,
    threads: int = 4,
) -> dict[str, object]:
    """Write ``<evidence_dir>/splice/`` and return the ``splice_summary.tsv`` row.

    Competitor mode has no molecule verdicts, so the row is ``not_assessed`` and nothing is computed.
    Without a counted virus-best molecule the row is ``not_evaluable`` with the five counts.
    ``evidence_manifest_sha256`` is empty: the evidence manifest is written after this check and lists the
    hash of ``splice_manifest.json`` in its ``outputs`` instead.
    """
    out = evidence_dir / SPLICE_DIR
    out.mkdir(exist_ok=True)
    for stale in ("splice.bam", "splice.bam.bai", "target_genome.fa"):
        (out / stale).unlink(
            missing_ok=True
        )  # an earlier run's alignment must not outlive a refusal
    summary: dict[str, object] = dict.fromkeys(SUMMARY_FIELDS, "")
    summary["note"] = NOTE
    manifest: dict[str, Any] = {
        "schema_version": "1",
        "check": "evid_02_splice",
        "unit": "docs/evidence_unit.md: counted virus-best molecules, all member reads",
        "constants": {
            "MIN_MINIMAP2_VERSION": list(MIN_MINIMAP2_VERSION),
            "MAX_INTRON_LENGTH": MAX_INTRON_LENGTH,
            "MIN_INTRON_LENGTH": MIN_INTRON_LENGTH,
            "MIN_ANCHOR": MIN_ANCHOR,
            "MAX_MISMATCH_FRACTION": MAX_MISMATCH_FRACTION,
            "MAX_CLIP_FRACTION": MAX_CLIP_FRACTION,
            "CANONICAL_MOTIFS": sorted(CANONICAL_MOTIFS),
        },
        "note": NOTE,
        "evidence_manifest_sha256": None,
    }
    verdicts, lineage, bam = (evidence_dir / n for n in (VERDICTS_FILE, LINEAGE_FILE, RAW_BAM))

    if competitor_mode or not verdicts.is_file():
        summary.update(
            status="not_assessed",
            reason="molecule_verdict_not_produced_in_competitor_mode"
            if competitor_mode
            else "molecule_verdicts_missing",
            annotation_status=ANNOT_NOT_EVALUABLE,
        )
        return _finish(out, summary, [], manifest)

    manifest["inputs"] = {
        n: sha256_file(p)
        for n, p in ((VERDICTS_FILE, verdicts), (LINEAGE_FILE, lineage), (RAW_BAM, bam))
    }
    counts, molecules = _molecule_categories(verdicts)
    summary.update({field: counts[field] for field in CATEGORY_FIELDS})
    if not molecules:
        summary.update(
            status="not_evaluable",
            reason="no_counted_virus_best_molecules",
            annotation_status=ANNOT_NOT_EVALUABLE,
        )
        return _finish(out, summary, [], manifest)
    if viral_fasta is None:
        raise SpliceError("--splice needs --viral-fasta")

    reads = _member_reads(lineage, molecules)
    sequences, contigs = _read_sequences(bam, reads)
    if not contigs:
        raise SpliceError(
            "Counted virus-best molecules have no viral alignment in the competitive BAM"
        )
    genome = _write_target_genome(viral_fasta, contigs, out / "target_genome.fa")
    version = require_splice_minimap2()
    manifest["inputs"]["viral_fasta"] = sha256_file(viral_fasta)
    manifest["inputs"]["annotation_gtfs"] = {str(p): sha256_file(p) for p in annotation_gtfs}
    manifest["minimap2"] = {
        "version": version,
        "options": list(MINIMAP2_SPLICE_OPTIONS),
        "argv": [
            "minimap2",
            *MINIMAP2_SPLICE_OPTIONS,
            "-t",
            "<threads>",
            "splice/target_genome.fa",
            "<member reads FASTA>",
        ],
        "target_contigs": sorted(genome),
    }

    with tempfile.TemporaryDirectory(dir=out) as tmp:
        reads_fasta = Path(tmp) / "member_reads.fa"
        reads_fasta.write_text(
            "".join(f">{n}\n{s}\n" for n, s in sorted(sequences.items())), encoding="utf-8"
        )
        _align(reads_fasta, out / "target_genome.fa", out / "splice.bam", threads)
    n_aligned, n_spliced, found, rejected = _tally(out / "splice.bam", reads)

    introns, evaluable = annotated_introns(annotation_gtfs, genome)
    status = (
        ANNOT_NOT_EVALUABLE
        if not evaluable
        else ANNOT_PARTIAL
        if evaluable != set(genome)
        else ANNOT_AVAILABLE
        if any(introns.values())
        else ANNOT_NO_INTRONS
    )
    rows: list[dict[str, object]] = []
    spliced_molecules: set[tuple[str, str]] = set()
    annotated_molecules: set[tuple[str, str]] = set()
    excluded_molecules: set[tuple[str, str]] = set()
    n_excluded = 0
    for (contig, start, end), support in sorted(found.items()):
        motif, canonical, strand = motif_of(genome[contig], start, end)
        transcripts = introns.get(contig, {}).get((start, end), [])
        matched = contig in evaluable
        # A junction counts toward spliced_molecules if its motif is canonical or it matches an annotated
        # intron exactly (EVID-DECISIONS-01); annotation matching has no tolerance.
        counted = canonical or bool(transcripts)
        rows.append(
            {
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "motif": motif,
                "motif_canonical": canonical,
                "supporting_molecules": len(support["molecules"]),
                "supporting_reads": len(support["reads"]),
                "annotated": (bool(transcripts) if matched else ""),
                "annotated_transcript": ",".join(transcripts),
                "counted_in_spliced_molecules": counted,
            }
        )
        if counted:
            spliced_molecules |= support["molecules"]
        else:
            excluded_molecules |= support["molecules"]
            n_excluded += 1
        if transcripts:
            annotated_molecules |= support["molecules"]
    summary.update(
        status="assessed",
        reads_realigned=len(sequences),
        reads_aligned=n_aligned,
        spliced_reads=n_spliced,
        spliced_molecules=len(spliced_molecules),
        spliced_molecules_annotated=len(annotated_molecules)
        if status != ANNOT_NOT_EVALUABLE
        else "",
        unique_junctions=len(rows),
        unique_junctions_excluded_noncanonical=n_excluded,
        molecules_excluded_noncanonical=len(excluded_molecules - spliced_molecules),
        junction_calls_rejected=sum(rejected.values()),
        annotation_status=status,
    )
    manifest["rejections"] = dict(sorted(rejected.items()))
    return _finish(out, summary, rows, manifest)
