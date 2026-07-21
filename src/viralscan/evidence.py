"""Read-level evidence extraction and validation for viral hits.

ViralScan calls viruses from pseudo-alignment (kallisto|bustools), which gives
no genome-coordinate alignment and no per-base identity — so a call's read-level
evidence cannot be inspected directly. This module closes that gap:

1. **Trace** the BUS records whose (barcode, UMI) were assigned to viral genes
   back to the reads that produced them.
2. **Extract** those reads (the cDNA mate) to FASTA/FASTQ.
3. **Re-align** them to the viral genome with ``minimap2`` -> sorted/indexed BAM,
   directly loadable in IGV alongside the viral FASTA/GTF.
4. **Quantify** evidence quality: per-virus genome coverage (breadth/depth/reads)
   and, optionally, BLAST identity of the reads to the target virus vs. host
   off-targets (the cross-homology false-positive check).

The trace/extract layer is pure and unit-tested; alignment/coverage/BLAST are
thin wrappers over ``minimap2`` / ``samtools`` / ``blast+`` (shelled out, list
form, no shell=True), matching the rest of the package.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Optional, cast

log = logging.getLogger("viralscan")

# Barcode/UMI geometry (cb_len, umi_len) per technology. Unlike the historical
# host_filter table this includes non-10x droplet chemistries so the trace works
# for Drop-seq etc. (the same gap as PLAN S1).
_TECH_GEOMETRY: dict[str, tuple[int, int]] = {
    "10xv1": (14, 10),
    "10xv2": (16, 10),
    "10xv3": (16, 12),
    "10xv3_5p": (16, 12),
    "dropseq": (12, 8),
}


def cb_umi_geometry(technology: str) -> tuple[int, int]:
    """Return ``(cb_len, umi_len)`` for *technology*.

    Accepts the named chemistries above (case-insensitive) or an explicit
    kallisto ``bc:umi:seq`` triplet like ``0,0,16:0,16,28:1,0,0`` from which the
    CB and UMI lengths are read directly. Raises ``ValueError`` for anything it
    cannot resolve — extracting reads with the wrong geometry silently yields
    nothing, so this fails loudly instead.
    """
    key = technology.strip().lower()
    if key in _TECH_GEOMETRY:
        return _TECH_GEOMETRY[key]
    # explicit "bc:umi:seq" with comma-separated (file,start,end) triplets
    if ":" in technology:
        try:
            bc, umi, _seq = technology.split(":")[:3]
            _, bcs, bce = (int(x) for x in bc.split(","))
            _, umis, umie = (int(x) for x in umi.split(","))
            return bce - bcs, umie - umis
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Cannot parse barcode geometry from -x {technology!r}: {exc}"
            ) from exc
    raise ValueError(
        f"Unknown technology {technology!r}; add it to _TECH_GEOMETRY or pass an "
        "explicit 'bc:umi:seq' geometry string."
    )


def viral_equivalence_classes(
    ec_map: dict[int, list[int]], viral_gene_indices: set[int]
) -> set[int]:
    """Return EC ids whose gene set includes at least one viral gene."""
    return {ec for ec, genes in ec_map.items() if any(g in viral_gene_indices for g in genes)}


def viral_assigned_keys(bus_text: Iterable[str], viral_ecs: set[int]) -> set[tuple[str, str]]:
    """Collect ``(barcode, umi)`` pairs whose EC is viral, from BUS text lines.

    *bus_text* yields ``bustools text`` output lines: ``barcode\\tumi\\tec\\tcount``.
    A pair is viral-assigned if any of its records maps to a viral EC.
    """
    keys: set[tuple[str, str]] = set()
    for line in bus_text:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        try:
            ec = int(parts[2])
        except ValueError:
            continue
        if ec in viral_ecs:
            keys.add((parts[0], parts[1]))
    return keys


def _open_maybe_gzip(path: str, mode: str = "rt") -> IO[str]:
    if str(path).endswith(".gz"):
        return cast(IO[str], gzip.open(path, mode))
    return open(path, mode)


@dataclass
class ExtractionStats:
    total_reads: int
    viral_reads: int


def extract_viral_reads(
    r1_path: str,
    r2_path: str,
    keys: set[tuple[str, str]],
    cb_len: int,
    umi_len: int,
    out_fasta: str,
    strip_suffix: str = "-1",
) -> ExtractionStats:
    """Write the cDNA (R2) reads whose (CB, UMI) from R1 is in *keys* to FASTA.

    The CB occupies bases ``0:cb_len`` of every R1 read and the UMI
    ``cb_len:cb_len+umi_len`` (10x/Drop-seq layout). Output is FASTA (one record
    per surviving read) named ``<CB>_<UMI>_<n>`` so reads remain traceable to
    their cell. Returns counts for a quick yield report.
    """
    bc_end = cb_len + umi_len
    total = 0
    kept = 0
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    with (
        _open_maybe_gzip(r1_path) as fq1,
        _open_maybe_gzip(r2_path) as fq2,
        open(out_fasta, "w") as out,
    ):
        while True:
            h1 = fq1.readline()
            s1 = fq1.readline()
            fq1.readline()
            fq1.readline()
            h2 = fq2.readline()
            s2 = fq2.readline()
            fq2.readline()
            fq2.readline()
            if not h1 or not h2:
                break
            total += 1
            seq1 = s1.rstrip("\n")
            cb = seq1[:cb_len]
            umi = seq1[cb_len:bc_end]
            # BUS barcodes carry no lane suffix; the counts matrix may. Match raw.
            if (cb, umi) in keys or (cb, umi.rstrip()) in keys:
                kept += 1
                out.write(f">{cb}_{umi}_{kept}\n{s2.rstrip()}\n")
    return ExtractionStats(total_reads=total, viral_reads=kept)


# ---------------------------------------------------------------------------
# Alignment / coverage / BLAST — thin tool wrappers
# ---------------------------------------------------------------------------


def _run(cmd: list[str], *, stdin: Optional[bytes] = None, capture: bool = False) -> bytes:
    """Run *cmd* (list form, no shell), surfacing the tool's stderr on failure.

    Returns stdout bytes when *capture* is True. On a non-zero exit raises
    ``RuntimeError`` that includes the tool name and its stderr — never silently
    swallows the error (the wrappers used to ``stderr=DEVNULL`` which hid the
    real cause behind a bare CalledProcessError).
    """
    proc = subprocess.run(
        cmd,
        input=stdin,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        err = (proc.stderr or b"").decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"{cmd[0]} failed (exit {proc.returncode}): {err or 'no stderr output'}")
    return proc.stdout if capture else b""


def align_reads_to_viral(reads_fasta: str, viral_fasta: str, out_bam: str, threads: int = 4) -> str:
    """minimap2 short-read align *reads_fasta* to *viral_fasta* -> sorted, indexed BAM."""
    out_bam = str(out_bam)
    Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
    sam = _run(
        ["minimap2", "-ax", "sr", "-t", str(threads), viral_fasta, reads_fasta], capture=True
    )
    _run(["samtools", "sort", "-@", str(threads), "-o", out_bam, "-"], stdin=sam)
    _run(["samtools", "index", out_bam])
    return out_bam


def _parse_coverage_output(text: str) -> list[dict[str, str]]:
    """Parse ``samtools coverage`` TSV text into dicts, keeping only covered references.

    Extracted from ``coverage_table`` so it can be unit-tested against synthetic
    tool output without requiring the ``samtools`` binary.
    """
    rows: list[dict[str, str]] = []
    header: list[str] = []
    for i, line in enumerate(text.splitlines()):
        cols = line.split("\t")
        if i == 0:
            header = [c.lstrip("#") for c in cols]
            continue
        rec = dict(zip(header, cols))
        if int(rec.get("numreads", 0) or 0) > 0:
            rows.append(rec)
    return rows


def coverage_table(bam: str) -> list[dict[str, str]]:
    """Per-reference coverage from ``samtools coverage`` (breadth, depth, #reads)."""
    stdout = _run(["samtools", "coverage", bam], capture=True).decode("utf-8", errors="replace")
    return _parse_coverage_output(stdout)


def _cigar_ref_span(cigar: str) -> int:
    """Reference bases consumed by a CIGAR string (M/D/N/=/X operations)."""
    if not cigar or cigar == "*":
        return 0
    span = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if num and ch in "MDN=X":
                span += int(num)
            num = ""
    return span


def _cb_umi(qname: str) -> Optional[tuple[str, str]]:
    """(CB, UMI) from an extracted-read name ``<CB>_<UMI>_<n>``, else None.

    Returns None if either the CB or UMI field is empty (a degenerate name like
    ``_UMI_1`` would otherwise yield an invalid empty SAM tag value)."""
    parts = qname.split("_")
    return (parts[0], parts[1]) if len(parts) >= 3 and parts[0] and parts[1] else None


def _parse_sam_read_starts(
    sam_text: str, *, dedup: str = "umi", strand_aware: bool = True, bin_size: int = 1
) -> list[dict[str, object]]:
    """Tally 5′ read-start positions per reference from ``samtools view`` text.

    Each primary alignment contributes its 5′ start: leftmost 0-based POS on the
    forward strand, or ``POS + reference_span − 1`` on the reverse strand when
    *strand_aware*. With ``dedup="umi"`` reads collapse to one per (CB, UMI,
    reference) — the scRNA-appropriate PCR-duplicate removal, read from the
    ``<CB>_<UMI>_<n>`` names ``extract_viral_reads`` writes; ``dedup="none"``
    keeps every read (use when dups were already removed, e.g. samtools markdup).

    Split from ``read_start_distribution`` so it is unit-testable without samtools.
    Returns rows sorted by (reference, position): {reference, position,
    n_read_starts, n_reads}.
    """
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")
    seen: set[tuple[object, str]] = set()
    hist: dict[tuple[str, int], int] = {}
    n_reads: dict[str, int] = {}
    for line in sam_text.splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        try:
            flag, pos1 = int(f[1]), int(f[3])  # SAM FLAG, POS (1-based)
        except ValueError:
            continue  # not a valid alignment record (e.g. a stray warning line)
        if flag & 0x904 or f[2] == "*":  # unmapped / secondary / supplementary / no ref
            continue
        qname, rname, pos0 = f[0], f[2], pos1 - 1
        # 5' end: leftmost POS on +, rightmost consumed ref base on - (strand-aware).
        reverse = strand_aware and bool(flag & 0x10)
        start = pos0 + max(_cigar_ref_span(f[5]) - 1, 0) if reverse else pos0
        if dedup == "umi":
            key = (_cb_umi(qname) or qname, rname)
            if key in seen:
                continue
            seen.add(key)
        n_reads[rname] = n_reads.get(rname, 0) + 1
        bin_pos = (start // bin_size) * bin_size
        hist[(rname, bin_pos)] = hist.get((rname, bin_pos), 0) + 1
    return [
        {"reference": rn, "position": p, "n_read_starts": c, "n_reads": n_reads[rn]}
        for (rn, p), c in sorted(hist.items())
    ]


def read_start_distribution(
    bam: str, *, dedup: str = "umi", strand_aware: bool = True, bin_size: int = 1
) -> list[dict[str, object]]:
    """Per-position 5′ read-start distribution along each viral reference.

    ``dedup``: ``umi`` (collapse per CB+UMI; default, scRNA-appropriate),
    ``markdup`` (``samtools markdup -r`` then tally survivors), or ``none``.
    Exposes 3′ bias, subgenomic-RNA junctions, and EVE/integration hotspots that
    the aggregate ``coverage_table`` cannot.
    """
    if dedup not in ("umi", "markdup", "none"):
        raise ValueError(f"dedup must be umi|markdup|none, got {dedup!r}")
    view_bam = bam
    if dedup == "markdup":
        marked = str(Path(bam).with_name(Path(bam).stem + ".markdup.bam"))
        _run(["samtools", "markdup", "-r", bam, marked])  # coordinate-sorted input assumed
        _run(["samtools", "index", marked])
        view_bam = marked
    sam = _run(["samtools", "view", view_bam], capture=True).decode("utf-8", errors="replace")
    return _parse_sam_read_starts(
        sam,
        dedup="none" if dedup == "markdup" else dedup,  # markdup already dropped dups
        strand_aware=strand_aware,
        bin_size=bin_size,
    )


def add_cell_tags_to_sam(sam_text: str) -> str:
    """Append ``CB:Z:<cb>`` and ``UB:Z:<umi>`` tags to each alignment record.

    The barcode/UMI are read from the ``<CB>_<UMI>_<n>`` read names that
    ``extract_viral_reads`` writes, so the viral BAM can be grouped by cell in a
    genome browser (IGV "Group by → tag → CB"). Header lines (``@...``) and
    records with an un-parseable name pass through unchanged. Split out for unit
    testing without samtools.
    """
    out: list[str] = []
    for line in sam_text.splitlines():
        fields = line.split("\t")
        # Header (@...) or a line without the 11 mandatory SAM fields: pass through
        # unchanged rather than append a tag after a non-optional field.
        if line.startswith("@") or len(fields) < 11:
            out.append(line)
            continue
        cbumi = _cb_umi(fields[0])
        out.append(f"{line}\tCB:Z:{cbumi[0]}\tUB:Z:{cbumi[1]}" if cbumi else line)
    return "\n".join(out) + "\n"


def write_tagged_bam(bam: str, out_bam: str) -> str:
    """Write a CB/UB-tagged, indexed copy of *bam* for per-cell IGV inspection."""
    sam = _run(["samtools", "view", "-h", bam], capture=True).decode("utf-8", errors="replace")
    _run(
        ["samtools", "view", "-b", "-o", str(out_bam), "-"],
        stdin=add_cell_tags_to_sam(sam).encode(),
    )
    _run(["samtools", "index", str(out_bam)])
    return str(out_bam)


def blast_identity(
    reads_fasta: str,
    viral_fasta: str,
    workdir: str,
    max_reads: int = 500,
    threads: int = 4,
) -> list[dict[str, str]]:
    """BLAST a sample of extracted reads against the viral reference.

    Builds a local BLAST db from *viral_fasta* and reports the best hit per read
    (subject, % identity, alignment length). This quantifies how well the
    viral-assigned reads actually match the target genome vs. an off-target —
    the cross-homology / false-positive check. Returns one dict per read hit.
    """
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    sample = work / "reads_sample.fasta"
    n_sampled = _head_fasta(reads_fasta, sample, max_reads)
    log.info("BLAST: sampling %d reads (cap %d) from %s", n_sampled, max_reads, reads_fasta)
    db = work / "viral_db"
    _run(["makeblastdb", "-in", viral_fasta, "-dbtype", "nucl", "-out", str(db)], capture=True)
    stdout = _run(
        [
            "blastn",
            "-query",
            str(sample),
            "-db",
            str(db),
            "-max_target_seqs",
            "1",
            "-num_threads",
            str(threads),
            "-outfmt",
            "6 qseqid sseqid pident length evalue",
        ],
        capture=True,
    ).decode("utf-8", errors="replace")
    return _parse_blast_output(stdout)


def _parse_blast_output(text: str) -> list[dict[str, str]]:
    """Parse ``blastn -outfmt '6 qseqid sseqid pident length evalue'`` text.

    Returns one dict per query (best hit only — first occurrence per read ID).
    Extracted from ``blast_identity`` so it can be unit-tested against synthetic
    tool output without requiring the ``blastn`` binary.
    """
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 4 and f[0] not in seen:  # best hit per read (first = top)
            seen.add(f[0])
            rows.append({"read": f[0], "subject": f[1], "pident": f[2], "length": f[3]})
    return rows


def _head_fasta(src: str, dst: Path, max_records: int) -> int:
    """Copy the first *max_records* FASTA records from *src* to *dst*; return the count."""
    n = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if n >= max_records:
                    break
                n += 1
            fout.write(line)
    return n


def have_tools(names: Iterable[str]) -> list[str]:
    """Return the subset of *names* that are NOT on PATH."""
    return [n for n in names if shutil.which(n) is None]
