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

import csv
import gzip
import hashlib
import heapq
import json
import logging
import math
import shutil
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Optional, cast

from viralscan.anellovirus import merged_name_map
from viralscan.virus_grouping import group_genes_by_virus

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

VIRUS_ALIASES: dict[str, str] = {
    "ebv": "Epstein-Barr virus",
    "hhv4": "Epstein-Barr virus",
    "hsv1": "Human herpesvirus 1",
    "hhv1": "Human herpesvirus 1",
    "hsv2": "Human herpesvirus 2",
    "hhv2": "Human herpesvirus 2",
    "hhv6a": "Human herpesvirus 6A",
    "hhv6b": "Human herpesvirus 6B",
    "kshv": "Kaposi sarcoma-associated herpesvirus",
    "hhv8": "Kaposi sarcoma-associated herpesvirus",
    "ttv": "Anelloviridae",
}


def resolve_viral_target(
    selector: str,
    viral_gene_ids: Iterable[str],
    *,
    detected_virus_names: Iterable[str] = (),
) -> tuple[str, list[str]]:
    """Resolve one exact accession/gene, alias, or detected canonical call.

    Substring matching is deliberately forbidden: a selector must equal a gene
    ID, a reference alias/prefix, a canonical virus label, or a detected label
    from ``viral_summary.tsv`` (case-insensitive).
    """
    query = selector.strip()
    if not query:
        raise ValueError("An exact --virus selector is required; refusing to trace every virus.")
    genes = list(dict.fromkeys(viral_gene_ids))
    by_lower_gene: dict[str, list[str]] = {}
    for gene in genes:
        by_lower_gene.setdefault(gene.casefold(), []).append(gene)
    exact_genes = by_lower_gene.get(query.casefold(), [])
    if exact_genes:
        if len(exact_genes) != 1:
            raise ValueError(f"Target {selector!r} matches duplicate gene IDs: {exact_genes}")
        return exact_genes[0], exact_genes

    name_map = merged_name_map()
    groups, _ = group_genes_by_virus(genes, name_map)
    canonical_by_lower = {name.casefold(): name for name in groups}
    canonical_by_lower.update(
        {str(name).casefold(): str(name) for name in detected_virus_names if str(name).strip()}
    )
    alias_to_name = {key.casefold(): value for key, value in name_map.items()}
    alias_to_name.update(VIRUS_ALIASES)
    canonical = alias_to_name.get(query.casefold()) or canonical_by_lower.get(query.casefold())
    if canonical and canonical in groups and groups[canonical]:
        return canonical, groups[canonical]

    choices = sorted(set(groups).union(detected_virus_names))
    raise ValueError(
        f"No exact viral target matches {selector!r}. Use an accession/gene ID, canonical "
        f"label, or detected call. Available calls include: {choices[:12]}"
    )


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


def _target_assignment(
    genes: list[int], target_genes: set[int], viral_genes: set[int], method: str
) -> tuple[str, float | None, str]:
    distinct = set(genes)
    target = distinct & target_genes
    host = distinct - viral_genes
    if len(distinct) == 1 and target:
        return "unique_target", 1.0, "candidate_unique"
    if host and target:
        weight = 0.0 if method == "host-conservative" else len(target) / len(distinct)
        return "host_virus_ambiguous", weight, "candidate_host_virus_ambiguous"
    if target:
        return (
            "virus_virus_ambiguous",
            len(target) / len(distinct),
            "candidate_virus_ambiguous",
        )
    return "not_target", 0.0, "not_detected"


def parse_flagged_target_bus(
    lines: Iterable[str],
    ec_map: dict[int, list[int]],
    target_gene_indices: set[int],
    viral_gene_indices: set[int],
    method: str,
) -> dict[int, dict[str, object]]:
    """Parse ``bustools text -f`` rows into exact read-number lineage."""
    result: dict[int, dict[str, object]] = {}
    for line_number, line in enumerate(lines, 1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            raise ValueError(f"Flagged BUS row {line_number} has fewer than five columns.")
        barcode, umi, ec_raw, _count, flag_raw = fields[:5]
        ec, flag = int(ec_raw), int(flag_raw)
        genes = sorted(set(ec_map.get(ec, [])))
        ambiguity, weight, tier = _target_assignment(
            genes, target_gene_indices, viral_gene_indices, method
        )
        if not set(genes).intersection(target_gene_indices):
            continue
        result[flag] = {
            "cb": barcode,
            "ub": umi,
            "ecs": str(ec),
            "compatible_genes": ",".join(str(gene) for gene in genes),
            "ambiguity_class": ambiguity,
            "assigned_weight": weight,
            "method": method,
            "evidence_tier": tier,
            "exclusion_reason": "",
        }
    return result


def _canonical_fastq_id(header: str) -> str:
    token = header.strip().split(maxsplit=1)[0].removeprefix("@")
    return token[:-2] if token.endswith(("/1", "/2")) else token


def extract_exact_reads_by_number(
    r1_path: str,
    r2_path: str,
    lineage_by_read_number: dict[int, dict[str, object]],
    out_fasta: str,
    lineage_tsv_gz: str,
) -> ExtractionStats:
    """Extract only exact FASTQ records named by kallisto's BUS flag column."""
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "read_number",
        "read_id",
        "cb",
        "ub",
        "ecs",
        "compatible_genes",
        "ambiguity_class",
        "assigned_weight",
        "method",
        "evidence_tier",
        "exclusion_reason",
    ]
    found: set[int] = set()
    total = 0
    with (
        _open_maybe_gzip(r1_path) as fq1,
        _open_maybe_gzip(r2_path) as fq2,
        open(out_fasta, "w") as fasta,
        gzip.open(lineage_tsv_gz, "wt", newline="") as lineage,
    ):
        writer = csv.DictWriter(lineage, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        while True:
            rec1 = [fq1.readline() for _ in range(4)]
            rec2 = [fq2.readline() for _ in range(4)]
            if not rec1[0] and not rec2[0]:
                break
            if not all(rec1) or not all(rec2):
                raise ValueError("Input FASTQs are truncated or have different record counts.")
            read_id = _canonical_fastq_id(rec1[0])
            if read_id != _canonical_fastq_id(rec2[0]):
                raise ValueError(f"FASTQ mate mismatch at read {total}: {read_id!r}")
            metadata = lineage_by_read_number.get(total)
            if metadata is not None:
                found.add(total)
                fasta.write(
                    f">{metadata['cb']}_{metadata['ub']}_{total}|{read_id}\n{rec2[1].rstrip()}\n"
                )
                writer.writerow({"read_number": total, "read_id": read_id, **metadata})
            total += 1
    missing = set(lineage_by_read_number).difference(found)
    if missing:
        raise ValueError(
            f"{len(missing)} BUS read numbers were absent from the FASTQs; exact lineage failed."
        )
    return ExtractionStats(total_reads=total, viral_reads=len(found))


def replay_exact_target_bus(
    *,
    index: str,
    technology: str,
    r1_path: str,
    r2_path: str,
    ec_file: str,
    transcripts_file: str,
    target_transcripts: Iterable[str],
    workdir: str,
    threads: int,
) -> Path:
    """Re-pseudoalign with read-number flags and capture only target transcripts."""
    missing = have_tools(["kallisto", "bustools"])
    if missing:
        raise RuntimeError(f"Exact read lineage requires: {', '.join(missing)}")
    work = Path(workdir)
    bus_dir = work / "lineage_bus"
    bus_dir.mkdir(parents=True, exist_ok=True)
    capture_list = work / "target_transcripts.txt"
    capture_list.write_text("\n".join(sorted(set(target_transcripts))) + "\n")
    _run(
        [
            "kallisto",
            "bus",
            "-i",
            index,
            "-o",
            str(bus_dir),
            "-x",
            technology,
            "-t",
            str(threads),
            "-n",
            r1_path,
            r2_path,
        ]
    )
    captured = work / "target.bus"
    _run(
        [
            "bustools",
            "capture",
            "-s",
            "-c",
            str(capture_list),
            "-e",
            ec_file,
            "-t",
            transcripts_file,
            "-o",
            str(captured),
            str(bus_dir / "output.bus"),
        ]
    )
    by_flag = work / "target.by_flag.bus"
    _run(["bustools", "sort", "--flags", "-t", str(threads), "-o", str(by_flag), str(captured)])
    flagged_text = work / "target.by_flag.bus.txt"
    _run(["bustools", "text", "-f", "-o", str(flagged_text), str(by_flag)])
    return flagged_text


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


def write_competitive_fasta(host_fasta: str, viral_fasta: str, output: str) -> str:
    """Combine full-host and exact-target FASTAs with unambiguous subject prefixes."""
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as target:
        for prefix, source in (("HOST", host_fasta), ("VIRUS", viral_fasta)):
            with _open_maybe_gzip(source) as handle:
                for line in handle:
                    if line.startswith(">"):
                        name, *description = line[1:].rstrip().split(maxsplit=1)
                        suffix = f" {description[0]}" if description else ""
                        target.write(f">{prefix}|{name}{suffix}\n")
                    else:
                        target.write(line)
    return str(out)


def write_igv_session(reference_fasta: str, bam_paths: Iterable[str], output: str) -> str:
    """Write a minimal portable IGV session referencing indexed BAM resources."""
    resources = "".join(f'<Resource path="{Path(path).resolve()}"/>' for path in bam_paths)
    Path(output).write_text(
        f'<Session genome="{Path(reference_fasta).resolve()}" version="8">'
        f"<Resources>{resources}</Resources></Session>\n"
    )
    return output


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


def _covered_intervals(positions: list[int]) -> str:
    if not positions:
        return ""
    intervals: list[str] = []
    start = previous = positions[0]
    for position in positions[1:]:
        if position != previous + 1:
            intervals.append(f"{start}-{previous}")
            start = position
        previous = position
    intervals.append(f"{start}-{previous}")
    return ",".join(intervals)


def _gini(values: list[int]) -> float:
    if not values or sum(values) == 0:
        return 0.0
    ordered = sorted(values)
    n = len(ordered)
    return (2 * sum((i + 1) * value for i, value in enumerate(ordered))) / (n * sum(ordered)) - (
        n + 1
    ) / n


def _alignment_qc_from_text(
    header_text: str, sam_text: str, depth_text: str
) -> list[dict[str, object]]:
    """Compute per-reference specificity, depth, and hotspot diagnostics."""
    lengths: dict[str, int] = {}
    for line in header_text.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" in fields and "LN" in fields:
            lengths[fields["SN"]] = int(fields["LN"])

    depths: dict[str, list[tuple[int, int]]] = {}
    for line in depth_text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            depths.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))

    stats: dict[str, dict[str, object]] = {}
    total_mapped = host_mapped = 0
    for line in sam_text.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & 0x904 or fields[2] == "*":
            continue
        reference = fields[2]
        total_mapped += 1
        host_mapped += int(reference.startswith("HOST|"))
        record = stats.setdefault(
            reference,
            {
                "reads": 0,
                "reverse": 0,
                "mapq": [],
                "identities": [],
                "starts": [],
                "cells": set(),
                "molecules": set(),
            },
        )
        record["reads"] = int(record["reads"]) + 1
        record["reverse"] = int(record["reverse"]) + int(bool(flag & 0x10))
        cast(list[int], record["mapq"]).append(int(fields[4]))
        pos = int(fields[3])
        cast(list[int], record["starts"]).append(pos)
        cbumi = _cb_umi(fields[0])
        if cbumi:
            cast(set[str], record["cells"]).add(cbumi[0])
            cast(set[tuple[str, str]], record["molecules"]).add(cbumi)
        nm = next((tag for tag in fields[11:] if tag.startswith("NM:i:")), None)
        aligned = _cigar_ref_span(fields[5])
        if nm and aligned:
            cast(list[float], record["identities"]).append(
                max(0.0, 100.0 * (aligned - int(nm.split(":")[-1])) / aligned)
            )

    rows: list[dict[str, object]] = []
    for reference, record in sorted(stats.items()):
        length = lengths.get(reference, 0)
        observed = sorted(depths.get(reference, []))
        depth_values = [depth for _pos, depth in observed]
        positions = [position for position, depth in observed if depth > 0]
        zeros = max(length - len(depth_values), 0)
        sorted_depth = sorted(depth_values)

        def depth_at(
            rank: int, zeros: int = zeros, sorted_depth: list[int] = sorted_depth
        ) -> float:
            return 0.0 if rank < zeros else float(sorted_depth[rank - zeros])

        if not length:
            median_depth = 0.0
        elif length % 2:
            median_depth = depth_at(length // 2)
        else:
            median_depth = (depth_at(length // 2 - 1) + depth_at(length // 2)) / 2
        starts = cast(list[int], record["starts"])
        start_counts: dict[int, int] = {}
        window_counts: dict[int, int] = {}
        for start in starts:
            start_counts[start] = start_counts.get(start, 0) + 1
            window = start // 50
            window_counts[window] = window_counts.get(window, 0) + 1
        counts = list(start_counts.values())
        entropy = -sum((count / len(starts)) * math.log2(count / len(starts)) for count in counts)
        reads = int(record["reads"])
        identities = cast(list[float], record["identities"])
        mapq = cast(list[int], record["mapq"])
        rows.append(
            {
                "reference": reference,
                "reference_class": "host" if reference.startswith("HOST|") else "virus",
                "reads": reads,
                "molecules": len(cast(set[tuple[str, str]], record["molecules"])),
                "cells": len(cast(set[str], record["cells"])),
                "breadth_1x": len(positions) / length if length else 0.0,
                "breadth_3x": sum(depth >= 3 for depth in depth_values) / length if length else 0.0,
                "breadth_10x": sum(depth >= 10 for depth in depth_values) / length
                if length
                else 0.0,
                "mean_depth": sum(depth_values) / length if length else 0.0,
                "median_depth": median_depth,
                "covered_intervals": _covered_intervals(positions),
                "read_start_entropy": entropy,
                "read_start_gini": _gini(counts),
                "max_50bp_window_fraction": max(window_counts.values()) / reads,
                "reverse_strand_fraction": int(record["reverse"]) / reads,
                "mean_mapping_quality": sum(mapq) / len(mapq),
                "mean_identity": sum(identities) / len(identities) if identities else "",
                "host_competitive_fraction": host_mapped / total_mapped if total_mapped else 0.0,
            }
        )
    return rows


def alignment_qc_table(bam: str) -> list[dict[str, object]]:
    """Compute alignment QC from samtools header, alignments, and covered depths."""
    header = _run(["samtools", "view", "-H", bam], capture=True).decode(errors="replace")
    sam = _run(["samtools", "view", bam], capture=True).decode(errors="replace")
    depth = _run(["samtools", "depth", bam], capture=True).decode(errors="replace")
    return _alignment_qc_from_text(header, sam, depth)


def _per_cell_qc_from_text(sam_text: str) -> list[dict[str, object]]:
    """Summarise primary competitive alignments per cell and reference class."""
    stats: dict[tuple[str, str], dict[str, object]] = {}
    for line in sam_text.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & 0x904 or fields[2] == "*":
            continue
        cbumi = _cb_umi(fields[0])
        if not cbumi:
            continue
        reference_class = "host" if fields[2].startswith("HOST|") else "virus"
        record = stats.setdefault(
            (cbumi[0], reference_class),
            {"reads": 0, "reverse": 0, "mapq": [], "identities": [], "molecules": set()},
        )
        record["reads"] = int(record["reads"]) + 1
        record["reverse"] = int(record["reverse"]) + int(bool(flag & 0x10))
        cast(list[int], record["mapq"]).append(int(fields[4]))
        cast(set[tuple[str, str]], record["molecules"]).add(cbumi)
        aligned = _cigar_ref_span(fields[5])
        nm = next((tag for tag in fields[11:] if tag.startswith("NM:i:")), None)
        if aligned and nm:
            cast(list[float], record["identities"]).append(
                max(0.0, 100.0 * (aligned - int(nm.split(":")[-1])) / aligned)
            )

    rows: list[dict[str, object]] = []
    for (cell, reference_class), record in sorted(stats.items()):
        reads = int(record["reads"])
        mapq = cast(list[int], record["mapq"])
        identities = cast(list[float], record["identities"])
        rows.append(
            {
                "cell_barcode": cell,
                "reference_class": reference_class,
                "reads": reads,
                "molecules": len(cast(set[tuple[str, str]], record["molecules"])),
                "reverse_strand_fraction": int(record["reverse"]) / reads,
                "mean_mapping_quality": sum(mapq) / len(mapq),
                "mean_identity": sum(identities) / len(identities) if identities else "",
            }
        )
    return rows


def per_cell_alignment_qc(bam: str) -> list[dict[str, object]]:
    """Return per-cell competitive alignment summaries for a BAM."""
    sam = _run(["samtools", "view", bam], capture=True).decode(errors="replace")
    return _per_cell_qc_from_text(sam)


def coverage_depth_points(bam: str) -> list[dict[str, object]]:
    """Return covered depth points for plotting without manufacturing zero depth."""
    text = _run(["samtools", "depth", "-a", bam], capture=True).decode(errors="replace")
    rows: list[dict[str, object]] = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            rows.append(
                {"reference": fields[0], "position": int(fields[1]), "depth": int(fields[2])}
            )
    return rows


def plot_coverage_comparison(raw_bam: str, dedup_bam: str, output: str) -> str:
    """Plot raw and deduplicated competitive depth tracks per viral reference."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layers = {
        "raw": coverage_depth_points(raw_bam),
        "deduplicated": coverage_depth_points(dedup_bam),
    }
    references = sorted(
        {
            str(row["reference"])
            for rows in layers.values()
            for row in rows
            if str(row["reference"]).startswith("VIRUS|")
        }
    )
    if not references:
        references = sorted({str(row["reference"]) for rows in layers.values() for row in rows})
    figure, axes = plt.subplots(
        max(len(references), 1), 1, figsize=(9, max(2.5, 2.5 * len(references))), squeeze=False
    )
    if not references:
        axes[0][0].text(0.5, 0.5, "No aligned coverage", ha="center", va="center")
        axes[0][0].set_axis_off()
    for axis, reference in zip(axes[:, 0], references):
        for layer, rows in layers.items():
            selected = [row for row in rows if row["reference"] == reference]
            axis.plot(
                [int(row["position"]) for row in selected],
                [int(row["depth"]) for row in selected],
                label=layer,
                linewidth=1,
            )
        axis.set(title=reference, xlabel="Position (bp)", ylabel="Depth")
        axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)
    return output


def interpretation_flags(
    qc_rows: Iterable[dict[str, object]],
    blast_rows: Iterable[dict[str, str]],
    lineage_rows: Iterable[dict[str, object]],
) -> list[dict[str, object]]:
    """Create transparent diagnostic flags; these are never biological conclusions."""
    qc = [row for row in qc_rows if row.get("reference_class") == "virus"]
    blast = list(blast_rows)
    lineage = list(lineage_rows)
    max_hotspot = max((float(row.get("max_50bp_window_fraction", 0)) for row in qc), default=0)
    max_breadth = max((float(row.get("breadth_1x", 0)) for row in qc), default=0)
    host_fraction = max((float(row.get("host_competitive_fraction", 0)) for row in qc), default=0)
    low_complexity_fraction = (
        sum(row.get("low_complexity") == "true" for row in blast) / len(blast) if blast else 0
    )
    host_preferred_fraction = (
        sum(float(row.get("viral_minus_host_bitscore") or 0) <= 0 for row in blast) / len(blast)
        if blast
        else 0
    )
    ambiguous_fraction = (
        sum(row.get("ambiguity_class") != "unique_target" for row in lineage) / len(lineage)
        if lineage
        else 0
    )
    metrics = [
        (
            "host_homology",
            max(host_fraction, host_preferred_fraction),
            0.25,
            "competitive host support",
        ),
        ("low_complexity", low_complexity_fraction, 0.25, "low-complexity sampled reads"),
        ("sibling_or_host_ambiguity", ambiguous_fraction, 0.25, "non-unique target lineage"),
        ("coverage_hotspot", max_hotspot, 0.50, "maximum 50-bp read-start window"),
        (
            "eve_or_integration_like_hotspot",
            max_hotspot if max_breadth < 0.10 else 0.0,
            0.50,
            "narrow breadth plus hotspot",
        ),
    ]
    rows = [
        {
            "flag": name,
            "status": "flagged" if value >= threshold else "not_flagged",
            "metric": value,
            "threshold": threshold,
            "basis": basis,
            "interpretation": "diagnostic_only",
        }
        for name, value, threshold, basis in metrics
    ]
    rows.extend(
        {
            "flag": name,
            "status": "not_assessed",
            "metric": "",
            "threshold": "",
            "basis": basis,
            "interpretation": "requires target-specific model",
        }
        for name, basis in (
            (
                "expected_3prime_chemistry_bias",
                "requires transcript orientation and chemistry model",
            ),
            ("subgenomic_rna_like_pattern", "requires virus-specific junction annotation"),
            ("contamination", "requires negative controls and ambient model"),
        )
    )
    return rows


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
            key = (_cb_umi(qname) or qname, rname, reverse, start, f[5])
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


def deduplicate_umi_sam(sam_text: str) -> str:
    """Keep one alignment per CB, UB, reference, strand, start, and CIGAR."""
    output: list[str] = []
    seen: set[tuple[object, ...]] = set()
    for line in sam_text.splitlines():
        if line.startswith("@"):
            output.append(line)
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        try:
            flag, pos = int(fields[1]), int(fields[3])
        except ValueError:
            continue
        if flag & 0x904 or fields[2] == "*":
            output.append(line)
            continue
        key = (
            _cb_umi(fields[0]) or fields[0],
            fields[2],
            bool(flag & 0x10),
            pos,
            fields[5],
        )
        if key not in seen:
            seen.add(key)
            output.append(line)
    return "\n".join(output) + "\n"


def deduplicate_bam(bam: str, out_bam: str, mode: str, threads: int = 4) -> str:
    """Create and index a separate deduplicated BAM while preserving the raw BAM."""
    if mode not in {"umi", "markdup", "none"}:
        raise ValueError(f"dedup mode must be umi|markdup|none, got {mode!r}")
    if mode == "none":
        shutil.copy2(bam, out_bam)
    elif mode == "markdup":
        _run(["samtools", "markdup", "-r", bam, out_bam])
    else:
        sam = _run(["samtools", "view", "-h", bam], capture=True).decode("utf-8", errors="replace")
        _run(
            ["samtools", "view", "-b", "-o", out_bam, "-"],
            stdin=deduplicate_umi_sam(sam).encode(),
        )
    _run(["samtools", "index", out_bam])
    return out_bam


def blast_identity(
    reads_fasta: str,
    viral_fasta: str,
    workdir: str,
    max_reads: int = 500,
    threads: int = 4,
    seed: int = 42,
    sampling_manifest: str | None = None,
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
    n_sampled, n_total = sample_fasta_deterministic(reads_fasta, sample, max_reads, seed)
    fraction = n_sampled / n_total if n_total else 0.0
    log.info(
        "BLAST: deterministically sampled %d/%d reads (seed %d) from %s",
        n_sampled,
        n_total,
        seed,
        reads_fasta,
    )
    if sampling_manifest:
        Path(sampling_manifest).write_text(
            json.dumps(
                {
                    "strategy": "lowest_sha256",
                    "seed": seed,
                    "sampled_reads": n_sampled,
                    "total_reads": n_total,
                    "sampling_fraction": fraction,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
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


def _parse_competitive_blast_output(text: str) -> list[dict[str, str]]:
    """Report the best host and viral hit per query plus viral-minus-host score."""
    grouped: dict[str, dict[str, list[str]]] = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 7:
            continue
        query, subject = fields[:2]
        category = "viral" if subject.startswith("VIRUS|") else "host"
        current = grouped.setdefault(query, {}).get(category)
        if current is None or float(fields[6]) > float(current[6]):
            grouped[query][category] = fields
    rows: list[dict[str, str]] = []
    for query, hits in sorted(grouped.items()):
        viral, host = hits.get("viral"), hits.get("host")
        viral_score = float(viral[6]) if viral else 0.0
        host_score = float(host[6]) if host else 0.0
        rows.append(
            {
                "read": query,
                "top_viral_hit": viral[1] if viral else "",
                "viral_identity": viral[2] if viral else "",
                "viral_query_coverage": viral[4] if viral else "",
                "viral_evalue": viral[5] if viral else "",
                "viral_bitscore": viral[6] if viral else "",
                "top_host_hit": host[1] if host else "",
                "host_identity": host[2] if host else "",
                "host_query_coverage": host[4] if host else "",
                "host_evalue": host[5] if host else "",
                "host_bitscore": host[6] if host else "",
                "viral_minus_host_bitscore": f"{viral_score - host_score:.6g}",
                "low_complexity": "not_assessed",
            }
        )
    return rows


def competitive_blast_identity(
    reads_fasta: str,
    competitive_fasta: str,
    workdir: str,
    *,
    max_reads: int = 500,
    threads: int = 4,
    seed: int = 42,
    sampling_manifest: str | None = None,
) -> list[dict[str, str]]:
    """BLAST a deterministic read sample against combined host and target virus."""
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    sample = work / "reads_sample.fasta"
    sampled, total = sample_fasta_deterministic(reads_fasta, sample, max_reads, seed)
    if sampling_manifest:
        Path(sampling_manifest).write_text(
            json.dumps(
                {
                    "strategy": "lowest_sha256",
                    "seed": seed,
                    "sampled_reads": sampled,
                    "total_reads": total,
                    "sampling_fraction": sampled / total if total else 0.0,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    db = work / "competitive_db"
    _run(
        ["makeblastdb", "-in", competitive_fasta, "-dbtype", "nucl", "-out", str(db)],
        capture=True,
    )
    stdout = _run(
        [
            "blastn",
            "-query",
            str(sample),
            "-db",
            str(db),
            "-max_target_seqs",
            "20",
            "-num_threads",
            str(threads),
            "-outfmt",
            "6 qseqid sseqid pident length qcovs evalue bitscore",
        ],
        capture=True,
    ).decode("utf-8", errors="replace")
    rows = _parse_competitive_blast_output(stdout)
    complexity: dict[str, str] = {}
    identifier = ""
    sequence: list[str] = []
    with sample.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if identifier:
                    complexity[identifier] = str(_is_low_complexity("".join(sequence))).lower()
                identifier = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if identifier:
        complexity[identifier] = str(_is_low_complexity("".join(sequence))).lower()
    for row in rows:
        row["low_complexity"] = complexity.get(row["read"], "not_assessed")
    return rows


def _is_low_complexity(sequence: str) -> bool:
    """Conservative sequence-composition flag for deterministic BLAST samples."""
    bases = [base for base in sequence.upper() if base in "ACGT"]
    if not bases:
        return True
    counts = [bases.count(base) for base in "ACGT"]
    fractions = [count / len(bases) for count in counts if count]
    entropy = -sum(fraction * math.log2(fraction) for fraction in fractions)
    return max(fractions) >= 0.80 or entropy < 1.20


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


def sample_fasta_deterministic(src: str, dst: Path, max_records: int, seed: int) -> tuple[int, int]:
    """Select records with the lowest seeded SHA-256 values, independent of order."""
    if max_records < 1:
        raise ValueError("max_records must be >= 1")
    heap: list[tuple[int, int, str]] = []
    total = 0

    def consider(record: str) -> None:
        nonlocal total
        if not record:
            return
        score = int.from_bytes(hashlib.sha256(f"{seed}\0{record}".encode()).digest()[:8], "big")
        item = (-score, total, record)
        total += 1
        if len(heap) < max_records:
            heapq.heappush(heap, item)
        elif score < -heap[0][0]:
            heapq.heapreplace(heap, item)

    current = ""
    with open(src) as handle:
        for line in handle:
            if line.startswith(">"):
                consider(current)
                current = line
            else:
                current += line
        consider(current)
    selected = sorted(heap, key=lambda item: (-item[0], item[1]))
    with dst.open("w") as handle:
        for _neg_score, _ordinal, record in selected:
            handle.write(record)
    return len(selected), total


def have_tools(names: Iterable[str]) -> list[str]:
    """Return the subset of *names* that are NOT on PATH."""
    return [n for n in names if shutil.which(n) is None]
