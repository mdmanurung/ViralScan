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
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any, Optional, cast

from viralscan.anello_align import (
    Alignment,
    has_reagent,
    is_complex_body,
    iter_fasta,
    parse_sam_line,
)
from viralscan.anellovirus import anello_name_map
from viralscan.chemistry import cb_umi_geometry  # noqa: F401  (re-exported)
from viralscan.validation import tool_path
from viralscan.virus_catalog import merged_name_map
from viralscan.virus_grouping import group_genes_by_virus
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

log = logging.getLogger("viralscan")


def to_int(value: object, default: int = 0) -> int:
    """``int(value)`` for the number-or-numeric-string cells of a result row."""
    return int(value) if isinstance(value, (int, float, str)) else default


def to_float(value: object, default: float = 0.0) -> float:
    """``float(value)`` for the number-or-numeric-string cells of a result row."""
    return float(value) if isinstance(value, (int, float, str)) else default


#: Selector aliases -> every display name that virus can resolve to, in
#: preference order. One name per alias used to fail: the maps produce legacy
#: ("Human herpesvirus 6b"), catalogue ("Human betaherpesvirus 6A") or RefSeq
#: ("Human gammaherpesvirus 8") labels, and none matched "Human herpesvirus 6B"
#: or "Kaposi sarcoma-associated herpesvirus", so `--virus hhv6b/kshv` raised
#: (SW-17). Legacy fallback: a Run with a Virus Identity table resolves through
#: the table instead (see :data:`SELECTOR_HANDLES`).
VIRUS_ALIASES: dict[str, tuple[str, ...]] = {
    "ebv": ("Epstein-Barr virus", "Human gammaherpesvirus 4"),
    "hhv4": ("Epstein-Barr virus", "Human gammaherpesvirus 4"),
    "hsv1": ("Human herpesvirus 1", "Human alphaherpesvirus 1"),
    "hhv1": ("Human herpesvirus 1", "Human alphaherpesvirus 1"),
    "hsv2": ("Human herpesvirus 2", "Human alphaherpesvirus 2"),
    "hhv2": ("Human herpesvirus 2", "Human alphaherpesvirus 2"),
    "hhv6a": ("Human herpesvirus 6A", "Human betaherpesvirus 6A", "Human herpesvirus 6"),
    "hhv6b": ("Human herpesvirus 6B", "Human betaherpesvirus 6B"),
    "kshv": (
        "Kaposi sarcoma-associated herpesvirus",
        "Human herpesvirus 8",
        "Human gammaherpesvirus 8",
    ),
    "hhv8": (
        "Kaposi sarcoma-associated herpesvirus",
        "Human herpesvirus 8",
        "Human gammaherpesvirus 8",
    ),
    "ttv": ("Anelloviridae",),
}

ANELLOVIRIDAE = "Anelloviridae"

#: Short selector handles -> the virus taxid they name, for Runs with an identity
#: table. A handle is a typing convenience for a taxid; no display name is
#: involved, so it cannot drift from the names the outputs print.
SELECTOR_HANDLES: dict[str, str] = {
    "ebv": "taxid:10376",
    "hhv4": "taxid:10376",
    "hsv1": "taxid:10298",
    "hhv1": "taxid:10298",
    "hsv2": "taxid:10310",
    "hhv2": "taxid:10310",
    "hhv6a": "taxid:32603",
    "hhv6b": "taxid:32604",
    "kshv": "taxid:37296",
    "hhv8": "taxid:37296",
}
#: Selector handles that name a whole family rather than one virus key.
FAMILY_HANDLES: dict[str, str] = {"ttv": ANELLOVIRIDAE}


def _resolve_by_identity(
    query: str, table: VirusIdentityTable, genes: list[str]
) -> tuple[str, list[str]] | None:
    """Resolve a selector through the Virus Identity table, or ``None`` if it names nothing.

    Tried in order, the first tier with a hit wins: virus key (``taxid:10376``,
    ``genus:Betatorquevirus``), bare taxid, short handle, display name
    (``virus_name``, which is the curated common name when there is one), NCBI
    organism or species name, then family. A tier that names more than one virus
    raises rather than guessing, except a family, which is a deliberate union.
    Only genes in ``genes`` (the Run's viral genes) are returned.
    """
    wanted = set(genes)
    viral = [g for g in table.genes if g.viral and g.gene_id in wanted]
    folded = query.casefold()
    handle = SELECTOR_HANDLES.get(folded)
    tiers: list[tuple[str, Callable[[GeneIdentity], bool]]] = [
        ("virus key", lambda g: g.virus_key.casefold() == folded),
        ("taxid", lambda g: bool(g.taxid) and g.taxid == query),
        ("handle", lambda g: handle is not None and g.virus_key == handle),
        ("name", lambda g: g.virus_name.casefold() == folded),
        (
            "organism",
            lambda g: folded in (g.organism.casefold(), g.species.casefold()) and bool(folded),
        ),
    ]
    family = FAMILY_HANDLES.get(folded, query)
    for tier, match in tiers:
        hits = [g for g in viral if match(g)]
        if not hits:
            continue
        names = sorted({g.virus_name for g in hits})
        if len({g.virus_key for g in hits}) > 1:
            raise ValueError(
                f"Selector {query!r} is ambiguous ({tier}): it names {len(names)} viruses, "
                f"{names[:6]}. Use a virus key such as 'taxid:<n>'."
            )
        return hits[0].virus_name, [g.gene_id for g in hits]
    hits = [g for g in viral if g.family.casefold() == family.casefold()]
    if hits:
        return hits[0].family, [g.gene_id for g in hits]
    return None


def _anellovirus_group_names() -> set[str]:
    """Every display name an Anelloviridae genome resolves to.

    The literal "Anelloviridae" group holds only genomes with no assigned genus,
    so the family selector must also cover each genus and the bundled RefSeq
    "Torque teno virus" label.
    """
    return set(anello_name_map().values()) | {ANELLOVIRIDAE, "Torque teno virus"}


def resolve_viral_target(
    selector: str,
    viral_gene_ids: Iterable[str],
    *,
    detected_virus_names: Iterable[str] = (),
    identity: VirusIdentityTable | None = None,
) -> tuple[str, list[str]]:
    """Resolve one exact accession/gene, alias, or detected canonical call.

    With ``identity`` (the Run's Virus Identity table) a virus selector resolves
    through the table by key, taxid, handle, name, organism or family; without
    it the legacy name maps and :data:`VIRUS_ALIASES` are used.

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

    if identity is not None:
        hit = _resolve_by_identity(query, identity, genes)
        if hit is not None:
            return hit
        choices = sorted(
            {g.virus_name for g in identity.genes if g.viral}.union(detected_virus_names)
        )
        raise ValueError(
            f"No exact viral target matches {selector!r}. Use an accession/gene ID, canonical "
            f"label, or detected call. Available calls include: {choices[:12]}"
        )

    name_map = merged_name_map()
    groups, _ = group_genes_by_virus(genes, name_map)
    canonical_by_lower = {name.casefold(): name for name in groups}
    canonical_by_lower.update(
        {str(name).casefold(): str(name) for name in detected_virus_names if str(name).strip()}
    )
    alias_to_name = {key.casefold(): value for key, value in name_map.items()}
    canonical = None
    for candidate in VIRUS_ALIASES.get(query.casefold(), ()):
        canonical = canonical_by_lower.get(candidate.casefold())
        if canonical or candidate == ANELLOVIRIDAE:
            canonical = canonical or ANELLOVIRIDAE
            break
    canonical = (
        canonical or alias_to_name.get(query.casefold()) or canonical_by_lower.get(query.casefold())
    )
    if canonical == ANELLOVIRIDAE:
        family = _anellovirus_group_names()
        family_genes = [gene for name in sorted(groups) if name in family for gene in groups[name]]
        if family_genes:
            return ANELLOVIRIDAE, family_genes
    if canonical and canonical in groups and groups[canonical]:
        return canonical, groups[canonical]

    choices = sorted(set(groups).union(detected_virus_names))
    raise ValueError(
        f"No exact viral target matches {selector!r}. Use an accession/gene ID, canonical "
        f"label, or detected call. Available calls include: {choices[:12]}"
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


def _open_maybe_gzip(path: str | Path, mode: str = "rt") -> IO[str]:
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


#: ViralScan strand vocabulary -> kallisto bus flag. ``None`` (no flag) means
#: "kallisto's default", which is version-dependent: the per-technology
#: forward default documented since kallisto v0.48.0 is dead code in v0.52.0
#: (the technology block sets ``opt.strand_specific`` before the finalizer can
#: apply it, so the run goes unstranded), so an unstranded primary run must be
#: replayed without a flag and a stranded one with its explicit flag.
_KALLISTO_STRAND_FLAGS = {
    "forward": "--fr-stranded",
    "reverse": "--rf-stranded",
    "unstranded": "--unstranded",
}


def replay_ec_path(workdir: str | Path) -> Path:
    """The ``matrix.ec`` that ``replay_exact_target_bus`` wrote beside its BUS file."""
    return Path(workdir) / "lineage_bus" / "matrix.ec"


def replay_exact_target_bus(
    *,
    index: str,
    technology: str,
    r1_path: str,
    r2_path: str,
    target_transcripts: Iterable[str],
    workdir: str,
    threads: int,
    whitelist: str | None = None,
    strand: str | None = None,
) -> Path:
    """Re-pseudoalign with read-number flags and capture only target transcripts.

    The replay must mirror the primary ``kb count`` quantification. Pass the
    primary run's on-list as ``whitelist`` so the same barcode correction runs
    before capture — without it, reads kb discarded via correction enter the
    evidence set. Pass the primary run's ``strand`` so a stranded run is not
    replayed with kallisto's (version-dependent) default.

    The capture and every later EC lookup must use the **replay's own**
    ``lineage_bus/matrix.ec``: kallisto numbers equivalence classes in the order it
    meets them, which differs between multithreaded runs, so the primary run's
    ``matrix.ec`` names different classes than the replay's BUS file does. Read it
    back with :func:`replay_ec_path`.
    """
    missing = have_tools(["kallisto", "bustools"])
    if missing:
        raise RuntimeError(f"Exact read lineage requires: {', '.join(missing)}")
    work = Path(workdir)
    bus_dir = work / "lineage_bus"
    bus_dir.mkdir(parents=True, exist_ok=True)
    capture_list = work / "target_transcripts.txt"
    capture_list.write_text("\n".join(sorted(set(target_transcripts))) + "\n")
    kallisto, bustools = tool_path("kallisto"), tool_path("bustools")
    if kallisto is None or bustools is None:  # have_tools() passed, so only a PATH race
        raise RuntimeError("Exact read lineage requires: kallisto, bustools")
    strand_flag = _KALLISTO_STRAND_FLAGS.get(strand) if strand else None
    if strand and strand_flag is None:
        raise ValueError(
            f"Unknown strand mode {strand!r}; expected one of {sorted(_KALLISTO_STRAND_FLAGS)}."
        )
    bus_cmd = [
        kallisto,
        "bus",
        "-i",
        index,
        "-o",
        str(bus_dir),
        "-x",
        technology,
        "-t",
        str(threads),
    ]
    if strand_flag:
        bus_cmd.append(strand_flag)
    bus_cmd += ["-n", r1_path, r2_path]
    _run(bus_cmd)
    quant_bus = bus_dir / "output.bus"
    if whitelist is not None:
        # Mirror kb count's chain: sort -> correct -> sort, then capture the
        # corrected records. bustools cannot read a gzipped on-list (it reads
        # the compressed bytes as barcodes), so decompress to the workdir.
        whitelist_path = Path(whitelist)
        if whitelist_path.suffix == ".gz":
            plain = work / "replay_whitelist.txt"
            with gzip.open(whitelist_path, "rt") as src, plain.open("w") as dst:
                dst.write(src.read())
            whitelist_path = plain
        sorted_raw = bus_dir / "output.sorted.bus"
        _run([bustools, "sort", "-t", str(threads), "-o", str(sorted_raw), str(quant_bus)])
        corrected = bus_dir / "output.corrected.bus"
        _run(
            [
                bustools,
                "correct",
                "-o",
                str(corrected),
                "-w",
                str(whitelist_path),
                str(sorted_raw),
            ]
        )
        quant_bus = bus_dir / "output.corrected.sorted.bus"
        _run([bustools, "sort", "-t", str(threads), "-o", str(quant_bus), str(corrected)])
    captured = work / "target.bus"
    _run(
        [
            bustools,
            "capture",
            "-s",
            "-c",
            str(capture_list),
            "-e",
            str(bus_dir / "matrix.ec"),
            "-t",
            str(bus_dir / "transcripts.txt"),
            "-o",
            str(captured),
            str(quant_bus),
        ]
    )
    by_flag = work / "target.by_flag.bus"
    _run([bustools, "sort", "--flags", "-t", str(threads), "-o", str(by_flag), str(captured)])
    flagged_text = work / "target.by_flag.bus.txt"
    _run([bustools, "text", "-f", "-o", str(flagged_text), str(by_flag)])
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
            seen: set[str] = set()
            with _open_maybe_gzip(source) as handle:
                for line in handle:
                    if line.startswith(">"):
                        name, *description = line[1:].rstrip().split(maxsplit=1)
                        if name in seen:
                            raise ValueError(
                                f"duplicate FASTA ID {name!r} in {source}; samtools faidx "
                                "would fail. Deduplicate the reference first."
                            )
                        seen.add(name)
                        suffix = f" {description[0]}" if description else ""
                        target.write(f">{prefix}|{name}{suffix}\n")
                    else:
                        target.write(line)
    return str(out)


# Explicit competitive references are opt-in; the HOST/VIRUS legacy format stays intact.
COMPETITOR_PREFIXES = {
    "target": "TARGET",
    "same_reporting_group": "SAMEGROUP",
    "related_virus": "RELATED",
    "host": "HOST",
    "decoy": "DECOY",
}
COMPETITOR_COLUMNS = (
    "reference_id",
    "source_role",
    "class",
    "virus_key",
    "reporting_group",
    "accession",
    "display_name",
    "sequence_sha256",
    "source",
    "relationship",
)

MOLECULE_COMPETITION_FIELDS = (
    "cell_barcode",
    "umi",
    "selected_reference_id",
    "selected_class",
    "n_candidate_alignments",
    "candidate_classes",
    "mapping_quality",
    "identity",
    "alignment_span",
    "ambiguity_class",
    "status",
    "candidate_search",
)
COMPETITOR_BLAST_FIELDS = (
    "read",
    "diagnostic_class",
    "top_class",
    "tie_classes",
    "top_tied_hits",
    "n_top_tied_hits",
    "search_status",
    "low_complexity",
    "reason",
    *(
        f"{key}"
        for cls in ("target", "same_group", "related", "host", "decoy")
        for key in (
            f"best_{cls}_hit",
            f"{cls}_winner_hits",
            f"{cls}_bitscore",
            f"{cls}_identity",
            f"{cls}_query_coverage",
        )
    ),
    "no_target_support",
    *(f"target_minus_{cls}_bitscore" for cls in ("same_group", "related", "host", "decoy")),
)


def reference_class(reference: str) -> str:
    """Class from an encoded reference identifier; legacy VIRUS and unprefixed ids stay virus."""
    prefix = reference.split("|", 1)[0]
    return next((cls for cls, encoded in COMPETITOR_PREFIXES.items() if prefix == encoded), "virus")


def _fasta_digests(path: str) -> dict[str, str]:
    """Hash DNA-IUPAC sequence bytes after uppercase/whitespace removal, streaming each record."""
    records: dict[str, str] = {}
    name = None
    digest = hashlib.sha256()
    length = 0
    with _open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    if not length:
                        raise ValueError(f"empty sequence {name!r} in {path}")
                    records[name] = digest.hexdigest()
                header = line[1:].split()
                if not header or "|" in header[0]:
                    raise ValueError(f"missing/reserved FASTA reference ID in {path}")
                name = header[0]
                if name in records:
                    raise ValueError(f"duplicate reference ID {name!r} in {path}")
                digest, length = hashlib.sha256(), 0
            else:
                sequence = "".join(line.split()).upper()
                if not sequence:
                    continue
                if name is None or set(sequence) - set("ACGTRYSWKMBDHVN"):
                    raise ValueError(f"invalid DNA-IUPAC sequence/header in {path}")
                digest.update(sequence.encode("ascii"))
                length += len(sequence)
    if name is not None:
        if not length:
            raise ValueError(f"empty sequence {name!r} in {path}")
        records[name] = digest.hexdigest()
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def validate_competitor_manifest(
    path: str | Path,
    *,
    target_fasta: str,
    host_fasta: str,
    competitor_fasta: str,
    identity: VirusIdentityTable | None,
    target_genes: Iterable[str],
) -> dict[str, dict[str, str]]:
    """Validate merged-reference TSV before any replay/output replacement.

    The ten required columns are :data:`COMPETITOR_COLUMNS`; source_role is
    target/host/competitor. Hashes normalize uppercase and all sequence whitespace
    (DNA IUPAC alphabet, gzip supported). Target accessions/keys must agree with
    the selected Run identity. same_reporting_group uses that exact reporting key;
    related_virus is prespecified explicitly and must have a different key/group.
    Returns validated rows indexed by the unique unencoded FASTA reference ID.
    """
    if identity is None:
        raise ValueError("competitor mode requires the run VirusIdentityTable")
    target_set = set(target_genes)
    selected = [g for g in identity.genes if g.gene_id in target_set and g.viral]
    keys = {g.virus_key for g in selected if g.virus_key}
    accessions = {g.genome_accession for g in selected if g.genome_accession}
    if not selected or not keys or not accessions:
        raise ValueError("selected target requires one resolved identity key and genome accessions")
    records: dict[str, tuple[str, str]] = {}
    for role, source in (
        ("target", target_fasta),
        ("host", host_fasta),
        ("competitor", competitor_fasta),
    ):
        for name, digest in _fasta_digests(source).items():
            if name in records:
                raise ValueError(f"duplicate merged FASTA reference ID {name!r}")
            records[name] = role, digest
    rows: dict[str, dict[str, str]] = {}
    sequence_classes: dict[str, str] = {}
    with _open_maybe_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if (
            reader.fieldnames is None
            or len(reader.fieldnames) != len(set(reader.fieldnames))
            or not set(COMPETITOR_COLUMNS) <= set(reader.fieldnames)
        ):
            raise ValueError(f"competitor manifest requires unique columns {COMPETITOR_COLUMNS}")
        for raw in reader:
            if None in raw or any(raw.get(c) is None for c in COMPETITOR_COLUMNS):
                raise ValueError("malformed competitor manifest row")
            row = {c: raw[c].strip() for c in COMPETITOR_COLUMNS}
            name, cls = row["reference_id"], row["class"]
            if name in rows or name not in records:
                raise ValueError(f"duplicate/unknown competitor reference ID {name!r}")
            role, digest = records[name]
            if cls not in COMPETITOR_PREFIXES or row["source_role"] != role:
                raise ValueError(f"invalid class/source role for {name!r}")
            allowed = {
                "target": {"target"},
                "host": {"host"},
                "competitor": {"same_reporting_group", "related_virus", "decoy"},
            }
            if cls not in allowed[role] or row["sequence_sha256"].lower() != digest:
                raise ValueError(f"class/sequence hash disagrees with FASTA source for {name!r}")
            if not all(row[c] for c in ("display_name", "source", "relationship", "accession")):
                raise ValueError(f"missing reference provenance for {name!r}")
            if cls in {"target", "same_reporting_group"} and (
                row["virus_key"] not in keys or row["reporting_group"] != row["virus_key"]
            ):
                raise ValueError(
                    f"target/same-group identity differs from selected run key for {name!r}"
                )
            if cls == "target" and not any(
                g.genome_accession == row["accession"] and g.virus_key == row["virus_key"]
                for g in selected
            ):
                raise ValueError(
                    f"target accession {row['accession']!r} is absent from selected run identity"
                )
            if cls == "related_virus" and (
                not row["virus_key"]
                or not row["reporting_group"]
                or row["virus_key"] in keys
                or row["reporting_group"] in keys
            ):
                raise ValueError(
                    f"related reference {name!r} requires a distinct explicit virus/reporting key"
                )
            if cls == "host" and (row["virus_key"] or row["reporting_group"]):
                raise ValueError(f"host reference {name!r} cannot carry viral grouping")
            if digest in sequence_classes and sequence_classes[digest] != cls:
                raise ValueError(
                    f"identical sequence has contradictory classes for {name!r}; declare an ambiguity-compatible reference identity first"
                )
            sequence_classes[digest] = cls
            rows[name] = row
    if set(rows) != set(records) or not any(row["class"] == "target" for row in rows.values()):
        raise ValueError(
            "competitor manifest must cover every merged FASTA record exactly once and include target"
        )
    return rows


def write_competitor_fasta(
    host_fasta: str,
    target_fasta: str,
    competitor_fasta: str,
    records: dict[str, dict[str, str]],
    output: str,
) -> str:
    """Write class-encoded merged FASTA from records returned by validate_competitor_manifest."""
    with open(output, "w") as target:
        for source in (host_fasta, target_fasta, competitor_fasta):
            with _open_maybe_gzip(source) as handle:
                for line in handle:
                    if line.startswith(">"):
                        name = line[1:].split()[0]
                        target.write(f">{COMPETITOR_PREFIXES[records[name]['class']]}|{name}\n")
                    else:
                        target.write("".join(line.split()).upper() + "\n")
    return output


def parse_competitor_blast_output(
    text: str, query_ids: Iterable[str], *, tie_delta: float = 0.0, search_complete: bool = True
) -> list[dict[str, object]]:
    """One deterministic diagnostic per sampled query, including no hits.

    Reduce HSPs by subject using highest bitscore, then identity/query coverage,
    alignment length and lexical row as deterministic tie breaks. qcovs is BLAST's
    per-subject query coverage; bitscores are per best HSP, never summed. Preserve
    every class winner and distinct subject within tie_delta of the global maximum.
    """
    if not math.isfinite(tie_delta) or tie_delta < 0:
        raise ValueError("blast_tie_delta must be finite and >= 0 bitscore units")
    queries = list(query_ids)
    if len(set(queries)) != len(queries):
        raise ValueError("sampled query IDs must be unique")
    hits: dict[str, dict[str, list[str]]] = {q: {} for q in queries}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 7 or fields[0] not in hits:
            raise ValueError("malformed/unknown query in competitive BLAST diagnostics")
        cls = reference_class(fields[1])
        if cls not in COMPETITOR_PREFIXES:
            raise ValueError(f"unencoded competitive BLAST subject {fields[1]!r}")
        numbers = [float(value) for value in fields[2:]]
        if not all(math.isfinite(v) for v in numbers) or numbers[-1] < 0:
            raise ValueError("competitive BLAST metrics must be finite and bitscore nonnegative")
        current = hits[fields[0]].get(fields[1])

        def rank(f: list[str]) -> tuple[object, ...]:
            return (float(f[6]), float(f[2]), float(f[4]), float(f[3]), tuple(f))

        if current is None or rank(fields) > rank(current):
            hits[fields[0]][fields[1]] = fields
    labels = {
        "target": "target_specific",
        "same_reporting_group": "within_group_support",
        "related_virus": "related_preferred",
        "host": "host_preferred",
        "decoy": "decoy_preferred",
    }
    pairs = {
        frozenset({"target", "same_reporting_group"}): "within_group_ambiguous",
        frozenset({"target", "related_virus"}): "related_virus_ambiguous",
        frozenset({"target", "host"}): "host_competitive",
    }
    rows = []
    for query in queries:
        subjects = hits[query]
        maximum = max((float(f[6]) for f in subjects.values()), default=None)
        tied = sorted(
            s
            for s, f in subjects.items()
            if maximum is not None and maximum - float(f[6]) <= tie_delta
        )
        classes = sorted({reference_class(s) for s in tied})
        diagnostic = (
            "no_hits"
            if not classes
            else (
                labels[classes[0]]
                if len(classes) == 1
                else pairs.get(frozenset(classes), "multi_class_ambiguous")
            )
        )
        row: dict[str, object] = dict(
            read=query,
            diagnostic_class=diagnostic if search_complete else "incomplete_search",
            top_class=classes[0] if len(classes) == 1 else None,
            tie_classes=";".join(classes),
            top_tied_hits=";".join(tied),
            n_top_tied_hits=len(tied),
            search_status="complete_subject_search" if search_complete else "incomplete",
            low_complexity="not_assessed",
            reason=None,
        )
        for cls in COMPETITOR_PREFIXES:
            short = {"same_reporting_group": "same_group", "related_virus": "related"}.get(cls, cls)
            candidates = sorted(s for s in subjects if reference_class(s) == cls)
            best = max((float(subjects[s][6]) for s in candidates), default=None)
            winners = [s for s in candidates if float(subjects[s][6]) == best]
            first = subjects[winners[0]] if winners else None
            row.update(
                {
                    f"best_{short}_hit": winners[0] if winners else None,
                    f"{short}_winner_hits": ";".join(winners),
                    f"{short}_bitscore": best,
                    f"{short}_identity": float(first[2]) if first else None,
                    f"{short}_query_coverage": float(first[4]) if first else None,
                }
            )
        row["no_target_support"] = row["target_bitscore"] is None
        for cls in ("same_group", "related", "host", "decoy"):
            t, c = row["target_bitscore"], row[f"{cls}_bitscore"]
            row[f"target_minus_{cls}_bitscore"] = (
                float(t) - float(c)
                if isinstance(t, (int, float)) and isinstance(c, (int, float))
                else None
            )
        rows.append(row)
    return rows


def competitor_summary(
    blast_rows: Iterable[dict[str, object]],
    lineage_rows: Iterable[dict[str, object]],
    *,
    target_id: str,
    manifest_sha256: str,
) -> dict[str, object]:
    """Read fractions use sampled queries only; molecule counts join exact read-number lineage.

    Molecules here are corrected (CB,UMI) candidates, without selecting a winning
    reference. No sampled-molecule classification is inferred from heterogeneous
    read evidence; reference-specific representatives are in molecule_competition.tsv.
    """
    blast, lineage = list(blast_rows), list(lineage_rows)
    by_number = {
        to_int(r["read_number"]): (str(r["cb"]), str(r["ub"]), str(r["read_id"])) for r in lineage
    }
    all_molecules = {(str(r["cb"]), str(r["ub"])) for r in lineage}
    sampled_numbers, sampled_molecules = set(), set()
    unresolved = 0
    for row in blast:
        read = str(row["read"])
        parts = read.split("_", 2)
        try:
            number = int(parts[2].split("|", 1)[0])
        except (IndexError, ValueError):
            unresolved += 1
            continue
        source = by_number.get(number)
        if source is None or read != f"{source[0]}_{source[1]}_{number}|{source[2]}":
            unresolved += 1
            continue
        sampled_numbers.add(number)
        sampled_molecules.add(source[:2])
    n = len(blast)
    out: dict[str, object] = dict(
        target_id=target_id,
        n_molecules=len(all_molecules),
        n_reads_total=len(lineage),
        n_reads_sampled_blast=n,
        n_sampled_molecules=len(sampled_molecules),
        n_unsampled_reads=len(lineage) - len(sampled_numbers),
        n_unsampled_molecules=len(all_molecules - sampled_molecules),
        n_unresolved_sampled_reads=unresolved,
        n_no_hit_reads=sum(r["diagnostic_class"] == "no_hits" for r in blast),
        n_failed_reads=sum(r["diagnostic_class"] in {"failed", "incomplete_search"} for r in blast),
        blast_status="assessed" if n else "not_assessed",
        blast_fraction_denominator="sampled_reads",
        molecule_unit="corrected_cb_umi_without_selected_reference",
        competitor_manifest_sha256=manifest_sha256,
    )
    for cls in (
        "target_specific",
        "within_group_ambiguous",
        "related_virus_ambiguous",
        "host_competitive",
        "multi_class_ambiguous",
        "within_group_support",
        "related_preferred",
        "host_preferred",
        "decoy_preferred",
        "no_hits",
    ):
        out[f"fraction_{cls}"] = sum(r["diagnostic_class"] == cls for r in blast) / n if n else None
    out["fraction_low_complexity"] = (
        sum(r.get("low_complexity") == "true" for r in blast) / n if n else None
    )
    for cls in ("related", "host"):
        values = sorted(
            to_float(r[f"target_minus_{cls}_bitscore"])
            for r in blast
            if r.get(f"target_minus_{cls}_bitscore") is not None
        )
        out[f"median_target_minus_{cls}_bitscore"] = (
            (values[(len(values) - 1) // 2] + values[len(values) // 2]) / 2 if values else None
        )
    if out["n_failed_reads"]:
        out["blast_status"] = "failed"
        for key in out:
            if key.startswith(("fraction_", "median_")):
                out[key] = None
    return out


def molecule_competition_rows(
    sam_text: str, lineage_rows: Iterable[dict[str, object]]
) -> list[dict[str, object]]:
    """Selected corrected (CB,UMI,reference) representatives and emitted candidate classes.

    Secondary records are candidate evidence only. Candidate counts describe emitted
    minimap2 records, never exhaustive search. Unaligned lineage molecules stay unresolved.
    """
    candidates: dict[tuple[str, str], list[list[str]]] = {}
    for line in sam_text.splitlines():
        fields = line.split("\t")
        if (
            len(fields) < 11
            or line.startswith("@")
            or int(fields[1]) & (EXCLUDE_FLAGS & ~0x100)
            or fields[2] == "*"
        ):
            continue
        key = _cb_umi(fields[0])
        if key:
            candidates.setdefault(key, []).append(fields)
    rows = []
    represented = set()
    for line in select_molecule_representatives(sam_text.splitlines())[0]:
        fields = line.split("\t")
        key = _cb_umi(fields[0])
        if key is None:
            continue
        represented.add(key)
        emitted = candidates.get(key, [])
        rows.append(
            dict(
                cell_barcode=key[0],
                umi=key[1],
                selected_reference_id=fields[2],
                selected_class=reference_class(fields[2]),
                n_candidate_alignments=len(emitted),
                candidate_classes=";".join(sorted({reference_class(f[2]) for f in emitted})),
                mapping_quality=int(fields[4]),
                identity=_identity_percent(fields),
                alignment_span=_cigar_ref_span(fields[5]),
                ambiguity_class="emitted_multiclass"
                if len({reference_class(f[2]) for f in emitted}) > 1
                else "candidate_search_not_exhaustive",
                status="selected",
                candidate_search="emitted_alignments_only",
            )
        )
    for key in sorted({(str(r["cb"]), str(r["ub"])) for r in lineage_rows} - represented):
        emitted = candidates.get(key, [])
        rows.append(
            dict(
                cell_barcode=key[0],
                umi=key[1],
                selected_reference_id=None,
                selected_class=None,
                n_candidate_alignments=len(emitted),
                candidate_classes=";".join(sorted({reference_class(f[2]) for f in emitted})),
                mapping_quality=None,
                identity=None,
                alignment_span=None,
                ambiguity_class="unresolved",
                status="unresolved",
                candidate_search="emitted_alignments_only",
            )
        )
    return sorted(
        rows, key=lambda r: (str(r["cell_barcode"]), str(r["umi"]), str(r["selected_reference_id"]))
    )


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
    stdout = _run(["samtools", "coverage", *_SAMTOOLS_COVERAGE_POLICY, bam], capture=True).decode(
        "utf-8", errors="replace"
    )
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


@dataclass
class _Tally:
    """Running per-reference (or per-cell) alignment counts for the QC tables."""

    reads: int = 0
    reverse: int = 0
    complex: int = 0
    reagent: int = 0
    mapq: list[int] = field(default_factory=list)
    identities: list[float] = field(default_factory=list)
    starts: list[int] = field(default_factory=list)
    cells: set[str] = field(default_factory=set)
    molecules: set[tuple[str, str]] = field(default_factory=set)


def _identity_percent(fields: list[str]) -> Optional[float]:
    """``100 * (aligned - NM) / aligned`` from a SAM record, or None without NM/CIGAR."""
    nm = next((tag for tag in fields[11:] if tag.startswith("NM:i:")), None)
    aligned = _cigar_ref_span(fields[5])
    if nm and aligned:
        return max(0.0, 100.0 * (aligned - int(nm.split(":")[-1])) / aligned)
    return None


def _alignment_qc_from_text(
    header_text: str, sam_text: str, depth_text: str
) -> list[dict[str, object]]:
    """Compute per-reference specificity, depth, and hotspot diagnostics."""
    lengths: dict[str, int] = {}
    for line in header_text.splitlines():
        if not line.startswith("@SQ"):
            continue
        sq = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
        if "SN" in sq and "LN" in sq:
            lengths[sq["SN"]] = int(sq["LN"])

    depths: dict[str, list[tuple[int, int]]] = {}
    for line in depth_text.splitlines():
        cols = line.split("\t")
        if len(cols) >= 3:
            depths.setdefault(cols[0], []).append((int(cols[1]), int(cols[2])))

    stats: dict[str, _Tally] = {}
    total_mapped = host_mapped = 0
    for line in sam_text.splitlines():
        fields = _primary_fields(line)
        if fields is None:
            continue
        flag = int(fields[1])
        reference = fields[2]
        total_mapped += 1
        host_mapped += int(reference.startswith("HOST|"))
        record = stats.setdefault(reference, _Tally())
        record.reads += 1
        record.reverse += int(bool(flag & 0x10))
        record.mapq.append(int(fields[4]))
        record.starts.append(int(fields[3]))
        cbumi = _cb_umi(fields[0])
        if cbumi:
            record.cells.add(cbumi[0])
            record.molecules.add(cbumi)
        # Labels, never filters. read_seq() restores sequencing orientation.
        aln = parse_sam_line(line)
        if aln is not None and aln.seq != "*":
            seq = aln.read_seq()
            record.complex += int(is_complex_body(seq, aligned=aln.query_span()))
            record.reagent += int(has_reagent(seq))
        identity = _identity_percent(fields)
        if identity is not None:
            record.identities.append(identity)

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
        starts = record.starts
        start_counts: dict[int, int] = {}
        window_counts: dict[int, int] = {}
        for start in starts:
            start_counts[start] = start_counts.get(start, 0) + 1
            window = start // 50
            window_counts[window] = window_counts.get(window, 0) + 1
        counts = list(start_counts.values())
        entropy = -sum((count / len(starts)) * math.log2(count / len(starts)) for count in counts)
        reads = record.reads
        identities = record.identities
        mapq = record.mapq
        rows.append(
            {
                "reference": reference,
                "reference_class": reference_class(reference),
                "reads": reads,
                "molecules": len(record.molecules),
                "cells": len(record.cells),
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
                "reverse_strand_fraction": record.reverse / reads,
                "mean_mapping_quality": sum(mapq) / len(mapq),
                "mean_identity": sum(identities) / len(identities) if identities else "",
                "host_competitive_fraction": host_mapped / total_mapped if total_mapped else 0.0,
                "complex_body_fraction": record.complex / reads,
                "reagent_fraction": record.reagent / reads,
            }
        )
    return rows


def alignment_qc_table(bam: str) -> list[dict[str, object]]:
    """Compute alignment QC from samtools header, alignments, and covered depths."""
    header = _run(["samtools", "view", "-H", bam], capture=True).decode(errors="replace")
    sam = _run(["samtools", "view", bam], capture=True).decode(errors="replace")
    depth = _run(["samtools", "depth", *_SAMTOOLS_DEPTH_POLICY, bam], capture=True).decode(
        errors="replace"
    )
    return _alignment_qc_from_text(header, sam, depth)


def _per_cell_qc_from_text(sam_text: str) -> list[dict[str, object]]:
    """Summarise primary competitive alignments per cell and reference class."""
    stats: dict[tuple[str, str], _Tally] = {}
    for line in sam_text.splitlines():
        fields = _primary_fields(line)
        if fields is None:
            continue
        flag = int(fields[1])
        cbumi = _cb_umi(fields[0])
        if not cbumi:
            continue
        category = reference_class(fields[2])
        record = stats.setdefault((cbumi[0], category), _Tally())
        record.reads += 1
        record.reverse += int(bool(flag & 0x10))
        record.mapq.append(int(fields[4]))
        record.molecules.add(cbumi)
        identity = _identity_percent(fields)
        if identity is not None:
            record.identities.append(identity)

    rows: list[dict[str, object]] = []
    for (cell, category), record in sorted(stats.items()):
        reads = record.reads
        mapq = record.mapq
        identities = record.identities
        rows.append(
            {
                "cell_barcode": cell,
                "reference_class": category,
                "reads": reads,
                "molecules": len(record.molecules),
                "reverse_strand_fraction": record.reverse / reads,
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
    text = _run(["samtools", "depth", *_SAMTOOLS_DEPTH_POLICY, bam], capture=True).decode(
        errors="replace"
    )
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
    explicit_classes = any(
        reference_class(str(row["reference"]))
        in {"target", "same_reporting_group", "related_virus", "decoy"}
        for rows in layers.values()
        for row in rows
    )
    references = sorted(
        {
            str(row["reference"])
            for rows in layers.values()
            for row in rows
            if explicit_classes or str(row["reference"]).startswith("VIRUS|")
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
                [to_int(row["position"]) for row in selected],
                [to_int(row["depth"]) for row in selected],
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
    qc = [row for row in qc_rows if row.get("reference_class") in {"virus", "target"}]
    blast = list(blast_rows)
    lineage = list(lineage_rows)
    max_hotspot = max((to_float(row.get("max_50bp_window_fraction")) for row in qc), default=0)
    max_breadth = max((to_float(row.get("breadth_1x")) for row in qc), default=0)
    host_fraction = max((to_float(row.get("host_competitive_fraction")) for row in qc), default=0)
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
    rows: list[dict[str, object]] = [
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
    assessment = {
        "host_homology": bool(qc or blast),
        "low_complexity": bool(blast),
        "sibling_or_host_ambiguity": bool(lineage),
        "coverage_hotspot": bool(qc),
        "eve_or_integration_like_hotspot": bool(qc),
    }
    for row in rows:
        if row["flag"] in assessment and not assessment[str(row["flag"])]:
            row.update(status="not_assessed", metric="")
    if any("diagnostic_class" in row for row in blast):
        for name, classes in (
            ("related_virus_competition", {"related_preferred", "related_virus_ambiguous"}),
            ("within_group_ambiguity", {"within_group_ambiguous"}),
        ):
            value = sum(r.get("diagnostic_class") in classes for r in blast) / len(blast)
            rows.append(
                dict(
                    flag=name,
                    status="flagged" if value >= 0.25 else "not_flagged",
                    metric=value,
                    threshold=0.25,
                    basis="sampled-read competitive BLAST",
                    interpretation="diagnostic_only",
                )
            )
        host_value = sum(
            r.get("diagnostic_class") in {"host_preferred", "host_competitive"} for r in blast
        ) / len(blast)
        host = next(r for r in rows if r["flag"] == "host_homology")
        host.update(
            metric=max(host_fraction, host_value),
            status="flagged" if max(host_fraction, host_value) >= 0.25 else "not_flagged",
        )
        if any(
            r.get("diagnostic_class") in {"no_hits", "incomplete_search", "failed"} for r in blast
        ):
            for row in rows:
                if (
                    row["flag"]
                    in {"host_homology", "related_virus_competition", "within_group_ambiguity"}
                    and row["status"] == "not_flagged"
                ):
                    row.update(
                        status="not_assessed",
                        basis="no-hit/incomplete/failed sampled-query evidence",
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


#: SAM FLAG bits excluded from every diagnostic layer (EVID-CORR-01): unmapped
#: 0x4, secondary 0x100, QC-fail 0x200, duplicate 0x400, supplementary 0x800.
#: samtools' own default omits supplementary, so ``coverage``/``depth`` are told
#: explicitly (see ``_SAMTOOLS_*_POLICY``) and the Python tallies use the same mask.
EXCLUDE_FLAGS = 0xF04

#: samtools policy matching the Python tallies: the FLAG mask above, no MAPQ or
#: base-quality filtering, deletions not counted as depth (``depth -J`` unset).
#: ``coverage`` replaces its default filter via ``--ff``; ``depth`` adds to its
#: default (UNMAP, SECONDARY, QCFAIL, DUP) via ``-G``, so only supplementary is named.
_SAMTOOLS_COVERAGE_POLICY = ["--ff", "0xF04", "-q", "0", "-Q", "0"]
_SAMTOOLS_DEPTH_POLICY = ["-G", "0x800", "-q", "0", "-Q", "0"]


def _primary_fields(line: str) -> Optional[list[str]]:
    """SAM fields of a primary mapped alignment, else None.

    Malformed records (fewer than 11 fields, non-integer FLAG/POS/MAPQ), records
    carrying any ``EXCLUDE_FLAGS`` bit and records without a reference return None.
    Every diagnostic layer filters through here so they cannot disagree.
    """
    if not line or line.startswith("@"):
        return None
    fields = line.split("\t")
    if len(fields) < 11:
        return None
    try:
        flag = int(fields[1])
        int(fields[3])
        int(fields[4])
    except ValueError:
        return None
    if flag & EXCLUDE_FLAGS or fields[2] == "*":
        return None
    return fields


def _int_tag(fields: list[str], name: str) -> Optional[int]:
    prefix = f"{name}:i:"
    for tag in fields[11:]:
        if tag.startswith(prefix):
            try:
                return int(tag[len(prefix) :])
            except ValueError:
                return None
    return None


def _representative_rank(fields: list[str], line: str) -> tuple[object, ...]:
    """Sort key, best first, for choosing a molecule's representative alignment.

    Highest MAPQ (255 = unavailable ranks below 0), highest AS, lowest NM (a
    missing AS/NM ranks after any present value), longest aligned query span
    (M/I/=/X bases; clips excluded), then lexical qname, POS, FLAG, CIGAR and the
    full record, so the choice never depends on input order.
    """
    mapq = int(fields[4])
    as_score, nm = _int_tag(fields, "AS"), _int_tag(fields, "NM")
    span = Alignment(fields[0], int(fields[1]), fields[2], int(fields[3]), fields[5], "*", {})
    return (
        1 if mapq == 255 else -mapq,
        as_score is None,
        -(as_score or 0),
        nm is None,
        nm or 0,
        -span.aligned_bases(),
        fields[0],
        int(fields[3]),
        int(fields[1]),
        fields[5],
        line,
    )


def select_molecule_representatives(
    lines: Iterable[str], ref_order: Optional[dict[str, int]] = None
) -> tuple[list[str], int]:
    """One deterministic representative alignment per corrected molecule.

    A molecule is a corrected (CB, UMI) on one reference, read from the
    ``<CB>_<UMI>_<n>`` names ``extract_viral_reads`` writes. Start, strand and
    CIGAR are deliberately NOT part of the key, so one molecule aligned at two
    starts survives once. Among a molecule's primary mapped alignments the best by
    ``_representative_rank`` wins. Records without a parseable CB/UMI have no
    molecule lineage: each (qname, reference) is kept read-level, never merged
    into a molecule, and counted in the returned *unresolved* total.

    Returns ``(records, unresolved)``. Records are sorted by (``ref_order`` rank,
    reference, POS, FLAG, qname, record), so the output is independent of input
    order and coordinate-sorted for BAM encoding. Header, malformed and
    ``EXCLUDE_FLAGS`` records are dropped. BAM deduplication and the direct
    read-start profile both call this function.
    """
    best: dict[tuple[str, ...], tuple[tuple[object, ...], list[str], str]] = {}
    for line in lines:
        fields = _primary_fields(line)
        if fields is None:
            continue
        cbumi = _cb_umi(fields[0])
        key = ("m", *cbumi, fields[2]) if cbumi else ("r", fields[0], fields[2])
        rank = _representative_rank(fields, line)
        if key not in best or rank < best[key][0]:
            best[key] = (rank, fields, line)
    order = ref_order or {}
    chosen = sorted(
        ((fields, line) for _rank, fields, line in best.values()),
        key=lambda fl: (
            order.get(fl[0][2], len(order)),
            fl[0][2],
            int(fl[0][3]),
            int(fl[0][1]),
            fl[0][0],
            fl[1],
        ),
    )
    unresolved = sum(1 for key in best if key[0] == "r")
    return [line for _fields, line in chosen], unresolved


def _parse_sam_read_starts(
    sam_text: str, *, dedup: str = "umi", strand_aware: bool = True, bin_size: int = 1
) -> list[dict[str, object]]:
    """Tally 5′ read-start positions per reference from ``samtools view`` text.

    Each primary alignment contributes its 5′ start: leftmost 0-based POS on the
    forward strand, or ``POS + reference_span − 1`` on the reverse strand when
    *strand_aware*. With ``dedup="umi"`` the alignments are first reduced to one
    representative per corrected (CB, UMI, reference) molecule by
    :func:`select_molecule_representatives` — the same choice the deduplicated BAM
    makes, independent of start, strand and CIGAR. ``dedup="none"`` keeps every
    primary alignment (use when dups were already removed: the deduplicated BAM or
    samtools markdup). ``n_reads`` counts the alignments tallied per reference.

    Split from ``read_start_distribution`` so it is unit-testable without samtools.
    Returns rows sorted by (reference, position): {reference, position,
    n_read_starts, n_reads}.
    """
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")
    lines: Iterable[str] = sam_text.splitlines()
    if dedup == "umi":
        lines = select_molecule_representatives(lines)[0]
    hist: dict[tuple[str, int], int] = {}
    n_reads: dict[str, int] = {}
    for line in lines:
        f = _primary_fields(line)
        if f is None:
            continue
        rname, pos0 = f[2], int(f[3]) - 1  # SAM POS is 1-based
        # 5' end: leftmost POS on +, rightmost consumed ref base on - (strand-aware).
        reverse = strand_aware and bool(int(f[1]) & 0x10)
        start = pos0 + max(_cigar_ref_span(f[5]) - 1, 0) if reverse else pos0
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
    """Header plus one representative alignment per corrected (CB, UMI, reference).

    Selection and tie-break are :func:`select_molecule_representatives`; records
    come out coordinate-sorted in ``@SQ`` order so the text can be encoded and
    indexed. Unmapped/secondary/supplementary/QC-fail/duplicate records and the
    non-representative candidates are not carried (they stay in the raw BAM).
    """
    lines = sam_text.splitlines()
    header = [ln for ln in lines if ln.startswith("@")]
    ref_order: dict[str, int] = {}
    for ln in header:
        if ln.startswith("@SQ"):
            tags = dict(t.split(":", 1) for t in ln.split("\t")[1:] if ":" in t)
            if "SN" in tags:
                ref_order.setdefault(tags["SN"], len(ref_order))
    kept, _unresolved = select_molecule_representatives(lines, ref_order)
    return "".join(f"{ln}\n" for ln in [*header, *kept])


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
    multi_class: bool = False,
    tie_delta: float = 0.0,
) -> list[dict[str, Any]]:
    """BLAST a deterministic read sample against combined host and target virus."""
    if not math.isfinite(tie_delta) or tie_delta < 0:
        raise ValueError("blast_tie_delta must be finite and >= 0")
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
    query_sequences = dict(iter_fasta(sample))
    if multi_class and len(query_sequences) != sampled:
        raise ValueError("sampled query IDs must be unique")
    n_subjects = 20
    if multi_class:
        with _open_maybe_gzip(competitive_fasta) as handle:
            n_subjects = sum(line.startswith(">") for line in handle)
    if sampled == 0 and multi_class:
        if sampling_manifest:
            record = json.loads(Path(sampling_manifest).read_text())
            record.update(query_ids=[], search_status="not_assessed", subject_count=n_subjects)
            Path(sampling_manifest).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        (work / "blast_raw_hits.tsv").write_text("")
        return []
    db = work / "competitive_db"
    failure = None
    try:
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
                str(max(1, n_subjects)),
                "-num_threads",
                str(threads),
                "-outfmt",
                "6 qseqid sseqid pident length qcovs evalue bitscore",
            ],
            capture=True,
        ).decode("utf-8", errors="replace")
    except RuntimeError as exc:
        if not multi_class:
            raise
        failure = str(exc)
        (work / "blast_failure.txt").write_text(failure + "\n")
        stdout = ""
    if multi_class:
        (work / "blast_raw_hits.tsv").write_text(stdout)
        rows = parse_competitor_blast_output(stdout, query_sequences, tie_delta=tie_delta)
        if failure is not None:
            for row in rows:
                row.update(diagnostic_class="failed", search_status="failed", reason=failure)
        if sampling_manifest:
            record = json.loads(Path(sampling_manifest).read_text())
            record.update(
                query_ids=list(query_sequences),
                subject_count=n_subjects,
                max_target_seqs=max(1, n_subjects),
                search_status="complete_subject_search" if failure is None else "failed",
                failure=failure,
                hsp_reduction="best_bitscore_then_identity_coverage_length_lexical",
                tie_delta=tie_delta,
            )
            Path(sampling_manifest).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    else:
        rows = [dict(row) for row in _parse_competitive_blast_output(stdout)]
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
        row["low_complexity"] = complexity.get(str(row["read"]), "not_assessed")
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
    """Return the subset of *names* that cannot be resolved (see ``tool_path``)."""
    return [n for n in names if tool_path(n) is None]
