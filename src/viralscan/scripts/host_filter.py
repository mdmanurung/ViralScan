"""
Optional host-subtraction pre-step for ViralScan.

Runs before ``kb_count`` when the user supplies ``--host-filter`` and
``--host-index``.  Two modes are supported:

starsolo
    STAR genome alignment (STARsolo barcode-aware mode).  Reads that do NOT
    align to the host genome are written as FASTQ by STAR via
    ``--outReadsUnmapped Fastx`` and collected as the filtered output.
    Requires a STAR genome directory built with ``STAR --runMode genomeGenerate``.

kallisto
    Pseudo-alignment of R2 (cDNA) against a host cDNA kallisto index.
    ``bustools`` converts the BUS file to text; the resulting (barcode, UMI)
    pairs that mapped to the host are used to filter the original FASTQ files
    in a single Python pass — any read pair whose (CB, UMI) was NOT seen in
    the host BUS is kept.
    Requires a kallisto index file (``.idx``) built from the host cDNA FASTA.

Output
    {output}host_filtered/R1.fastq.gz   — barcode+UMI read (R1 in 10x convention)
    {output}host_filtered/R2.fastq.gz   — cDNA read (R2 in 10x convention)

These paths are pre-registered in config["kb_r1"] / config["kb_r2"] by
``createconfig.py``, so the downstream ``kb_count`` rule consumes them
transparently with no further changes.

CB/UMI geometry is resolved via ``viralscan.evidence.cb_umi_geometry``, which
handles 10x v1/v2/v3, Drop-seq, and explicit kallisto ``bc:umi:seq`` triplets,
and raises ``ValueError`` for unknown chemistries instead of silently
mis-slicing barcodes (fixes PLAN S1/S6).
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from viralscan.evidence import _open_maybe_gzip, cb_umi_geometry
from viralscan.utils import load_config, setup_script_logging

log = setup_script_logging()


# ── Helper: gzip-copy a plain-text file to a .gz destination ─────────────────
def _gzip_file(src: Path, dst: str) -> None:
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


# ── FASTQ pair filter (pure — testable without Snakemake) ────────────────────
def filter_fastq_pairs(
    r1_path: str,
    r2_path: str,
    out_r1: str,
    out_r2: str,
    cb_len: int,
    umi_len: int,
    host_mapped: set[tuple[str, str]],
) -> tuple[int, int]:
    """Write read pairs whose (CB, UMI) is NOT in *host_mapped* to *out_r1*/*out_r2*.

    The CB occupies bases ``0:cb_len`` and the UMI ``cb_len:cb_len+umi_len``
    of every R1 read (10x/Drop-seq layout).  Returns ``(kept, total)`` counts.
    """
    bc_end = cb_len + umi_len
    kept = 0
    total = 0

    with (
        _open_maybe_gzip(r1_path) as fq1,
        _open_maybe_gzip(r2_path) as fq2,
        gzip.open(out_r1, "wt") as out1,
        gzip.open(out_r2, "wt") as out2,
    ):
        while True:
            lines1 = [fq1.readline() for _ in range(4)]
            lines2 = [fq2.readline() for _ in range(4)]
            if not lines1[0] or not lines2[0]:  # EOF on either file
                break
            if not lines1[3] or not lines2[3]:
                # Mid-record truncation: partial record at end of file.
                raise ValueError(
                    f"Truncated FASTQ: {r1_path!r} or {r2_path!r} ends mid-record"
                )
            total += 1
            seq1 = lines1[1].rstrip()
            cb = seq1[:cb_len]
            umi = seq1[cb_len:bc_end]
            if (cb, umi) not in host_mapped:
                out1.writelines(lines1)
                out2.writelines(lines2)
                kept += 1

    return kept, total


# ── STARsolo mode ─────────────────────────────────────────────────────────────
def _starsolo_filter(
    r1: str,
    r2: str,
    host_index: str,
    technology: str,
    whitelist: Optional[str],
    out_dir: Path,
    filtered_r1: str,
    filtered_r2: str,
    n_threads: int,
) -> None:
    """Run STARsolo with ``--outReadsUnmapped Fastx``.

    STAR writes unmapped mates to:
      Unmapped.out.mate1  — the *first* file given to --readFilesIn (= R2/cDNA in 10x)
      Unmapped.out.mate2  — the *second* file (= R1/barcode+UMI in 10x)

    In the 10x convention we pass R2 first and R1 second so that STAR's
    soloType CB_UMI_Simple knows which read carries the barcode.
    """
    star_tmp = out_dir / "star_tmp"
    star_tmp.mkdir(exist_ok=True)

    cb_len, umi_len = cb_umi_geometry(technology)

    read_files_cmd = "zcat" if r1.endswith(".gz") or r2.endswith(".gz") else "-"

    cmd = [
        "STAR",
        "--runThreadN",
        str(n_threads),
        "--genomeDir",
        host_index,
        # 10x convention: cDNA read (R2) first, barcode+UMI read (R1) second
        "--readFilesIn",
        r2,
        r1,
        "--readFilesCommand",
        read_files_cmd,
        "--soloType",
        "CB_UMI_Simple",
        "--soloCBstart",
        "1",
        "--soloCBlen",
        str(cb_len),
        "--soloUMIstart",
        str(cb_len + 1),
        "--soloUMIlen",
        str(umi_len),
        "--outSAMtype",
        "None",
        "--outReadsUnmapped",
        "Fastx",
        "--outFileNamePrefix",
        str(star_tmp) + os.sep,
    ]
    if whitelist:
        cmd += ["--soloCBwhitelist", whitelist]
    else:
        # Without a whitelist STARsolo accepts any barcode; pass "None" (STAR literal)
        cmd += ["--soloCBwhitelist", "None"]

    log.info("Running STARsolo host filter...")
    subprocess.run(cmd, check=True)

    # Mate1 = cDNA (→ R2), Mate2 = barcode+UMI (→ R1)
    unmapped_cdna = star_tmp / "Unmapped.out.mate1"
    unmapped_bc = star_tmp / "Unmapped.out.mate2"

    log.info("Compressing unmapped reads → %s", out_dir)
    _gzip_file(unmapped_bc, filtered_r1)  # barcode+UMI → R1
    _gzip_file(unmapped_cdna, filtered_r2)  # cDNA       → R2

    n_r1 = sum(1 for _ in gzip.open(filtered_r1, "rt")) // 4
    log.info("STARsolo host filter complete: %d unmapped read pairs retained.", n_r1)


# ── kallisto mode ─────────────────────────────────────────────────────────────
def _kallisto_filter(
    r1: str,
    r2: str,
    host_index: str,
    technology: str,
    out_dir: Path,
    filtered_r1: str,
    filtered_r2: str,
    n_threads: int,
) -> None:
    """Pseudo-align R1+R2 against a host cDNA kallisto index; keep only read
    pairs whose (barcode, UMI) was NOT seen in the host BUS file.

    Steps
    -----
    1. ``kallisto bus``  — pseudo-align; produces output.bus (unmapped reads are
       simply absent from the BUS file).
    2. ``bustools sort`` — sort for downstream text conversion.
    3. ``bustools text`` — convert sorted BUS to tab-delimited text
       (columns: barcode  umi  EC  count).
    4. Python pass       — build a set of host-mapped (CB, UMI) pairs, then
       scan the original FASTQs and keep pairs not in that set.
    """
    bus_dir = out_dir / "kb_host"
    bus_dir.mkdir(exist_ok=True)

    log.info("Running kallisto bus against host index...")
    subprocess.run(
        ["kallisto", "bus", "-i", host_index, "-o", str(bus_dir), "-x", technology, r1, r2],
        check=True,
    )

    sorted_bus = str(bus_dir / "sorted.bus")
    subprocess.run(
        ["bustools", "sort", "-t", str(n_threads), "-o", sorted_bus, str(bus_dir / "output.bus")],
        check=True,
    )

    bus_text = str(bus_dir / "mapped.txt")
    subprocess.run(["bustools", "text", "-o", bus_text, sorted_bus], check=True)

    host_mapped: set[tuple[str, str]] = set()
    with open(bus_text) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                host_mapped.add((parts[0], parts[1]))

    log.info(
        "kallisto host BUS: %d unique (barcode, UMI) pairs mapped to host.",
        len(host_mapped),
    )

    cb_len, umi_len = cb_umi_geometry(technology)
    kept, total = filter_fastq_pairs(r1, r2, filtered_r1, filtered_r2, cb_len, umi_len, host_mapped)
    pct = 100.0 * kept / total if total else 0.0
    log.info(
        "kallisto host filter complete: kept %d / %d read pairs (%.1f%% passed host filter).",
        kept,
        total,
        pct,
    )


# ── Entry point ───────────────────────────────────────────────────────────────
def main(config: dict, n_threads: int, done_path: str) -> None:
    output = config["output"]
    aligner = config.get("host_filter_aligner") or "starsolo"
    host_index = config["host_index"]
    r1 = config["sample1"]
    r2 = config["sample2"]
    technology = config.get("technology", "10xv3")
    whitelist = config.get("whitelist") or None

    out_dir = Path(output) / "host_filtered"
    out_dir.mkdir(parents=True, exist_ok=True)
    filtered_r1 = str(out_dir / "R1.fastq.gz")
    filtered_r2 = str(out_dir / "R2.fastq.gz")

    log.info("Host pre-subtraction: aligner=%s, host_index=%s", aligner, host_index)

    if aligner == "starsolo":
        _starsolo_filter(r1, r2, host_index, technology, whitelist, out_dir, filtered_r1, filtered_r2, n_threads)
    elif aligner == "kallisto":
        _kallisto_filter(r1, r2, host_index, technology, out_dir, filtered_r1, filtered_r2, n_threads)
    else:
        raise ValueError(
            f"Unknown host_filter_aligner: {aligner!r}. Choose 'starsolo' or 'kallisto'."
        )

    Path(done_path).touch()


# ── Snakemake wiring (only runs under snakemake) ─────────────────────────────
if "snakemake" in globals():
    _config = load_config(snakemake.params.configfile)  # noqa: F821
    main(
        config=_config,
        n_threads=snakemake.threads,  # noqa: F821
        done_path=str(snakemake.output.done),  # noqa: F821
    )
