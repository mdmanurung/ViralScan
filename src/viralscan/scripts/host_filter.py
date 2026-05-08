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
"""

import gzip
import os
import shutil
import subprocess
from pathlib import Path

from viralscan.utils import load_config, setup_script_logging

log = setup_script_logging()

# ── Snakemake bindings ────────────────────────────────────────────────────────
configfile = snakemake.params.configfile  # noqa: F821  (snakemake magic global)
config = load_config(configfile)

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

# ── Technology parameters ─────────────────────────────────────────────────────
# (barcode_length, umi_length) for the most common 10x Chromium versions.
# Used by the kallisto mode to parse CB+UMI from R1.
_TECH_PARAMS: dict[str, tuple[int, int]] = {
    "10xv1": (14, 10),
    "10xv2": (16, 10),
    "10xv3": (16, 12),
    "10xv3_5p": (16, 12),
}


def _cb_umi_lengths() -> tuple[int, int]:
    """Return (cb_len, umi_len) for the configured technology."""
    return _TECH_PARAMS.get(technology, (16, 12))


# ── Helper: gzip-copy a plain-text file to a .gz destination ─────────────────
def _gzip_file(src: Path, dst: str) -> None:
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


# ── STARsolo mode ─────────────────────────────────────────────────────────────
def _starsolo_filter() -> None:
    """
    Run STARsolo with ``--outReadsUnmapped Fastx``.

    STAR writes unmapped mates to:
      Unmapped.out.mate1  — the *first* file given to --readFilesIn (= R2/cDNA in 10x)
      Unmapped.out.mate2  — the *second* file (= R1/barcode+UMI in 10x)

    In the 10x convention we pass R2 first and R1 second so that STAR's
    soloType CB_UMI_Simple knows which read carries the barcode.
    """
    star_tmp = out_dir / "star_tmp"
    star_tmp.mkdir(exist_ok=True)

    cb_len, umi_len = _cb_umi_lengths()

    read_files_cmd = "zcat" if r1.endswith(".gz") or r2.endswith(".gz") else "-"

    cmd = [
        "STAR",
        "--runThreadN", str(snakemake.threads),  # noqa: F821
        "--genomeDir", host_index,
        # 10x convention: cDNA read (R2) first, barcode+UMI read (R1) second
        "--readFilesIn", r2, r1,
        "--readFilesCommand", read_files_cmd,
        "--soloType", "CB_UMI_Simple",
        "--soloCBstart", "1",
        "--soloCBlen", str(cb_len),
        "--soloUMIstart", str(cb_len + 1),
        "--soloUMIlen", str(umi_len),
        "--outSAMtype", "None",
        "--outReadsUnmapped", "Fastx",
        "--outFileNamePrefix", str(star_tmp) + os.sep,
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
    _gzip_file(unmapped_bc, filtered_r1)   # barcode+UMI → R1
    _gzip_file(unmapped_cdna, filtered_r2)  # cDNA       → R2

    n_r1 = sum(1 for _ in gzip.open(filtered_r1, "rt")) // 4
    log.info("STARsolo host filter complete: %d unmapped read pairs retained.", n_r1)


# ── kallisto mode ─────────────────────────────────────────────────────────────
def _kallisto_filter() -> None:
    """
    Pseudo-align R1+R2 against a host cDNA kallisto index; keep only read
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

    # Step 1: pseudo-align
    log.info("Running kallisto bus against host index...")
    subprocess.run(
        [
            "kallisto", "bus",
            "-i", host_index,
            "-o", str(bus_dir),
            "-x", technology,
            r1, r2,
        ],
        check=True,
    )

    # Step 2: sort BUS file
    sorted_bus = str(bus_dir / "sorted.bus")
    subprocess.run(
        [
            "bustools", "sort",
            "-t", str(snakemake.threads),  # noqa: F821
            "-o", sorted_bus,
            str(bus_dir / "output.bus"),
        ],
        check=True,
    )

    # Step 3: convert to text
    bus_text = str(bus_dir / "mapped.txt")
    subprocess.run(
        ["bustools", "text", "-o", bus_text, sorted_bus],
        check=True,
    )

    # Step 4a: build set of host-mapped (CB, UMI) pairs
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

    # Step 4b: filter original FASTQs
    _filter_fastq_pairs(host_mapped)


def _filter_fastq_pairs(host_mapped: set[tuple[str, str]]) -> None:
    """
    Write read pairs whose (CB, UMI) is NOT in *host_mapped* to the
    filtered output FASTQs.

    The CB and UMI are extracted from the first ``cb_len + umi_len`` bases
    of every R1 sequence (10x Chromium layout: CB occupies bases 1..cb_len,
    UMI occupies bases cb_len+1..cb_len+umi_len).
    """
    cb_len, umi_len = _cb_umi_lengths()
    bc_end = cb_len + umi_len

    def _open_fq(path: str):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

    kept = 0
    total = 0

    with (
        _open_fq(r1) as fq1,
        _open_fq(r2) as fq2,
        gzip.open(filtered_r1, "wt") as out1,
        gzip.open(filtered_r2, "wt") as out2,
    ):
        while True:
            lines1 = [fq1.readline() for _ in range(4)]
            lines2 = [fq2.readline() for _ in range(4)]
            if not lines1[0]:  # EOF
                break
            total += 1
            seq1 = lines1[1].rstrip()
            cb = seq1[:cb_len]
            umi = seq1[cb_len:bc_end]
            if (cb, umi) not in host_mapped:
                out1.writelines(lines1)
                out2.writelines(lines2)
                kept += 1

    pct = 100.0 * kept / total if total else 0.0
    log.info(
        "kallisto host filter complete: kept %d / %d read pairs (%.1f%% passed host filter).",
        kept,
        total,
        pct,
    )


# ── Entry point ───────────────────────────────────────────────────────────────
def main() -> None:
    log.info("Host pre-subtraction: aligner=%s, host_index=%s", aligner, host_index)

    if aligner == "starsolo":
        _starsolo_filter()
    elif aligner == "kallisto":
        _kallisto_filter()
    else:
        raise ValueError(f"Unknown host_filter_aligner: {aligner!r}. Choose 'starsolo' or 'kallisto'.")

    # Signal completion to Snakemake
    Path(snakemake.output.done).touch()  # noqa: F821


main()
