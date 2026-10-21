"""STARsolo alignment of host-unmapped reads against one host + viral-panel index.

The host filter drops reads that STAR aligns to the host, so a coverage map built
from the kallisto/competitive BAMs never shows host competition for the reads that
survive. This module aligns the *host-unmapped* read pairs to a single STAR index
that holds the host genome (contigs ``HOST|<name>``) and the full viral panel
(``VIRUS|<accession>``) together, so the coordinate-sorted BAM carries CB/UB and the
prefixes ``plot_coverage_maps`` / ``is_viral_reference`` expect (PLAN STAR-HOSTVIRUS-01).

Python standard library only at import time; ``evidence`` is imported lazily.
Command line: ``python -m viralscan.star_hostvirus {reference,index,align} ...``.
"""

import argparse
import json
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from viralscan.chemistry import cb_umi_geometry

#: ``--genomeSAindexNbases 14`` is STAR's value for a 3 Gb genome (the host filter
#: index uses it too). ``--genomeChrBinNbits 15``: with ~2,500 contigs STAR's
#: default of 18 pads each ~5 kb viral contig to 262 kb (~0.6 Gb wasted); 2^15 bins
#: cost 77 Mb. The suffix array is not sparse (host filter used ``SAsparseD 3``)
#: so seed search on divergent viral reads is as sensitive as STAR allows.
GENOME_GENERATE_ARGS: tuple[str, ...] = (
    "--genomeChrBinNbits",
    "15",
    "--sjdbOverhang",
    "100",
)

#: Alignment settings. Host-unmapped reads are the input, so the aim is to keep every
#: placement a read has: ``--outFilterMultimapNmax 100`` (host filter used 20) with
#: ``--outSAMmultNmax 100`` writes all of them (HI/NH tell them apart) and
#: ``--winAnchorMultimapNmax 200`` lets seeds find that many loci; reads over 100 loci
#: are "mapped to too many loci" and absent from the BAM. Mismatch/score filters are the
#: relaxed anello_align values (8 % of read length, 80 % matched) so a divergent strain
#: still aligns, and, being applied against the host too, lets host compete.
#: No ``--alignIntronMax``: it must stay open or spliced host reads would be lost to the
#: viral contigs only because host introns were forbidden.
ALIGN_ARGS: tuple[str, ...] = (
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
    "CR",
    "CY",
    "UR",
    "UY",
    "CB",
    "UB",
    "--outSAMtype",
    "BAM",
    "SortedByCoordinate",
    "--soloStrand",
    "Unstranded",
    "--soloFeatures",
    "GeneFull",
)


def _prefixed_gtf(source: str | Path, prefix: str, target: TextIO) -> None:
    with open(source, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                target.write(line)
            else:
                target.write(f"{prefix}|{line}")


def write_reference(
    host_fasta: str | Path,
    viral_fasta: str | Path,
    host_gtf: str | Path,
    viral_gtf: str | Path,
    out_dir: str | Path,
) -> tuple[Path, Path]:
    """Write ``hostvirus.fa`` and ``hostvirus.gtf`` with ``HOST|`` / ``VIRUS|`` seqids."""
    from viralscan.evidence import write_competitive_fasta

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fasta = out / "hostvirus.fa"
    gtf = out / "hostvirus.gtf"
    write_competitive_fasta(str(host_fasta), str(viral_fasta), str(fasta))
    with open(gtf, "w", encoding="utf-8") as target:
        _prefixed_gtf(host_gtf, "HOST", target)
        _prefixed_gtf(viral_gtf, "VIRUS", target)
    return fasta, gtf


def genome_generate_cmd(
    star: str,
    index_dir: str | Path,
    fasta: str | Path,
    gtf: str | Path,
    threads: int,
    ram_bytes: int,
    sa_index_nbases: int = 14,
) -> list[str]:
    """STAR genomeGenerate command (pure)."""
    return [
        star, "--runMode", "genomeGenerate", "--runThreadN", str(threads),
        "--genomeDir", str(index_dir), "--genomeFastaFiles", str(fasta),
        "--sjdbGTFfile", str(gtf), "--limitGenomeGenerateRAM", str(ram_bytes),
        "--genomeSAindexNbases", str(sa_index_nbases),
        *GENOME_GENERATE_ARGS,
    ]  # fmt: skip


def build_star_hostvirus_index(
    fasta: str | Path,
    gtf: str | Path,
    index_dir: str | Path,
    threads: int = 8,
    ram_bytes: int = 85_000_000_000,
    star: str = "STAR",
    sa_index_nbases: int = 14,
) -> Path:
    """Run STAR genomeGenerate for the combined reference; returns *index_dir*."""
    index = Path(index_dir)
    index.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        genome_generate_cmd(star, index, fasta, gtf, threads, ram_bytes, sa_index_nbases),
        check=True,
    )
    return index


def align_cmd(
    star: str,
    index: str | Path,
    r1: str | Path,
    r2: str | Path,
    technology: str,
    prefix: str,
    threads: int,
    whitelist: str | Path | None = None,
    bam_sort_ram: int = 30_000_000_000,
) -> list[str]:
    """STARsolo command (pure). cDNA (R2) first, barcode+UMI (R1) second.

    ``--soloBarcodeReadLength 0``: R1 can be longer than barcode+UMI (5' libraries ship
    150 bp). ``whitelist=None`` passes STAR's literal ``None`` like the host filter did
    for these libraries, so CB is the raw barcode (no correction).
    """
    cb_len, umi_len = cb_umi_geometry(technology)
    gz = str(r1).endswith(".gz") or str(r2).endswith(".gz")
    return [
        star, "--runThreadN", str(threads), "--genomeDir", str(index),
        "--readFilesIn", str(r2), str(r1),
        "--readFilesCommand", "zcat" if gz else "-",
        "--outFileNamePrefix", prefix,
        "--soloType", "CB_UMI_Simple",
        "--soloCBstart", "1", "--soloCBlen", str(cb_len),
        "--soloUMIstart", str(cb_len + 1), "--soloUMIlen", str(umi_len),
        "--soloBarcodeReadLength", "0",
        "--soloCBwhitelist", str(whitelist) if whitelist else "None",
        "--limitBAMsortRAM", str(bam_sort_ram),
        *ALIGN_ARGS,
    ]  # fmt: skip


def align_host_unmapped(
    r1: str | Path,
    r2: str | Path,
    index: str | Path,
    out_dir: str | Path,
    technology: str,
    whitelist: str | Path | None = None,
    threads: int = 8,
    star: str = "STAR",
    samtools: str = "samtools",
) -> Path:
    """Align host-unmapped pairs; returns the indexed ``Aligned.sortedByCoord.out.bam``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = str(out) + "/"  # STAR concatenates; the trailing slash is the directory form
    cmd = align_cmd(star, index, r1, r2, technology, prefix, threads, whitelist)
    subprocess.run(cmd, check=True)
    bam = out / "Aligned.sortedByCoord.out.bam"
    subprocess.run([samtools, "index", str(bam)], check=True)
    return bam


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python -m viralscan.star_hostvirus")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reference")
    for a in ("host-fasta", "viral-fasta", "host-gtf", "viral-gtf", "out-dir"):
        r.add_argument(f"--{a}", required=True)
    i = sub.add_parser("index")
    for a in ("fasta", "gtf", "index-dir"):
        i.add_argument(f"--{a}", required=True)
    i.add_argument("--threads", type=int, default=8)
    i.add_argument("--ram-bytes", type=int, default=85_000_000_000)
    a_ = sub.add_parser("align")
    for a in ("r1", "r2", "index", "out-dir"):
        a_.add_argument(f"--{a}", required=True)
    a_.add_argument(
        "--chem-json", help="chem/<lib>.json; its chemistry.chemistry is the technology"
    )
    a_.add_argument("--technology")
    a_.add_argument("--whitelist")
    a_.add_argument("--threads", type=int, default=8)
    args = p.parse_args(argv)
    if args.cmd == "reference":
        fa, gtf = write_reference(
            args.host_fasta, args.viral_fasta, args.host_gtf, args.viral_gtf, args.out_dir
        )
        print(fa, gtf)
    elif args.cmd == "index":
        build_star_hostvirus_index(
            args.fasta, args.gtf, args.index_dir, args.threads, args.ram_bytes
        )
    else:
        tech = args.technology
        if args.chem_json:
            tech = json.loads(Path(args.chem_json).read_text())["chemistry"]["chemistry"]
        if not tech:
            p.error("give --technology or --chem-json")
        if shutil.which("STAR") is None:
            sys.exit("STAR not on PATH")
        print(align_host_unmapped(
            args.r1, args.r2, args.index, args.out_dir, tech, args.whitelist, args.threads
        ))  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
