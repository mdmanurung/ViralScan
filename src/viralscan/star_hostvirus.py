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
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TextIO

from viralscan.chemistry import cb_umi_geometry, onlist_path
from viralscan.star_run import plain_whitelist, run_star

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
    "--soloFeatures",
    "GeneFull",
)


#: viralscan strand -> ``--soloStrand``. STAR's own help: ``Forward`` = "read strand same as the
#: original RNA molecule", ``Reverse`` = "opposite to the original RNA molecule". The cDNA read
#: (R2, STAR mate 1) is what is aligned, so viralscan ``forward`` (R2 sense, 10x 3', kb
#: ``--strand forward``) is ``Forward`` and ``reverse`` (R2 antisense, 10x 5' kits; forward there
#: pseudoaligns only 6.5-8.9 % of reads, F-020) is ``Reverse``.
SOLO_STRAND: dict[str, str] = {
    "forward": "Forward",
    "reverse": "Reverse",
    "unstranded": "Unstranded",
}


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
    strand: str = "unstranded",
) -> list[str]:
    """STARsolo command (pure). cDNA (R2) first, barcode+UMI (R1) second.

    ``--soloBarcodeReadLength 0``: R1 can be longer than barcode+UMI (5' libraries ship
    150 bp). ``whitelist=None`` passes STAR's literal ``None`` (CB is the raw barcode, no
    correction); ``align_host_unmapped`` resolves the chemistry's on-list before calling
    this. ``strand`` is a viralscan strand, see ``SOLO_STRAND``.
    """
    if strand not in SOLO_STRAND:
        raise ValueError(f"strand must be one of {sorted(SOLO_STRAND)}, got {strand!r}")
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
        "--soloStrand", SOLO_STRAND[strand],
        *ALIGN_ARGS,
    ]  # fmt: skip


def resolve_whitelist(
    technology: str, whitelist: str | Path | None, work_dir: str | Path
) -> str | None:
    """The STAR whitelist: an explicit one wins, else kb's on-list (gunzipped), else None.

    None (STAR's literal ``None``) only for chemistries with no list (Drop-seq, explicit
    geometry strings), as ``chemistry.kb_whitelist_arg`` does (SW-21).
    """
    if whitelist:
        return plain_whitelist(str(whitelist), Path(work_dir))
    try:
        found = onlist_path(technology)
    except ValueError as exc:
        if type(exc) is not ValueError:  # ChemistryError (no ngs_tools) must not look list-less
            raise
        return None
    return plain_whitelist(found, Path(work_dir))


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
    strand: str = "unstranded",
    prune: bool = False,
    bam_sort_ram: int = 30_000_000_000,
) -> Path:
    """Align host-unmapped pairs; returns the indexed ``Aligned.sortedByCoord.out.bam``.

    With ``prune=True`` the per-virus summary is written to ``star_support_by_virus.tsv`` in
    *out_dir* first, then the full BAM and ``_STARtmp`` are replaced by a ``VIRUS|``-only BAM,
    whose path is returned (``prune_to_viral``).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = str(out) + "/"  # STAR concatenates; the trailing slash is the directory form
    wl = resolve_whitelist(technology, whitelist, out)
    cmd = align_cmd(star, index, r1, r2, technology, prefix, threads, wl, bam_sort_ram, strand)
    run_star(cmd, out, cb_len=cb_umi_geometry(technology)[0])
    bam = out / "Aligned.sortedByCoord.out.bam"
    subprocess.run([samtools, "index", str(bam)], check=True)
    if prune:
        write_summary(summarise_by_virus(bam, samtools), out / "star_support_by_virus.tsv")
        return prune_to_viral(bam, samtools)
    return bam


# ── per-virus summary ───────────────────────────────────────────────────────────────
#: One row per VIRUS| contig (= viral accession, the key ``detection`` uses for kb's targets).
STAR_SUMMARY_COLUMNS: tuple[str, ...] = (
    "accession",
    "star_reads_primary",
    "star_reads_any",
    "star_molecules",
    "star_host_competing",
    "star_other_virus",
)


def _sam_lines(samtools: str, *args: str) -> Iterator[list[str]]:
    proc = subprocess.Popen(
        [samtools, "view", *args], stdout=subprocess.PIPE, text=True, encoding="utf-8"
    )
    assert proc.stdout is not None
    with proc.stdout:
        for line in proc.stdout:
            yield line.rstrip("\n").split("\t")
    if proc.wait():
        raise subprocess.CalledProcessError(proc.returncode, proc.args)


def _viral_bed(bam: Path, samtools: str, bed: Path) -> int:
    """Write a BED of the VIRUS| contigs with mapped reads (idxstats); returns their count."""
    stats = subprocess.run(
        [samtools, "idxstats", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    rows = [ln.split("\t") for ln in stats.splitlines()]
    viral = [(r[0], r[1]) for r in rows if r[0].startswith("VIRUS|") and int(r[2]) > 0]
    bed.write_text("".join(f"{c}\t0\t{n}\n" for c, n in viral), encoding="utf-8")
    return len(viral)


def _tag(fields: list[str], key: str) -> str | None:
    for t in fields[11:]:
        if t.startswith(key):
            return t[len(key) :]
    return None


def summarise_by_virus(bam: str | Path, samtools: str = "samtools") -> list[dict[str, object]]:
    """Per-``VIRUS|`` contig support of a STAR host+virus BAM (rows keyed by STAR_SUMMARY_COLUMNS).

    Unit is the read (QNAME; the STAR BAM is single-end, only the cDNA mate is aligned). STAR
    writes one record per locus (NH records, same QNAME). For a read with a record on contig
    ``c`` (any record, secondary included, ``star_reads_any``):

    * ``star_reads_primary``: the read's primary record (flag 0x100 clear) is on ``c``;
    * ``star_molecules``: distinct Solo-assigned (CB, UB) among those reads. STARsolo sets UB
      only for reads it counts to a GeneFull feature (Unique mode, ``--soloStrand`` respected),
      so intergenic, wrong-strand and multimapping reads have ``-`` and are dropped: this is
      not read support (use ``star_reads_any``) and shifts with ``--strand`` and the GTF;
    * ``star_host_competing``: the read also has a ``HOST|`` record;
    * ``star_other_virus``: the read also has a record on another ``VIRUS|`` contig.

    Two passes: (1) viral contigs only (BED from idxstats, uses the index) collects the reads;
    (2) for those with NH > 1 only, ``samtools view -N`` over the whole BAM fetches the other
    loci. Memory is one small record per viral read plus the contig set of each viral
    multimapper; host-only reads are never held. Rows are sorted by accession.
    """
    bam = Path(bam)
    reads: dict[str, dict[str, object]] = {}
    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "viral.bed"
        if not _viral_bed(bam, samtools, bed):
            return []
        multi: set[str] = set()
        for f in _sam_lines(samtools, "-L", str(bed), str(bam)):
            r = reads.setdefault(f[0], {"contigs": set(), "primary": None, "cb": None, "ub": None})
            r["contigs"].add(f[2])  # type: ignore[attr-defined]
            if not int(f[1]) & 0x100:
                r["primary"] = f[2]
            r["cb"], r["ub"] = _tag(f, "CB:Z:"), _tag(f, "UB:Z:")
            if int(_tag(f, "NH:i:") or 1) > 1:
                multi.add(f[0])
        if multi:
            names = Path(tmp) / "names.txt"
            names.write_text("".join(f"{q}\n" for q in sorted(multi)), encoding="utf-8")
            for f in _sam_lines(samtools, "-N", str(names), str(bam)):
                if f[0] in multi and f[2] != "*":
                    r = reads[f[0]]
                    r["contigs"].add(f[2])  # type: ignore[attr-defined]
                    if not int(f[1]) & 0x100:
                        r["primary"] = f[2]
    acc: dict[str, dict[str, object]] = {}
    for r in reads.values():
        contigs: set[str] = r["contigs"]  # type: ignore[assignment]
        viral = {c for c in contigs if c.startswith("VIRUS|")}
        host = any(c.startswith("HOST|") for c in contigs)
        for c in viral:
            row = acc.setdefault(
                c,
                {k: 0 for k in STAR_SUMMARY_COLUMNS[1:]} | {"_mol": set()},
            )
            row["star_reads_any"] += 1  # type: ignore[operator]
            row["star_reads_primary"] += r["primary"] == c  # type: ignore[operator]
            row["star_host_competing"] += host  # type: ignore[operator]
            row["star_other_virus"] += len(viral) > 1  # type: ignore[operator]
            if r["cb"] not in (None, "-") and r["ub"] not in (None, "-"):
                row["_mol"].add((r["cb"], r["ub"]))  # type: ignore[attr-defined]
    out = []
    for c in sorted(acc):
        row = acc[c]
        row["star_molecules"] = len(row.pop("_mol"))  # type: ignore[arg-type]
        out.append(
            {"accession": c.removeprefix("VIRUS|")} | {k: row[k] for k in STAR_SUMMARY_COLUMNS[1:]}
        )
    return out


def write_summary(rows: Sequence[dict[str, object]], path: str | Path) -> Path:
    """Write ``summarise_by_virus`` rows as a TSV with the STAR_SUMMARY_COLUMNS header."""
    path = Path(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=STAR_SUMMARY_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return path


def prune_to_viral(bam: str | Path, samtools: str = "samtools") -> Path:
    """Keep a ``VIRUS|``-only BAM (+ index) next to *bam*; delete *bam*, its index and ``_STARtmp``.

    Call after ``summarise_by_virus``: the HOST| partners of viral multimappers go with the
    full BAM, so competition categories can no longer be rebuilt from the pruned file.
    """
    bam = Path(bam)
    viral = bam.with_name("Aligned.viral.bam")
    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "viral.bed"
        _viral_bed(bam, samtools, bed)
        subprocess.run(
            [samtools, "view", "-b", "-L", str(bed), "-o", str(viral), str(bam)], check=True
        )
    subprocess.run([samtools, "index", str(viral)], check=True)
    bam.unlink()
    bam.with_name(bam.name + ".bai").unlink(missing_ok=True)
    shutil.rmtree(bam.parent / "_STARtmp", ignore_errors=True)
    return viral


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
    a_.add_argument("--strand", choices=sorted(SOLO_STRAND), default="unstranded")
    a_.add_argument(
        "--prune", action="store_true", help="write the per-virus TSV, keep only a VIRUS| BAM"
    )
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
            args.r1, args.r2, args.index, args.out_dir, tech, args.whitelist, args.threads,
            strand=args.strand, prune=args.prune,
        ))  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
