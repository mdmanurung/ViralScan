"""Snakemake rule ``anello_align``: STARsolo on host-unmapped reads (PLAN ANDET-09).

Runs only when the host filter ran and an ``anello_star/`` index sits next to
the kb index (the Snakefile gates it). Writes:

* ``results/anello_alignment_by_accession.tsv`` — per-accession evidence;
* ``results/anello_alignment_by_virus.tsv`` — virus-level columns that
  ``detection`` merges into ``viral_summary.tsv``.

The BAM and STARsolo matrices stay in ``anello_align/`` for read review.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from pathlib import Path

from viralscan.anello_align import (
    ACCESSION_TSV,
    SUMMARY_COLUMNS,
    accession_metrics,
    align_cmd,
    iter_fasta,
    parse_sam_line,
    virus_molecules,
    virus_summary,
    write_accession_tsv,
)
from viralscan.chemistry import onlist_path
from viralscan.runconfig import RunConfig
from viralscan.scripts.host_filter import _plain_whitelist, starsolo_barcode_args

log = logging.getLogger(__name__)

VIRUS_TSV = "results/anello_alignment_by_virus.tsv"


def _acc_to_virus(identity_tsv: Path) -> dict[str, str]:
    """genome_accession -> virus_name for Anelloviridae rows of virus_identity.tsv."""
    out = {}
    with open(identity_tsv, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("family") == "Anelloviridae" and row.get("genome_accession"):
                out[row["genome_accession"]] = row["virus_name"]
    return out


def _alignments(bam: Path):
    """Stream mapped records of *bam* through ``samtools view``."""
    with subprocess.Popen(
        ["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            a = parse_sam_line(line)
            if a is not None:
                yield a
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, ["samtools", "view", str(bam)])


def _onlist(technology: str):
    """kb's bundled on-list for *technology*, or None (geometry strings, no list)."""
    try:
        return onlist_path(technology)
    except ValueError:
        return None


def _star_version() -> str:
    return subprocess.run(["STAR", "--version"], capture_output=True, text=True, check=True).stdout.strip()


def main(config: RunConfig, identity_tsv: str, n_threads: int, done_path: str) -> None:
    out = Path(config.output)
    work = out / "anello_align"
    work.mkdir(parents=True, exist_ok=True)
    index = config.anello_index
    assert index, "anello_align runs only when an anello_star index was resolved"

    # Correct barcodes against the list kb uses, so molecules line up with the
    # kallisto matrix. With no list (SW-21 chemistries) STARsolo keeps raw CBs.
    whitelist = _plain_whitelist(config.whitelist or _onlist(config.technology), work)
    cmd = align_cmd(
        "STAR",
        index,
        config.kb_r2,  # cDNA
        config.kb_r1,  # barcode + UMI
        starsolo_barcode_args(config.technology, whitelist),
        f"{work}/",
        n_threads,
    )
    log.info("anello_align: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    bam = work / "Aligned.sortedByCoord.out.bam"
    if not bam.is_file():
        raise FileNotFoundError(f"STAR wrote no BAM: {bam}")
    subprocess.run(["samtools", "index", str(bam)], check=True)

    lengths = {name: len(seq) for name, seq in iter_fasta(Path(index) / "anello.fa")}
    acc_to_virus = _acc_to_virus(Path(identity_tsv))
    acc_rows = accession_metrics(_alignments(bam), lengths)
    molecules = virus_molecules(_alignments(bam), acc_to_virus)
    meta = {
        "star_version": _star_version(),
        "anello_index": str(index),
        "n_contigs": len(lengths),
        "barcode_list": whitelist or "none (uncorrected)",
    }
    write_accession_tsv(acc_rows, out / ACCESSION_TSV, meta)

    summary = virus_summary(acc_rows, molecules, acc_to_virus)
    with open(out / VIRUS_TSV, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["virus_name", *SUMMARY_COLUMNS[2:]])
        for virus, cols in summary.items():
            writer.writerow([virus, *(cols[c] for c in SUMMARY_COLUMNS[2:])])
    log.info(
        "anello_align: %d accessions, %d viruses with alignment evidence",
        len(acc_rows),
        len(summary),
    )
    # Two-pass mode leaves its pass-1 genome behind; it is large and unused.
    shutil.rmtree(work / "_STARgenome", ignore_errors=True)
    Path(done_path).touch()


if "snakemake" in globals():
    main(
        config=RunConfig.from_yaml(snakemake.params.configfile),  # noqa: F821
        identity_tsv=str(snakemake.input.virus_identity),  # noqa: F821
        n_threads=snakemake.threads,  # noqa: F821
        done_path=str(snakemake.output.done),  # noqa: F821
    )
