"""Live kallisto/bustools exact read-number lineage smoke test."""

from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from argparse import Namespace
from pathlib import Path

import pytest

from viralscan.evidence import (
    extract_exact_reads_by_number,
    have_tools,
    parse_flagged_target_bus,
    replay_exact_target_bus,
)
from viralscan.runconfig import RunConfig
from viralscan.scripts.evidence_run import run_evidence


@pytest.mark.integration
def test_kallisto_read_number_replay_extracts_only_target(tmp_path: Path) -> None:
    missing = have_tools(["kallisto", "bustools"])
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

    fixture = Path(__file__).parents[1] / "data" / "evidence_tiny"
    index = tmp_path / "index.idx"
    subprocess.run(
        ["kallisto", "index", "-i", str(index), str(fixture / "transcripts.fasta")],
        check=True,
    )

    initial = tmp_path / "initial"
    subprocess.run(
        [
            "kallisto",
            "bus",
            "-i",
            str(index),
            "-o",
            str(initial),
            "-x",
            "10xv3",
            "-n",
            str(fixture / "R1.fastq"),
            str(fixture / "R2.fastq"),
        ],
        check=True,
    )
    flagged = replay_exact_target_bus(
        index=str(index),
        technology="10xv3",
        r1_path=str(fixture / "R1.fastq"),
        r2_path=str(fixture / "R2.fastq"),
        ec_file=str(initial / "matrix.ec"),
        transcripts_file=str(initial / "transcripts.txt"),
        target_transcripts=["viral_tx"],
        workdir=str(tmp_path / "replay"),
        threads=1,
    )
    with flagged.open() as handle:
        lineage = parse_flagged_target_bus(
            handle,
            {0: [0], 1: [1]},
            target_gene_indices={0},
            viral_gene_indices={0},
            method="host-conservative",
        )
    assert set(lineage) == {0}

    fasta = tmp_path / "viral_reads.fasta"
    lineage_path = tmp_path / "read_lineage.tsv.gz"
    stats = extract_exact_reads_by_number(
        str(fixture / "R1.fastq"),
        str(fixture / "R2.fastq"),
        lineage,
        str(fasta),
        str(lineage_path),
    )
    assert stats.total_reads == 2 and stats.viral_reads == 1
    assert "viral_read" in fasta.read_text() and "host_read" not in fasta.read_text()
    with gzip.open(lineage_path, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["read_id"] for row in rows] == ["viral_read"]

    # Build the minimal completed-run layout and exercise the documented CLI
    # service through exact replay, competitive alignment/BLAST, dedup, QC, and IGV.
    extra_missing = have_tools(["minimap2", "samtools", "blastn", "makeblastdb"])
    if extra_missing:
        pytest.skip(f"Evidence-chain binaries not on PATH: {', '.join(extra_missing)}")
    run_dir = tmp_path / "run"
    kb_dir = run_dir / "kb-python"
    counts = kb_dir / "counts_unfiltered"
    log_dir = run_dir / "log"
    results = run_dir / "results"
    for directory in (counts, log_dir, results):
        directory.mkdir(parents=True)
    for name in ("matrix.ec", "transcripts.txt"):
        shutil.copy2(initial / name, kb_dir / name)
    subprocess.run(
        [
            "bustools",
            "text",
            "-o",
            str(kb_dir / "output.resolved.sorted.bus.txt"),
            str(initial / "output.bus"),
        ],
        check=True,
    )
    (counts / "cells_x_genes.genes.txt").write_text("VIRUS_TARGET\nHOST_GENE\n")
    (log_dir / "analysis.txt").write_text("VIRUS_TARGET\n")
    (results / "viral_summary.tsv").write_text(
        "virus_name\tviral_molecules_total_est\nVIRUS_TARGET\t1\n"
    )
    RunConfig(
        output=f"{run_dir}/",
        index=str(index),
        transcripts=str(fixture / "t2g.txt"),
        sample1=str(fixture / "R1.fastq"),
        sample2=str(fixture / "R2.fastq"),
        technology="10xv3",
        multimap_method="host-conservative",
    ).to_yaml(run_dir / "config.yaml")
    evidence_out = tmp_path / "evidence"
    run_evidence(
        Namespace(
            run_dir=str(run_dir),
            output=str(evidence_out),
            virus="VIRUS_TARGET",
            viral_fasta=str(fixture / "viral.fasta"),
            host_fasta=str(fixture / "host.fasta"),
            blast=True,
            read_start_profile=True,
            cell_tags=True,
            dedup="umi",
            bin_size=1,
            sampling_seed=11,
            cores=1,
            verbose=False,
            quiet=True,
        )
    )
    for required in (
        "read_lineage.tsv.gz",
        "evidence_manifest.json",
        "competitive_reads.raw.bam",
        "competitive_reads.raw.bam.bai",
        "competitive_reads.umi_dedup.bam",
        "coverage.raw.tsv",
        "coverage.deduplicated.tsv",
        "alignment_qc.tsv",
        "per_cell_alignment_qc.tsv",
        "coverage.raw_vs_deduplicated.png",
        "interpretation_flags.tsv",
        "blast_identity.tsv",
        "blast_sampling.json",
        "read_start_profile.tsv",
        "viralscan_evidence.igv.xml",
    ):
        assert (evidence_out / required).is_file(), required
