"""ANDET-09: real STARsolo run of the anellovirus branch on the tiny fixture.

Checks the two things the pure tests cannot: that ALIGN_ARGS is a command
STAR 2.7.11b accepts, and that STARsolo writes CB/UB tags with a
one-gene-per-contig GTF, which is what molecule counting rests on.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from viralscan import anello_align as aa
from viralscan.evidence import have_tools
from viralscan.scripts import anello_align as rule
from viralscan.scripts.host_filter import starsolo_barcode_args

pytestmark = pytest.mark.integration

FIXTURE = Path(__file__).parents[1] / "data" / "evidence_tiny"


def test_starsolo_branch_tags_viral_read_and_drops_host(tmp_path: Path) -> None:
    missing = have_tools(["STAR", "samtools"])
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
    index = tmp_path / "anello_star"
    fasta, gtf, n = aa.write_star_reference(FIXTURE / "viral.fasta", {"viral_tx"}, index)
    assert n == 1
    cmd = aa.genome_generate_cmd("STAR", index, fasta, gtf, 1)
    # An 80-nt genome needs a tiny suffix-array prefix; the shipped value is for ~6 Mb.
    cmd[cmd.index("--genomeSAindexNbases") + 1] = "3"
    subprocess.run(cmd, check=True, capture_output=True, cwd=index)
    assert aa.resolve_index(str(tmp_path / "panel.idx")) == str(index.resolve())

    # The shared fixture's UMI is a homopolymer (CCCCCCCCCCCC), which STARsolo
    # filters, leaving CB:Z:-. Use a normal UMI here; the barcode is unchanged.
    r1 = tmp_path / "R1.fastq"
    r1.write_text(
        (FIXTURE / "R1.fastq")
        .read_text()
        .replace("AAACCCAAGAAACACTCCCCCCCCCCCC", "AAACCCAAGAAACACTACGTACGTACGT")
    )
    onlist = tmp_path / "onlist.txt"
    onlist.write_text("AAACCCAAGAAACACT\nTTTTTTTTTTTTTTTT\n")
    for whitelist in (str(onlist), None):
        work = tmp_path / f"out_{whitelist is None}"
        work.mkdir()
        align = aa.align_cmd(
            "STAR",
            str(index),
            str(FIXTURE / "R2.fastq"),
            str(r1),
            starsolo_barcode_args("10xv3", whitelist),
            f"{work}/",
            1,
        )
        subprocess.run(align, check=True, capture_output=True)
        bam = work / "Aligned.sortedByCoord.out.bam"

        records = list(rule._alignments(bam))
        assert {a.qname.split("/")[0] for a in records} == {"viral_read"}
        a = records[0]
        assert a.rname == "viral_tx" and a.nh == 1
        assert a.tags.get("CB") == "AAACCCAAGAAACACT"
        assert a.tags.get("UB") == "ACGTACGTACGT"

        rows = aa.accession_metrics(rule._alignments(bam), {"viral_tx": 80})
        assert rows[0]["reads"] == 1 and rows[0]["median_identity"] == 1.0
        mol = aa.virus_molecules(rule._alignments(bam), {"viral_tx": "Testvirus"})
        assert mol == {"Testvirus": {"molecules": 1, "cells": 1}}
