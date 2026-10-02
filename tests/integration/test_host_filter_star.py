"""SW-07: real STAR host filter on the tiny fixture keeps removed-fragment lineage."""

from __future__ import annotations

import csv
import gzip
import subprocess
from pathlib import Path

import pytest

from viralscan.evidence import have_tools
from viralscan.scripts import host_filter

pytestmark = pytest.mark.integration

FIXTURE = Path(__file__).parents[1] / "data" / "evidence_tiny"


def test_star_host_filter_lineage_lists_removed_fragments(tmp_path: Path) -> None:
    missing = have_tools(["STAR"])
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
    genome = tmp_path / "genome"
    genome.mkdir()
    # STARsolo needs gene annotation in the index (geneInfo.tab).
    gtf = tmp_path / "host.gtf"
    gtf.write_text(
        'host_tx\tt\texon\t1\t80\t.\t+\t.\tgene_id "HOST_GENE"; transcript_id "host_tx";\n'
    )
    subprocess.run(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            str(genome),
            "--genomeFastaFiles",
            str(FIXTURE / "host.fasta"),
            "--sjdbGTFfile",
            str(gtf),
            "--genomeSAindexNbases",
            "3",
            "--outFileNamePrefix",
            str(tmp_path / "gen_"),
        ],
        check=True,
        capture_output=True,
    )
    out = tmp_path / "out"
    out.mkdir()
    r1, r2 = out / "R1.fastq.gz", out / "R2.fastq.gz"
    host_filter._starsolo_filter(
        str(FIXTURE / "R1.fastq"),
        str(FIXTURE / "R2.fastq"),
        str(genome),
        "10xv3",
        None,
        out,
        str(r1),
        str(r2),
        1,
    )
    with gzip.open(out / "fragment_lineage.tsv.gz", "rt") as fh:
        rows = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["host_read"]["filter_decision"] == "removed"
    assert rows["host_read"]["reason"] == "host_mapped"
    assert rows["viral_read"]["filter_decision"] == "retained"

    truth = tmp_path / "truth.tsv"
    truth.write_text("read_id\tlabel\nviral_read\tviral\nhost_read\thost\n")
    got = host_filter.lost_truth_counts(str(truth), str(out / "fragment_lineage.tsv.gz"))
    assert got["d15_removed_fragments"] == 0 and got["d15_truth_fragments"] == 1
