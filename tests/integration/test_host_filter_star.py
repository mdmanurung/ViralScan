"""SW-07: real STAR host filter on the tiny fixture keeps removed-fragment lineage."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

from viralscan.evidence import have_tools
from viralscan.scripts import host_filter

pytestmark = pytest.mark.integration

FIXTURE = Path(__file__).parents[1] / "data" / "evidence_tiny"


def _star_genome(tmp_path: Path) -> Path:
    """STAR index of the fixture host (skips without STAR)."""
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
    return genome


@pytest.mark.parametrize("param_set", sorted(host_filter.STAR_FILTER_PARAM_SETS))
def test_star_host_filter_lineage_lists_removed_fragments(tmp_path: Path, param_set: str) -> None:
    """HOST-05b: the spiked viral read survives both STAR parameter sets (d15 loss 0 on the fixture)."""
    genome = _star_genome(tmp_path)
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
        param_set=param_set,
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
    assert got["d15_loss_fraction"] == 0 and got["not_in_lineage"] == 0


def test_host_loss_on_a_real_starsolo_run_and_its_primary_tree_guard(tmp_path: Path) -> None:
    """REVIEW-FIX-E end to end: host-loss reads a real --host-filter starsolo run, writes under
    posthoc/, and refuses an --output elsewhere inside the run tree, --force included."""
    from tests.integration.test_tiny_end_to_end import _run_tiny

    genome = _star_genome(tmp_path)
    out = _run_tiny(
        tmp_path, "--host-filter", "starsolo", "--host-index", str(genome), "--tchem", "off"
    )
    (sample,) = [d for d in out.iterdir() if (d / "config.yaml").is_file()]
    assert (sample / "host_filtered" / "fragment_lineage.tsv.gz").is_file()
    truth = tmp_path / "truth.tsv"
    truth.write_text("read_id\tlabel\nviral_read\tviral\nhost_read\thost\n")

    def host_loss(*extra: str) -> subprocess.CompletedProcess:
        cmd = ["host-loss", "--run-dir", str(sample), "--truth", str(truth), *extra]
        return subprocess.run(
            [sys.executable, "-m", "viralscan.menu", *cmd], capture_output=True, text=True
        )

    done = host_loss()
    assert done.returncode == 0, done.stderr
    payload = json.loads((sample / "posthoc" / "host_loss" / "host_loss.json").read_text())
    assert payload["d15_truth_fragments"] == 1 and payload["d15_removed_fragments"] == 0

    for target in (sample / "results" / "host_loss", out / "host_loss"):
        for force in ((), ("--force",)):
            refused = host_loss("--output", str(target), *force)
            assert refused.returncode != 0
            assert "is inside the ViralScan run" in refused.stderr
            assert not target.exists()
