"""M-01: the strand a sample was counted with is the strand its config.yaml records and replays use."""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan.menu import _write_run_config
from viralscan.runconfig import RunConfig


def _cfg(tmp_path: Path, strand: str | None) -> dict[str, object]:
    return {
        "output": str(tmp_path / "s") + "/",
        "index": "/ref/index.idx",
        "transcripts": "/ref/t2g.txt",
        "sample1": "R1.fastq.gz",
        "sample2": "R2.fastq.gz",
        "gtf": None,
        "fasta": None,
        "visual": True,
        "f1": None,
        "reference": False,
        "umap": False,
        "technology": "10xv3",
        "whitelist": None,
        "multimapping": True,
        "strand": strand,
    }


@pytest.mark.parametrize("strand", [None, "forward", "reverse", "unstranded"])
def test_config_yaml_round_trips_the_strand_the_count_used(tmp_path: Path, strand) -> None:
    """``evidence`` re-reads ``config.yaml``; the Snakefile's ``--strand`` comes from the same file."""
    path = _write_run_config(RunConfig.from_snakemake_config(_cfg(tmp_path, strand)))
    assert RunConfig.from_yaml(path).strand == strand
