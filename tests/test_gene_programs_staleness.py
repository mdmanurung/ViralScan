"""AUD-13 A-15: compare regenerated marker partitions, independent of group labels."""

import csv
import os
from collections import defaultdict
from pathlib import Path

import pytest

from extras import build_gene_programs as generator
from viralscan.scripts.ncbi_fetch import DEFAULT_CACHE_DIR


def partitions(rows):
    groups = defaultdict(set)
    for row in rows:
        groups[(row["virus"], row["overlap_group"])].add(row["gene_id_bundled"])
    return {(virus, frozenset(members)) for (virus, _), members in groups.items()}


def test_partition_comparison_accepts_renumbering_but_detects_splits():
    rows = [{"virus": "V", "overlap_group": "g1", "gene_id_bundled": gene} for gene in ("A", "B")]
    relabelled = [{**row, "overlap_group": "g9"} for row in rows]
    split = [relabelled[0], {**relabelled[1], "overlap_group": "g10"}]
    assert partitions(rows) == partitions(relabelled)
    assert partitions(rows) != partitions(split)


@pytest.mark.research
def test_shipped_partitions_and_independence_match_full_regeneration(tmp_path):
    cache = Path(os.environ.get("VIRALSCAN_TEST_NCBI_CACHE", DEFAULT_CACHE_DIR))
    if not cache.is_dir() or not any(cache.glob("*/*.fasta")):
        pytest.skip("Full cached reference FASTA/GenBank inputs are required for regeneration")
    missing = [
        name
        for name in generator.VIRUS_GTF.values()
        if not (Path(generator.DATA_DIR) / name).is_file()
    ]
    if missing:
        pytest.skip(f"Local reference GTF inputs are missing: {missing}")
    generated = tmp_path / "regenerated.tsv"
    assert generator.main(["--out", str(generated), "--cache-dir", str(cache)]) == 0
    with generated.open() as handle:
        actual = list(csv.DictReader(handle, delimiter="\t"))
    with open(generator.DEFAULT_OUT) as handle:
        expected = list(csv.DictReader(handle, delimiter="\t"))
    assert partitions(actual) == partitions(expected), "Shipped marker partition is stale"

    def flags(rows):
        return {(r["virus"], r["gene_id_bundled"]): r["non_overlapping"] for r in rows}

    assert flags(actual) == flags(expected), "Shipped marker independence flags are stale"
