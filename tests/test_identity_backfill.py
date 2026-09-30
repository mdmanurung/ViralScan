"""Backfill of ``results/virus_identity.tsv`` for Runs made before it existed (PLAN MECH-A 4c)."""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests.test_virus_identity import ANELLO_GENUS, CATALOGUE, _combined_t2g, _write_t2g
from viralscan import virus_identity
from viralscan.virus_identity import (
    TABLE_FILENAME,
    VirusIdentityTable,
    backfill_identity_table,
    build_identity_table,
    write_identity_table,
)

GENE_SET = {"NC_007605.1_gene1", "NC_007605.1_gene2", "NC_009334_gene1", "MY_VIRUS_g1"}


@pytest.fixture(autouse=True)
def _packaged_catalogue(monkeypatch):
    # Keep the test independent of the packaged catalogue.
    real = build_identity_table

    def patched(t2g, gtf, **kw):
        kw.setdefault("catalogue_rows", CATALOGUE)
        kw.setdefault("anello_genus", ANELLO_GENUS)
        return real(t2g, gtf, **kw)

    monkeypatch.setattr(virus_identity, "build_identity_table", patched)


def _old_run(root: Path, t2g: Path, gene_ids=GENE_SET) -> SimpleNamespace:
    """A Run directory as written before 35940ec: ``analysis.txt``, no table."""
    (root / "log").mkdir(parents=True)
    (root / "results").mkdir()
    gene_list = root / "log" / "analysis.txt"
    gene_list.write_text("".join(f"{g}\n" for g in sorted(gene_ids)))
    os.utime(gene_list, ns=(1_000_000_000_000_000_000, 1_000_000_000_000_000_000))
    return SimpleNamespace(transcripts=str(t2g), index="", output=str(root) + os.sep)


class TestBackfillIdentityTable:
    def test_writes_the_table_from_the_runs_own_gene_list(self, tmp_path):
        t2g = _combined_t2g(tmp_path)
        config = _old_run(tmp_path / "run", t2g)
        out = backfill_identity_table(config)
        assert out == tmp_path / "run" / "results" / TABLE_FILENAME
        viral = set(VirusIdentityTable.read_tsv(out).viral_gene_ids())
        # MY_VIRUS_g1 is uncatalogued: viral only because analysis.txt lists it.
        assert {"NC_007605.1_gene1", "MY_VIRUS_g1"} <= viral
        assert not viral & {"ENSG1.1", "MYSTERY_g1"}

    def test_matches_the_table_the_analysis_rule_writes(self, tmp_path):
        # Backfill parity: same inputs give a byte-identical table.
        t2g = _combined_t2g(tmp_path)
        backfilled = backfill_identity_table(_old_run(tmp_path / "old", t2g))
        fresh_dir = tmp_path / "fresh"
        (fresh_dir / "results").mkdir(parents=True)
        fresh = write_identity_table(
            SimpleNamespace(transcripts=str(t2g), index="", output=str(fresh_dir)), GENE_SET
        )
        assert backfilled.read_bytes() == fresh.read_bytes()

    def test_takes_the_gene_lists_mtime_so_nothing_downstream_turns_stale(self, tmp_path):
        config = _old_run(tmp_path / "run", _combined_t2g(tmp_path))
        out = backfill_identity_table(config)
        gene_list = tmp_path / "run" / "log" / "analysis.txt"
        assert out.stat().st_mtime_ns == gene_list.stat().st_mtime_ns

    def test_leaves_an_existing_table_alone(self, tmp_path):
        config = _old_run(tmp_path / "run", _combined_t2g(tmp_path))
        existing = tmp_path / "run" / "results" / TABLE_FILENAME
        existing.write_text("sentinel\n")
        assert backfill_identity_table(config) is None
        assert existing.read_text() == "sentinel\n"

    def test_skips_a_run_without_a_gene_list(self, tmp_path):
        (tmp_path / "run").mkdir()
        config = SimpleNamespace(
            transcripts=str(_combined_t2g(tmp_path)), index="", output=str(tmp_path / "run")
        )
        assert backfill_identity_table(config) is None
        assert not (tmp_path / "run" / "results" / TABLE_FILENAME).exists()

    def test_a_run_with_no_viral_gene_raises(self, tmp_path):
        host_only = _write_t2g(tmp_path / "t2g.txt", ["ENST1.1\tENSG1.1\t\t\tENST1.1\t1\t400\t+"])
        config = _old_run(tmp_path / "run", host_only, gene_ids={"NOT_INDEXED"})
        with pytest.raises(ValueError):
            backfill_identity_table(config)
