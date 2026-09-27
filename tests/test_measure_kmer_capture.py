"""Tests for scripts/measure_kmer_capture.py (PLAN CAT-11).

The tool's job is to say whether a reference panel would actually see a strain
that is not in it. These tests pin the arithmetic on inputs whose answers are
known by hand, so a refactor cannot quietly change what "coverage" means.
"""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "measure_kmer_capture",
    Path(__file__).resolve().parents[1] / "scripts" / "measure_kmer_capture.py",
)
assert _SPEC and _SPEC.loader
mkc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(mkc)


def _write(path: Path, records: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
    return path


@pytest.fixture
def genomes() -> dict[str, str]:
    rng = random.Random(1234)
    return {
        "A.1": "".join(rng.choice("ACGT") for _ in range(1500)),
        "B.1": "".join(rng.choice("ACGT") for _ in range(1500)),
    }


class TestKmers:
    def test_counts_distinct_windows(self) -> None:
        assert len(mkc.kmers("ACGT" * 10, k=4)) == 4

    def test_skips_ambiguity_codes(self) -> None:
        """An N cannot match a concrete k-mer in an index, so it must not count."""
        assert mkc.kmers("ACGTN", k=4) == {"ACGT"}

    def test_sequence_shorter_than_k_has_no_kmers(self) -> None:
        assert mkc.kmers("ACG", k=31) == set()


class TestSelfCoverage:
    def test_genome_against_itself_is_total(self, tmp_path, genomes) -> None:
        fa = _write(tmp_path / "one.fa", {"A.1": genomes["A.1"]})
        rows, summary = mkc.measure(fa, fa)
        assert rows[0]["coverage"] == 1.0
        assert summary["median_coverage"] == 1.0

    def test_unrelated_panel_shares_nothing(self, tmp_path, genomes) -> None:
        panel = _write(tmp_path / "panel.fa", {"B.1": genomes["B.1"]})
        pop = _write(tmp_path / "pop.fa", {"A.1": genomes["A.1"]})
        rows, summary = mkc.measure(panel, pop)
        assert rows[0]["coverage"] == 0.0
        assert rows[0]["p_fragment"] == 0.0
        assert summary["zero_coverage_fraction"] == 1.0


class TestLeaveOneOut:
    """The headline number: what a panel sees of a strain it does not contain."""

    def test_own_genome_does_not_count_towards_its_own_capture(self, tmp_path, genomes) -> None:
        panel = _write(tmp_path / "panel.fa", genomes)
        pop = _write(tmp_path / "pop.fa", {"A.1": genomes["A.1"]})
        rows, _ = mkc.measure(panel, pop)
        assert rows[0]["coverage"] == 1.0, "A.1 is in the panel"
        assert rows[0]["leave_one_out"] < 0.01, "but must not cover itself"

    def test_a_duplicate_genome_keeps_coverage(self, tmp_path, genomes) -> None:
        """Removing one copy of a duplicated genome removes no k-mer."""
        panel = _write(tmp_path / "panel.fa", {"A.1": genomes["A.1"], "COPY.1": genomes["A.1"]})
        pop = _write(tmp_path / "pop.fa", {"A.1": genomes["A.1"]})
        rows, _ = mkc.measure(panel, pop)
        assert rows[0]["leave_one_out"] == 1.0
        assert rows[0]["p_fragment"] == 1.0

    def test_a_near_relative_gives_partial_capture(self, tmp_path, genomes) -> None:
        """A relative sharing a conserved block covers that block and no more."""
        target = genomes["A.1"]
        relative = target[:500] + genomes["B.1"][500:]
        panel = _write(tmp_path / "panel.fa", {"A.1": target, "REL.1": relative})
        pop = _write(tmp_path / "pop.fa", {"A.1": target})
        rows, _ = mkc.measure(panel, pop)
        loo = float(rows[0]["leave_one_out"])
        assert 0.0 < loo < 1.0
        # ~470 of ~1470 shared 31-mers fall in the conserved 500 bp prefix.
        assert 0.25 < loo < 0.40


class TestFragmentCapture:
    def test_reports_the_fraction_of_reads_that_would_align(self, tmp_path) -> None:
        """p_fragment is about reads, not k-mers: a conserved half captures ~half."""
        rng = random.Random(7)
        shared = "".join(rng.choice("ACGT") for _ in range(1000))
        unique = "".join(rng.choice("ACGT") for _ in range(1000))
        target = shared + unique
        panel = _write(tmp_path / "panel.fa", {"REL.1": shared})
        pop = _write(tmp_path / "pop.fa", {"T.1": target})
        rows, _ = mkc.measure(panel, pop, read_length=90)
        assert 0.4 < float(rows[0]["p_fragment"]) < 0.6

    def test_genome_shorter_than_one_read_is_still_scored(self, tmp_path) -> None:
        rng = random.Random(11)
        short = "".join(rng.choice("ACGT") for _ in range(60))
        fa = _write(tmp_path / "s.fa", {"S.1": short})
        rows, _ = mkc.measure(fa, fa, read_length=90)
        assert rows[0]["coverage"] == 1.0


class TestGrouping:
    def test_per_group_summary(self, tmp_path, genomes) -> None:
        groups_tsv = tmp_path / "groups.tsv"
        groups_tsv.write_text("accession\tgenus\nA.1\tAlpha\nB.1\tBeta\n")
        groups = mkc.load_groups(groups_tsv, "accession", "genus")
        assert groups["A.1"] == "Alpha"
        assert groups["A"] == "Alpha", "bare accession must resolve too"
        panel = _write(tmp_path / "panel.fa", genomes)
        pop = _write(tmp_path / "pop.fa", genomes)
        _, summary = mkc.measure(panel, pop, groups=groups)
        assert set(summary["groups"]) == {"Alpha", "Beta"}
        assert summary["groups"]["Alpha"]["genomes"] == 1


def test_self_check_entrypoint_passes() -> None:
    assert mkc.self_check() == 0
