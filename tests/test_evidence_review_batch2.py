"""EVID-REVIEW-02: the --dedup mode and the 5' read-start metrics."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan.evidence import _alignment_qc_from_text, interpretation_flags
from viralscan.scripts import evidence_run


def test_dedup_mode_sets_the_profile_bam_only(tmp_path: Path) -> None:
    made: dict[str, str] = {}

    def fake_dedup(bam: str, out_bam: str, mode: str, threads: int = 4) -> str:
        made[mode] = out_bam
        return out_bam

    with patch.object(evidence_run, "deduplicate_bam", fake_dedup):
        profile, umi = evidence_run.dedup_layers("raw.bam", tmp_path, "none", 1)
    assert profile == made["none"]
    assert umi == made["umi"]
    assert profile != umi
    with patch.object(evidence_run, "deduplicate_bam", fake_dedup):
        assert evidence_run.dedup_layers("raw.bam", tmp_path, "umi", 1)[0] == made["umi"]


def test_hotspot_window_uses_the_five_prime_start() -> None:
    header = "@SQ\tSN:VIRUS|v\tLN:1000\n"
    # Reverse read POS 1 spanning 100 bp has its 5' end at 100; leftmost POS would be 1.
    sam = (
        "CB0_U0_0\t16\tVIRUS|v\t1\t60\t100M\t*\t0\t0\t*\t*\tNM:i:0\n"
        "CB1_U1_1\t0\tVIRUS|v\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
        "CB2_U2_2\t0\tVIRUS|v\t140\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
    )
    (row,) = _alignment_qc_from_text(header, sam, "")
    assert row["max_50bp_window_fraction"] == pytest.approx(1.0)


@pytest.mark.parametrize(
    "method, counted",
    [("host-conservative", False), ("equal", True), ("unique-weighted", True), ("em-cell", True)],
)
def test_host_virus_ambiguous_reads_are_counted_as_the_primary_matrix_counts_them(
    method: str, counted: bool
) -> None:
    # multimapping.py gives these reads a share above zero for every method but host-conservative
    # (pseudocount > 0), so the evidence replay's counted flag agrees with the primary matrix.
    from viralscan.evidence import _target_assignment

    _, weight, _ = _target_assignment([1, 9], {1}, {1}, method)
    assert (weight > 0) is counted


def _flags(*args, **kwargs) -> dict[str, dict[str, object]]:
    return {r["flag"]: r for r in interpretation_flags(*args, **kwargs)}


HOST_ROW = {"reference_class": "host", "reads": 40, "host_competitive_fraction": 1.0}


def test_all_host_best_reads_host_best_all_not_not_assessed() -> None:
    row = _flags([HOST_ROW], [], [])["host_homology"]
    assert row["status"] == "host_best_all"
    assert row["metric"] == 1.0
    # BLAST rows do not turn it back into not_flagged
    blast = [{"viral_minus_host_bitscore": "50", "low_complexity": "false"}]
    assert _flags([HOST_ROW], blast, [])["host_homology"]["status"] == "host_best_all"
    competitor = [{"diagnostic_class": "target_specific", "low_complexity": "false"}]
    assert _flags([HOST_ROW], competitor, [])["host_homology"]["status"] == "host_best_all"


def test_host_homology_stays_not_assessed_without_molecules_or_alignments() -> None:
    assert _flags([], [], [])["host_homology"]["status"] == "not_assessed"
    assert _flags([{**HOST_ROW, "reads": 0}], [], [])["host_homology"]["status"] == "not_assessed"


def test_host_best_all_ignored_when_a_target_row_exists() -> None:
    target = {"reference_class": "virus", "reads": 5, "host_competitive_fraction": 0.1}
    assert _flags([HOST_ROW, target], [], [])["host_homology"]["status"] == "not_flagged"


def test_blast_fractions_divide_by_all_queried_reads_and_report_counts() -> None:
    # 2 low-complexity hits out of 10 queried reads: 0.2, not 2 / 4 = 0.5
    blast = [{"viral_minus_host_bitscore": "-5", "low_complexity": "true"}] * 2 + [
        {"viral_minus_host_bitscore": "9", "low_complexity": "false"}
    ] * 2
    flags = _flags([], blast, [], n_queried=10)
    assert flags["low_complexity"]["metric"] == pytest.approx(0.2)
    assert flags["low_complexity"]["status"] == "not_flagged"
    assert (flags["low_complexity"]["n_queried"], flags["low_complexity"]["n_hit"]) == (10, 4)
    assert flags["host_homology"]["metric"] == pytest.approx(0.2)
    assert (flags["host_homology"]["n_queried"], flags["host_homology"]["n_hit"]) == (10, 4)
    # without n_queried the denominator is the rows given
    assert _flags([], blast, [])["low_complexity"]["metric"] == pytest.approx(0.5)
    assert _flags([], [], [])["low_complexity"]["n_queried"] == ""


def test_competitor_blast_fractions_count_no_hit_rows_and_carry_counts() -> None:
    blast = [{"diagnostic_class": "host_preferred", "low_complexity": "false"}] * 3 + [
        {"diagnostic_class": "no_hits", "low_complexity": "false"}
    ] * 5
    flags = _flags([], blast, [])
    assert flags["host_homology"]["metric"] == pytest.approx(3 / 8)
    assert (flags["host_homology"]["n_queried"], flags["host_homology"]["n_hit"]) == (8, 3)
    assert flags["related_virus_competition"]["n_queried"] == 8


def test_ambiguous_fraction_keeps_weight_zero_reads() -> None:
    lineage = [
        {"ambiguity_class": "unique_target", "assigned_weight": 1.0},
        {"ambiguity_class": "host_virus_ambiguous", "assigned_weight": 0.0},
    ]
    row = _flags([], [], lineage)["sibling_or_host_ambiguity"]
    assert row["metric"] == pytest.approx(0.5)
    assert row["status"] == "flagged"
