"""EVID-REVIEW-02: the --dedup mode and the 5' read-start metrics."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan.evidence import _alignment_qc_from_text
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
