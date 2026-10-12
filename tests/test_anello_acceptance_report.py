"""Candidate rows must not become detections in the acceptance report (AUD-21)."""

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "anello_acceptance_report",
    Path(__file__).parent.parent / "scripts" / "anello_acceptance_report.py",
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


@pytest.mark.parametrize("with_detected", [False, True])
def test_arm_summary_filters_candidates_and_accepts_legacy_rows(tmp_path, with_detected):
    results = tmp_path / "sample" / "results"
    results.mkdir(parents=True)
    columns = (
        "virus_name\talignment_status\tviral_molecules_total_est\t"
        "alignment_molecules_unique\tdetection_source\talignment_homopolymer_fraction"
    )
    rows = [
        "candidate\tmeasured\t100\t50\tkallisto\t0.0",
        "called\tmeasured\t4\t3\tkallisto\t0.0",
        "alignment-only\tmeasured\t0\t2\talignment_only\t0.0",
    ]
    if with_detected:
        columns += "\tdetected"
        rows = [row + "\t" + value for row, value in zip(rows, [" False ", "True", ""])]
    (results / "viral_summary.tsv").write_text(columns + "\n" + "\n".join(rows) + "\n")
    summary = mod.arm_summary(tmp_path)
    assert summary["kallisto_anellovirus_molecules"] == (4 if with_detected else 104)
    assert summary["alignment_anellovirus_molecules"] == (5 if with_detected else 55)
    assert len(summary["rows"]) == (2 if with_detected else 3)
