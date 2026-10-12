"""DSR-16: round 2 vs round 1 comparison on tiny synthetic round directories."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "dsr_round_compare", Path(__file__).resolve().parents[1] / "scripts" / "dsr_round_compare.py"
)
cmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cmp)

HEAD = "virus_name\tviral_molecules_total_est\tinfected_called\tn_called_cells\tpct_infected_called"


def _run(root: Path, ds: str, arm: str, sample: str, calls, detected_col=True, final=True):
    """calls: [(virus, molecules, infected, cells, pct, detected)]"""
    res = root / "runs" / ds / arm / sample / sample / "results"
    res.mkdir(parents=True)
    head = HEAD + ("\tdetected" if detected_col else "")
    lines = [head] + [
        f"{v}\t{m}\t{i}\t{c}\t{p}" + (f"\t{d}" if detected_col else "")
        for v, m, i, c, p, d in calls
    ]
    (res / "viral_summary.tsv").write_text("\n".join(lines) + "\n")
    if arm == "twostep_v2" and final:
        (root / "runs" / ds / arm / sample / "run_complete.json").write_text("{}")


def _tsv(prefix: Path):
    with Path(f"{prefix}.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_twostep_baseline_flags_and_detected_filter(tmp_path):
    r1, r2 = tmp_path / "r1", tmp_path / "r2"
    # round 1: twostep (pre-DEF-05) must be ignored, twostep_v2 is the baseline
    _run(r1, "ebv", "twostep", "S1", [("EBV", 5, 5, 10, 50.0, "True")])
    _run(
        r1,
        "ebv",
        "twostep_v2",
        "S1",
        [("EBV", 9, 8, 10, 80.0, "True"), ("HSV", 2, 2, 10, 20.0, "True")],
    )
    _run(r1, "ebv", "combined_off", "S1", [("EBV", 9, 8, 10, 80.0, "True")])
    _run(
        r1, "ebv", "combined_off", "GONE", [("EBV", 1, 1, 10, 10.0, "True")]
    )  # dropped from round 2
    _run(
        r2,
        "ebv",
        "twostep",
        "S1",
        [("EBV", 10, 9, 10, 90.0, "True"), ("HPV", 3, 3, 10, 30.0, "True")],
    )
    # a candidate-only row (detected False) is not a call
    _run(
        r2,
        "ebv",
        "combined_off",
        "S1",
        [("EBV", 9, 8, 10, 80.0, "True"), ("VARV", 1, 1, 10, 10.0, "False")],
    )
    _run(r1, "sfl_tonsil", "combined_off", "x223", [("EBV", 4, 4, 40, 10.0, "True")])
    _run(r2, "sfl_tonsil", "combined_off", "x223", [("EBV", 4, 4, 40, 10.0, "True")])

    assert cmp.main([str(r1), str(r2), "--out", str(tmp_path / "out" / "cmp")]) == 0
    rows = {
        (r["dataset"], r["sample"], r["arm"], r["virus"]): r for r in _tsv(tmp_path / "out" / "cmp")
    }

    assert rows[("ebv", "S1", "twostep", "EBV")]["status"] == "both"
    assert (
        rows[("ebv", "S1", "twostep", "EBV")]["r1_viral_molecules_total_est"] == "9"
    )  # twostep_v2, not v1's 5
    assert rows[("ebv", "S1", "twostep", "EBV")]["delta_pct_infected_called"] == "10.0000"
    assert rows[("ebv", "S1", "twostep", "HSV")]["status"] == "disappeared"
    assert rows[("ebv", "S1", "twostep", "HPV")]["status"] == "appeared"
    assert ("ebv", "S1", "combined_off", "VARV") not in rows
    assert rows[("ebv", "GONE", "combined_off", "EBV")]["status"] == "run_only_round1"
    assert rows[("sfl_tonsil", "x223", "combined_off", "EBV")]["group"] == "institutional"

    md = (tmp_path / "out" / "cmp.md").read_text()
    assert "## Public datasets" in md and "## Institutional datasets" in md
    assert "| ebv | S1 | twostep | HPV | appeared |" in md
    assert "VARV" not in md


def test_unfinished_twostep_v2_is_not_a_baseline_and_legacy_summaries_counted(tmp_path):
    r1, r2 = tmp_path / "r1", tmp_path / "r2"
    _run(
        r1,
        "ebv",
        "twostep_v2",
        "S1",
        [("EBV", 9, 8, 10, 80.0, "")],
        detected_col=False,
        final=False,
    )
    _run(r1, "ebv", "combined_off", "S1", [("EBV", 9, 8, 10, 80.0, "")], detected_col=False)
    _run(r2, "ebv", "twostep", "S1", [("EBV", 9, 8, 10, 80.0, "True")])
    assert cmp.main([str(r1), str(r2), "--out", str(tmp_path / "cmp")]) == 0
    rows = {(r["arm"], r["virus"]): r for r in _tsv(tmp_path / "cmp")}
    assert rows[("twostep", "EBV")]["status"] == "run_only_round2"
    assert rows[("combined_off", "EBV")]["status"] == "run_only_round1"
    assert "without a `detected` column" in (tmp_path / "cmp.md").read_text()
