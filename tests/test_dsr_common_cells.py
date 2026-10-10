"""DSR-16: every arm is scored over the combined_off called cells."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "dsr_common_cells", Path(__file__).resolve().parents[1] / "scripts" / "dsr_common_cells.py"
)
dsr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dsr)


def _sample(root: Path, arm: str, called: list[str], per_cell: list[tuple[str, str, float]]):
    res = root / "runs" / "ebv" / arm / "S1" / "S1" / "results"
    res.mkdir(parents=True)
    (res / "called_cells.tsv").write_text("barcode\n" + "".join(f"{b}\n" for b in called))
    (res / "per_cell_viral.tsv").write_text(
        "barcode\tvirus_name\tviral_molecules_total_est\n"
        + "".join(f"{b}\t{v}\t{n}\n" for b, v, n in per_cell)
    )


def test_arms_share_the_combined_off_denominator(tmp_path):
    _sample(
        tmp_path,
        "combined_off",
        ["A", "B", "C", "D"],
        [("A", "EBV", 3.0), ("B", "EBV", 1.0), ("Z", "EBV", 9.0)],
    )
    # twostep called only the infected cells (2 of 2 = 100 %); over the reference it is 2/4.
    _sample(
        tmp_path, "twostep", ["A", "B"], [("A", "EBV", 4.0), ("B", "EBV", 2.0), ("Q", "EBV", 5.0)]
    )
    got = {(r["arm"], r["virus"]): r for r in dsr.rows(tmp_path)}
    assert got[("twostep", "EBV")]["n_reference_cells"] == 4
    assert got[("twostep", "EBV")]["pct_infected_reference"] == "50.0000"
    assert got[("twostep", "EBV")]["molecules_in_reference"] == "6.0"  # Q is outside the set
    assert got[("combined_off", "EBV")]["infected_in_reference"] == 2  # Z is outside the set


def test_vendor_cells_replace_the_emptydrops_set(tmp_path):
    _sample(tmp_path, "combined_off", ["A", "B", "C", "D"], [("A", "EBV", 3.0), ("Z", "EBV", 9.0)])
    bc = tmp_path / "barcodes.tsv"
    bc.write_text("A-1\nZ-1\n")
    got = list(dsr.rows(tmp_path, {"ebv/S1": bc}))
    assert got[0]["n_reference_cells"] == 2 and got[0]["infected_in_reference"] == 2


def test_twostep_v2_counts_only_once_finished(tmp_path):
    _sample(tmp_path, "combined_off", ["A", "B"], [("A", "EBV", 1.0)])
    _sample(tmp_path, "twostep", ["A"], [("A", "EBV", 5.0)])
    _sample(tmp_path, "twostep_v2", ["A"], [("A", "EBV", 7.0)])  # stale copy of v1 until the marker
    arms = {r["arm"] for r in dsr.rows(tmp_path)}
    assert arms == {"combined_off", "twostep"}
    (tmp_path / "runs" / "ebv" / "twostep_v2" / "S1" / "run_complete.json").write_text("{}")
    assert {r["arm"] for r in dsr.rows(tmp_path)} == {"combined_off", "twostep_v2"}


def test_twostep_calls_need_a_viral_best_verdict(tmp_path):
    """F-028: HPV77 is called by both arms on host reads; only a viral_best verdict lets it through."""
    _sample(tmp_path, "combined_off", ["A", "B"], [("A", "EBV", 1.0), ("A", "HPV 77", 1.0)])
    _sample(
        tmp_path,
        "twostep",
        ["A", "B"],
        [("A", "EBV", 5.0), ("A", "HPV 77", 2.0), ("B", "HPV 16", 2.0)],
    )

    def status(verdicts):
        return {
            (r["arm"], r["virus"]): r["evidence_status"]
            for r in dsr.rows(tmp_path, verdicts=verdicts)
        }

    got = status({})
    assert got[("combined_off", "EBV")] == "not_gated"
    assert (
        got[("twostep", "HPV 77")] == "needs_evidence"
    )  # called by combined_off too, still unchecked
    got = status(
        {"ebv__S1__combined_off__HPV_77": "host_best", "ebv__S1__twostep__HPV_16": "viral_best"}
    )
    assert (
        got[("twostep", "HPV 77")] == "rejected_host_best"
    )  # falls back to combined_off's verdict
    assert got[("twostep", "HPV 16")] == "verified"
