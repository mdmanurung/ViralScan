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
    _sample(tmp_path, "combined_off", ["A", "B", "C", "D"],
            [("A", "EBV", 3.0), ("B", "EBV", 1.0), ("Z", "EBV", 9.0)])
    # twostep called only the infected cells (2 of 2 = 100 %); over the reference it is 2/4.
    _sample(tmp_path, "twostep", ["A", "B"], [("A", "EBV", 4.0), ("B", "EBV", 2.0), ("Q", "EBV", 5.0)])
    got = {(r["arm"], r["virus"]): r for r in dsr.rows(tmp_path)}
    assert got[("twostep", "EBV")]["n_reference_cells"] == 4
    assert got[("twostep", "EBV")]["pct_infected_reference"] == "50.0000"
    assert got[("twostep", "EBV")]["molecules_in_reference"] == "6.0"  # Q is outside the set
    assert got[("combined_off", "EBV")]["infected_in_reference"] == 2  # Z is outside the set
