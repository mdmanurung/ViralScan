"""Q17: cells the host+virus call keeps but the host-only anchor drops, ranked by viral load."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "dsr_host_shutoff_diagnostic",
    Path(__file__).resolve().parents[1] / "scripts" / "dsr_host_shutoff_diagnostic.py",
)
diag = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(diag)


def _lib(root: Path, anchors: Path, anchor: list[str], called: list[str], per_cell, ds="hsv1"):
    res = root / "runs" / ds / "combined_off" / "S1" / "S1" / "results"
    res.mkdir(parents=True)
    (res / "called_cells.tsv").write_text("barcode\n" + "".join(f"{b}\n" for b in called))
    (res / "per_cell_viral.tsv").write_text(
        "barcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction\tis_called_cell\n"
        + "".join(f"{b}\t{v}\t{n}\t1000\t0.1\tTrue\n" for b, v, n in per_cell)
    )
    anchors.mkdir(exist_ok=True)
    (anchors / f"{ds}__S1.tsv").write_text("barcode\n" + "".join(f"{b}\n" for b in anchor))


def test_dropped_infected_cells_are_ranked_and_flagged(tmp_path):
    _lib(
        tmp_path, tmp_path / "anc", anchor=["A", "B", "C"], called=["A", "D", "E", "F"],
        per_cell=[("A", "Human herpesvirus 1", 5), ("D", "Human herpesvirus 1", 150),
                  ("E", "Human herpesvirus 1", 12), ("F", "Human herpesvirus 1", 1),
                  ("D", "Human herpesvirus 2", 99), ("B", "Human herpesvirus 1", 30)],
    )  # fmt: skip
    ((row, top),) = list(diag.diagnose(tmp_path, tmp_path / "anc"))
    assert (row["n_both"], row["n_combined_only"], row["n_anchor_only"]) == (1, 3, 2)
    assert (
        row["combined_only_with_virus"],
        row["combined_only_ge10"],
        row["combined_only_ge100"],
    ) == (3, 2, 1)
    assert row["combined_only_max"] == "150.0" and row["combined_only_molecules"] == "163.0"
    assert row["pct_combined_called_molecules_dropped"] == f"{100 * 163 / 168:.2f}"
    assert (row["anchor_only_with_virus"], row["anchor_only_ge10"]) == (1, 1)  # B
    assert row["flag"] == "HIGH_LOAD_CELLS_LOST" and row["barcode_format"] == "len=1;suffix=no"
    assert [t["barcode"] for t in top] == [
        "D",
        "E",
        "F",
    ]  # by load, HHV-2 of D is not the expected virus


def test_barcode_format_mismatch_is_reported_and_cli_writes_tables(tmp_path):
    _lib(tmp_path, tmp_path / "anc", ["AAAA"], ["AAAA-1"], [("AAAA-1", "Human herpesvirus 1", 1)])
    assert (
        diag.main(
            [str(tmp_path), "--anchors", str(tmp_path / "anc"), "--out-dir", str(tmp_path / "o")]
        )
        == 0
    )
    with (tmp_path / "o" / "host_shutoff_diagnostic.tsv").open() as fh:
        (row,) = list(csv.DictReader(fh, delimiter="\t"))
    assert row["barcode_format"].startswith("MISMATCH") and row["flag"] == ""
