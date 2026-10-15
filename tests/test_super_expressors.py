"""scripts/super_expressors.py: >=10 viral UMIs per cell (Lareau 2023), denominator = the cell set."""

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "super_expressors", Path(__file__).resolve().parents[1] / "scripts" / "super_expressors.py"
)
se = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(se)

HEADER = "barcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction\tis_called_cell\n"


def _write(tmp_path, rows):
    p = tmp_path / "per_cell_viral.tsv"
    p.write_text(HEADER + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return p


def test_threshold_denominator_and_share(tmp_path):
    p = _write(
        tmp_path,
        [
            ("A", "HHV6B", 10, 1000, 0.01, "True"),
            ("B", "HHV6B", 9, 1000, 0.009, "True"),
            ("C", "HHV6B", 90, 100, 0.9, "False"),  # outside the called set
        ],
    )
    (row,) = se.summarise(p, min_umi=10)  # denominator = is_called_cell (A, B)
    assert (row["n_cells"], row["n_super_expressors"], row["pct_cells_super"]) == (2, 1, 50.0)
    assert row["share_molecules_in_super"] == round(10 / 19, 4)
    assert row["max_viral_per_10k"] == 100.0
    cells = tmp_path / "cells.tsv"
    cells.write_text("A\nB\nC\nD\n")  # frozen set: C counts, D has no viral UMI
    (row,) = se.summarise(p, se.read_cells(cells), min_umi=10)
    assert (row["n_cells"], row["n_super_expressors"], row["n_ge_3"]) == (4, 2, 3)
    assert row["n_cells_viral_umi"] == 3
