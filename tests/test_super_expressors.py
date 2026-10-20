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
    (row,) = se.summarise(p, {"A", "B"}, min_umi=10)  # explicit complete cell set
    assert (row["n_cells"], row["n_super_expressors"], row["pct_cells_super"]) == (2, 1, 50.0)
    assert row["share_molecules_in_super"] == round(10 / 19, 4)
    assert row["max_viral_per_10k"] == 100.0
    cells = tmp_path / "cells.tsv"
    cells.write_text("A\nB\nC\nD\n")  # frozen set: C counts, D has no viral UMI
    (row,) = se.summarise(p, se.read_cells(cells), min_umi=10)
    assert (row["n_cells"], row["n_super_expressors"], row["n_ge_3"]) == (4, 2, 3)
    assert row["n_cells_viral_umi"] == 3


def test_complete_universe_legacy_and_fractional(tmp_path, caplog):
    rows = [(f"C{i}", "V", 10 if i < 2 else 9.6, 100, 0.1, True) for i in range(50)]
    p = _write(tmp_path, rows)
    cells = {f"C{i}" for i in range(1000)}
    (row,) = se.summarise(p, cells)
    assert (row["n_cells"], row["n_super_expressors"], row["pct_cells_super"]) == (1000, 2, 0.2)
    (legacy,) = se.summarise(p, denominator="viral_positive_called")
    assert legacy["pct_cells_super"] == 4.0
    assert legacy["denominator_kind"] == "viral_positive_called"
    assert "not prevalence" in caplog.text


def test_auto_resolution_missing_empty_and_suffix(tmp_path):
    import pytest

    results = tmp_path / "results"
    results.mkdir()
    p = _write(results, [("A", "V", 10, 100, 0.1, True)])
    with pytest.raises(SystemExit, match="called_cells.tsv"):
        se.summarise(p)
    (results / "called_cells.tsv").write_text("barcode\nA\nB\n")
    (row,) = se.summarise(p)
    assert row["n_cells"] == 2
    with pytest.raises(SystemExit, match="suffixes"):
        se.summarise(p, {"A-1", "B-1"})
    p.write_text(HEADER)
    out = tmp_path / "summary.tsv"
    assert se.main(["--per-cell-viral", str(p), "--out", str(out)]) == 0
    assert out.read_text().strip().split("\t") == list(se.SUMMARY_FIELDS)
