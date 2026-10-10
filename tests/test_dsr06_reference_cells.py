"""CELLS-01: one frozen emptyDrops call per sample, used for every arm and every cell-level table."""

from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1] / "scripts"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _ROOT / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


common, ref = _load("dsr_common_cells"), _load("dsr06_reference_cells")


def _run(
    rnd: Path,
    arm: str,
    sample: str,
    inner: str,
    called: list[str],
    per_cell: list[tuple[str, str, float]],
):
    res = rnd / "runs" / "ds" / arm / sample / inner / "results"
    res.mkdir(parents=True)
    (res / "called_cells.tsv").write_text("barcode\n" + "".join(f"{b}\n" for b in called))
    (res / "per_cell_viral.tsv").write_text(
        "barcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction\n"
        + "".join(f"{b}\t{v}\t{n}\t100\t0.01\n" for b, v, n in per_cell)
    )
    return res.parent


def test_freeze_is_idempotent_and_refuses_a_changed_call(tmp_path):
    inner = _run(tmp_path, "combined_off", "S1", "S1", ["B", "A"], [("A", "EBV", 2.0)])
    (inner / "config.yaml").write_text("emptydrops_fdr: 0.01\nemptydrops_lower: 100\n")
    rows = ref.freeze(tmp_path)
    frozen = tmp_path / "reference_cells" / "ds__S1.tsv"
    assert frozen.read_text() == "barcode\nA\nB\n"  # sorted
    assert (
        rows[0]["n_cells"] == 2
        and rows[0]["emptydrops_fdr"] == "0.01"
        and rows[0]["emptydrops_lower"] == "100"
    )
    assert (
        ref.freeze(tmp_path)[0]["frozen_sha256"] == rows[0]["frozen_sha256"]
    )  # same call: no change

    (inner / "results" / "called_cells.tsv").write_text("barcode\nA\nB\nC\n")  # the run was redone
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        ref.freeze(tmp_path)
    assert ref.freeze(tmp_path, force=True)[0]["n_cells"] == 3


def test_every_arm_is_scored_over_the_frozen_cells(tmp_path):
    _run(
        tmp_path,
        "combined_off",
        "S1",
        "S1",
        ["A", "B", "C", "D"],
        [("A", "EBV", 3.0), ("Z", "EBV", 9.0)],
    )
    _run(
        tmp_path,
        "twostep_v2",
        "S1",
        "S1",
        ["A"],
        [("A", "EBV", 5.0), ("B", "EBV", 1.0), ("Q", "EBV", 7.0)],
    )
    (tmp_path / "runs" / "ds" / "twostep_v2" / "S1" / "run_complete.json").write_text("{}")
    ref.freeze(tmp_path)
    # the arm's own called set is one cell, but the frozen reference has four
    got = {(r["arm"], r["virus"]): r for r in common.rows(tmp_path)}
    assert (
        got["twostep_v2", "EBV"]["n_reference_cells"] == 4
        and got["twostep_v2", "EBV"]["infected_in_reference"] == 2
    )
    assert got["combined_off", "EBV"]["infected_in_reference"] == 1  # Z is not a called cell

    # the frozen list wins over a later change of the source call; an explicit vendor set still wins over both
    (
        tmp_path / "runs" / "ds" / "combined_off" / "S1" / "S1" / "results" / "called_cells.tsv"
    ).write_text("barcode\nZ\n")
    assert next(common.rows(tmp_path))["n_reference_cells"] == 4
    vendor = tmp_path / "vendor.tsv"
    vendor.write_text("Z-1\n")
    assert next(common.rows(tmp_path, {"ds/S1": vendor}))["n_reference_cells"] == 1


def test_cell_level_table_holds_only_frozen_cells_of_final_arms(tmp_path):
    _run(
        tmp_path,
        "combined_off",
        "S1",
        "S1",
        ["A", "B"],
        [("A", "EBV", 3.0), ("B", "EBV", 0.0), ("Z", "EBV", 9.0)],
    )
    _run(
        tmp_path, "twostep_v2", "S1", "S1", ["A"], [("A", "EBV", 5.0)]
    )  # no run_complete.json: still v1's copy
    ref.freeze(tmp_path)
    assert ref.export_cell_level(tmp_path) == {("ds", "S1"): 1}
    with gzip.open(tmp_path / "cell_level" / "ds__S1.tsv.gz", "rt") as fh:
        got = list(csv.DictReader(fh, delimiter="\t"))
    assert [(r["arm"], r["barcode"]) for r in got] == [
        ("combined_off", "A")
    ]  # B has 0 molecules, Z is not a cell


def test_run_root_name_differs_from_the_inner_directory(tmp_path):
    """x223 holds LUM-SJ-x223: it must be kept, and named by its run root."""
    _run(tmp_path, "combined_off", "x223", "LUM-SJ-x223", ["A", "B"], [("A", "EBV", 4.0)])
    ref.freeze(tmp_path)
    assert (tmp_path / "reference_cells" / "ds__x223.tsv").is_file()
    assert [(r["sample"], r["arm"], r["infected_in_reference"]) for r in common.rows(tmp_path)] == [
        ("x223", "combined_off", 1)
    ]
