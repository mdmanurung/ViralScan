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


def _host_run(rnd: Path, host_called: list[str]):
    """A twostep arm: host-only STARsolo call + its emptyDrops table (knee / inflection per row)."""
    inner = _run(rnd, "twostep", "S1", "S1", ["X"], [("A", "EBV", 5.0)])
    (inner / "results" / "host_called_cells.tsv").write_text(
        "barcode\n" + "".join(f"{b}\n" for b in host_called)
    )
    raw = inner / "host_filtered" / "star_tmp" / "Solo.out" / "GeneFull" / "raw"
    raw.mkdir(parents=True)
    (raw / "emptydrops_cells.tsv").write_text(
        "barcode\ttotal\tFDR\tis_cell\tknee\tinflection\n"
        + "".join(f"{b}\t9\t0\tTRUE\t500.5\t60.25\n" for b in host_called)
        + "Q\t1\tNA\tFALSE\t500.5\t60.25\n"
    )
    (inner / "config.yaml").write_text("emptydrops_fdr: 0.01\nemptydrops_seed: 100\n")
    meta = rnd / "runs" / "ds" / "twostep" / "S1.meta"
    meta.mkdir()
    (meta / "code_sha.txt").write_text("abc123\n")
    return inner


def test_host_called_source_into_an_out_dir_leaves_the_round_untouched(tmp_path):
    rnd, out = tmp_path / "round", tmp_path / "anchors"
    _run(rnd, "combined_off", "S1", "S1", ["A", "B", "C"], [("A", "EBV", 3.0), ("H", "EBV", 8.0)])
    _host_run(rnd, ["H", "A"])
    rows = ref.freeze(rnd, source="host_called", out_dir=out)
    assert (out / "ds__S1.tsv").read_text() == "barcode\nA\nH\n"
    assert not (rnd / "reference_cells").exists()
    r = rows[0]
    assert (r["n_cells"], r["method"]) == (2, "emptydrops-host-only-starsolo-genefull")
    assert (r["barcodes_tested"], r["knee"], r["inflection"]) == (3, "500.5", "60.25")
    assert (r["emptydrops_fdr"], r["emptydrops_seed"], r["code_sha"]) == ("0.01", "100", "abc123")
    assert r["source"].endswith("twostep/S1/S1/results/host_called_cells.tsv")
    assert (out / "index.tsv").read_text().splitlines()[1].startswith("ds\tS1\t2\t")
    assert r["status"] == "frozen"
    assert (
        ref.freeze(rnd, source="host_called", out_dir=out, provisional=("ds",))[0]["status"]
        == "provisional"
    )
    # idempotent; a changed host call is refused
    assert (
        ref.freeze(rnd, source="host_called", out_dir=out)[0]["frozen_sha256"] == r["frozen_sha256"]
    )
    (rnd / "runs/ds/twostep/S1/S1/results/host_called_cells.tsv").write_text("barcode\nH\n")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        ref.freeze(rnd, source="host_called", out_dir=out)

    # CLI: cell_level goes under out-dir and is scored over the anchors (A and H, not B / C)
    (rnd / "runs/ds/twostep/S1/S1/results/host_called_cells.tsv").write_text("barcode\nH\nA\n")
    assert ref.main([str(rnd), "--source", "host_called", "--out-dir", str(out)]) == 0
    assert not (rnd / "cell_level").exists()
    with gzip.open(out / "cell_level" / "ds__S1.tsv.gz", "rt") as fh:
        got = {(r["arm"], r["barcode"]) for r in csv.DictReader(fh, delimiter="\t")}
    assert got == {("combined_off", "A"), ("combined_off", "H"), ("twostep", "A")}
