"""DOSSIER-02: scripts/dsr05_dossier.py statuses and claim ladder on a synthetic round."""

from __future__ import annotations

import csv
import importlib.util
import json
import os
from pathlib import Path

SPEC_COLS = [
    "item_id",
    "section",
    "scope",
    "arms",
    "path_pattern",
    "check",
    "gates",
    "stale_if_older_than",
    "producer",
    "gap_note",
]
ROLE_COLS = [
    "dataset",
    "sample",
    "role",
    "expected",
    "vendor_barcodes",
    "cell_labels",
    "findings",
    "acceptance",
]

_spec = importlib.util.spec_from_file_location(
    "dsr05_dossier", Path(__file__).resolve().parents[1] / "scripts" / "dsr05_dossier.py"
)
dos = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dos)

SPEC = [
    (
        "A10",
        "A",
        "arm",
        "combined_off,twostep_v2",
        "{root}/run_complete.json",
        "exists",
        "C0",
        "",
        "",
        "",
    ),
    (
        "D01",
        "D",
        "arm",
        "combined_off",
        "{inner}/results/viral_summary.tsv",
        "exists",
        "C1",
        "",
        "",
        "",
    ),
    (
        "E07",
        "E",
        "call",
        "",
        "{ev}/molecule_verdict_summary.tsv",
        "exists",
        "C2",
        "",
        "",
        "needs VERDICT-01",
    ),
    ("E09", "E", "call", "", "verdicts.tsv", "verdict_row", "C2", "", "", ""),
    ("H01", "H", "sample", "", "{cell_labels}", "labels", "C3", "", "", ""),
    (
        "G01",
        "G",
        "sample",
        "",
        "{inner}/results/positive_control.json",
        "json_ne:status=not-configured",
        "N1",
        "",
        "",
        "",
    ),
]
ROLES = [
    ("ds", "S1", "positive", "", "", "", "", "x"),
    ("ds", "S2", "negative", "", "", "", "", "x"),
    ("ds", "S3", "n_a", "", "", "", "", "x"),
    ("newds", "-", "not_started", "", "", "", "", "blocked: access"),
]


def _write_tsv(path: Path, header: list[str], rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _run(round_dir: Path, sample: str, arm: str = "combined_off", molecules: str = "5") -> Path:
    root = round_dir / "runs" / "ds" / arm / sample
    (root / sample / "results").mkdir(parents=True)
    (root / "run_complete.json").write_text("{}")
    (root / sample / "results" / "positive_control.json").write_text(
        json.dumps({"status": "not-configured"})
    )
    _write_tsv(
        root / sample / "results" / "viral_summary.tsv",
        ["virus_name", "viral_molecules_total_est"],
        [("Epstein-Barr virus", molecules)],
    )
    return root


def _setup(tmp_path: Path, monkeypatch) -> Path:
    spec, roles = tmp_path / "spec.tsv", tmp_path / "roles.tsv"
    _write_tsv(spec, SPEC_COLS, SPEC)
    _write_tsv(roles, ROLE_COLS, ROLES)
    monkeypatch.setattr(dos, "SPEC", spec)
    monkeypatch.setattr(dos, "ROLES", roles)
    return tmp_path / "round"


def _status(index: list[dict], item: str, unit_prefix: str, arm: str | None = None) -> list[str]:
    return [
        r["status"]
        for r in index
        if r["item_id"] == item
        and r["unit"].startswith(unit_prefix)
        and (arm is None or r["arm"] == arm)
    ]


def test_statuses_and_ladder(tmp_path, monkeypatch):
    rnd = _setup(tmp_path, monkeypatch)
    _run(rnd, "S1")
    (rnd / "status").mkdir()
    _write_tsv(
        rnd / "status" / "n_a.tsv",
        ["dataset", "sample", "arm", "state", "reason"],
        [("ds", "S3", "all", "n/a", "not GEX")],
    )
    index, ctx = dos.evaluate(rnd, tmp_path / "out", derived=False)

    # the twostep_v2 arm has no directory: missing; a directory without run_complete.json: failed
    assert _status(index, "A10", "ds/S1", "combined_off") == ["ok"]
    assert _status(index, "A10", "ds/S1", "twostep_v2") == ["missing"]
    (rnd / "runs" / "ds" / "twostep_v2" / "S1").mkdir(parents=True)
    index, _ = dos.evaluate(rnd, tmp_path / "out", derived=False)
    assert _status(index, "A10", "ds/S1", "twostep_v2") == ["failed"]
    assert _status(index, "A10", "ds/S3", "combined_off") == ["n_a"]
    assert _status(index, "-", "newds") == ["not_started"]

    # S1 has one call and no molecule verdicts: C1 at most; the negative control is not certified
    reached, blockers, top = dos.rung_table(index, dos._tsv(dos.SPEC), "ds", "S1", ctx["calls"])
    assert top == "C0" or top == "C1"
    assert reached["C0"] and not reached["C2"]
    assert any(b.startswith("E07 missing") for b in blockers["C2"])
    assert not reached["N1"] and any(
        b.startswith("G01") for b in blockers["N1"]
    )  # positive_control not-configured
    assert dos.rung_table(index, dos._tsv(dos.SPEC), "ds", "S3", ctx["calls"])[2] == "n_a"


def test_c2_needs_the_molecule_summary_and_a_verdict_row_and_a_call(tmp_path, monkeypatch):
    rnd = _setup(tmp_path, monkeypatch)
    _run(rnd, "S1")
    _run(rnd, "S2", molecules="1")  # below the 3-molecule rule: no call, so C2 would be vacuous
    spec = dos._tsv(dos.SPEC)
    ev = rnd / "evidence_v2" / "ds__S1__combined_off__Epstein_Barr_virus"
    ev.mkdir(parents=True)
    (ev / "molecule_verdict_summary.tsv").write_text("scope\tverdict\tmolecules\n")
    index, ctx = dos.evaluate(rnd, tmp_path / "out", derived=False)
    # the summary exists but the roll-up verdicts.tsv has no row yet
    assert _status(index, "E07", "ds/S1") == ["ok"] and _status(index, "E09", "ds/S1") == [
        "missing"
    ]
    assert not dos.rung_table(index, spec, "ds", "S1", ctx["calls"])[0]["C2"]

    _write_tsv(
        rnd / "evidence_v2" / "verdicts.tsv", ["evidence_dir", "verdict"], [(ev.name, "viral_best")]
    )
    index, ctx = dos.evaluate(rnd, tmp_path / "out", derived=False)
    reached, _, top = dos.rung_table(index, spec, "ds", "S1", ctx["calls"])
    assert reached["C2"] and top == "C2"

    reached, blockers, top = dos.rung_table(index, spec, "ds", "S2", ctx["calls"])
    assert not reached["C2"] and "no combined_off call to read-validate" in blockers["C2"]


def test_stale_when_older_than_its_input(tmp_path, monkeypatch):
    rnd = _setup(tmp_path, monkeypatch)
    root = _run(rnd, "S1")
    rows = [
        (
            "E01",
            "E",
            "call",
            "",
            "{ev}/interpretation_flags.tsv",
            "exists",
            "",
            "{arm_run_complete}",
            "",
            "",
        )
    ]
    _write_tsv(dos.SPEC, SPEC_COLS, rows)
    ev = rnd / "evidence" / "ds__S1__combined_off__Epstein_Barr_virus"
    ev.mkdir(parents=True)
    flags = ev / "interpretation_flags.tsv"
    flags.write_text("flag\tstatus\n")
    os.utime(flags, (1, 1))  # the evidence predates the run
    os.utime(root / "run_complete.json", (1000, 1000))
    index, _ = dos.evaluate(rnd, tmp_path / "out", derived=False)
    assert _status(index, "E01", "ds/S1") == ["stale"]


def test_x223_style_inner_directory_is_found(tmp_path):
    root = tmp_path / "runs" / "d" / "combined_off" / "x223"
    (root / "LUM-SJ-x223" / "log").mkdir(parents=True)
    assert dos.inner_name(root, "x223") == "LUM-SJ-x223"
    assert dos.slug("Human papillomavirus 77") == "Human_papillomavirus_77"
