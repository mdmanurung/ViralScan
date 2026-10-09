"""scripts/dsr02_enumerate_calls.py keeps >= 3-molecule calls and every anellovirus call."""

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "dsr02_enumerate_calls", Path(__file__).parent.parent / "scripts" / "dsr02_enumerate_calls.py"
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def _summary(tmp_path, dataset, rows):
    d = tmp_path / "runs" / dataset / "combined_off" / "S1" / "S1" / "results"
    d.mkdir(parents=True)
    body = "virus_name\tviral_molecules_total_est\n" + "".join(f"{n}\t{m}\n" for n, m in rows)
    (d / "viral_summary.tsv").write_text(body)


def test_threshold_anellovirus_and_expected(tmp_path):
    _summary(
        tmp_path,
        "ebv",
        [
            ("Epstein-Barr virus", 749314),
            ("Human betaherpesvirus 6B", 5),
            ("Human betaherpesvirus 6A", 2),
            ("Torque teno virus 1", 1),
            ("Torque teno virus 2", 0),
        ],
    )
    got = {r["virus"]: r for r in mod.calls(tmp_path)}
    assert set(got) == {"Epstein-Barr virus", "Human betaherpesvirus 6B", "Torque teno virus 1"}
    assert got["Epstein-Barr virus"]["expected"] == "true"
    assert got["Human betaherpesvirus 6B"]["expected"] == "false"
    assert got["Torque teno virus 1"]["anellovirus"] == "true"
    assert got["Epstein-Barr virus"]["sample"] == "S1"


def test_twostep_v2_replaces_v1_only_once_finished(tmp_path):
    def arm(name, rows):
        d = tmp_path / "runs" / "ebv" / name / "S1" / "S1" / "results"
        d.mkdir(parents=True)
        (d / "viral_summary.tsv").write_text(
            "virus_name\tviral_molecules_total_est\n" + "".join(f"{n}\t{m}\n" for n, m in rows)
        )

    arm("twostep", [("Epstein-Barr virus", 5)])
    arm("twostep_v2", [("Human papillomavirus 77", 9)])  # stale copy until the marker exists
    assert {(r["arm"], r["virus"]) for r in mod.calls(tmp_path)} == {
        ("twostep", "Epstein-Barr virus")
    }
    (tmp_path / "runs" / "ebv" / "twostep_v2" / "S1" / "run_complete.json").write_text("{}")
    assert {(r["arm"], r["virus"]) for r in mod.calls(tmp_path)} == {
        ("twostep_v2", "Human papillomavirus 77")
    }
