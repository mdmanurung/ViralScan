"""Parity compare and report logic with fake per-environment results (REL-09, REL-10)."""

import json

import pytest

from scripts.release_parity import compare, explained_reason, main

H1, H2 = "1" * 64, "2" * 64


def _result(name, files=None, molecules=None, ok=True):
    return {
        "name": name,
        "kind": "conda",
        "digest": f"sha256:{name}",
        "versions": {"kallisto": "0.51.1"},
        "commands": [{"command": "viralscan run", "returncode": 0}],
        "validate_run_ok": ok,
        "molecules": molecules or {"viral_molecules_total_est": 3.0, "viruses_detected": ["V"]},
        "files": files or {"s/results/viral_summary.tsv": H1, "s/config.yaml": H1},
    }


def test_identical_environments_pass_and_report_lists_inputs():
    report = compare([_result("conda"), _result("docker"), _result("apptainer")])
    assert report["ok"] and report["problems"] == []
    assert report["reference"] == "conda"
    assert [e["name"] for e in report["environments"]] == ["conda", "docker", "apptainer"]
    assert report["environments"][1]["digest"] == "sha256:docker"
    assert report["deterministic_files"] == {"s/results/viral_summary.tsv": H1}


def test_explained_file_may_differ_but_is_listed_with_its_reason():
    files = {"s/results/viral_summary.tsv": H1, "s/config.yaml": H2}
    report = compare([_result("conda"), _result("docker", files)])
    assert report["ok"]
    entry = report["explained_nondeterministic_files"]["s/config.yaml"]
    assert entry["sha256"] == {"conda": H1, "docker": H2} and entry["reason"]


def test_unexplained_hash_difference_fails():
    files = {"s/results/viral_summary.tsv": H2, "s/config.yaml": H1}
    report = compare([_result("conda"), _result("docker", files)])
    assert not report["ok"]
    assert any("viral_summary.tsv" in p and "hashes differ" in p for p in report["problems"])


def test_file_missing_from_one_environment_fails():
    report = compare([_result("conda"), _result("docker", {"s/config.yaml": H1})])
    assert any("missing in ['docker']" in p for p in report["problems"])


def test_molecule_difference_and_failed_validation_fail():
    other = {"viral_molecules_total_est": 4.0, "viruses_detected": ["V"]}
    report = compare([_result("conda"), _result("docker", molecules=other, ok=False)])
    assert any("molecules" in p for p in report["problems"])
    assert any("validate-run" in p for p in report["problems"])


def test_one_environment_or_duplicate_names_are_refused():
    with pytest.raises(ValueError, match="at least two"):
        compare([_result("conda")])
    with pytest.raises(ValueError, match="unique"):
        compare([_result("conda"), _result("conda")])


@pytest.mark.parametrize(
    "path,explained",
    [
        ("sample/results/viral_summary.tsv", False),
        ("sample/results/count_audit.tsv", False),
        ("sample/kb-python/output.resolved.sorted.bus", False),
        ("sample/config.yaml", True),
        ("README.txt", True),
        ("sample/results/chemistry_sanity.json", True),
        ("sample/results/artifact_inventory.json", True),
        ("run_manifest.json", True),
        ("sample/log/kb.done", True),
        ("sample/kb-python/counts_unfiltered/adata_multimap.h5ad", True),
        ("sample/.snakemake/metadata/x", True),
        ("sample/results/my_config.yaml", False),
    ],
)
def test_explained_list_is_exact(path, explained):
    assert (explained_reason(path) is not None) is explained


def test_compare_cli_writes_report_and_sets_exit_code(tmp_path, monkeypatch):
    for name, files in (("a", None), ("b", {"s/results/viral_summary.tsv": H2})):
        (tmp_path / f"result_{name}.json").write_text(json.dumps(_result(name, files)))
    out = tmp_path / "analysis" / "release_parity" / "parity_report.json"
    argv = [
        "release_parity",
        "compare",
        str(tmp_path / "result_a.json"),
        str(tmp_path / "result_b.json"),
        "--output",
        str(out),
    ]
    monkeypatch.setattr("sys.argv", argv)
    assert main() == 1
    assert json.loads(out.read_text())["ok"] is False
