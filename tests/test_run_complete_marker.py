"""SW-06: run_complete.json completion marker."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tests.test_validation import _valid_run
from viralscan.run_safety import (
    RUN_COMPLETE,
    RUN_MANIFEST,
    build_run_manifest,
    clear_run_complete,
    prepare_output_directory,
    restamp_run_complete,
    write_run_complete,
)
from viralscan.validation import packaged_schema_resource, validate_json_schema, validate_run


def _codes(report: dict) -> set[str]:
    return {i["code"] for i in report["issues"]}


def _with_artifact(run: Path) -> Path:
    sample = run / "sample"
    (sample / "results").mkdir()
    (sample / "results" / "viral_summary.tsv").write_text("virus\nA\n")
    write_run_complete(run)
    return sample


def test_marker_records_fingerprint_and_artifact_hashes(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    _with_artifact(run)
    marker = json.loads((run / RUN_COMPLETE).read_text())
    manifest = json.loads((run / RUN_MANIFEST).read_text())
    assert marker["run_fingerprint"] == manifest["run_fingerprint"]
    assert marker["samples"] == ["sample"]
    assert "sample/results/viral_summary.tsv" in marker["artifacts"]
    assert "sample/kb-python/counts_unfiltered/adata_multimap.h5ad" in marker["artifacts"]
    assert not validate_json_schema(marker, packaged_schema_resource("run_complete.schema.json"))
    assert validate_run(run)["ok"]


def test_new_run_missing_marker_is_error(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    clear_run_complete(run)
    assert "missing_completion_marker" in _codes(validate_run(run))


def test_old_run_missing_marker_is_warning_only(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    clear_run_complete(run)
    manifest = json.loads((run / RUN_MANIFEST).read_text())
    del manifest["completion_marker"]
    (run / RUN_MANIFEST).write_text(json.dumps(manifest))
    report = validate_run(run)
    assert report["ok"]
    assert [i["level"] for i in report["issues"]] == ["warning"]


def test_tampered_artifact_and_wrong_fingerprint_are_errors(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    sample = _with_artifact(run)
    (sample / "results" / "viral_summary.tsv").write_text("tampered\n")
    assert "completion_artifact_mismatch" in _codes(validate_run(run))
    marker = json.loads((run / RUN_COMPLETE).read_text())
    marker["run_fingerprint"] = "0" * 64
    (run / RUN_COMPLETE).write_text(json.dumps(marker))
    assert "completion_fingerprint_mismatch" in _codes(validate_run(run))


def test_prepare_resume_removes_marker(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    manifest = json.loads((run / RUN_MANIFEST).read_text())
    assert (run / RUN_COMPLETE).is_file()
    prepare_output_directory(run, manifest, resume=True, overwrite=False, yes=False)
    assert not (run / RUN_COMPLETE).exists()


def test_restamp_only_when_marker_existed(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    sample = _with_artifact(run)
    (sample / "results" / "viral_summary.tsv").write_text("changed\n")
    assert restamp_run_complete(run)
    assert validate_run(run)["ok"]
    clear_run_complete(run)
    assert not restamp_run_complete(run)
    assert not (run / RUN_COMPLETE).exists()


def test_new_manifest_declares_marker() -> None:
    args = argparse.Namespace(multimap_method="equal", cell_calling="none", called_cells_file=None)
    assert build_run_manifest(args)["completion_marker"] is True


def test_new_manifest_stamps_provisional_defaults_outside_the_fingerprint() -> None:
    """DEF-09: dev builds say their defaults are provisional, without breaking --resume."""
    args = argparse.Namespace(multimap_method="equal", cell_calling="none", called_cells_file=None)
    manifest = build_run_manifest(args)
    assert manifest["defaults_status"] == "provisional"
    hashed = {
        k: v
        for k, v in manifest.items()
        if k not in {"run_fingerprint", "completion_marker", "defaults_status", "software_identity"}
    }
    canonical = json.dumps(hashed, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(canonical).hexdigest() == manifest["run_fingerprint"]
