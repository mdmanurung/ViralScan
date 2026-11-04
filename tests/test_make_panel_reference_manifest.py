"""scripts/make_panel_reference_manifest.py (PLAN REF-07 / REF-04): manifest + audit from a build dir."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import compare_reference_reproducibility as cmp  # noqa: E402
import make_panel_reference_manifest as make  # noqa: E402


def _build_dir(root: Path, extra: str = "") -> Path:
    build = root / "build"
    build.mkdir()
    viral = ">V1.1\nACGTCAGTACGTACGGTCAT\n>V2.1\nTTGACCAGTTGACAGGCATC\n" + extra
    (build / "viral.fa").write_text(viral)
    (build / "combined.fa").write_text(">ENST1\n" + "ACGGT" * 12 + "\n" + viral)
    for name in ("viral_panel.gtf", "combined.gtf", "cdna.fa", "panel.t2g", "panel.idx"):
        (build / name).write_text(f"{name}\n")
    return build


def _run(build: Path, out: Path, monkeypatch) -> dict:
    monkeypatch.setattr(
        sys,
        "argv",
        ["x", "--build-dir", str(build), "--out-dir", str(out), "--cache-dir", str(build)],
    )
    make.main()
    return json.loads((out / "reference_manifest.json").read_text())


def test_manifest_rows_and_audit_compare_clean(tmp_path, monkeypatch):
    build = _build_dir(tmp_path)
    first = _run(build, tmp_path / "a", monkeypatch)
    _run(build, tmp_path / "b", monkeypatch)
    assert [row["accession_version"] for row in first["sequences"]] == ["V1.1", "V2.1"]
    assert first["sequences"][0]["source_licence"]["status"] == "unreviewed"
    assert (
        first["sequences"][0]["retrieval_status"] == "unknown_legacy_cache"
    )  # no receipt, no date
    assert "panel.idx" not in json.dumps(first["content_files"])  # binary index is not claimed
    assert first["build_receipt"]["index"]["sha256"]
    audits = [
        json.loads((tmp_path / n / "reference_reproducibility.json").read_text()) for n in "ab"
    ]
    assert cmp.compare_reproducibility(*audits) == []


def test_a_changed_build_file_is_reported_and_an_existing_out_dir_is_refused(tmp_path, monkeypatch):
    build = _build_dir(tmp_path)
    _run(build, tmp_path / "a", monkeypatch)
    (build / "panel.t2g").write_text("changed\n")
    _run(build, tmp_path / "b", monkeypatch)
    audits = [
        json.loads((tmp_path / n / "reference_reproducibility.json").read_text()) for n in "ab"
    ]
    assert any(
        line.startswith("content_files.t2g") for line in cmp.compare_reproducibility(*audits)
    )
    with pytest.raises(FileExistsError):
        _run(build, tmp_path / "a", monkeypatch)
