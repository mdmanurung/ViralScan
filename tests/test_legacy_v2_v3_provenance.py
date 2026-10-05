"""Tests for the ignored legacy-v2/v3 provenance freezer."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from scripts import freeze_legacy_v2_v3_provenance

pytestmark = pytest.mark.research


def _executable(path: Path, output: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!/bin/sh\nprintf '%s\\n' '{output}'\n", encoding="utf-8")
    path.chmod(0o755)


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    v2_prefix = tmp_path / "v2"
    site_packages = v2_prefix / "lib/python3.12/site-packages"
    viralscan = site_packages / "viralscan"
    viralscan.mkdir(parents=True)
    (viralscan / "__init__.py").write_text("__version__ = '2.2.0'\n", encoding="utf-8")
    (viralscan / "menu.py").write_text("def main(): return 0\n", encoding="utf-8")
    (viralscan / "Snakefile").write_text("rule all:\n    input: []\n", encoding="utf-8")
    viral_metadata = site_packages / "viralscan-2.2.0.dist-info/METADATA"
    viral_metadata.parent.mkdir()
    viral_metadata.write_text(
        "Metadata-Version: 2.1\nName: ViralScan\nVersion: 2.2.0\n",
        encoding="utf-8",
    )

    kb_python = site_packages / "kb_python"
    kb_python.mkdir()
    (kb_python / "__init__.py").write_text("", encoding="utf-8")
    kb_metadata = site_packages / "kb_python-0.29.5.dist-info/METADATA"
    kb_metadata.parent.mkdir()
    kb_metadata.write_text(
        "Metadata-Version: 2.1\nName: kb_python\nVersion: 0.29.5\n",
        encoding="utf-8",
    )
    _executable(v2_prefix / "bin/python", "Python 3.12.11")
    _executable(
        kb_python / "bins/linux/kallisto/kallisto",
        "kallisto, version 0.51.1",
    )
    _executable(
        kb_python / "bins/linux/bustools/bustools",
        "bustools, version 0.45.1",
    )

    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "index_serratus.idx").write_bytes(b"index")
    (reference / "t2g_serratus.txt").write_text("tx\tgene\n", encoding="utf-8")
    (reference / "transcriptome.fa").write_text(">tx\nACGT\n", encoding="utf-8")

    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / ".gitignore").write_text("benchmark_runs/\n", encoding="utf-8")
    (repo / "tracked.txt").write_text("before\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    subprocess.run(["git", "-C", str(repo), "add", ".gitignore", "tracked.txt"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "fixture"], check=True)
    (repo / "tracked.txt").write_text("after\n", encoding="utf-8")
    return v2_prefix, reference, repo


def test_freeze_emits_complete_json_and_tsv_records(tmp_path: Path) -> None:
    v2_prefix, reference, repo = _fixture(tmp_path)
    output_dir = repo / "benchmark_runs/provenance"

    assert (
        freeze_legacy_v2_v3_provenance.main(
            [
                "--v2-prefix",
                str(v2_prefix),
                "--reference-dir",
                str(reference),
                "--repo-root",
                str(repo),
                "--output-dir",
                str(output_dir),
            ]
        )
        == 0
    )

    payload = json.loads((output_dir / "provenance.json").read_text(encoding="utf-8"))
    assert payload["schema_version"] == "1.0.0"
    assert payload["v2"]["packages"] == {"ViralScan": "2.2.0", "kb-python": "0.29.5"}
    assert payload["v2"]["tools"] == {"kallisto": "0.51.1", "bustools": "0.45.1"}
    assert len(payload["v2"]["viralscan_code_tree_sha256"]) == 64
    assert set(payload["reference"]["files"]) == {
        "index_serratus.idx",
        "t2g_serratus.txt",
        "transcriptome.fa",
    }
    assert len(payload["repository"]["head"]) == 40
    assert len(payload["repository"]["diff_sha256"]) == 64
    assert len(payload["files_tsv_sha256"]) == 64
    assert (output_dir / "provenance_files.tsv").is_file()


def test_freeze_fails_closed_before_writing_when_a_required_path_is_missing(
    tmp_path: Path,
) -> None:
    v2_prefix, reference, repo = _fixture(tmp_path)
    (reference / "transcriptome.fa").unlink()
    output_dir = repo / "benchmark_runs/provenance"

    with pytest.raises(SystemExit, match="missing required reference file transcriptome.fa"):
        freeze_legacy_v2_v3_provenance.main(
            [
                "--v2-prefix",
                str(v2_prefix),
                "--reference-dir",
                str(reference),
                "--repo-root",
                str(repo),
                "--output-dir",
                str(output_dir),
            ]
        )

    assert not output_dir.exists()


def test_repository_diff_hash_includes_untracked_file_content(tmp_path: Path) -> None:
    v2_prefix, reference, repo = _fixture(tmp_path)
    before, _ = freeze_legacy_v2_v3_provenance.collect_provenance(
        v2_prefix=v2_prefix,
        reference_dir=reference,
        repo_root=repo,
    )

    (repo / "new_source.py").write_text("value = 1\n", encoding="utf-8")
    after, _ = freeze_legacy_v2_v3_provenance.collect_provenance(
        v2_prefix=v2_prefix,
        reference_dir=reference,
        repo_root=repo,
    )

    assert before["repository"]["diff_sha256"] != after["repository"]["diff_sha256"]


def test_site_packages_compatibility_symlink_is_not_treated_as_ambiguity(
    tmp_path: Path,
) -> None:
    v2_prefix, reference, repo = _fixture(tmp_path)
    (v2_prefix / "lib/python3.1").symlink_to("python3.12")

    payload, _ = freeze_legacy_v2_v3_provenance.collect_provenance(
        v2_prefix=v2_prefix,
        reference_dir=reference,
        repo_root=repo,
    )

    assert payload["v2"]["packages"]["ViralScan"] == "2.2.0"
