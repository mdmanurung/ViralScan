from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from viralscan.run_safety import RunSafetyError, build_run_manifest, prepare_output_directory


def _args(tmp_path: Path, **changes):
    values = {
        "sample1": str(tmp_path / "R1.fastq"),
        "sample2": str(tmp_path / "R2.fastq"),
        "output": str(tmp_path / "out"),
        "multimap_method": "host-conservative",
        "cell_calling": "emptydrops",
        "called_cells_file": None,
        "resume": False,
        "overwrite": False,
        "yes": False,
        "verbose": False,
        "quiet": False,
    }
    values.update(changes)
    Path(values["sample1"]).write_text("r1\n")
    Path(values["sample2"]).write_text("r2\n")
    return argparse.Namespace(**values)


def test_new_output_writes_atomic_manifest(tmp_path: Path) -> None:
    args = _args(tmp_path)
    manifest = build_run_manifest(args)
    out = Path(args.output)
    assert (
        prepare_output_directory(out, manifest, resume=False, overwrite=False, yes=False) == "new"
    )
    assert (
        json.loads((out / "run_manifest.json").read_text())["run_fingerprint"]
        == manifest["run_fingerprint"]
    )
    assert not (out / ".run_manifest.json.tmp").exists()


def test_nonempty_output_refused_without_mode_even_with_yes(tmp_path: Path) -> None:
    args = _args(tmp_path)
    out = Path(args.output)
    out.mkdir()
    marker = out / "keep.txt"
    marker.write_text("keep")
    with pytest.raises(RunSafetyError, match="non-empty"):
        prepare_output_directory(
            out, build_run_manifest(args), resume=False, overwrite=False, yes=True
        )
    assert marker.read_text() == "keep"


def test_resume_requires_identical_fingerprint(tmp_path: Path) -> None:
    args = _args(tmp_path)
    out = Path(args.output)
    manifest = build_run_manifest(args)
    prepare_output_directory(out, manifest, resume=False, overwrite=False, yes=False)
    (out / "partial.txt").write_text("partial")
    assert (
        prepare_output_directory(out, manifest, resume=True, overwrite=False, yes=False) == "resume"
    )

    changed = build_run_manifest(_args(tmp_path, multimap_method="equal"))
    with pytest.raises(RunSafetyError, match="fingerprint"):
        prepare_output_directory(out, changed, resume=True, overwrite=False, yes=False)


def test_overwrite_is_explicit_and_scoped(tmp_path: Path) -> None:
    args = _args(tmp_path)
    out = Path(args.output)
    out.mkdir()
    (out / "old").write_text("old")
    outside = tmp_path / "outside"
    outside.write_text("safe")
    manifest = build_run_manifest(args)
    assert (
        prepare_output_directory(out, manifest, resume=False, overwrite=True, yes=True)
        == "overwrite"
    )
    assert not (out / "old").exists()
    assert outside.read_text() == "safe"


def test_overwrite_cancel_is_non_destructive(tmp_path: Path) -> None:
    args = _args(tmp_path)
    out = Path(args.output)
    out.mkdir()
    marker = out / "keep"
    marker.write_text("safe")
    with pytest.raises(RunSafetyError, match="cancelled"):
        prepare_output_directory(
            out,
            build_run_manifest(args),
            resume=False,
            overwrite=True,
            yes=False,
            confirm=lambda _prompt: "no",
        )
    assert marker.read_text() == "safe"
