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


def test_manifest_records_strand_only_when_set(tmp_path: Path) -> None:
    assert build_run_manifest(_args(tmp_path, strand="reverse"))["options"]["strand"] == "reverse"
    assert "strand" not in build_run_manifest(_args(tmp_path, strand=None))["options"]


def test_manifest_omits_unset_star_params(tmp_path: Path) -> None:
    """Manifests written before --host-filter-star-params still resume when it is unset."""
    unset = build_run_manifest(_args(tmp_path, host_filter_star_params=None))
    assert "host_filter_star_params" not in unset["options"]
    assert "read_filter" not in build_run_manifest(_args(tmp_path, read_filter=None))["options"]
    chosen = build_run_manifest(_args(tmp_path, host_filter_star_params="star-default"))
    assert chosen["options"]["host_filter_star_params"] == "star-default"


def test_resume_old_manifest_matches_unset_strand(tmp_path: Path) -> None:
    """A manifest written before --strand existed resumes when strand is unset."""
    old_args = _args(tmp_path)  # no strand attribute at all, as before this change
    out = Path(old_args.output)
    prepare_output_directory(
        out, build_run_manifest(old_args), resume=False, overwrite=False, yes=False
    )
    (out / "partial.txt").write_text("partial")
    new = build_run_manifest(_args(tmp_path, strand=None))
    assert prepare_output_directory(out, new, resume=True, overwrite=False, yes=False) == "resume"


def test_resume_old_manifest_refuses_explicit_strand(tmp_path: Path) -> None:
    """Old counts were made with kb's default strand: never reuse them for --strand."""
    old_args = _args(tmp_path)
    out = Path(old_args.output)
    prepare_output_directory(
        out, build_run_manifest(old_args), resume=False, overwrite=False, yes=False
    )
    (out / "partial.txt").write_text("partial")
    new = build_run_manifest(_args(tmp_path, strand="reverse"))
    with pytest.raises(RunSafetyError, match="--strand differs"):
        prepare_output_directory(out, new, resume=True, overwrite=False, yes=False)


def test_resume_refuses_changed_strand_both_set(tmp_path: Path) -> None:
    args = _args(tmp_path, strand="reverse")
    out = Path(args.output)
    prepare_output_directory(
        out, build_run_manifest(args), resume=False, overwrite=False, yes=False
    )
    (out / "partial.txt").write_text("partial")
    with pytest.raises(RunSafetyError, match="--strand differs"):
        prepare_output_directory(
            out,
            build_run_manifest(_args(tmp_path, strand="unstranded")),
            resume=True,
            overwrite=False,
            yes=False,
        )
    same = build_run_manifest(_args(tmp_path, strand="reverse"))
    assert prepare_output_directory(out, same, resume=True, overwrite=False, yes=False) == "resume"


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


def test_force_technology_does_not_change_the_fingerprint(tmp_path: Path) -> None:
    """The override bypasses a check; the resolved -x is what is fingerprinted."""
    plain = build_run_manifest(_args(tmp_path, technology="10xv3"))
    forced = build_run_manifest(_args(tmp_path, technology="10xv3", force_technology=True))
    assert plain["run_fingerprint"] == forced["run_fingerprint"]
    assert plain["options"]["technology"] == "10xv3"
