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


def test_manifest_omits_unset_positive_control_scope(tmp_path: Path) -> None:
    """Existing runs keep their fingerprint: an unset scope/key adds no manifest key."""
    absent = build_run_manifest(_args(tmp_path))
    unset = build_run_manifest(
        _args(tmp_path, positive_control_scope=None, positive_control_virus_key=None)
    )
    assert "positive_control_scope" not in unset["options"]
    assert "positive_control_virus_key" not in unset["options"]
    assert unset["run_fingerprint"] == absent["run_fingerprint"]


def test_positive_control_scope_changes_the_fingerprint_when_set(tmp_path: Path) -> None:
    base = build_run_manifest(_args(tmp_path))["run_fingerprint"]
    scoped = build_run_manifest(
        _args(
            tmp_path,
            positive_control_scope="exact_sequence",
            positive_control_virus_key="Torque teno virus",
        )
    )
    assert scoped["options"]["positive_control_scope"] == "exact_sequence"
    assert scoped["run_fingerprint"] != base
    other = build_run_manifest(
        _args(tmp_path, positive_control_scope="exact_sequence", positive_control_virus_key="X")
    )
    assert other["run_fingerprint"] != scoped["run_fingerprint"]


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


# Fingerprint of a plain run at aca8998 (before C-13), with __version__ pinned.
_PLAIN_RUN_FINGERPRINT = "3199343e8d011316ed404647f812115debcfd904f1f80146f6db197d4656c2c9"


def test_plain_run_fingerprint_is_unchanged_by_c13(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Runs that use none of the newly hashed inputs must keep resuming."""
    import viralscan.run_safety as rs

    monkeypatch.setattr(rs, "__version__", "0.0.test")
    monkeypatch.chdir(tmp_path)
    for name, text in [("R1.fastq", "r1\n"), ("R2.fastq", "r2\n"), ("idx.idx", "x\n")]:
        Path(name).write_text(text)
    Path("t2g.txt").write_text("a\tb\n")
    args = argparse.Namespace(
        sample1="R1.fastq", sample2="R2.fastq", index="idx.idx", transcripts="t2g.txt",
        gtf=None, fasta=None, f1=None, whitelist=None, called_cells_file=None,
        host_index=None, host_h5ad=None, cell_types=None, host_fasta=None, viral_fasta=None,
        positive_control_receipt=None, cores=4, multimap_method="equal",
        cell_calling="emptydrops", technology="10xv3", output="o", resume=False,
        overwrite=False, yes=False, verbose=False, quiet=False, ncbi_accession=None,
        ncbi_api_key="k", ncbi_email="e", strand=None, umap=True,
    )  # fmt: skip
    assert rs.build_run_manifest(args)["run_fingerprint"] == _PLAIN_RUN_FINGERPRINT


def test_input_paths_are_resolved_and_outside_the_fingerprint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Same bytes and option strings from two directories: paths differ, fingerprint not."""
    manifests = []
    for name in ("a", "b"):
        where = tmp_path / name
        where.mkdir()
        monkeypatch.chdir(where)
        manifests.append(build_run_manifest(_args(where, sample1="R1.fastq", sample2="R2.fastq")))
    first, second = manifests
    assert first["input_paths"]["sample1:0"] == str((tmp_path / "a" / "R1.fastq").resolve())
    assert first["input_paths"] != second["input_paths"]
    assert first["run_fingerprint"] == second["run_fingerprint"]


def test_resume_mismatch_names_the_differing_options(tmp_path: Path) -> None:
    """C-12: say which options changed (cores stays in the fingerprint)."""
    old = _args(tmp_path, cores=4)
    out = Path(old.output)
    prepare_output_directory(out, build_run_manifest(old), resume=False, overwrite=False, yes=False)
    (out / "partial.txt").write_text("partial")
    new = build_run_manifest(_args(tmp_path, cores=12))
    with pytest.raises(RunSafetyError, match=r"Options that differ: --cores\."):
        prepare_output_directory(out, new, resume=True, overwrite=False, yes=False)


def test_host_index_directory_and_extra_inputs_are_hashed(tmp_path: Path) -> None:
    """C-13: a STAR index directory and the extra file inputs enter the fingerprint."""
    idx = tmp_path / "star"
    idx.mkdir()
    (idx / "SA").write_text("one")
    extra = tmp_path / "host.h5ad"
    extra.write_text("h5")
    args = _args(tmp_path, host_index=str(idx), host_h5ad=str(extra))
    base = build_run_manifest(args)
    assert {"host_index:0", "host_h5ad:0"} <= set(base["input_fingerprints"])
    (idx / "SA").write_text("two")
    changed_dir = build_run_manifest(args)
    assert changed_dir["run_fingerprint"] != base["run_fingerprint"]
    assert changed_dir["reference_fingerprint"] != base["reference_fingerprint"]
    extra.write_text("h5-changed")
    assert build_run_manifest(args)["run_fingerprint"] != changed_dir["run_fingerprint"]
