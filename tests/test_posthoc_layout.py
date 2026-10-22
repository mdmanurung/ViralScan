"""CLI-UNIFY-02: post-hoc verbs write under <run-dir>/posthoc/; layout_version 1 and 2 both read."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.test_rerun_multimap import TestRerunLeavesSourceUntouched as _Rerun
from viralscan import menu
from viralscan.outputs import LAYOUT_VERSION, V3Outputs, is_posthoc, posthoc_dir
from viralscan.validation import validate_run


def _tree_hashes(root: Path, *, skip_posthoc: bool = True) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file() and not (skip_posthoc and is_posthoc(p, root))
    }


def _rerun_args(source: Path, method: str, output: str | None = None) -> argparse.Namespace:
    return argparse.Namespace(
        run_dir=str(source),
        output=output,
        multimap_method=method,
        multimap_em_max_iter=None,
        multimap_em_tol=None,
        cores=1,
        verbose=False,
        quiet=True,
    )


def _rerun(args: argparse.Namespace) -> None:
    with (
        patch.object(menu, "_check_required_tools"),
        patch.object(menu, "_check_cell_caller_tools"),
        patch.object(menu.subprocess, "run"),
    ):
        menu._run_rerun_multimap(args)


def test_layout_version_is_2_and_posthoc_path(tmp_path: Path) -> None:
    assert LAYOUT_VERSION == 2
    assert V3Outputs(tmp_path).posthoc_dir == tmp_path / "posthoc"
    assert posthoc_dir(tmp_path, "evidence", "EBV") == tmp_path / "posthoc" / "evidence" / "EBV"


def test_rerun_multimap_defaults_into_posthoc_and_leaves_run_untouched(tmp_path: Path) -> None:
    source = _Rerun._source_run(tmp_path)
    before = _tree_hashes(source)
    _rerun(_rerun_args(source, "host-conservative"))
    out = source / "posthoc" / "multimap_host_conservative"
    assert (out / "SAMPLE" / "config.yaml").is_file()
    assert _tree_hashes(source) == before, "post-hoc verb changed the run's own files"

    # a second rerun does not copy the first one's posthoc/ into itself
    _rerun(_rerun_args(source, "equal"))
    assert not (source / "posthoc" / "multimap_equal" / "posthoc").exists()
    assert _tree_hashes(source) == before


def test_rerun_multimap_still_refuses_non_empty_output(tmp_path: Path) -> None:
    source = _Rerun._source_run(tmp_path)
    _rerun(_rerun_args(source, "equal"))
    with pytest.raises(SystemExit):
        _rerun(_rerun_args(source, "equal"))


def test_rerun_multimap_output_overrides_and_must_stay_out_of_run(tmp_path: Path) -> None:
    source = _Rerun._source_run(tmp_path)
    _rerun(_rerun_args(source, "equal", str(tmp_path / "elsewhere")))
    assert (tmp_path / "elsewhere" / "SAMPLE" / "config.yaml").is_file()
    assert not (source / "posthoc").exists()
    with pytest.raises(SystemExit):  # inside the run but not under posthoc/
        _rerun(_rerun_args(source, "equal", str(source / "other")))


@pytest.mark.parametrize(
    ("version", "unsupported"), [(1, False), (2, False), (None, False), (3, True)]
)
def test_validate_run_layout_versions(tmp_path: Path, version, unsupported) -> None:
    manifest = {"schema_version": "3.0.0"}
    if version is not None:
        manifest["layout_version"] = version
    (tmp_path / "run_manifest.json").write_text(json.dumps(manifest))
    codes = {i["code"] for i in validate_run(tmp_path, verify_inputs=False)["issues"]}
    assert ("unsupported_layout_version" in codes) is unsupported


def test_validate_run_ignores_posthoc_contents(tmp_path: Path) -> None:
    (tmp_path / "run_manifest.json").write_text(json.dumps({"schema_version": "3.0.0"}))
    bad = tmp_path / "posthoc" / "evidence" / "X"
    bad.mkdir(parents=True)
    (bad / "evidence_manifest.json").write_text("{not json")
    (bad / "adata_multimap.h5ad").write_text("junk")
    report = validate_run(tmp_path, verify_inputs=False)
    assert report["h5ad_files"] == []
    assert all("posthoc" not in i["path"] for i in report["issues"])


def test_rerun_programs_writes_under_posthoc_and_leaves_run_untouched(tmp_path: Path) -> None:
    from viralscan.runconfig import RunConfig
    from viralscan.scripts import gene_programs as gp

    sample = tmp_path / "SAMPLE"
    (sample / "log").mkdir(parents=True)
    (sample / "log" / "detection.done").touch()
    h5ad = sample / "kb-python" / "counts_unfiltered" / "adata_multimap.h5ad"
    h5ad.parent.mkdir(parents=True)
    h5ad.write_bytes(b"x")
    RunConfig(output=str(sample) + "/", cell_calling="none").to_yaml(sample / "config.yaml")
    before = _tree_hashes(sample)
    args = argparse.Namespace(
        run_dir=str(sample),
        output=None,
        programme_min_breadth=None,
        programme_latent_min_breadth=None,
        programme_min_umi=None,
        verbose=False,
        quiet=True,
    )
    with patch.object(gp, "main") as main:
        menu._run_rerun_programs(args)
    assert main.call_args.args[2] == ""
    assert (
        Path(main.call_args.kwargs["results_root"]) == (sample / "posthoc" / "programs").resolve()
    )
    assert _tree_hashes(sample) == before


@pytest.mark.parametrize(
    "samples",
    [["/d/posthoc_R1.fastq.gz"], ["/d/S1_R1.fastq.gz", "/d/S1_L002_R1.fastq.gz"]],
    ids=["reserved", "duplicate"],
)
def test_run_rejects_reserved_or_duplicate_sample_id(samples: list[str]) -> None:
    with pytest.raises(SystemExit):
        menu._check_sample_ids(samples)
    menu._check_sample_ids(["/d/S1_R1.fastq.gz", "/d/S2_R1.fastq.gz"])
