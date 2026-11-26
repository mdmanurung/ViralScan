"""TOOLCHECK-01: one PATH preflight and one STARsolo barcode-argument builder."""

from pathlib import Path

import pytest

from viralscan import star_hostvirus, star_run, validation
from viralscan.scripts import host_filter


def test_missing_tools_keeps_order_and_reports_only_absent(monkeypatch) -> None:
    monkeypatch.setattr(
        validation.shutil, "which", lambda tool: None if tool in ("STAR", "samtools") else "/x"
    )
    assert validation.missing_tools(("samtools", "minimap2", "STAR")) == ["samtools", "STAR"]
    assert validation.missing_tools(()) == []


@pytest.mark.parametrize("whitelist", [None, "", "/x/wl.txt", Path("/x/wl.txt")])
def test_builder_whitelist_forms(whitelist) -> None:
    args = star_run.starsolo_barcode_args("10xv3", whitelist)
    expected = str(whitelist) if whitelist else "None"
    assert args[-2:] == ["--soloCBwhitelist", expected]
    assert all(isinstance(a, str) for a in args)


def test_builder_argument_order_is_pinned() -> None:
    assert star_run.starsolo_barcode_args("10xv2", None) == [
        "--soloType", "CB_UMI_Simple",
        "--soloCBstart", "1", "--soloCBlen", "16",
        "--soloUMIstart", "17", "--soloUMIlen", "10",
        "--soloBarcodeReadLength", "0",
        "--soloCBwhitelist", "None",
    ]  # fmt: skip


def test_host_filter_and_hostvirus_use_the_one_builder() -> None:
    assert host_filter.starsolo_barcode_args is star_run.starsolo_barcode_args
    cmd = star_hostvirus.align_cmd(
        "STAR", "/i", "r1.fq", "r2.fq", "10xv3", "p/", 2, whitelist="/x/wl.txt"
    )
    block = star_run.starsolo_barcode_args("10xv3", "/x/wl.txt")
    start = cmd.index("--soloType")
    assert cmd[start : start + len(block)] == block
