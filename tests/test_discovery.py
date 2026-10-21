"""DISC-02: `viralscan discover` refuses unsafe inputs and never touches the source run."""

import hashlib
from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan import discovery
from viralscan.discovery import DiscoveryError


def _run(tmp_path: Path, samples=("S1",), complete: bool = True) -> Path:
    root = tmp_path / "run"
    for sample in samples:
        (root / sample / "results").mkdir(parents=True)
        (root / sample / "config.yaml").write_text("output: x\n")
        (root / sample / "results" / "viral_summary.tsv").write_text("virus\tdetected\nEBV\tTrue\n")
    (root / "run_manifest.json").write_text("{}\n")
    if complete:
        (root / "run_complete.json").write_text("{}\n")
    return root


def _db(tmp_path: Path, empty: str | None = None) -> Path:
    db = tmp_path / "k2db"
    db.mkdir()
    for name in discovery.KRAKEN2_DB_FILES:
        (db / name).write_bytes(b"" if name == empty else b"x")
    return db


def _inventory(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_resolves_single_sample(tmp_path):
    root = _run(tmp_path)
    assert discovery.resolve_sample_dir(root, None) == root / "S1"


@pytest.mark.parametrize(
    "samples,sample,complete,match",
    [
        (("S1",), None, False, "not a completed"),
        (("S1", "S2"), None, True, "several samples"),
        (("S1",), "S9", True, "not in run"),
    ],
)
def test_refuses_bad_run(tmp_path, samples, sample, complete, match):
    root = _run(tmp_path, samples, complete)
    with pytest.raises(DiscoveryError, match=match):
        discovery.resolve_sample_dir(root, sample)


def test_refuses_output_inside_run_or_non_empty(tmp_path):
    root = _run(tmp_path)
    for inside in (root, root / "S1" / "discovery"):
        with pytest.raises(DiscoveryError, match="outside --run-dir"):
            discovery.check_output_dir(inside, root)
    busy = tmp_path / "busy"
    busy.mkdir()
    (busy / "old.tsv").write_text("x")
    with pytest.raises(DiscoveryError, match="empty or absent"):
        discovery.check_output_dir(busy, root)
    # A leftover staging directory alone does not count as results.
    stale = tmp_path / "stale"
    (stale / discovery.PARTIAL).mkdir(parents=True)
    discovery.check_output_dir(stale, root)


def test_refuses_incomplete_db(tmp_path):
    with pytest.raises(DiscoveryError, match="does not exist"):
        discovery.validate_kraken2_db(tmp_path / "absent")
    with pytest.raises(DiscoveryError, match="hash.k2d"):
        discovery.validate_kraken2_db(_db(tmp_path, empty="hash.k2d"))


def test_refuses_missing_tools():
    with patch("viralscan.discovery.shutil.which", return_value=None):
        with pytest.raises(DiscoveryError, match="kraken2, k2"):
            discovery.check_tools()


def test_stage_and_publish(tmp_path):
    out = tmp_path / "out"
    partial = discovery.stage(out)
    (partial / "leftover.txt").write_text("from an interrupted run")
    partial = discovery.stage(out)  # restart discards the leftover
    (partial / "qc.json").write_text("{}")
    discovery.publish(out)
    assert sorted(p.name for p in out.iterdir()) == ["qc.json"]


def test_cli_leaves_source_run_untouched(tmp_path):
    """G2: the source inventory is byte-identical and the skeleton publishes nothing."""
    from viralscan.menu import main

    root = _run(tmp_path)
    before = _inventory(root)
    out = tmp_path / "disc"
    argv = [
        "viralscan",
        "discover",
        "--run-dir",
        str(root),
        "-o",
        str(out),
        "--kraken2-db",
        str(_db(tmp_path)),
    ]
    with (
        patch("sys.argv", argv),
        patch("viralscan.discovery.shutil.which", return_value="/bin/true"),
        pytest.raises(SystemExit) as exit_info,
    ):
        main()
    assert exit_info.value.code == 1
    assert _inventory(root) == before
    assert not out.exists()


def test_run_parser_gains_no_discovery_options():
    """G1: discovery adds nothing to the primary run's options (and so to its fingerprint)."""
    from viralscan.menu import build_parser

    with patch("sys.argv", ["viralscan"]):
        args = build_parser().parse_args([])
    assert not {"kraken2_db", "sample", "threads"} & set(vars(args))
