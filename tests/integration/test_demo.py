"""``viralscan demo`` end to end (PLAN DEMO-01): the real run, ``check run`` and the pinned digests.

Run with::

    PYTHONPATH=src python -m pytest -m integration tests/integration/test_demo.py

Needs ``kallisto``, ``bustools``, ``kb`` and ``snakemake`` on PATH; skipped otherwise.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from viralscan.evidence import have_tools

from ._golden import load_golden

pytestmark = pytest.mark.integration


def _demo(*args: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "viralscan", "demo", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )


def test_demo_passes_end_to_end(tmp_path) -> None:
    missing = have_tools(["kallisto", "bustools", "kb", "snakemake"])
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

    done = _demo("--keep", cwd=tmp_path)  # default output: ./viralscan_demo
    assert done.returncode == 0, done.stdout + done.stderr
    out = tmp_path / "viralscan_demo"
    for name in load_golden():
        assert f"PASS      {name}" in done.stdout
    assert "PASS      check run" in done.stdout
    report = json.loads((out / "demo_report.json").read_text())
    assert report["ok"] is True and report["check_run"]["ok"] is True
    assert set(report["files"].values()) == {"PASS"}
    assert (out / "run" / "run_manifest.json").is_file()

    again = _demo(cwd=tmp_path)  # the freshness guard: refuse an existing non-empty output
    assert again.returncode == 2
    assert "already holds files" in again.stderr
