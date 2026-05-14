"""Integration smoke tests for CLI entrypoints.

These tests are intentionally lightweight and opt-in:
- marked with ``integration`` so they are excluded from default unit test runs.
- network-dependent tests are additionally marked with ``network``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.integration
class TestSmoke:
    def test_cli_help(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "viralscan", "--help"],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONPATH": "src"},
        )
        assert result.returncode == 0, result.stderr
        assert "ViralScan" in result.stdout


@pytest.mark.integration
@pytest.mark.network
def test_build_ref_no_kb(tmp_path: Path) -> None:
    ncbi_email = os.environ.get("NCBI_EMAIL")
    if not ncbi_email:
        pytest.skip("NCBI_EMAIL is required for live NCBI integration tests")

    out_dir = tmp_path / "viralscan_test_ref"
    cmd = [
        sys.executable,
        "-m",
        "viralscan.menu",
        "build-ref",
        "--no-kb-ref",
        "--host",
        "human",
        "--virus-accessions",
        "NC_045512.2",
        "--output",
        str(out_dir),
        "--ncbi-email",
        ncbi_email,
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONPATH": "src"},
    )

    assert result.returncode == 0, result.stderr
    assert (out_dir / "combined.fa").exists()
    assert (out_dir / "combined.gtf").exists()
