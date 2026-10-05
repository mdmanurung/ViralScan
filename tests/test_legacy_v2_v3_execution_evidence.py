"""Tests for frozen execution evidence and deterministic pilot comparison."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from scripts import freeze_legacy_v2_v3_execution

pytestmark = pytest.mark.research


def _executable(path: Path, body: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path


def test_freeze_records_exact_sources_environment_versions_tools_and_doctor(
    tmp_path: Path,
) -> None:
    wheel = tmp_path / "dist/viralscan.whl"
    wheel.parent.mkdir()
    wheel.write_bytes(b"wheel")
    script = tmp_path / "scripts/benchmark.py"
    script.parent.mkdir()
    script.write_text("print('benchmark')\n", encoding="utf-8")

    prefix = tmp_path / "env"
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin/python").symlink_to(sys.executable)
    (prefix / "pyvenv.cfg").write_text(
        "include-system-site-packages = false\n",
        encoding="utf-8",
    )
    doctor = _executable(
        prefix / "bin/viralscan",
        "printf '%s\\n' "
        """'{"profile":"full","ok":true,"python":{"numpy":true},"tools":{},"schemas":{}}'"""
        "\n",
    )
    bustools = _executable(
        tmp_path / "tools/bustools",
        "printf '%s\\n' 'bustools, version 0.45.1'\n",
    )
    output = tmp_path / "evidence.json"

    assert (
        freeze_legacy_v2_v3_execution.main(
            [
                "--wheel",
                str(wheel),
                "--script",
                str(script),
                "--environment-prefix",
                str(prefix),
                "--tool",
                f"bustools={bustools}",
                "--package",
                "pip",
                "--doctor-executable",
                str(doctor),
                "--doctor-profile",
                "full",
                "--output",
                str(output),
            ]
        )
        == 0
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "1.0.0"
    assert payload["sources"]["wheel"]["sha256"]
    assert payload["sources"]["scripts"][0]["sha256"]
    assert payload["environment"]["python"]["is_symlink"] is True
    assert payload["environment"]["packages"]["pip"]
    assert payload["environment"]["tools"]["bustools"]["version"] == "0.45.1"
    assert payload["environment"]["tools"]["bustools"]["sha256"]
    assert payload["doctor"]["return_code"] == 0
    assert payload["doctor"]["report"]["ok"] is True
    assert payload["evidence_sha256"]
