"""tchem is pinned to v0.1.0 and `check env` reports it (PLAN TCHEM-SOURCE-01, FIRSTRUN-01)."""

import json
import os
import re
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import install_tchem
from viralscan import tchem
from viralscan.validation import doctor_report

ROOT = Path(__file__).resolve().parents[1]
LINUX = "x86_64-unknown-linux-musl"
HEX64 = re.compile(r"\b[0-9a-f]{64}\b")


def _fake_tchem(tmp_path: Path, banner: str) -> Path:
    exe = tmp_path / "tchem"
    exe.write_text(f"#!/bin/sh\necho '{banner}'\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return exe


@pytest.mark.parametrize("name", ["Dockerfile", "Singularity.pixi.def"])
def test_build_recipes_pin_the_release_asset_and_its_digest(name: str) -> None:
    text = (ROOT / name).read_text()
    asset = f"tchem-v0.1.0-{LINUX}.tar.gz"
    assert f"releases/download/v0.1.0/{asset}" in text
    # the line that checks the archive carries the pinned 64-hex digest (the base image pin is another one)
    check = next(line for line in text.splitlines() if "sha256sum -c" in line)
    assert HEX64.findall(check) == [install_tchem.SHA256[LINUX]]
    assert "sha256sum -c" in text


def test_release_singularity_template_checks_the_version_without_a_post_section() -> None:
    text = (ROOT / "Singularity.def").read_text()
    assert "tchem 0.1.0" in text and "%post" not in text


def test_install_script_pins_four_assets_with_64_hex_digests() -> None:
    assert install_tchem.VERSION == tchem.PINNED_VERSION == "0.1.0"
    assert set(install_tchem.SHA256) == {
        "x86_64-unknown-linux-musl",
        "aarch64-unknown-linux-musl",
        "x86_64-apple-darwin",
        "aarch64-apple-darwin",
    }
    assert all(HEX64.fullmatch(d) for d in install_tchem.SHA256.values())
    docs = (ROOT / "docs" / "installation.md").read_text()
    assert all(d in docs for d in install_tchem.SHA256.values())


def test_doctor_reports_a_matching_tchem(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(tchem.ENV_BIN, str(_fake_tchem(tmp_path, "tchem 0.1.0")))
    block = doctor_report("pip")["tchem"]
    assert block["version"] == "0.1.0" and block["expected"] == "0.1.0"
    assert block["warning"] is None


def test_doctor_warns_but_stays_ok_for_another_version(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(tchem.ENV_BIN, str(_fake_tchem(tmp_path, "tchem 0.2.0")))
    report = doctor_report("pip")
    assert report["tchem"]["version"] == "0.2.0"
    assert "0.1.0" in report["tchem"]["warning"]
    assert report["ok"] is True


def test_doctor_warns_when_tchem_is_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(tchem.ENV_BIN, raising=False)
    monkeypatch.setenv("PATH", "/nonexistent")
    block = doctor_report("pip")["tchem"]
    assert block["path"] is None and "not installed" in block["warning"]


def test_check_env_prints_the_tchem_block_and_the_warning(tmp_path: Path) -> None:
    env = {
        **os.environ,
        "PYTHONPATH": str(ROOT / "src"),
        tchem.ENV_BIN: str(_fake_tchem(tmp_path, "tchem 0.2.0")),
    }
    proc = subprocess.run(
        [sys.executable, "-m", "viralscan", "check", "env", "--profile", "pip"],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["tchem"]["version"] == "0.2.0"
    assert "warning: tchem 0.2.0" in proc.stderr
