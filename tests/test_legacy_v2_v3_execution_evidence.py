"""Freeze regression follows the current environment and executable source contract."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts import freeze_legacy_v2_v3_execution as freeze

pytestmark = pytest.mark.research


@pytest.fixture
def execution_inputs(tmp_path, monkeypatch):
    prefix = tmp_path / "env"
    (prefix / "bin").mkdir(parents=True)
    for name in ("python", "viralscan", "bustools"):
        (prefix / "bin" / name).write_text(f"frozen {name}\n")
    source = tmp_path / "source"
    source.mkdir()
    (source / "benchmark.py").write_text("print('benchmark')\n")
    (source / "__pycache__").mkdir()
    (source / "__pycache__/benchmark.pyc").write_bytes(b"excluded cache")
    spec = tmp_path / "environment.yml"
    spec.write_text("name: frozen\n")
    conda = tmp_path / "conda"
    conda.write_text("frozen conda\n")
    doctor = {"profile": "full", "ok": True, "tools": {"bustools": str(prefix / "bin/bustools")}}
    commands = []

    def run(command, *, env=None):
        commands.append((command, env))
        if command[0] == str(conda):
            return "@EXPLICIT\nhttps://example.test/frozen.conda#012345\n"
        if command[0] == str(prefix / "bin/viralscan"):
            return json.dumps(doctor)
        if command[0] == str(prefix / "bin/python"):
            return json.dumps({"ViralScan": "3.0.0.dev1", "kb-python": "0.30.2"})
        raise AssertionError(f"unexpected executable: {command}")

    monkeypatch.setattr(freeze, "_run", run)
    return prefix, source, spec, conda, doctor, commands


def test_freeze_records_exact_sources_environment_versions_tools_and_doctor(
    tmp_path: Path, execution_inputs
):
    prefix, source, spec, conda, doctor, commands = execution_inputs
    output = tmp_path / "evidence"
    assert (
        freeze.main(
            [
                "--env-prefix",
                str(prefix),
                "--source-dir",
                str(source),
                "--resolved-spec",
                str(spec),
                "--conda-executable",
                str(conda),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    payload = json.loads((output / "execution_provenance.json").read_text())
    assert payload["sources"] == {
        "benchmark.py": {
            "size_bytes": (source / "benchmark.py").stat().st_size,
            "sha256": hashlib.sha256((source / "benchmark.py").read_bytes()).hexdigest(),
        }
    }
    assert payload["package_versions"] == {"ViralScan": "3.0.0.dev1", "kb-python": "0.30.2"}
    assert payload["tools"]["bustools"]["path_relative_to_env"] == "bin/bustools"
    assert (
        payload["tools"]["bustools"]["sha256"]
        == hashlib.sha256((prefix / "bin/bustools").read_bytes()).hexdigest()
    )
    assert json.loads((output / "doctor_full.json").read_text()) == doctor
    for field, filename in [
        ("doctor_sha256", "doctor_full.json"),
        ("conda_explicit_sha256", "conda_explicit.txt"),
    ]:
        assert payload[field] == hashlib.sha256((output / filename).read_bytes()).hexdigest()
    assert commands[1][0] == [str(conda), "list", "--explicit", "--md5", "--prefix", str(prefix)]
    assert commands[0][1]["PATH"].split(":")[0] == str(prefix / "bin")
    assert not (output / "execution_provenance.json.tmp").exists()


@pytest.mark.parametrize("failure", ["doctor", "explicit", "source"])
def test_incomplete_execution_cannot_produce_a_frozen_receipt(
    tmp_path, monkeypatch, execution_inputs, failure
):
    prefix, source, spec, conda, doctor, _ = execution_inputs
    if failure == "doctor":
        doctor["ok"] = False
    elif failure == "explicit":
        original = freeze._run
        monkeypatch.setattr(
            freeze,
            "_run",
            lambda cmd, **kw: "incomplete" if cmd[0] == str(conda) else original(cmd, **kw),
        )
    else:
        (source / "benchmark.py").unlink()
    output = tmp_path / "evidence"
    with pytest.raises(freeze.ExecutionFreezeError):
        freeze.freeze_execution(
            env_prefix=prefix,
            source_dir=source,
            resolved_spec=spec,
            conda_executable=conda,
            output_dir=output,
        )
    assert not output.exists()
