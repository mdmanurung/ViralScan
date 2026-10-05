from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import freeze_legacy_v2_v3_execution as freeze

pytestmark = pytest.mark.research


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    env = tmp_path / "env"
    (env / "bin").mkdir(parents=True)
    for name in ("python", "viralscan", "bustools"):
        (env / "bin" / name).write_bytes(name.encode())
    source = tmp_path / "source"
    source.mkdir()
    (source / "benchmark.py").write_text("value = 1\n")
    (source / "__pycache__").mkdir()
    (source / "__pycache__" / "benchmark.cpython-311.pyc").write_bytes(b"mutable")
    spec = tmp_path / "environment.yml"
    spec.write_text("dependencies: [python=3.11]\n")
    conda = tmp_path / "conda"
    conda.write_bytes(b"conda")
    return env, source, spec, conda


def test_freeze_execution_records_green_doctor_sources_and_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    env, source, spec, conda = _fixture(tmp_path)

    def fake_run(command: list[str], *, env=None) -> str:
        del env
        if "doctor" in command:
            return json.dumps(
                {
                    "ok": True,
                    "tools": {"bustools": str(command[0]).replace("viralscan", "bustools")},
                }
            )
        if "list" in command:
            return "@EXPLICIT\nhttps://example.invalid/package.conda#abc\n"
        return json.dumps(
            {
                "ViralScan": "3.0.0.dev0",
                "kb-python": "0.28.2",
                "anndata": "0.12.19",
                "numpy": "2.4.6",
                "scipy": "1.17.1",
                "pandas": "2.3.3",
            }
        )

    monkeypatch.setattr(freeze, "_run", fake_run)
    output = tmp_path / "out"
    payload = freeze.freeze_execution(
        env_prefix=env,
        source_dir=source,
        resolved_spec=spec,
        conda_executable=conda,
        output_dir=output,
    )

    assert payload["package_versions"]["ViralScan"] == "3.0.0.dev0"
    assert payload["sources"]["benchmark.py"]["sha256"]
    assert not any("__pycache__" in path for path in payload["sources"])
    assert not any(path.endswith((".pyc", ".pyo")) for path in payload["sources"])
    assert payload["tools"]["bustools"]["path_relative_to_env"] == "bin/bustools"
    assert (output / "doctor_full.json").is_file()
    assert (output / "conda_explicit.txt").is_file()
    assert (output / "execution_provenance.json").is_file()


def test_freeze_execution_rejects_failed_doctor_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    env, source, spec, conda = _fixture(tmp_path)
    monkeypatch.setattr(
        freeze,
        "_run",
        lambda *_args, **_kwargs: json.dumps({"ok": False, "tools": {}}),
    )
    output = tmp_path / "out"

    with pytest.raises(freeze.ExecutionFreezeError, match="doctor"):
        freeze.freeze_execution(
            env_prefix=env,
            source_dir=source,
            resolved_spec=spec,
            conda_executable=conda,
            output_dir=output,
        )

    assert not output.exists()
