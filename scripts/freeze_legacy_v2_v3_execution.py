#!/usr/bin/env python3
"""Freeze the executable source and environment used by the legacy diagnostic."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from collections.abc import Sequence
from pathlib import Path


class ExecutionFreezeError(RuntimeError):
    """Raised when the execution snapshot is incomplete."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_file(path: Path) -> Path:
    if not path.is_file():
        raise ExecutionFreezeError(f"missing required file: {path}")
    return path.resolve()


def _run(command: list[str], *, env: dict[str, str] | None = None) -> str:
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ExecutionFreezeError(f"command failed: {command!r}: {exc}") from exc
    return result.stdout


def freeze_execution(
    *,
    env_prefix: Path,
    source_dir: Path,
    resolved_spec: Path,
    conda_executable: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Write a fail-closed environment/source identity below the ignored run root."""

    env_prefix = env_prefix.resolve()
    source_dir = source_dir.resolve()
    resolved_spec = _required_file(resolved_spec)
    conda_executable = _required_file(conda_executable)
    output_dir = output_dir.resolve()
    python = _required_file(env_prefix / "bin/python")
    viralscan = _required_file(env_prefix / "bin/viralscan")

    runtime_env = os.environ.copy()
    runtime_env["PATH"] = f"{env_prefix / 'bin'}:{runtime_env.get('PATH', '')}"
    doctor = json.loads(
        _run(
            [str(viralscan), "doctor", "--profile", "full", "--json"],
            env=runtime_env,
        )
    )
    if doctor.get("ok") is not True:
        raise ExecutionFreezeError("full-profile doctor did not pass")

    explicit = _run(
        [
            str(conda_executable),
            "list",
            "--explicit",
            "--md5",
            "--prefix",
            str(env_prefix),
        ]
    )
    if "@EXPLICIT" not in explicit:
        raise ExecutionFreezeError("Conda explicit export is incomplete")

    package_versions = json.loads(
        _run(
            [
                str(python),
                "-c",
                (
                    "import importlib.metadata as m,json;"
                    "print(json.dumps({n:m.version(n) for n in "
                    "['ViralScan','kb-python','anndata','numpy','scipy','pandas']}))"
                ),
            ],
            env=runtime_env,
        )
    )
    source_paths = sorted(
        path
        for path in source_dir.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and path.suffix not in {".pyc", ".pyo"}
    )
    if not source_paths:
        raise ExecutionFreezeError("source snapshot is empty")
    sources = {
        path.relative_to(source_dir).as_posix(): {
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in source_paths
    }
    tools = {}
    for name, value in sorted(doctor["tools"].items()):
        path = _required_file(Path(value))
        tools[name] = {
            "path_relative_to_env": path.relative_to(env_prefix).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    doctor_path = output_dir / "doctor_full.json"
    explicit_path = output_dir / "conda_explicit.txt"
    doctor_path.write_text(json.dumps(doctor, indent=2, sort_keys=True) + "\n")
    explicit_path.write_text(explicit)
    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "environment_prefix": str(env_prefix),
        "resolved_spec": {
            "sha256": _sha256(resolved_spec),
            "size_bytes": resolved_spec.stat().st_size,
        },
        "doctor_sha256": _sha256(doctor_path),
        "conda_explicit_sha256": _sha256(explicit_path),
        "package_versions": package_versions,
        "sources": sources,
        "tools": tools,
    }
    staging = output_dir / "execution_provenance.json.tmp"
    staging.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    staging.replace(output_dir / "execution_provenance.json")
    return payload


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--resolved-spec", type=Path, required=True)
    parser.add_argument("--conda-executable", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    freeze_execution(
        env_prefix=args.env_prefix,
        source_dir=args.source_dir,
        resolved_spec=args.resolved_spec,
        conda_executable=args.conda_executable,
        output_dir=args.output_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
