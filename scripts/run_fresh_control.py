#!/usr/bin/env python3
"""Run one frozen legacy-diagnostic control through one ViralScan stack."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import subprocess
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any


class FreshControlError(RuntimeError):
    """Raised when a fresh control run violates its frozen execution contract."""


def ensure_fresh_output(output: Path) -> None:
    """Fail before execution rather than resume or overwrite a prior attempt."""

    if output.exists():
        raise FreshControlError(f"fresh control output already exists: {output}")


def ensure_fresh_attempt(*paths: Path) -> None:
    """Refuse to replace status or log evidence from an earlier attempt."""

    for path in paths:
        if path.exists():
            raise FreshControlError(f"attempt evidence already exists: {path}")


def verify_frozen_fastq(
    path: Path,
    *,
    expected_bytes: int,
    expected_sha256: str,
) -> dict[str, object]:
    """Verify runtime size against an identity established by a full-stream audit."""

    if not path.is_file():
        raise FreshControlError(f"missing required fresh-run input: {path}")
    if expected_bytes <= 0:
        raise FreshControlError(f"invalid expected stored byte count: {path}")
    if path.stat().st_size != expected_bytes:
        raise FreshControlError(f"stored byte count drifted: {path}")
    if re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None:
        raise FreshControlError(f"invalid audited storage SHA-256: {path}")
    return {
        "path": str(path.resolve()),
        "storage_bytes": expected_bytes,
        "audited_storage_sha256": expected_sha256,
        "runtime_size_verified": True,
    }


def resolve_v2_output_root(output: Path) -> tuple[Path, list[str]]:
    """Locate the directory the legacy 2.2.0 workflow actually wrote.

    The legacy CLI nests its result tree one level beneath the requested output
    directory, named for the sample. Resolve that nesting instead of assuming a
    flat layout, and keep the resolution non-vacuous: neither a missing root nor
    an ambiguous one may pass as success.
    """

    if (output / "config.yaml").is_file():
        return output, []
    if not output.is_dir():
        return output, [f"missing v2 output directory: {output}"]
    candidates = sorted(
        child for child in output.iterdir() if child.is_dir() and (child / "config.yaml").is_file()
    )
    if len(candidates) == 1:
        return candidates[0], []
    if not candidates:
        return output, [f"no config.yaml at {output} or any immediate child directory"]
    named = ", ".join(candidate.name for candidate in candidates)
    return output, [f"ambiguous v2 output roots beneath {output}: {named}"]


def _v2_artifact_errors(output: Path) -> tuple[Path, list[str]]:
    required = (
        "config.yaml",
        "summary.txt",
        "kb-python/run_info.json",
        "kb-python/output.bus",
        "kb-python/counts_unfiltered/adata.h5ad",
        "kb-python/counts_unfiltered/adata_multimap.h5ad",
    )
    root, resolution_errors = resolve_v2_output_root(output)
    if resolution_errors:
        return root, resolution_errors
    return root, [
        relative
        for relative in required
        if not (root / relative).is_file() or (root / relative).stat().st_size == 0
    ]


def build_command(
    *,
    stack: str,
    viralscan: Path,
    output: Path,
    read1: Path,
    read2: Path,
    index: Path,
    t2g: Path,
    whitelist: Path,
    cores: int,
) -> list[str]:
    """Build an explicit argv for the frozen v2.2.0 or v3 control policy."""

    common = [
        "-t",
        str(t2g),
        "-i",
        str(index),
        "-o",
        f"{output}/",
        "-s1",
        str(read1),
        "-s2",
        str(read2),
        "-c",
        str(cores),
        "-x",
        "10xv2",
        "-w",
        str(whitelist),
    ]
    if stack == "v2":
        return [
            str(viralscan),
            *common,
            "-mm",
            "True",
            "-v",
            "True",
        ]
    if stack == "v3":
        return [
            str(viralscan),
            *common,
            "--multimapping",
            "--multimap-method",
            "host-conservative",
            "--cell-calling",
            "auto",
            "--visual",
        ]
    raise FreshControlError(f"unsupported stack: {stack}")


def run_control(
    *,
    sample_id: str,
    stack: str,
    viralscan: Path,
    output: Path,
    read1: Path,
    read2: Path,
    index: Path,
    t2g: Path,
    whitelist: Path,
    cores: int,
    status_path: Path,
    stdout_path: Path,
    stderr_path: Path,
    attempt_id: str,
    viralscan_cache: Path | None = None,
    read1_storage_bytes: int | None = None,
    read2_storage_bytes: int | None = None,
    read1_storage_sha256: str | None = None,
    read2_storage_sha256: str | None = None,
) -> dict[str, Any]:
    """Execute one fresh attempt and retain terminal status on command failure."""

    ensure_fresh_output(output)
    ensure_fresh_attempt(status_path, stdout_path, stderr_path)
    for path in (viralscan, index, t2g, whitelist):
        if not path.is_file():
            raise FreshControlError(f"missing required fresh-run input: {path}")
    identity_values = (
        read1_storage_bytes,
        read2_storage_bytes,
        read1_storage_sha256,
        read2_storage_sha256,
    )
    if any(value is not None for value in identity_values) and any(
        value is None for value in identity_values
    ):
        raise FreshControlError("incomplete frozen FASTQ identity")
    if all(value is not None for value in identity_values):
        assert read1_storage_bytes is not None
        assert read2_storage_bytes is not None
        assert read1_storage_sha256 is not None
        assert read2_storage_sha256 is not None
        input_storage = {
            "read1": verify_frozen_fastq(
                read1,
                expected_bytes=read1_storage_bytes,
                expected_sha256=read1_storage_sha256,
            ),
            "read2": verify_frozen_fastq(
                read2,
                expected_bytes=read2_storage_bytes,
                expected_sha256=read2_storage_sha256,
            ),
        }
    else:
        input_storage = None
        for path in (read1, read2):
            if not path.is_file():
                raise FreshControlError(f"missing required fresh-run input: {path}")
    command = build_command(
        stack=stack,
        viralscan=viralscan,
        output=output,
        read1=read1,
        read2=read2,
        index=index,
        t2g=t2g,
        whitelist=whitelist,
        cores=cores,
    )
    command_sha256 = hashlib.sha256(json.dumps(command, separators=(",", ":")).encode()).hexdigest()
    status_path.parent.mkdir(parents=True, exist_ok=True)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    runtime_env = os.environ.copy()
    runtime_env["PATH"] = f"{viralscan.parent}:{runtime_env.get('PATH', '')}"
    if viralscan_cache is not None:
        cache_manifest = viralscan_cache / "data" / "manifest.json"
        if not cache_manifest.is_file():
            raise FreshControlError(f"missing viral-data cache manifest: {cache_manifest}")
        runtime_env["VIRALSCAN_CACHE"] = str(viralscan_cache.resolve())
    started = time.monotonic()
    validation_command: list[str] | None = None
    validation_exit_code: int | None = None
    artifact_validation_errors: list[str] = []
    with (
        stdout_path.open("wb") as stdout,
        stderr_path.open("wb") as stderr,
    ):
        completed = subprocess.run(
            command,
            check=False,
            stdout=stdout,
            stderr=stderr,
            env=runtime_env,
        )
        if stack == "v3" and completed.returncode == 0:
            validation_command = [
                str(viralscan),
                "validate-run",
                str(output),
            ]
            validation_exit_code = subprocess.run(
                validation_command,
                check=False,
                stdout=stdout,
                stderr=stderr,
                env=runtime_env,
            ).returncode
    v2_output_root: str | None = None
    if stack == "v2" and completed.returncode == 0:
        resolved_root, artifact_validation_errors = _v2_artifact_errors(output)
        v2_output_root = str(resolved_root.resolve())
    exit_code = (
        validation_exit_code if validation_exit_code not in (None, 0) else completed.returncode
    )
    if exit_code == 0 and artifact_validation_errors:
        exit_code = 65
    if completed.returncode != 0:
        stage = "workflow"
    elif artifact_validation_errors or validation_exit_code not in (None, 0):
        stage = "artifact_validation"
    else:
        stage = "complete"
    payload: dict[str, Any] = {
        "schema_version": "1.1.0",
        "sample_id": sample_id,
        "run_id": sample_id,
        "stack": stack,
        "stage": stage,
        "attempt_id": attempt_id,
        "status": "success" if exit_code == 0 else "failed",
        "exit_code": exit_code,
        "workflow_exit_code": completed.returncode,
        "validation_exit_code": validation_exit_code,
        "artifact_validation_errors": artifact_validation_errors,
        "v2_output_root": v2_output_root,
        "elapsed_seconds": time.monotonic() - started,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "command": command,
        "command_sha256": command_sha256,
        "scientific_parameter_hash": command_sha256,
        "validation_command": validation_command,
        "input_storage": input_storage,
        "output": str(output.resolve()),
        "stdout": str(stdout_path.resolve()),
        "stderr": str(stderr_path.resolve()),
        "stderr_path": str(stderr_path.resolve()),
    }
    staging = status_path.with_name(f"{status_path.name}.tmp")
    staging.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    staging.replace(status_path)
    return payload


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stack", choices=("v2", "v3"), required=True)
    parser.add_argument("--viralscan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--read1-storage-bytes", type=int, required=True)
    parser.add_argument("--read2-storage-bytes", type=int, required=True)
    parser.add_argument("--read1-storage-sha256", required=True)
    parser.add_argument("--read2-storage-sha256", required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--t2g", type=Path, required=True)
    parser.add_argument("--whitelist", type=Path, required=True)
    parser.add_argument("--cores", type=int, required=True)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--stdout", type=Path, required=True)
    parser.add_argument("--stderr", type=Path, required=True)
    parser.add_argument(
        "--attempt-id",
        required=True,
        help="frozen attempt identifier retained in the status record",
    )
    parser.add_argument(
        "--viralscan-cache",
        type=Path,
        default=None,
        help="pinned viral-annotation cache root exported as VIRALSCAN_CACHE",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    payload = run_control(
        sample_id=args.sample_id,
        stack=args.stack,
        viralscan=args.viralscan,
        output=args.output,
        read1=args.read1,
        read2=args.read2,
        index=args.index,
        t2g=args.t2g,
        whitelist=args.whitelist,
        cores=args.cores,
        status_path=args.status,
        stdout_path=args.stdout,
        stderr_path=args.stderr,
        attempt_id=args.attempt_id,
        viralscan_cache=args.viralscan_cache,
        read1_storage_bytes=args.read1_storage_bytes,
        read2_storage_bytes=args.read2_storage_bytes,
        read1_storage_sha256=args.read1_storage_sha256,
        read2_storage_sha256=args.read2_storage_sha256,
    )
    return int(payload["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
