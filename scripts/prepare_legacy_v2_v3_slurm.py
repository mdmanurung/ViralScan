#!/usr/bin/env python3
"""Prepare a bounded, tiered SLURM plan without submitting any jobs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

FOUR_GIB = 4 * 1024**3
EXPECTED_ROWS = 44
MAX_CONCURRENCY = 4
MANIFEST_FIELDS = (
    "array_index",
    "run_id",
    "logical_id",
    "technical_repeat_group",
    "chemistry",
    "sample_dir",
    "output_bus_path",
    "output_bus_bytes",
    "whitelist_path",
    "t2g_path",
    "source_root",
    "raw_manifest_path",
    "raw_manifest_sha256",
    "resource_tier",
    "cpus",
    "memory_gib",
    "walltime",
    "scratch_required_bytes",
    "output_path",
    "stdout_path",
    "stderr_path",
    "frozen_python_path",
    "frozen_python_sha256",
    "frozen_benchmark_path",
    "frozen_benchmark_sha256",
    "frozen_compare_path",
    "frozen_compare_sha256",
    "frozen_helper_path",
    "frozen_helper_sha256",
    "frozen_source_dir",
    "execution_provenance_path",
    "execution_provenance_sha256",
    "comparison_argv_json",
)


class SlurmPreparationError(RuntimeError):
    """Raised when a safe frozen execution plan cannot be generated."""


@dataclass(frozen=True)
class PreparedOutputs:
    manifest: Path
    submission_plan: Path
    wrapper: Path


def _read_tsv(path: Path, required_fields: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise SlurmPreparationError(f"missing manifest: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = sorted(required_fields - fields)
        if missing:
            raise SlurmPreparationError(f"{path}: missing fields {missing}")
        return list(reader)


def _required_absolute_file(value: str, label: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise SlurmPreparationError(f"{label} must be an absolute path: {value!r}")
    if not path.is_file():
        raise SlurmPreparationError(f"missing required {label}: {path}")
    return path.resolve()


def _required_absolute_directory(path: Path, label: str) -> Path:
    if not path.is_absolute():
        raise SlurmPreparationError(f"{label} must be an absolute path: {path}")
    if not path.is_dir():
        raise SlurmPreparationError(f"missing required {label}: {path}")
    return path.resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_execution_provenance(path: Path) -> dict[str, object]:
    try:
        provenance = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SlurmPreparationError(
            f"invalid execution provenance JSON: {path}"
        ) from exc
    if not isinstance(provenance, dict):
        raise SlurmPreparationError("execution provenance must be a JSON object")
    if not isinstance(provenance.get("environment_prefix"), str):
        raise SlurmPreparationError("execution provenance lacks environment_prefix")
    if not isinstance(provenance.get("sources"), dict):
        raise SlurmPreparationError("execution provenance lacks sources")
    return provenance


def _verify_frozen_source(
    *,
    path: Path,
    source_dir: Path,
    provenance: dict[str, object],
    label: str,
    expected_sha256: str | None = None,
) -> str:
    try:
        relative = path.relative_to(source_dir).as_posix()
    except ValueError as exc:
        raise SlurmPreparationError(
            f"frozen {label} is outside frozen source directory: {path}"
        ) from exc
    sources = provenance["sources"]
    assert isinstance(sources, dict)
    record = sources.get(relative)
    if not isinstance(record, dict):
        raise SlurmPreparationError(
            f"frozen {label} is absent from execution provenance sources: {relative}"
        )
    recorded_size = record.get("size_bytes")
    recorded_sha256 = record.get("sha256")
    actual_sha256 = _sha256(path)
    if (
        not isinstance(recorded_size, int)
        or recorded_size != path.stat().st_size
        or not isinstance(recorded_sha256, str)
        or recorded_sha256 != actual_sha256
        or (expected_sha256 is not None and expected_sha256 != actual_sha256)
    ):
        raise SlurmPreparationError(
            f"frozen {label} does not match execution provenance: {relative}"
        )
    return actual_sha256


def _verify_all_provenance_sources(
    *, source_dir: Path, provenance: dict[str, object]
) -> None:
    sources = provenance["sources"]
    assert isinstance(sources, dict)
    for relative, record in sources.items():
        if not isinstance(relative, str) or not isinstance(record, dict):
            raise SlurmPreparationError("invalid execution provenance source record")
        path = (source_dir / relative).resolve()
        try:
            path.relative_to(source_dir)
        except ValueError as exc:
            raise SlurmPreparationError(
                f"execution provenance source escapes frozen source directory: {relative}"
            ) from exc
        if not path.is_file():
            raise SlurmPreparationError(
                f"missing frozen provenance source: {relative}"
            )
        recorded_size = record.get("size_bytes")
        recorded_sha256 = record.get("sha256")
        if (
            not isinstance(recorded_size, int)
            or recorded_size != path.stat().st_size
            or not isinstance(recorded_sha256, str)
            or recorded_sha256 != _sha256(path)
        ):
            raise SlurmPreparationError(
                f"frozen provenance source does not match execution provenance: "
                f"{relative}"
            )


def _resolve_frozen_execution(
    *,
    frozen_python: Path,
    frozen_benchmark: Path,
    frozen_compare: Path,
    frozen_helper: Path,
    frozen_source_dir: Path,
    execution_provenance: Path,
) -> dict[str, str]:
    python = _required_absolute_file(str(frozen_python), "frozen Python")
    benchmark = _required_absolute_file(str(frozen_benchmark), "frozen benchmark")
    compare = _required_absolute_file(str(frozen_compare), "frozen compare wrapper")
    helper = _required_absolute_file(str(frozen_helper), "frozen helper")
    source_dir = _required_absolute_directory(
        frozen_source_dir, "frozen source directory"
    )
    provenance_path = _required_absolute_file(
        str(execution_provenance), "execution provenance"
    )
    provenance = _read_execution_provenance(provenance_path)
    environment_prefix = Path(str(provenance["environment_prefix"]))
    if not environment_prefix.is_absolute():
        raise SlurmPreparationError(
            "execution provenance environment_prefix must be absolute"
        )
    expected_python = (environment_prefix / "bin/python").resolve()
    if python != expected_python:
        raise SlurmPreparationError(
            f"frozen Python is not provenance environment Python: {python}"
        )
    benchmark_sha256 = _verify_frozen_source(
        path=benchmark,
        source_dir=source_dir,
        provenance=provenance,
        label="benchmark",
    )
    compare_sha256 = _verify_frozen_source(
        path=compare,
        source_dir=source_dir,
        provenance=provenance,
        label="compare wrapper",
    )
    helper_sha256 = _verify_frozen_source(
        path=helper,
        source_dir=source_dir,
        provenance=provenance,
        label="helper",
    )
    _verify_all_provenance_sources(
        source_dir=source_dir,
        provenance=provenance,
    )
    return {
        "frozen_python_path": str(python),
        "frozen_python_sha256": _sha256(python),
        "frozen_benchmark_path": str(benchmark),
        "frozen_benchmark_sha256": benchmark_sha256,
        "frozen_compare_path": str(compare),
        "frozen_compare_sha256": compare_sha256,
        "frozen_helper_path": str(helper),
        "frozen_helper_sha256": helper_sha256,
        "frozen_source_dir": str(source_dir),
        "execution_provenance_path": str(provenance_path),
        "execution_provenance_sha256": _sha256(provenance_path),
    }


def _safe_run_id(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9._-]+", value):
        raise SlurmPreparationError(f"unsafe run_id: {value!r}")
    return value


def _tier(bus_bytes: int) -> tuple[str, int, str]:
    if bus_bytes < FOUR_GIB:
        return "bus-below-4-gib", 32, "06:00:00"
    return "bus-at-least-4-gib", 128, "24:00:00"


def _write_tsv(path: Path, rows: list[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _array_spec(indices: list[int]) -> str:
    if not indices:
        raise SlurmPreparationError("resource tier has no array indices")
    return ",".join(str(index) for index in indices) + f"%{MAX_CONCURRENCY}"


def _assert_ignored(repo_root: Path, run_root: Path) -> None:
    resolved_repo = repo_root.resolve()
    resolved_run = run_root.resolve()
    try:
        relative = resolved_run.relative_to(resolved_repo)
    except ValueError as exc:
        raise SlurmPreparationError("run root must be inside the repository") from exc
    result = subprocess.run(
        ["git", "-C", str(resolved_repo), "check-ignore", "-q", str(relative)],
        check=False,
    )
    if result.returncode != 0:
        raise SlurmPreparationError("run root must be ignored by Git")


def prepare(
    *,
    raw_manifest: Path,
    cohort_manifest: Path,
    repo_root: Path,
    run_root: Path,
    frozen_python: Path,
    frozen_benchmark: Path,
    frozen_compare: Path,
    frozen_helper: Path,
    frozen_source_dir: Path,
    execution_provenance: Path,
    source_root: Path,
) -> PreparedOutputs:
    """Generate the 44-row execution manifest, tier plan, and non-submitting wrapper."""

    repo_root = repo_root.resolve()
    run_root = run_root.resolve()
    _assert_ignored(repo_root, run_root)
    raw_manifest = _required_absolute_file(str(raw_manifest), "raw manifest")
    source_root = _required_absolute_directory(source_root, "archive source root")
    frozen = _resolve_frozen_execution(
        frozen_python=frozen_python,
        frozen_benchmark=frozen_benchmark,
        frozen_compare=frozen_compare,
        frozen_helper=frozen_helper,
        frozen_source_dir=frozen_source_dir,
        execution_provenance=execution_provenance,
    )
    raw_rows = _read_tsv(
        raw_manifest,
        {
            "run_id",
            "logical_id",
            "technical_repeat_group",
            "chemistry",
            "output_bus_path",
            "output_bus_bytes",
            "whitelist_path",
            "t2g_path",
        },
    )
    cohort_rows = _read_tsv(
        cohort_manifest,
        {"run_id", "logical_id", "technical_repeat_group", "chemistry"},
    )
    if len(raw_rows) != EXPECTED_ROWS or len(cohort_rows) != EXPECTED_ROWS:
        raise SlurmPreparationError(
            f"expected {EXPECTED_ROWS} raw and cohort rows, found "
            f"{len(raw_rows)} and {len(cohort_rows)}"
        )
    raw_ids = [row["run_id"] for row in raw_rows]
    cohort_ids = [row["run_id"] for row in cohort_rows]
    if len(set(raw_ids)) != EXPECTED_ROWS or set(raw_ids) != set(cohort_ids):
        raise SlurmPreparationError("raw and cohort manifests do not contain the same 44 run IDs")

    cohort_by_id = {row["run_id"]: row for row in cohort_rows}
    rows: list[dict[str, object]] = []
    tier_indices: dict[str, list[int]] = {
        "bus-below-4-gib": [],
        "bus-at-least-4-gib": [],
    }
    for index, raw in enumerate(sorted(raw_rows, key=lambda row: row["run_id"])):
        run_id = _safe_run_id(raw["run_id"])
        cohort = cohort_by_id[run_id]
        for field in ("logical_id", "technical_repeat_group", "chemistry"):
            if raw[field] != cohort[field]:
                raise SlurmPreparationError(f"{run_id}: raw/cohort mismatch for {field}")
        bus = _required_absolute_file(raw["output_bus_path"], f"{run_id} BUS")
        whitelist = _required_absolute_file(raw["whitelist_path"], f"{run_id} whitelist")
        t2g = _required_absolute_file(raw["t2g_path"], f"{run_id} t2g")
        try:
            bus_bytes = int(raw["output_bus_bytes"])
        except ValueError as exc:
            raise SlurmPreparationError(f"{run_id}: invalid output_bus_bytes") from exc
        if bus_bytes < 0 or bus.stat().st_size != bus_bytes:
            raise SlurmPreparationError(f"{run_id}: BUS size does not match raw manifest")

        tier, memory_gib, walltime = _tier(bus_bytes)
        tier_indices[tier].append(index)
        output_path = run_root / "rows" / run_id
        stdout_path = run_root / "logs" / f"{tier}-%A_%a.out"
        stderr_path = run_root / "logs" / f"{tier}-%A_%a.err"
        command = [
            frozen["frozen_python_path"],
            frozen["frozen_compare_path"],
            "run-row",
            "--raw-manifest",
            str(raw_manifest),
            "--source-root",
            str(source_root),
            "--run-id",
            run_id,
            "--output-root",
            str(run_root / "rows"),
            "--threads",
            "8",
        ]
        rows.append(
            {
                "array_index": index,
                "run_id": run_id,
                "logical_id": raw["logical_id"],
                "technical_repeat_group": raw["technical_repeat_group"],
                "chemistry": raw["chemistry"],
                "sample_dir": str(bus.parent.parent),
                "output_bus_path": str(bus),
                "output_bus_bytes": bus_bytes,
                "whitelist_path": str(whitelist),
                "t2g_path": str(t2g),
                "source_root": str(source_root),
                "raw_manifest_path": str(raw_manifest),
                "raw_manifest_sha256": _sha256(raw_manifest),
                "resource_tier": tier,
                "cpus": 8,
                "memory_gib": memory_gib,
                "walltime": walltime,
                "scratch_required_bytes": bus_bytes * 3,
                "output_path": str(output_path),
                "stdout_path": str(stdout_path),
                "stderr_path": str(stderr_path),
                **frozen,
                "comparison_argv_json": json.dumps(command, separators=(",", ":")),
            }
        )

    manifest = run_root / "slurm_manifest.tsv"
    wrapper = run_root / "run_array_row.sh"
    submission_plan = run_root / "submission_plan.tsv"
    (run_root / "logs").mkdir(parents=True, exist_ok=True)
    _write_tsv(manifest, rows, MANIFEST_FIELDS)
    wrapper.parent.mkdir(parents=True, exist_ok=True)
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        ': "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"\n'
        f'exec "{frozen["frozen_python_path"]}" '
        f'"{frozen["frozen_helper_path"]}" run-row --manifest "{manifest}" '
        '--index "$SLURM_ARRAY_TASK_ID"\n',
        encoding="utf-8",
    )
    wrapper.chmod(0o755)

    plan_rows = []
    for tier, memory_gib, walltime in (
        ("bus-below-4-gib", 32, "06:00:00"),
        ("bus-at-least-4-gib", 128, "24:00:00"),
    ):
        indices = tier_indices[tier]
        if not indices:
            continue
        stdout_path = str(run_root / "logs" / f"{tier}-%A_%a.out")
        stderr_path = str(run_root / "logs" / f"{tier}-%A_%a.err")
        array_spec = _array_spec(indices)
        sbatch_argv = [
            "sbatch",
            f"--array={array_spec}",
            "--cpus-per-task=8",
            f"--mem={memory_gib}G",
            f"--time={walltime}",
            f"--output={stdout_path}",
            f"--error={stderr_path}",
            str(wrapper),
        ]
        plan_rows.append(
            {
                "resource_tier": tier,
                "array_spec": array_spec,
                "max_concurrency": MAX_CONCURRENCY,
                "cpus": 8,
                "memory_gib": memory_gib,
                "walltime": walltime,
                "stdout_path": stdout_path,
                "stderr_path": stderr_path,
                "wrapper_path": str(wrapper),
                "submit_sequentially": "true",
                "sbatch_argv_json": json.dumps(sbatch_argv, separators=(",", ":")),
            }
        )
    _write_tsv(
        submission_plan,
        plan_rows,
        (
            "resource_tier",
            "array_spec",
            "max_concurrency",
            "cpus",
            "memory_gib",
            "walltime",
            "stdout_path",
            "stderr_path",
            "wrapper_path",
            "submit_sequentially",
            "sbatch_argv_json",
        ),
    )
    return PreparedOutputs(
        manifest=manifest,
        submission_plan=submission_plan,
        wrapper=wrapper,
    )


def preflight_tmpdir(tmpdir: Path, required_bytes: int) -> int:
    """Fail before execution unless TMPDIR is writable and has 3x BUS capacity."""

    if not tmpdir.is_dir():
        raise SlurmPreparationError(f"TMPDIR is not a directory: {tmpdir}")
    if not os.access(tmpdir, os.W_OK):
        raise SlurmPreparationError(f"TMPDIR is not writable: {tmpdir}")
    available = shutil.disk_usage(tmpdir).free
    if available < required_bytes:
        raise SlurmPreparationError(
            f"TMPDIR has {available} bytes free; {required_bytes} bytes required"
        )
    return available


def _verify_row_frozen_execution(row: dict[str, str]) -> list[str]:
    provenance_path = _required_absolute_file(
        row["execution_provenance_path"], "execution provenance"
    )
    if _sha256(provenance_path) != row["execution_provenance_sha256"]:
        raise SlurmPreparationError(
            "execution provenance does not match the prepared manifest"
        )
    frozen = _resolve_frozen_execution(
        frozen_python=Path(row["frozen_python_path"]),
        frozen_benchmark=Path(row["frozen_benchmark_path"]),
        frozen_compare=Path(row["frozen_compare_path"]),
        frozen_helper=Path(row["frozen_helper_path"]),
        frozen_source_dir=Path(row["frozen_source_dir"]),
        execution_provenance=provenance_path,
    )
    for field, value in frozen.items():
        if row[field] != value:
            raise SlurmPreparationError(
                f"{field} does not match the prepared frozen execution inputs"
            )
    raw_manifest = _required_absolute_file(row["raw_manifest_path"], "raw manifest")
    if _sha256(raw_manifest) != row["raw_manifest_sha256"]:
        raise SlurmPreparationError(
            "raw manifest does not match the prepared execution inputs"
        )
    source_root = _required_absolute_directory(
        Path(row["source_root"]), "archive source root"
    )
    try:
        command = json.loads(row["comparison_argv_json"])
    except json.JSONDecodeError as exc:
        raise SlurmPreparationError("invalid comparison command JSON") from exc
    expected_command = [
        row["frozen_python_path"],
        row["frozen_compare_path"],
        "run-row",
        "--raw-manifest",
        str(raw_manifest),
        "--source-root",
        str(source_root),
        "--run-id",
        row["run_id"],
        "--output-root",
        str(Path(row["output_path"]).parent),
        "--threads",
        "8",
    ]
    if command != expected_command:
        raise SlurmPreparationError(
            "comparison command does not match the frozen row execution contract"
        )
    return command


def run_row(*, manifest: Path, index: int) -> int:
    rows = _read_tsv(manifest, set(MANIFEST_FIELDS))
    selected = [row for row in rows if int(row["array_index"]) == index]
    if len(selected) != 1:
        raise SlurmPreparationError(f"array index {index} does not select exactly one row")
    row = selected[0]
    command = _verify_row_frozen_execution(row)
    tmpdir_value = os.environ.get("TMPDIR")
    if not tmpdir_value:
        raise SlurmPreparationError("TMPDIR is required")
    tmpdir = Path(tmpdir_value).resolve()
    preflight_tmpdir(tmpdir, int(row["scratch_required_bytes"]))
    numba_cache_dir = tmpdir / "viralscan-numba-cache" / _safe_run_id(row["run_id"])
    numba_cache_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    frozen_bin = str(Path(row["frozen_python_path"]).parent)
    current_path = environment.get("PATH")
    environment["PATH"] = (
        f"{frozen_bin}{os.pathsep}{current_path}" if current_path else frozen_bin
    )
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["NUMBA_CACHE_DIR"] = str(numba_cache_dir)
    current_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        f"{row['frozen_source_dir']}{os.pathsep}{current_pythonpath}"
        if current_pythonpath
        else row["frozen_source_dir"]
    )
    result = subprocess.run(command, check=False, env=environment)
    return result.returncode


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--raw-manifest", type=Path, required=True)
    prepare_parser.add_argument("--cohort-manifest", type=Path, required=True)
    prepare_parser.add_argument("--repo-root", type=Path, required=True)
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--frozen-python", type=Path, required=True)
    prepare_parser.add_argument("--frozen-benchmark", type=Path, required=True)
    prepare_parser.add_argument("--frozen-compare", type=Path, required=True)
    prepare_parser.add_argument("--frozen-helper", type=Path, required=True)
    prepare_parser.add_argument("--frozen-source-dir", type=Path, required=True)
    prepare_parser.add_argument("--execution-provenance", type=Path, required=True)
    prepare_parser.add_argument("--source-root", type=Path, required=True)
    run_parser = subparsers.add_parser("run-row")
    run_parser.add_argument("--manifest", type=Path, required=True)
    run_parser.add_argument("--index", type=int, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "prepare":
            outputs = prepare(
                raw_manifest=args.raw_manifest,
                cohort_manifest=args.cohort_manifest,
                repo_root=args.repo_root,
                run_root=args.run_root,
                frozen_python=args.frozen_python,
                frozen_benchmark=args.frozen_benchmark,
                frozen_compare=args.frozen_compare,
                frozen_helper=args.frozen_helper,
                frozen_source_dir=args.frozen_source_dir,
                execution_provenance=args.execution_provenance,
                source_root=args.source_root,
            )
            print(outputs.manifest)
            print(outputs.submission_plan)
            print(outputs.wrapper)
            return 0
        if args.command == "run-row":
            return run_row(manifest=args.manifest, index=args.index)
    except SlurmPreparationError as exc:
        raise SystemExit(f"SLURM preparation failed: {exc}") from exc
    raise AssertionError(f"unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
