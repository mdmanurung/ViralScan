"""Tests for preparing, but never submitting, the 44-row legacy SLURM plan."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import prepare_legacy_v2_v3_slurm

pytestmark = pytest.mark.research


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_inputs(
    tmp_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / ".gitignore").write_text("benchmark_runs/\n", encoding="utf-8")
    source_dir = tmp_path / "source-dist"
    benchmark = source_dir / "scripts/benchmark_v3_multimap.py"
    benchmark.parent.mkdir(parents=True)
    benchmark.write_text("raise SystemExit(0)\n", encoding="utf-8")
    compare = source_dir / "scripts/compare_legacy_v2_v3.py"
    compare.write_text("raise SystemExit(0)\n", encoding="utf-8")
    helper = source_dir / "scripts/prepare_legacy_v2_v3_slurm.py"
    helper.write_text("raise SystemExit(0)\n", encoding="utf-8")
    wheel = source_dir / "dist/viralscan-frozen.whl"
    wheel.parent.mkdir()
    wheel.write_bytes(b"frozen-wheel")
    env_prefix = tmp_path / "frozen-env"
    frozen_python = env_prefix / "bin/python"
    frozen_python.parent.mkdir(parents=True)
    frozen_python.write_bytes(b"frozen-python")
    frozen_python.chmod(0o755)
    provenance = tmp_path / "execution_provenance.json"
    provenance.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "environment_prefix": str(env_prefix.resolve()),
                "sources": {
                    "scripts/benchmark_v3_multimap.py": {
                        "size_bytes": benchmark.stat().st_size,
                        "sha256": _sha256(benchmark),
                    },
                    "scripts/compare_legacy_v2_v3.py": {
                        "size_bytes": compare.stat().st_size,
                        "sha256": _sha256(compare),
                    },
                    "scripts/prepare_legacy_v2_v3_slurm.py": {
                        "size_bytes": helper.stat().st_size,
                        "sha256": _sha256(helper),
                    },
                    "dist/viralscan-frozen.whl": {
                        "size_bytes": wheel.stat().st_size,
                        "sha256": _sha256(wheel),
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    subprocess.run(["git", "init", "-q", str(repo)], check=True)

    t2g = tmp_path / "reference/t2g.txt"
    t2g.parent.mkdir()
    t2g.write_text("tx\tgene\n", encoding="utf-8")
    whitelist = tmp_path / "reference/whitelist.txt"
    whitelist.write_text("AAAA\n", encoding="utf-8")
    raw_rows = []
    cohort_rows = []
    for index in range(44):
        run_id = f"run-{index:02d}"
        sample_dir = tmp_path / "source" / run_id
        bus = sample_dir / "kb-python/output.bus"
        bus.parent.mkdir(parents=True)
        bus.touch()
        with bus.open("r+b") as handle:
            handle.truncate(4 * 1024**3 if index == 43 else 1024 + index)
        row = {
            "run_id": run_id,
            "logical_id": f"logical-{index:02d}",
            "technical_repeat_group": "",
            "chemistry": "10xv3",
            "output_bus_path": str(bus),
            "output_bus_bytes": str(bus.stat().st_size),
            "whitelist_path": str(whitelist),
            "t2g_path": str(t2g),
        }
        raw_rows.append(row)
        cohort_rows.append(
            {
                key: row[key]
                for key in (
                    "run_id",
                    "logical_id",
                    "technical_repeat_group",
                    "chemistry",
                )
            }
        )

    raw_manifest = tmp_path / "raw.tsv"
    cohort_manifest = tmp_path / "cohort.tsv"
    with raw_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=raw_rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(raw_rows)
    with cohort_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cohort_rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(cohort_rows)
    return (
        repo,
        raw_manifest,
        cohort_manifest,
        frozen_python,
        benchmark,
        helper,
        source_dir,
        provenance,
    )


def test_prepare_writes_frozen_44_row_tiers_and_non_submitting_wrapper(
    tmp_path: Path,
) -> None:
    (
        repo,
        raw_manifest,
        cohort_manifest,
        frozen_python,
        benchmark,
        helper,
        source_dir,
        provenance,
    ) = _write_inputs(tmp_path)
    run_root = repo / "benchmark_runs/legacy-v2-v3"

    outputs = prepare_legacy_v2_v3_slurm.prepare(
        raw_manifest=raw_manifest,
        cohort_manifest=cohort_manifest,
        repo_root=repo,
        run_root=run_root,
        frozen_python=frozen_python,
        frozen_benchmark=benchmark,
        frozen_compare=source_dir / "scripts/compare_legacy_v2_v3.py",
        frozen_helper=helper,
        frozen_source_dir=source_dir,
        execution_provenance=provenance,
        source_root=tmp_path / "source",
    )

    with outputs.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with outputs.submission_plan.open(newline="", encoding="utf-8") as handle:
        plans = list(csv.DictReader(handle, delimiter="\t"))
    wrapper = outputs.wrapper.read_text(encoding="utf-8")

    assert len(rows) == 44
    assert [row["array_index"] for row in rows] == [str(index) for index in range(44)]
    assert rows[0]["resource_tier"] == "bus-below-4-gib"
    assert rows[0]["memory_gib"] == "32"
    assert rows[0]["walltime"] == "06:00:00"
    assert rows[43]["resource_tier"] == "bus-at-least-4-gib"
    assert rows[43]["memory_gib"] == "128"
    assert rows[43]["walltime"] == "24:00:00"
    assert int(rows[43]["scratch_required_bytes"]) == 3 * 4 * 1024**3
    assert all(Path(row["stdout_path"]).is_absolute() for row in rows)
    assert all(str(run_root.resolve()) in row["stdout_path"] for row in rows)
    assert (run_root / "logs").is_dir()

    assert {plan["max_concurrency"] for plan in plans} == {"4"}
    assert all(plan["array_spec"].endswith("%4") for plan in plans)
    assert all(plan["submit_sequentially"] == "true" for plan in plans)
    assert "sbatch" not in wrapper
    assert "SLURM_ARRAY_TASK_ID" in wrapper
    assert "run-row" in wrapper
    subprocess.run(["bash", "-n", str(outputs.wrapper)], check=True)

    command = json.loads(rows[0]["comparison_argv_json"])
    assert command == [
        str(frozen_python.resolve()),
        str((source_dir / "scripts/compare_legacy_v2_v3.py").resolve()),
        "run-row",
        "--raw-manifest",
        str(raw_manifest.resolve()),
        "--source-root",
        str((tmp_path / "source").resolve()),
        "--run-id",
        "run-00",
        "--output-root",
        str((run_root / "rows").resolve()),
        "--threads",
        "8",
    ]
    assert rows[0]["output_path"] == str((run_root / "rows/run-00").resolve())
    assert rows[0]["raw_manifest_path"] == str(raw_manifest.resolve())
    assert rows[0]["raw_manifest_sha256"] == _sha256(raw_manifest)
    assert rows[0]["source_root"] == str((tmp_path / "source").resolve())
    assert rows[0]["frozen_compare_path"] == str(
        (source_dir / "scripts/compare_legacy_v2_v3.py").resolve()
    )
    assert rows[0]["frozen_helper_path"] == str(helper.resolve())
    assert rows[0]["frozen_source_dir"] == str(source_dir.resolve())
    assert rows[0]["execution_provenance_path"] == str(provenance.resolve())
    assert f'exec "{frozen_python.resolve()}" "{helper.resolve()}" run-row' in wrapper


def test_tmpdir_preflight_enforces_the_exact_three_times_requirement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    required = 3 * 4096
    monkeypatch.setattr(
        prepare_legacy_v2_v3_slurm.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=required),
    )
    assert prepare_legacy_v2_v3_slurm.preflight_tmpdir(tmp_path, required) == required

    monkeypatch.setattr(
        prepare_legacy_v2_v3_slurm.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=required - 1),
    )
    with pytest.raises(
        prepare_legacy_v2_v3_slurm.SlurmPreparationError,
        match="bytes required",
    ):
        prepare_legacy_v2_v3_slurm.preflight_tmpdir(tmp_path, required)


def test_prepare_fails_before_writing_when_manifest_is_not_44_rows(
    tmp_path: Path,
) -> None:
    (
        repo,
        raw_manifest,
        cohort_manifest,
        frozen_python,
        benchmark,
        helper,
        source_dir,
        provenance,
    ) = _write_inputs(tmp_path)
    lines = raw_manifest.read_text(encoding="utf-8").splitlines()
    raw_manifest.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
    run_root = repo / "benchmark_runs/legacy-v2-v3"

    with pytest.raises(
        prepare_legacy_v2_v3_slurm.SlurmPreparationError,
        match="expected 44 raw and cohort rows",
    ):
        prepare_legacy_v2_v3_slurm.prepare(
            raw_manifest=raw_manifest,
            cohort_manifest=cohort_manifest,
            repo_root=repo,
            run_root=run_root,
            frozen_python=frozen_python,
            frozen_benchmark=benchmark,
            frozen_compare=source_dir / "scripts/compare_legacy_v2_v3.py",
            frozen_helper=helper,
            frozen_source_dir=source_dir,
            execution_provenance=provenance,
            source_root=tmp_path / "source",
        )

    assert not run_root.exists()


def test_run_row_uses_frozen_runtime_and_tmpdir_numba_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (
        repo,
        raw_manifest,
        cohort_manifest,
        frozen_python,
        benchmark,
        helper,
        source_dir,
        provenance,
    ) = _write_inputs(tmp_path)
    outputs = prepare_legacy_v2_v3_slurm.prepare(
        raw_manifest=raw_manifest,
        cohort_manifest=cohort_manifest,
        repo_root=repo,
        run_root=repo / "benchmark_runs/legacy-v2-v3",
        frozen_python=frozen_python,
        frozen_benchmark=benchmark,
        frozen_compare=source_dir / "scripts/compare_legacy_v2_v3.py",
        frozen_helper=helper,
        frozen_source_dir=source_dir,
        execution_provenance=provenance,
        source_root=tmp_path / "source",
    )
    captured: dict[str, object] = {}
    row_output = repo / "benchmark_runs/legacy-v2-v3/rows/run-00"

    def fake_run(
        command: list[str],
        *,
        check: bool,
        env: dict[str, str],
    ) -> SimpleNamespace:
        assert not row_output.exists()
        captured.update(command=command, check=check, env=env)
        return SimpleNamespace(returncode=0)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("TMPDIR", str(scratch))
    monkeypatch.setattr(prepare_legacy_v2_v3_slurm.subprocess, "run", fake_run)

    assert prepare_legacy_v2_v3_slurm.run_row(manifest=outputs.manifest, index=0) == 0

    assert captured["command"][:3] == [
        str(frozen_python.resolve()),
        str((source_dir / "scripts/compare_legacy_v2_v3.py").resolve()),
        "run-row",
    ]
    environment = captured["env"]
    assert isinstance(environment, dict)
    assert environment["PATH"].split(":")[0] == str(frozen_python.resolve().parent)
    assert environment["PYTHONDONTWRITEBYTECODE"] == "1"
    assert environment["PYTHONPATH"].split(":")[0] == str(source_dir.resolve())
    numba_cache = Path(environment["NUMBA_CACHE_DIR"])
    assert numba_cache.is_dir()
    assert numba_cache.is_relative_to(scratch.resolve())
    assert not row_output.exists()


@pytest.mark.parametrize(
    "frozen_input_name",
    ["benchmark", "compare wrapper", "helper", "other", "raw manifest"],
)
def test_run_row_rejects_mutated_frozen_source_before_output_creation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    frozen_input_name: str,
) -> None:
    (
        repo,
        raw_manifest,
        cohort_manifest,
        frozen_python,
        benchmark,
        helper,
        source_dir,
        provenance,
    ) = _write_inputs(tmp_path)
    outputs = prepare_legacy_v2_v3_slurm.prepare(
        raw_manifest=raw_manifest,
        cohort_manifest=cohort_manifest,
        repo_root=repo,
        run_root=repo / "benchmark_runs/legacy-v2-v3",
        frozen_python=frozen_python,
        frozen_benchmark=benchmark,
        frozen_compare=source_dir / "scripts/compare_legacy_v2_v3.py",
        frozen_helper=helper,
        frozen_source_dir=source_dir,
        execution_provenance=provenance,
        source_root=tmp_path / "source",
    )
    frozen_input = {
        "benchmark": benchmark,
        "compare wrapper": source_dir / "scripts/compare_legacy_v2_v3.py",
        "helper": helper,
        "other": source_dir / "dist/viralscan-frozen.whl",
        "raw manifest": raw_manifest,
    }[frozen_input_name]
    frozen_input.write_text("mutated\n", encoding="utf-8")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("TMPDIR", str(scratch))

    def unexpected_run(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("benchmark must not run after provenance failure")

    monkeypatch.setattr(prepare_legacy_v2_v3_slurm.subprocess, "run", unexpected_run)

    error_pattern = (
        f"frozen {frozen_input_name} does not match execution provenance"
        if frozen_input_name in {"benchmark", "compare wrapper", "helper"}
        else (
            "frozen provenance source does not match execution provenance"
            if frozen_input_name == "other"
            else "raw manifest does not match the prepared execution inputs"
        )
    )
    with pytest.raises(
        prepare_legacy_v2_v3_slurm.SlurmPreparationError,
        match=error_pattern,
    ):
        prepare_legacy_v2_v3_slurm.run_row(manifest=outputs.manifest, index=0)

    assert not (repo / "benchmark_runs/legacy-v2-v3/rows/run-00").exists()
