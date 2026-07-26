import csv
import hashlib
from pathlib import Path

import pytest

from scripts.prepare_fresh_controls import (
    FreshControlPreparationError,
    prepare_tasks,
)


def _make_viral_cache(tmp_path: Path) -> tuple[Path, str]:
    """Create a pinned viral-annotation cache and return its root and manifest digest."""

    cache = tmp_path / "viralscan_cache"
    manifest = cache / "data" / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text('{"files": {}}\n', encoding="utf-8")
    return cache, hashlib.sha256(manifest.read_bytes()).hexdigest()


def test_prepare_fresh_controls_pairs_both_stacks_on_identical_inputs(
    tmp_path: Path,
) -> None:
    sample_ids = (
        "SRR12682296",
        "SRR12682297",
        "SRR12682298",
        "SRR6825024",
        "SRR6825025",
    )
    raw_manifest = tmp_path / "control_inputs.raw.tsv"
    with raw_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "sample_id",
                "read1_path",
                "read2_path",
                "read1_storage_bytes",
                "read2_storage_bytes",
                "read1_storage_sha256",
                "read2_storage_sha256",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        for sample_id in sample_ids:
            read1 = tmp_path / f"{sample_id}_1.fastq"
            read2 = tmp_path / f"{sample_id}_2.fastq"
            large = sample_id.startswith("SRR68")
            read1_bytes = 70 * 1024**3 if large else 10
            read2_bytes = 10
            with read1.open("wb") as handle:
                handle.truncate(read1_bytes)
            with read2.open("wb") as handle:
                handle.truncate(read2_bytes)
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "read1_path": read1,
                    "read2_path": read2,
                    "read1_storage_bytes": read1_bytes,
                    "read2_storage_bytes": read2_bytes,
                    "read1_storage_sha256": "1" * 64,
                    "read2_storage_sha256": "2" * 64,
                }
            )
    paths = {}
    for name in ("v2_viralscan", "v3_viralscan", "index", "t2g", "whitelist"):
        paths[name] = tmp_path / name
        paths[name].write_text(name, encoding="utf-8")

    cache, cache_sha256 = _make_viral_cache(tmp_path)

    tasks = prepare_tasks(
        raw_manifest=raw_manifest,
        output_root=tmp_path / "fresh",
        task_manifest=tmp_path / "packet/tasks.tsv",
        cores=8,
        attempt_id="attempt2",
        viralscan_cache=cache,
        viralscan_cache_manifest_sha256=cache_sha256,
        **paths,
    )

    assert len(tasks) == 10
    assert sum(task["tier"] == "small" for task in tasks) == 6
    assert sum(task["tier"] == "large" for task in tasks) == 4
    assert len((tmp_path / "packet/tasks.small.tsv").read_text().splitlines()) == 7
    assert len((tmp_path / "packet/tasks.large.tsv").read_text().splitlines()) == 5
    for sample_id in sample_ids:
        paired = [task for task in tasks if task["sample_id"] == sample_id]
        assert {task["stack"] for task in paired} == {"v2", "v3"}
        assert len({task["read1_path"] for task in paired}) == 1
        assert len({task["read2_path"] for task in paired}) == 1
        assert len({task["read1_storage_bytes"] for task in paired}) == 1
        assert len({task["read2_storage_bytes"] for task in paired}) == 1
        assert len({task["read1_storage_sha256"] for task in paired}) == 1
        assert len({task["read2_storage_sha256"] for task in paired}) == 1


def test_prepare_fresh_controls_rejects_storage_size_drift(tmp_path: Path) -> None:
    raw_manifest = tmp_path / "control_inputs.raw.tsv"
    with raw_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "sample_id",
                "read1_path",
                "read2_path",
                "read1_storage_bytes",
                "read2_storage_bytes",
                "read1_storage_sha256",
                "read2_storage_sha256",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        for sample_id in (
            "SRR12682296",
            "SRR12682297",
            "SRR12682298",
            "SRR6825024",
            "SRR6825025",
        ):
            read1 = tmp_path / f"{sample_id}_1.fastq"
            read2 = tmp_path / f"{sample_id}_2.fastq"
            read1.write_bytes(b"r1")
            read2.write_bytes(b"r2")
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "read1_path": read1,
                    "read2_path": read2,
                    "read1_storage_bytes": 999 if sample_id == "SRR12682296" else 2,
                    "read2_storage_bytes": 2,
                    "read1_storage_sha256": "1" * 64,
                    "read2_storage_sha256": "2" * 64,
                }
            )
    paths = {}
    for name in ("v2_viralscan", "v3_viralscan", "index", "t2g", "whitelist"):
        paths[name] = tmp_path / name
        paths[name].write_text(name, encoding="utf-8")

    cache, cache_sha256 = _make_viral_cache(tmp_path)

    with pytest.raises(FreshControlPreparationError, match="stored byte count drifted"):
        prepare_tasks(
            raw_manifest=raw_manifest,
            output_root=tmp_path / "fresh",
            task_manifest=tmp_path / "packet/tasks.tsv",
            cores=8,
            attempt_id="attempt2",
            viralscan_cache=cache,
            viralscan_cache_manifest_sha256=cache_sha256,
            **paths,
        )


def _write_five_control_manifest(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    """Write a valid five-sample raw manifest and the frozen executable inputs."""

    raw_manifest = tmp_path / "control_inputs.raw.tsv"
    with raw_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "sample_id",
                "read1_path",
                "read2_path",
                "read1_storage_bytes",
                "read2_storage_bytes",
                "read1_storage_sha256",
                "read2_storage_sha256",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        for sample_id in (
            "SRR12682296",
            "SRR12682297",
            "SRR12682298",
            "SRR6825024",
            "SRR6825025",
        ):
            read1 = tmp_path / f"{sample_id}_1.fastq"
            read2 = tmp_path / f"{sample_id}_2.fastq"
            large = sample_id.startswith("SRR68")
            read1_bytes = 70 * 1024**3 if large else 10
            with read1.open("wb") as stream:
                stream.truncate(read1_bytes)
            with read2.open("wb") as stream:
                stream.truncate(10)
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "read1_path": read1,
                    "read2_path": read2,
                    "read1_storage_bytes": read1_bytes,
                    "read2_storage_bytes": 10,
                    "read1_storage_sha256": "1" * 64,
                    "read2_storage_sha256": "2" * 64,
                }
            )
    paths = {}
    for name in ("v2_viralscan", "v3_viralscan", "index", "t2g", "whitelist"):
        paths[name] = tmp_path / name
        paths[name].write_text(name, encoding="utf-8")
    return raw_manifest, paths


def test_prepare_fresh_controls_emits_only_the_selected_stack(tmp_path: Path) -> None:
    """Attempt 3 reruns v3 alone; the four succeeded v2 rows are not re-executed."""

    raw_manifest, paths = _write_five_control_manifest(tmp_path)
    cache, cache_sha256 = _make_viral_cache(tmp_path)

    tasks = prepare_tasks(
        raw_manifest=raw_manifest,
        output_root=tmp_path / "fresh",
        task_manifest=tmp_path / "packet/tasks.tsv",
        cores=8,
        attempt_id="attempt3",
        viralscan_cache=cache,
        viralscan_cache_manifest_sha256=cache_sha256,
        stacks=("v3",),
        **paths,
    )

    assert len(tasks) == 5
    assert {task["stack"] for task in tasks} == {"v3"}
    assert all(task["attempt_id"] == "attempt3" for task in tasks)
    assert all(task["viralscan_cache_manifest_sha256"] == cache_sha256 for task in tasks)
    assert len((tmp_path / "packet/tasks.small.tsv").read_text().splitlines()) == 4
    assert len((tmp_path / "packet/tasks.large.tsv").read_text().splitlines()) == 3


def test_prepare_fresh_controls_emits_only_the_selected_samples(tmp_path: Path) -> None:
    raw_manifest, paths = _write_five_control_manifest(tmp_path)
    cache, cache_sha256 = _make_viral_cache(tmp_path)

    tasks = prepare_tasks(
        raw_manifest=raw_manifest,
        output_root=tmp_path / "fresh",
        task_manifest=tmp_path / "packet/tasks.tsv",
        cores=8,
        attempt_id="attempt3",
        viralscan_cache=cache,
        viralscan_cache_manifest_sha256=cache_sha256,
        stacks=("v3",),
        sample_ids=("SRR6825024",),
        **paths,
    )

    assert len(tasks) == 1
    assert tasks[0]["task_id"] == "v3__SRR6825024"


def test_prepare_fresh_controls_rejects_a_drifted_viral_cache_manifest(tmp_path: Path) -> None:
    raw_manifest, paths = _write_five_control_manifest(tmp_path)
    cache, _ = _make_viral_cache(tmp_path)

    with pytest.raises(FreshControlPreparationError, match="cache manifest drifted"):
        prepare_tasks(
            raw_manifest=raw_manifest,
            output_root=tmp_path / "fresh",
            task_manifest=tmp_path / "packet/tasks.tsv",
            cores=8,
            attempt_id="attempt3",
            viralscan_cache=cache,
            viralscan_cache_manifest_sha256="a" * 64,
            **paths,
        )


def test_prepare_fresh_controls_rejects_a_missing_viral_cache(tmp_path: Path) -> None:
    raw_manifest, paths = _write_five_control_manifest(tmp_path)

    with pytest.raises(FreshControlPreparationError, match="missing viral-data cache manifest"):
        prepare_tasks(
            raw_manifest=raw_manifest,
            output_root=tmp_path / "fresh",
            task_manifest=tmp_path / "packet/tasks.tsv",
            cores=8,
            attempt_id="attempt3",
            viralscan_cache=tmp_path / "absent_cache",
            viralscan_cache_manifest_sha256="a" * 64,
            **paths,
        )


def test_prepare_fresh_controls_rejects_an_unknown_sample_id(tmp_path: Path) -> None:
    raw_manifest, paths = _write_five_control_manifest(tmp_path)
    cache, cache_sha256 = _make_viral_cache(tmp_path)

    with pytest.raises(FreshControlPreparationError, match="unknown sample ids"):
        prepare_tasks(
            raw_manifest=raw_manifest,
            output_root=tmp_path / "fresh",
            task_manifest=tmp_path / "packet/tasks.tsv",
            cores=8,
            attempt_id="attempt3",
            viralscan_cache=cache,
            viralscan_cache_manifest_sha256=cache_sha256,
            sample_ids=("SRR0000000",),
            **paths,
        )
