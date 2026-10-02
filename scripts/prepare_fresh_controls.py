#!/usr/bin/env python3
"""Prepare frozen, paired v2.2.0/v3 fresh-control SLURM task rows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections.abc import Sequence
from pathlib import Path


class FreshControlPreparationError(RuntimeError):
    """Raised when fresh-control tasks cannot be prepared without drift."""


EXPECTED_IDS = frozenset(
    {
        "SRR12682296",
        "SRR12682297",
        "SRR12682298",
        "SRR6825024",
        "SRR6825025",
    }
)
TASK_FIELDS = (
    "task_id",
    "sample_id",
    "stack",
    "tier",
    "attempt_id",
    "viralscan_cache_path",
    "viralscan_cache_manifest_sha256",
    "gtf_path",
    "gtf_sha256",
    "read1_path",
    "read2_path",
    "read1_storage_bytes",
    "read2_storage_bytes",
    "read1_storage_sha256",
    "read2_storage_sha256",
    "viralscan_path",
    "output_path",
    "index_path",
    "t2g_path",
    "whitelist_path",
    "cores",
    "status_path",
    "stdout_path",
    "stderr_path",
)
LARGE_INPUT_BYTES = 64 * 1024**3
# Protocol 1.1.0 tier bus-at-least-4-gib-highmem: assigned only to task ids that
# already ran out of memory at the 128 GiB tier (never by input size).
TIERS = ("small", "large", "highmem")
EMPTY = "-"  # tab-separated shell `read` collapses empty fields, so use a sentinel
_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gtf_gene_ids(path: Path) -> set[str]:
    """Return every ``gene_id`` attribute in a GTF (same rule as the v3 analysis step)."""
    ids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 9 and (match := _GENE_ID_RE.search(cols[8])):
                ids.add(match.group(1))
    return ids


def build_panel_gtf(v2_data_dir: Path, dest: Path) -> str:
    """Concatenate the v2 arm's packaged GTFs (sorted filenames) into ``dest``; return sha256."""
    files = sorted(v2_data_dir.glob("*.gtf"), key=lambda p: p.name)
    if not files:
        raise FreshControlPreparationError(f"no .gtf files under {v2_data_dir}")
    if dest.exists():
        raise FreshControlPreparationError(f"panel GTF already exists: {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("wb") as out:
        for file in files:
            data = file.read_bytes()
            out.write(data if data.endswith(b"\n") else data + b"\n")
    return _sha256(dest)


def gtf_t2g_parity(gtf: Path, t2g: Path) -> dict[str, object]:
    """Compare GTF gene ids with the frozen index t2g gene column (column 2).

    Passes only when the GTF set is a subset of the t2g gene set, i.e. every viral
    gene the v3 run would call viral exists in the frozen index. t2g genes absent
    from the GTF are recorded but do not fail the check.
    """
    gtf_genes = gtf_gene_ids(gtf)
    with t2g.open(encoding="utf-8") as handle:
        t2g_genes = {cols[1] for line in handle if len(cols := line.rstrip("\n").split("\t")) > 1}
    gtf_only = gtf_genes - t2g_genes
    # Assumption: host genes carry Ensembl ids (ENS*); every other t2g gene is treated as viral.
    t2g_viral = {g for g in t2g_genes if not g.startswith("ENS")}
    return {
        "gtf_sha256": _sha256(gtf),
        "t2g_sha256": _sha256(t2g),
        "gtf_genes": len(gtf_genes),
        "t2g_genes": len(t2g_genes),
        "t2g_non_ensembl_genes": len(t2g_viral),
        "intersection": len(gtf_genes & t2g_genes),
        "gtf_only": len(gtf_only),
        "gtf_only_genes": sorted(gtf_only),
        "t2g_non_ensembl_not_in_gtf": len(t2g_viral - gtf_genes),
        "t2g_non_ensembl_not_in_gtf_genes": sorted(t2g_viral - gtf_genes),
        "rule": "gtf_genes subset of t2g_genes (gtf_only == 0)",
        "passed": not gtf_only,
    }


def _write_tasks(path: Path, tasks: list[dict[str, object]]) -> None:
    staging = path.with_name(f"{path.name}.tmp")
    with staging.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TASK_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(tasks)
    staging.replace(path)


def prepare_tasks(
    *,
    raw_manifest: Path,
    output_root: Path,
    task_manifest: Path,
    v2_viralscan: Path,
    v3_viralscan: Path,
    index: Path,
    t2g: Path,
    whitelist: Path,
    cores: int,
    attempt_id: str,
    viralscan_cache: Path | None = None,
    viralscan_cache_manifest_sha256: str | None = None,
    gtf: Path | None = None,
    gtf_sha256: str | None = None,
    highmem_task_ids: Sequence[str] = (),
    stacks: Sequence[str] = ("v2", "v3"),
    sample_ids: Sequence[str] | None = None,
) -> list[dict[str, object]]:
    """Write paired task rows with identical per-sample inputs across both stacks.

    ``stacks`` and ``sample_ids`` narrow the emitted rows. Both default to the
    full ten-row matrix, so an unfiltered call is unchanged.
    """

    if cores <= 0:
        raise FreshControlPreparationError("cores must be positive")
    if not attempt_id:
        raise FreshControlPreparationError("attempt_id must be non-empty")
    unknown_stacks = sorted(set(stacks) - {"v2", "v3"})
    if unknown_stacks:
        raise FreshControlPreparationError(f"unsupported stacks: {', '.join(unknown_stacks)}")
    if not stacks:
        raise FreshControlPreparationError("at least one stack is required")
    selected_ids = EXPECTED_IDS if sample_ids is None else frozenset(sample_ids)
    unknown_ids = sorted(selected_ids - EXPECTED_IDS)
    if unknown_ids:
        raise FreshControlPreparationError(f"unknown sample ids: {', '.join(unknown_ids)}")
    if not selected_ids:
        raise FreshControlPreparationError("at least one sample id is required")
    if viralscan_cache is not None:
        if (
            viralscan_cache_manifest_sha256 is None
            or re.fullmatch(r"[0-9a-f]{64}", viralscan_cache_manifest_sha256) is None
        ):
            raise FreshControlPreparationError("invalid viral-data cache manifest SHA-256")
        cache_manifest = viralscan_cache / "data" / "manifest.json"
        if not cache_manifest.is_file():
            raise FreshControlPreparationError(
                f"missing viral-data cache manifest: {cache_manifest}"
            )
        if _sha256(cache_manifest) != viralscan_cache_manifest_sha256:
            raise FreshControlPreparationError(
                f"viral-data cache manifest drifted: {cache_manifest}"
            )
    elif viralscan_cache_manifest_sha256 is not None:
        raise FreshControlPreparationError("cache manifest SHA-256 given without a cache")
    if gtf is not None:
        if gtf_sha256 is None or re.fullmatch(r"[0-9a-f]{64}", gtf_sha256) is None:
            raise FreshControlPreparationError("invalid explicit GTF SHA-256")
        if not gtf.is_file() or _sha256(gtf) != gtf_sha256:
            raise FreshControlPreparationError(f"explicit GTF missing or drifted: {gtf}")
    elif "v3" in stacks and viralscan_cache is None:
        raise FreshControlPreparationError(
            "v3 rows need either --viralscan-cache or an explicit GTF"
        )
    known_tasks = {f"{st}__{sid}" for st in ("v2", "v3") for sid in EXPECTED_IDS}
    unknown_highmem = sorted(set(highmem_task_ids) - known_tasks)
    if unknown_highmem:
        raise FreshControlPreparationError(
            f"unknown highmem task ids: {', '.join(unknown_highmem)}"
        )
    for path in (
        raw_manifest,
        v2_viralscan,
        v3_viralscan,
        index,
        t2g,
        whitelist,
    ):
        if not path.is_file():
            raise FreshControlPreparationError(f"missing required file: {path}")
    manifest_paths = [
        task_manifest,
        task_manifest.with_name(f"{task_manifest.stem}.small.tsv"),
        task_manifest.with_name(f"{task_manifest.stem}.large.tsv"),
        task_manifest.with_name(f"{task_manifest.stem}.highmem.tsv"),
    ]
    if any(path.exists() for path in manifest_paths):
        raise FreshControlPreparationError(
            f"fresh task manifest already exists below: {task_manifest.parent}"
        )
    with raw_manifest.open(newline="", encoding="utf-8") as handle:
        input_rows = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["sample_id"]: row for row in input_rows}
    if len(input_rows) != len(EXPECTED_IDS) or set(by_id) != EXPECTED_IDS:
        raise FreshControlPreparationError(
            "control input manifest must contain exactly the frozen five samples"
        )

    packet_root = task_manifest.parent.resolve()
    tasks: list[dict[str, object]] = []
    for sample_id in sorted(selected_ids):
        source = by_id[sample_id]
        read1 = Path(source["read1_path"])
        read2 = Path(source["read2_path"])
        if not read1.is_file() or not read2.is_file():
            raise FreshControlPreparationError(f"missing frozen FASTQ pair for {sample_id}")
        storage_bytes = {mate: int(source[f"read{mate}_storage_bytes"]) for mate in (1, 2)}
        storage_sha256 = {mate: source[f"read{mate}_storage_sha256"] for mate in (1, 2)}
        for mate, path in ((1, read1), (2, read2)):
            if path.stat().st_size != storage_bytes[mate]:
                raise FreshControlPreparationError(f"stored byte count drifted: {path}")
            if not re.fullmatch(r"[0-9a-f]{64}", storage_sha256[mate]):
                raise FreshControlPreparationError(
                    f"invalid storage SHA-256 for {sample_id} read {mate}"
                )
            # Deliberately size-only here. The recorded digest is re-derived from
            # the bytes in verify_frozen_fastq, immediately before the run reads
            # the file — the point where a post-audit swap actually matters.
            # Hashing again at manifest-build time would re-read a quarter of a
            # terabyte to close no additional window.
        total_bytes = sum(storage_bytes.values())
        tier = "large" if total_bytes >= LARGE_INPUT_BYTES else "small"
        for stack, executable in (
            ("v2", v2_viralscan),
            ("v3", v3_viralscan),
        ):
            if stack not in stacks:
                continue
            task_id = f"{stack}__{sample_id}"
            row_tier = "highmem" if task_id in highmem_task_ids else tier
            output = (output_root / stack / sample_id).resolve()
            if output.exists():
                raise FreshControlPreparationError(f"fresh output already exists: {output}")
            tasks.append(
                {
                    "task_id": task_id,
                    "sample_id": sample_id,
                    "stack": stack,
                    "tier": row_tier,
                    "attempt_id": attempt_id,
                    "viralscan_cache_path": (
                        str(viralscan_cache.resolve()) if viralscan_cache else EMPTY
                    ),
                    "viralscan_cache_manifest_sha256": viralscan_cache_manifest_sha256 or EMPTY,
                    # Only the v3 arm gets a GTF; v2 reads its own packaged data.
                    "gtf_path": str(gtf.resolve()) if gtf and stack == "v3" else EMPTY,
                    "gtf_sha256": gtf_sha256 if gtf and stack == "v3" else EMPTY,
                    "read1_path": str(read1.resolve()),
                    "read2_path": str(read2.resolve()),
                    "read1_storage_bytes": storage_bytes[1],
                    "read2_storage_bytes": storage_bytes[2],
                    "read1_storage_sha256": storage_sha256[1],
                    "read2_storage_sha256": storage_sha256[2],
                    "viralscan_path": str(executable.resolve()),
                    "output_path": str(output),
                    "index_path": str(index.resolve()),
                    "t2g_path": str(t2g.resolve()),
                    "whitelist_path": str(whitelist.resolve()),
                    "cores": cores,
                    "status_path": str(packet_root / "status" / f"{task_id}.json"),
                    "stdout_path": str(packet_root / "logs" / f"{task_id}.out"),
                    "stderr_path": str(packet_root / "logs" / f"{task_id}.err"),
                }
            )

    task_manifest.parent.mkdir(parents=True, exist_ok=True)
    _write_tasks(task_manifest, tasks)
    for tier, path in zip(TIERS, manifest_paths[1:], strict=True):
        _write_tasks(
            path,
            [task for task in tasks if task["tier"] == tier],
        )
    return tasks


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--v2-viralscan", type=Path, required=True)
    parser.add_argument("--v3-viralscan", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--t2g", type=Path, required=True)
    parser.add_argument("--whitelist", type=Path, required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--attempt-id", required=True)
    parser.add_argument("--viralscan-cache", type=Path, default=None)
    parser.add_argument("--viralscan-cache-manifest-sha256", default=None)
    parser.add_argument(
        "--v2-data-dir",
        type=Path,
        default=None,
        help="v2.2.0 packaged data dir; its sorted GTFs become the explicit v3 -gtf "
        "(written to <packet>/reference/v2_panel.gtf with a t2g parity report)",
    )
    parser.add_argument(
        "--highmem-task",
        action="append",
        default=[],
        help="task id (e.g. v2__SRR6825024) to place in the 384 GiB highmem tier",
    )
    parser.add_argument(
        "--stack",
        action="append",
        choices=("v2", "v3"),
        dest="stacks",
        help="restrict emitted rows to this stack; repeatable, defaults to both",
    )
    parser.add_argument(
        "--sample-id",
        action="append",
        dest="sample_ids",
        help="restrict emitted rows to this sample; repeatable, defaults to all five",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    gtf = gtf_sha = None
    if args.v2_data_dir is not None:
        reference = args.task_manifest.parent / "reference"
        gtf = reference / "v2_panel.gtf"
        gtf_sha = build_panel_gtf(args.v2_data_dir, gtf)
        parity = gtf_t2g_parity(gtf, args.t2g)
        (reference / "gtf_t2g_parity.json").write_text(
            json.dumps(parity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(parity, indent=2, sort_keys=True))
        if not parity["passed"]:
            raise FreshControlPreparationError(
                "GTF/t2g parity FAILED; packet not prepared (see reference/gtf_t2g_parity.json)"
            )
    prepare_tasks(
        raw_manifest=args.raw_manifest,
        output_root=args.output_root,
        task_manifest=args.task_manifest,
        v2_viralscan=args.v2_viralscan,
        v3_viralscan=args.v3_viralscan,
        index=args.index,
        t2g=args.t2g,
        whitelist=args.whitelist,
        cores=args.cores,
        attempt_id=args.attempt_id,
        viralscan_cache=args.viralscan_cache,
        viralscan_cache_manifest_sha256=args.viralscan_cache_manifest_sha256,
        gtf=gtf,
        gtf_sha256=gtf_sha,
        highmem_task_ids=tuple(args.highmem_task),
        stacks=tuple(args.stacks) if args.stacks else ("v2", "v3"),
        sample_ids=tuple(args.sample_ids) if args.sample_ids else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
