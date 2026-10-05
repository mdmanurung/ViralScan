#!/usr/bin/env python3
"""Freeze read-only v2, reference, and current-repository provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import re
import subprocess
from collections.abc import Sequence
from email.parser import Parser
from pathlib import Path


class ProvenanceError(RuntimeError):
    """Raised when required provenance cannot be frozen safely."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise ProvenanceError(f"missing required {label}: {path}")
    return path


def _one(paths: list[Path], label: str) -> Path:
    if len(paths) != 1:
        raise ProvenanceError(f"expected exactly one {label}, found {len(paths)}")
    return paths[0]


def _distribution(site_packages: Path, pattern: str, expected_name: str) -> tuple[str, Path]:
    metadata = _one(sorted(site_packages.glob(f"{pattern}.dist-info/METADATA")), expected_name)
    parsed = Parser().parsestr(metadata.read_text(encoding="utf-8"))
    observed_name = re.sub(r"[-_.]+", "-", parsed.get("Name", "").lower())
    normalized_expected = re.sub(r"[-_.]+", "-", expected_name.lower())
    if observed_name != normalized_expected:
        raise ProvenanceError(f"unexpected package name in {metadata}")
    version = parsed.get("Version")
    if not version:
        raise ProvenanceError(f"missing package version in {metadata}")
    return version, metadata


def _run(command: list[str], label: str, *, cwd: Path | None = None) -> str:
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ProvenanceError(f"failed to capture {label}: {exc}") from exc
    output = (result.stdout or result.stderr).strip()
    if not output:
        raise ProvenanceError(f"empty version output for {label}")
    return output


def _version(output: str, label: str) -> str:
    match = re.search(r"\bversion[,\s]+([0-9][0-9A-Za-z.+_-]*)", output, re.IGNORECASE)
    if match is None:
        raise ProvenanceError(f"could not parse {label} version from: {output!r}")
    return match.group(1)


def _file_row(scope: str, root: Path, path: Path) -> dict[str, object]:
    return {
        "scope": scope,
        "relative_path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _tree_sha256(rows: list[dict[str, object]]) -> str:
    digest = hashlib.sha256()
    for row in sorted(rows, key=lambda item: str(item["relative_path"])):
        digest.update((f"{row['relative_path']}\0{row['size_bytes']}\0{row['sha256']}\n").encode())
    return digest.hexdigest()


def _git(repo_root: Path, *args: str, binary: bool = False) -> str | bytes:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), *args],
            check=True,
            capture_output=True,
            text=not binary,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ProvenanceError(f"failed git {' '.join(args)}: {exc}") from exc
    return result.stdout


def collect_provenance(
    *,
    v2_prefix: Path,
    reference_dir: Path,
    repo_root: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Collect provenance without writing to any input tree."""

    for path, label in (
        (v2_prefix, "v2 prefix"),
        (reference_dir, "reference directory"),
        (repo_root, "repository root"),
    ):
        if not path.is_dir():
            raise ProvenanceError(f"missing required {label}: {path}")

    site_packages = _one(
        sorted(
            {
                path.resolve()
                for path in (v2_prefix / "lib").glob("python*/site-packages")
                if (path / "viralscan").is_dir() and (path / "kb_python").is_dir()
            }
        ),
        "v2 site-packages directory",
    )
    viralscan_dir = site_packages / "viralscan"
    kb_python_dir = site_packages / "kb_python"
    viralscan_version, viralscan_metadata = _distribution(site_packages, "viralscan-*", "ViralScan")
    kb_version, kb_metadata = _distribution(site_packages, "kb_python-*", "kb-python")

    python_bin = _required_file(v2_prefix / "bin/python", "v2 Python")
    os_name = platform.system().lower()
    kallisto = _required_file(
        kb_python_dir / f"bins/{os_name}/kallisto/kallisto",
        "embedded kallisto",
    )
    bustools = _required_file(
        kb_python_dir / f"bins/{os_name}/bustools/bustools",
        "embedded bustools",
    )
    python_version = _run([str(python_bin), "--version"], "v2 Python")
    kallisto_version = _version(
        _run([str(kallisto), "version"], "embedded kallisto"),
        "kallisto",
    )
    bustools_version = _version(
        _run([str(bustools), "version"], "embedded bustools"),
        "bustools",
    )

    code_paths = sorted(
        path
        for path in viralscan_dir.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and (path.suffix in {".py", ".R", ".j2"} or path.name == "Snakefile")
    )
    if not code_paths:
        raise ProvenanceError(f"no installed ViralScan code files found in {viralscan_dir}")

    file_rows = [_file_row("v2-viralscan-code", viralscan_dir, path) for path in code_paths]
    metadata_rows = [
        _file_row("v2-package-metadata", site_packages, viralscan_metadata),
        _file_row("v2-package-metadata", site_packages, kb_metadata),
    ]
    binary_rows = [
        _file_row("v2-embedded-binary", site_packages, kallisto),
        _file_row("v2-embedded-binary", site_packages, bustools),
    ]
    file_rows.extend(metadata_rows)
    file_rows.extend(binary_rows)

    reference_names = (
        "index_serratus.idx",
        "t2g_serratus.txt",
        "transcriptome.fa",
    )
    reference_rows = [
        _file_row(
            "original-reference",
            reference_dir,
            _required_file(reference_dir / name, f"reference file {name}"),
        )
        for name in reference_names
    ]
    file_rows.extend(reference_rows)

    head = str(_git(repo_root, "rev-parse", "HEAD")).strip()
    if not re.fullmatch(r"[0-9a-f]{40}", head):
        raise ProvenanceError(f"invalid repository HEAD: {head!r}")
    diff = _git(repo_root, "diff", "--binary", "HEAD", "--", binary=True)
    assert isinstance(diff, bytes)
    untracked_output = _git(
        repo_root,
        "ls-files",
        "--others",
        "--exclude-standard",
        "-z",
        binary=True,
    )
    assert isinstance(untracked_output, bytes)
    untracked_paths = sorted(
        path.decode("utf-8", errors="strict") for path in untracked_output.split(b"\0") if path
    )
    combined_diff = hashlib.sha256()
    combined_diff.update(diff)
    combined_diff.update(b"\0UNTRACKED\0")
    untracked_records = []
    for relative_path in untracked_paths:
        path = _required_file(
            repo_root / relative_path, f"untracked repository file {relative_path}"
        )
        record = {
            "relative_path": relative_path,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        untracked_records.append(record)
        combined_diff.update(
            (f"{record['relative_path']}\0{record['size_bytes']}\0{record['sha256']}\n").encode()
        )
    status = str(_git(repo_root, "status", "--short", "--untracked-files=all"))

    reference_files = {
        str(row["relative_path"]): {
            "size_bytes": row["size_bytes"],
            "sha256": row["sha256"],
        }
        for row in reference_rows
    }
    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "v2": {
            "prefix": str(v2_prefix.resolve()),
            "site_packages": str(site_packages.resolve()),
            "python_version": python_version.removeprefix("Python ").strip(),
            "packages": {
                "ViralScan": viralscan_version,
                "kb-python": kb_version,
            },
            "tools": {
                "kallisto": kallisto_version,
                "bustools": bustools_version,
            },
            "viralscan_code_tree_sha256": _tree_sha256(
                [row for row in file_rows if row["scope"] == "v2-viralscan-code"]
            ),
        },
        "reference": {
            "root": str(reference_dir.resolve()),
            "files": reference_files,
            "tree_sha256": _tree_sha256(reference_rows),
        },
        "repository": {
            "root": str(repo_root.resolve()),
            "head": head,
            "diff_sha256": combined_diff.hexdigest(),
            "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
            "untracked_files": untracked_records,
            "status_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
            "dirty": bool(status),
        },
    }
    return payload, sorted(
        file_rows,
        key=lambda row: (str(row["scope"]), str(row["relative_path"])),
    )


def write_provenance(
    *,
    payload: dict[str, object],
    file_rows: list[dict[str, object]],
    output_dir: Path,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=("scope", "relative_path", "size_bytes", "sha256"),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(file_rows)
    tsv_text = buffer.getvalue()
    payload = dict(payload)
    payload["files_tsv_sha256"] = hashlib.sha256(tsv_text.encode("utf-8")).hexdigest()

    tsv_path = output_dir / "provenance_files.tsv"
    json_path = output_dir / "provenance.json"
    tsv_path.write_text(tsv_text, encoding="utf-8")
    json_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return json_path, tsv_path


def freeze_provenance(
    *,
    v2_prefix: Path,
    reference_dir: Path,
    repo_root: Path,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Collect and write provenance below the repository's ignored run tree."""

    resolved_repo = repo_root.resolve()
    resolved_output = output_dir.resolve()
    try:
        relative_output = resolved_output.relative_to(resolved_repo)
    except ValueError as exc:
        raise ProvenanceError("output directory must be inside the repository") from exc
    if not relative_output.parts or relative_output.parts[0] != "benchmark_runs":
        raise ProvenanceError("output directory must be below benchmark_runs/")
    ignored = subprocess.run(
        ["git", "-C", str(resolved_repo), "check-ignore", "-q", str(relative_output)],
        check=False,
    )
    if ignored.returncode != 0:
        raise ProvenanceError("output directory is not ignored by Git")

    payload, rows = collect_provenance(
        v2_prefix=v2_prefix,
        reference_dir=reference_dir,
        repo_root=repo_root,
    )
    return write_provenance(payload=payload, file_rows=rows, output_dir=output_dir)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v2-prefix", required=True, type=Path)
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        json_path, tsv_path = freeze_provenance(
            v2_prefix=args.v2_prefix,
            reference_dir=args.reference_dir,
            repo_root=args.repo_root,
            output_dir=args.output_dir,
        )
    except ProvenanceError as exc:
        raise SystemExit(f"provenance freeze failed: {exc}") from exc
    print(json_path)
    print(tsv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
