"""V3 run fingerprints and non-destructive output-directory handling."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any, Callable

from viralscan import __version__

RUN_MANIFEST = "run_manifest.json"


class RunSafetyError(RuntimeError):
    """Raised when output reuse would be unsafe or irreproducible."""


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def build_run_manifest(args: Any) -> dict[str, Any]:
    """Build a canonical manifest from scientific options and input bytes."""
    path_fields = (
        "sample1",
        "sample2",
        "index",
        "transcripts",
        "gtf",
        "fasta",
        "f1",
        "whitelist",
        "called_cells_file",
        "host_index",
    )
    input_fingerprints: dict[str, str] = {}
    for field in path_fields:
        raw = getattr(args, field, None)
        if not raw:
            continue
        for i, value in enumerate(str(raw).split(",")):
            path = Path(value).expanduser().resolve()
            if path.is_file():
                input_fingerprints[f"{field}:{i}"] = sha256_file(path)

    excluded = {
        "output",
        "yes",
        "resume",
        "overwrite",
        "verbose",
        "quiet",
        "_subcommand",
    }
    options = {
        key: value
        for key, value in sorted(vars(args).items())
        if key not in excluded and isinstance(value, (str, int, float, bool, type(None)))
    }
    reference_hashes = {
        key: value
        for key, value in input_fingerprints.items()
        if key.split(":", 1)[0] in {"index", "transcripts", "gtf", "fasta", "f1", "host_index"}
    }
    reference_fingerprint = hashlib.sha256(
        json.dumps(reference_hashes, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    payload = {
        "schema_version": "3.0.0",
        "viralscan_version": __version__,
        "quantification_unit": "bustools-resolved-cb-umi-molecule",
        "allocation_method": getattr(args, "multimap_method", None),
        "reference_fingerprint": reference_fingerprint,
        "cell_calling": {
            "method": getattr(args, "cell_calling", None),
            "external_file": getattr(args, "called_cells_file", None),
        },
        "input_fingerprints": input_fingerprints,
        "options": options,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["run_fingerprint"] = hashlib.sha256(canonical).hexdigest()
    return payload


def _write_manifest_atomic(output_dir: Path, manifest: dict[str, Any]) -> None:
    target = output_dir / RUN_MANIFEST
    staging = output_dir / f".{RUN_MANIFEST}.tmp"
    staging.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staging.replace(target)


def prepare_output_directory(
    output_dir: Path,
    manifest: dict[str, Any],
    *,
    resume: bool,
    overwrite: bool,
    yes: bool,
    confirm: Callable[[str], str] = input,
) -> str:
    """Create, resume, or explicitly replace one resolved output directory."""
    output_dir = output_dir.resolve()
    if resume and overwrite:
        raise RunSafetyError("--resume and --overwrite are mutually exclusive.")

    nonempty = output_dir.is_dir() and any(output_dir.iterdir())
    if not nonempty:
        output_dir.mkdir(parents=True, exist_ok=True)
        _write_manifest_atomic(output_dir, manifest)
        return "new"

    if resume:
        manifest_path = output_dir / RUN_MANIFEST
        if not manifest_path.is_file():
            raise RunSafetyError("Cannot resume: run_manifest.json is missing.")
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous.get("run_fingerprint") != manifest.get("run_fingerprint"):
            raise RunSafetyError("Cannot resume: run fingerprint does not match this invocation.")
        return "resume"

    if not overwrite:
        raise RunSafetyError(
            "Output directory is non-empty. Use --resume for an identical run or "
            "--overwrite to replace it."
        )
    if not yes:
        answer = (
            confirm(f"Overwrite all contents of {output_dir}? Type 'yes' to continue: ")
            .strip()
            .lower()
        )
        if answer != "yes":
            raise RunSafetyError("Overwrite cancelled; output directory was not changed.")

    for child in output_dir.iterdir():
        if child.is_dir() and not child.is_symlink():
            shutil.rmtree(child)
        else:
            child.unlink()
    _write_manifest_atomic(output_dir, manifest)
    return "overwrite"
