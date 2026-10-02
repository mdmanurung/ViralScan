"""V3 run fingerprints and non-destructive output-directory handling."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from viralscan import __version__

RUN_MANIFEST = "run_manifest.json"
RUN_COMPLETE = "run_complete.json"
#: Per-sample artifacts whose sha256 the completion marker records (each only if
#: present). Deliberately narrow: the headline tables and the count matrix, not
#: logs, plots, or every intermediate file.
MARKER_ARTIFACTS = (
    "results/viral_summary.tsv",
    "results/virus_identity.tsv",
    "results/multimap_evidence.tsv",
    "kb-python/counts_unfiltered/adata_multimap.h5ad",
)


class RunSafetyError(RuntimeError):
    """Raised when output reuse would be unsafe or irreproducible."""


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
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
    # Options added after v3.0 manifests were first written. An unset (None) value
    # is omitted, so a manifest that predates the option still matches an
    # invocation that does not use it; any non-None value changes the fingerprint
    # and refuses --resume against such a manifest (old counts were made with
    # kb's default, not this value).
    omit_when_unset = {"strand"}
    options = {
        key: value
        for key, value in sorted(vars(args).items())
        if key not in excluded
        and isinstance(value, (str, int, float, bool, type(None)))
        and not (key in omit_when_unset and value is None)
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
    # Added after hashing so --resume still matches manifests written without it.
    payload["completion_marker"] = True
    return payload


def build_run_complete(run_root: Path) -> dict[str, Any]:
    """Describe a finished run: its fingerprint and the hashes of MARKER_ARTIFACTS."""
    run_root = Path(run_root)
    manifest = json.loads((run_root / RUN_MANIFEST).read_text(encoding="utf-8"))
    samples = sorted(p.parent.name for p in run_root.glob("*/config.yaml"))
    artifacts = {
        f"{sample}/{rel}": sha256_file(run_root / sample / rel)
        for sample in samples
        for rel in MARKER_ARTIFACTS
        if (run_root / sample / rel).is_file()
    }
    return {
        "schema_version": "3.0.0",
        "run_fingerprint": manifest.get("run_fingerprint"),
        "viralscan_version": __version__,
        "samples": samples,
        "artifacts": artifacts,
        "completed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def write_run_complete(run_root: Path) -> None:
    """Atomically (re)write ``run_complete.json`` at the run root."""
    run_root = Path(run_root)
    staging = run_root / f".{RUN_COMPLETE}.tmp"
    staging.write_text(
        json.dumps(build_run_complete(run_root), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    staging.replace(run_root / RUN_COMPLETE)


def clear_run_complete(run_root: Path) -> None:
    (Path(run_root) / RUN_COMPLETE).unlink(missing_ok=True)


def restamp_run_complete(run_root: Path) -> bool:
    """Recompute the marker after an in-place mutation; no-op if the run had none."""
    if not (Path(run_root) / RUN_COMPLETE).is_file():
        return False
    write_run_complete(run_root)
    return True


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
            old_strand = (previous.get("options") or {}).get("strand")
            new_strand = (manifest.get("options") or {}).get("strand")
            reason = ""
            if old_strand != new_strand:
                reason = (
                    f" (--strand differs: previous run used {old_strand or 'the kb default'!r}, "
                    f"this invocation uses {new_strand or 'the kb default'!r}; "
                    "existing counts are not reusable)"
                )
            raise RunSafetyError(
                f"Cannot resume: run fingerprint does not match this invocation{reason}."
            )
        # A resumed run is in progress again: the old marker no longer vouches
        # for it. (New and overwrite starts have no marker: empty dir / wiped.)
        clear_run_complete(output_dir)
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
