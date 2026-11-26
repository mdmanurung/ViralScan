"""File hashing and run-manifest writes, kept below every module that records provenance.

``run_safety`` re-exports these names. They live here so ``virus_identity`` and other core
modules can import them at module level without a cycle through ``run_safety`` (LAYERS-01).
This module imports only the standard library.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

RUN_MANIFEST = "run_manifest.json"


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_path(path: str | Path) -> str:
    """sha256 of a file, or of a directory (sorted relative paths + file hashes)."""
    path = Path(path)
    if not path.is_dir():
        return sha256_file(path)
    digest = hashlib.sha256()
    for child in sorted(p for p in path.rglob("*") if p.is_file()):
        digest.update(child.relative_to(path).as_posix().encode() + b"\0")
        digest.update(bytes.fromhex(sha256_file(child)))
    return digest.hexdigest()


def write_manifest_atomic(output_dir: Path, manifest: dict[str, Any]) -> None:
    target = output_dir / RUN_MANIFEST
    staging = output_dir / f".{RUN_MANIFEST}.tmp"
    staging.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staging.replace(target)


def record_manifest_block(output_dir: Path, key: str, sample: str, block: dict[str, Any]) -> None:
    """Add a sample's ``key`` block to the manifest, outside the run fingerprint."""
    output_dir = Path(output_dir)
    manifest = json.loads((output_dir / RUN_MANIFEST).read_text(encoding="utf-8"))
    manifest.setdefault(key, {})[sample] = block
    write_manifest_atomic(output_dir, manifest)
