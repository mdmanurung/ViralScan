#!/usr/bin/env python3
"""Freeze and checksum the executable fresh-control runner packet."""

from __future__ import annotations

import argparse
import hashlib
from collections.abc import Sequence
from pathlib import Path


class FreshPacketFreezeError(RuntimeError):
    """Raised when the runner packet cannot be frozen without replacement."""


PACKET_INPUTS = (
    "run_task.sh",
    "tasks.tsv",
    "tasks.small.tsv",
    "tasks.large.tsv",
    "control_inputs.raw.tsv",
)
FROZEN_RUNNER = "source/run_fresh_control.py"
FROZEN_EVIDENCE = "source/run_fresh_control_evidence.py"
# Hashed when present (highmem tier manifest; explicit-GTF packets).
OPTIONAL_INPUTS = (
    "tasks.highmem.tsv",
    "reference/v2_panel.gtf",
    "reference/gtf_t2g_parity.json",
    "reference/panel_provenance.json",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def freeze_packet(
    *, packet_root: Path, runner_source: Path, evidence_source: Path | None = None
) -> Path:
    """Copy the runner and write a checksum manifest for every packet input."""

    packet_root = packet_root.resolve()
    runner_source = runner_source.resolve()
    manifest = packet_root / "packet.sha256"
    frozen_runner = packet_root / FROZEN_RUNNER
    if manifest.exists() or frozen_runner.exists():
        raise FreshPacketFreezeError(f"fresh packet snapshot already exists: {packet_root}")
    required = [packet_root / relative for relative in PACKET_INPUTS]
    for path in (*required, runner_source):
        if not path.is_file():
            raise FreshPacketFreezeError(f"missing packet input: {path}")

    frozen_runner.parent.mkdir(parents=True, exist_ok=False)
    runner_staging = frozen_runner.with_name(f"{frozen_runner.name}.tmp")
    runner_staging.write_bytes(runner_source.read_bytes())
    runner_staging.chmod(0o444)
    runner_staging.replace(frozen_runner)

    frozen_paths = [
        *required,
        *(packet_root / rel for rel in OPTIONAL_INPUTS if (packet_root / rel).is_file()),
        frozen_runner,
    ]
    if evidence_source is not None:
        evidence_source = evidence_source.resolve()
        if not evidence_source.is_file():
            raise FreshPacketFreezeError(f"missing packet input: {evidence_source}")
        frozen_evidence = packet_root / FROZEN_EVIDENCE
        frozen_evidence.write_bytes(evidence_source.read_bytes())
        frozen_evidence.chmod(0o444)
        frozen_paths.append(frozen_evidence)
    lines = [
        f"{_sha256(path)}  {path.relative_to(packet_root).as_posix()}\n"
        for path in sorted(frozen_paths)
    ]
    staging = manifest.with_name(f"{manifest.name}.tmp")
    staging.write_text("".join(lines), encoding="utf-8")
    staging.chmod(0o444)
    staging.replace(manifest)
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-root", type=Path, required=True)
    parser.add_argument("--runner-source", type=Path, required=True)
    parser.add_argument("--evidence-source", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    freeze_packet(
        packet_root=args.packet_root,
        runner_source=args.runner_source,
        evidence_source=args.evidence_source,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
