#!/usr/bin/env python3
"""Re-validate completed fresh-control outputs without rewriting attempt evidence.

Attempt 2 recorded four legacy v2.2.0 rows as failed with exit code 65 because
``_v2_artifact_errors`` checked a flat output layout while the legacy CLI nests
its result tree one level beneath the sample identifier. Those workflows had in
fact exited zero and written complete trees.

This tool re-runs the corrected artifact check against an already-terminal
output tree and writes a *superseding* record that cites, but never replaces,
the original status file. Original attempt evidence stays immutable.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from scripts.run_fresh_control import _v2_artifact_errors

REQUIRED_FAILURE_FIELDS = (
    "run_id",
    "stage",
    "command",
    "exit_code",
    "stderr_path",
    "attempt_id",
    "scientific_parameter_hash",
)


class RevalidationError(RuntimeError):
    """Raised when a re-validation would overwrite or fabricate evidence."""


def revalidate(
    *,
    status_path: Path,
    record_path: Path,
    attempt_id: str,
) -> dict[str, Any]:
    """Re-check one terminal v2 output tree and write a superseding record."""

    if not status_path.is_file():
        raise RevalidationError(f"missing original status record: {status_path}")
    if record_path.exists():
        raise RevalidationError(f"revalidation record already exists: {record_path}")
    original = json.loads(status_path.read_text(encoding="utf-8"))
    if original.get("stack") != "v2":
        raise RevalidationError(f"only v2 rows are revalidated here: {status_path}")
    if original.get("workflow_exit_code") != 0:
        raise RevalidationError(
            f"workflow did not succeed; nothing to revalidate: {status_path}"
        )
    output = Path(original["output"])
    resolved_root, errors = _v2_artifact_errors(output)
    exit_code = 0 if not errors else 65
    payload: dict[str, Any] = {
        "schema_version": "1.1.0",
        "record_type": "revalidation",
        "supersedes": str(status_path.resolve()),
        "sample_id": original["sample_id"],
        "run_id": original["sample_id"],
        "stack": "v2",
        "stage": "artifact_validation",
        "attempt_id": attempt_id,
        "status": "success" if exit_code == 0 else "failed",
        "exit_code": exit_code,
        "workflow_exit_code": original["workflow_exit_code"],
        "original_exit_code": original["exit_code"],
        "original_artifact_validation_errors": original.get("artifact_validation_errors", []),
        "artifact_validation_errors": errors,
        "v2_output_root": str(resolved_root.resolve()),
        "command": original["command"],
        "command_sha256": original["command_sha256"],
        "scientific_parameter_hash": original["command_sha256"],
        "output": original["output"],
        "stdout": original["stdout"],
        "stderr": original["stderr"],
        "stderr_path": original.get("stderr_path", original["stderr"]),
    }
    missing = [field for field in REQUIRED_FAILURE_FIELDS if field not in payload]
    if missing:
        raise RevalidationError(f"record omits required fields: {', '.join(missing)}")
    record_path.parent.mkdir(parents=True, exist_ok=True)
    staging = record_path.with_name(f"{record_path.name}.tmp")
    staging.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staging.replace(record_path)
    return payload


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--attempt-id", default="attempt2-revalidated")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    payload = revalidate(
        status_path=args.status,
        record_path=args.record,
        attempt_id=args.attempt_id,
    )
    return int(payload["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
