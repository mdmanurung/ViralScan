#!/usr/bin/env python3
"""Post-run exact-read evidence for one fresh v3 control (LVC-12 checklist).

Runs ``viralscan evidence`` (exact CB-UMI read replay and extraction) for the
expected target and for the largest non-target candidate in
``results/viral_summary.tsv``. Read extraction only: no alignment FASTA is given,
so no new reference enters the packet. Every step's outcome, including failures,
is retained in ``<status-dir>/<task_id>.evidence.json``; nothing is overwritten.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any


class FreshEvidenceError(RuntimeError):
    """Raised when evidence cannot be attempted without replacing prior evidence."""


def sample_root(output: Path) -> Path:
    """Directory holding config.yaml; the CLI nests results one level below ``-o``."""

    if (output / "config.yaml").is_file():
        return output
    children = sorted(c for c in output.iterdir() if (c / "config.yaml").is_file())
    if len(children) != 1:
        raise FreshEvidenceError(f"no unique config.yaml root beneath {output}")
    return children[0]


def largest_non_target(summary: Path, target: str, aliases: Sequence[str] = ()) -> str | None:
    """Return the highest-molecule virus whose name is not the target (or an alias)."""

    excluded = {target.casefold(), *(a.casefold() for a in aliases)}
    with summary.open(newline="", encoding="utf-8") as handle:
        rows = [
            (float(row["viral_molecules_total_est"]), row["virus_name"])
            for row in csv.DictReader(handle, delimiter="\t")
            if row["virus_name"].casefold() not in excluded
        ]
    positive = [row for row in rows if row[0] > 0]
    return max(positive)[1] if positive else None


def _resolve_target_label(run_dir: Path, target: str) -> str:
    """Canonical label of the target via the installed viralscan resolver, else as given."""

    try:
        from viralscan.evidence import resolve_viral_target

        ids = [line.strip() for line in (run_dir / "log" / "analysis.txt").open()]
        return resolve_viral_target(target, ids)[0]
    except Exception:  # noqa: BLE001 - resolver absence/failure must not stop the step
        return target


def run_evidence_steps(
    *,
    task_id: str,
    viralscan: Path,
    run_dir: Path,
    expected_target: str,
    status_dir: Path,
    cores: int,
) -> dict[str, Any]:
    record_path = status_dir / f"{task_id}.evidence.json"
    if record_path.exists():
        raise FreshEvidenceError(f"evidence record already exists: {record_path}")
    run_dir = sample_root(run_dir)
    steps: list[dict[str, Any]] = []
    summary = run_dir / "results" / "viral_summary.tsv"
    candidates: list[tuple[str, str | None]] = [("expected_target", expected_target)]
    if summary.is_file():
        label = _resolve_target_label(run_dir, expected_target)
        candidates.append(("largest_non_target", largest_non_target(summary, label, [label])))
    else:
        candidates.append(("largest_non_target", None))
    for role, virus in candidates:
        step: dict[str, Any] = {"role": role, "virus": virus}
        if virus is None:
            step.update(status="skipped", reason="no positive non-target candidate or no summary")
            steps.append(step)
            continue
        out = run_dir / "evidence" / role
        command = [
            str(viralscan),
            "evidence",
            "--run-dir",
            str(run_dir),
            "--output",
            str(out),
            "--virus",
            virus,
            "--cores",
            str(cores),
        ]
        if out.exists():
            step.update(status="failed", exit_code=None, reason=f"output exists: {out}")
        else:
            code = subprocess.run(command, check=False).returncode
            step.update(status="success" if code == 0 else "failed", exit_code=code)
        step["command"] = command
        steps.append(step)
    payload = {
        "schema_version": "1.0.0",
        "record_type": "exact_read_evidence",
        "task_id": task_id,
        "run_dir": str(run_dir),
        "mode": "read-extraction-only",
        "steps": steps,
        "status": "success" if all(s["status"] != "failed" for s in steps) else "failed",
    }
    status_dir.mkdir(parents=True, exist_ok=True)
    record_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--viralscan", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--expected-target", required=True)
    parser.add_argument("--status-dir", type=Path, required=True)
    parser.add_argument("--cores", type=int, default=8)
    args = parser.parse_args(argv)
    payload = run_evidence_steps(
        task_id=args.task_id,
        viralscan=args.viralscan,
        run_dir=args.run_dir,
        expected_target=args.expected_target,
        status_dir=args.status_dir,
        cores=args.cores,
    )
    return 0 if payload["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
