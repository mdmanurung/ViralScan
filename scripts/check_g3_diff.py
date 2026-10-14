#!/usr/bin/env python3
"""Mechanical diff check for the G3 signing commit (SCI-05 r15 N6).

The G3 commit differs from the reviewed draft by construction: it fills factor levels and the generation
seed, flips six section statuses with their digests, records the freeze, and refreshes ``pre_freeze_runs``.
This check fails if ANYTHING else changed between the reviewed revision and the signing commit: a protocol
path outside ``ALLOWED_PROTOCOL_PATHS``, a removed blocker other than the three G3 resolves, a file outside
``ALLOWED_FILES``, or a schema, script, source or test edit. It does not sign or freeze anything.

    python scripts/check_g3_diff.py <reviewed-sha> <g3-commit-sha>
"""

from __future__ import annotations

import re
import subprocess
import sys
from collections.abc import Sequence
from typing import Any

import yaml

PROTOCOL = "analysis/v3_validation/protocol.yaml"
SECTIONS = (
    "partitions|calibration|defaults_selection|workflow_matrix|failure_and_deviation_reporting"
)
RESOLVED_BLOCKERS = "partitions_metrics|independent_review_round_2|independent_review"
ALLOWED_PROTOCOL_PATHS = [
    r"factors\[[^\]]+\]\.(levels|status|pending_reason)",
    r"seeds\.generation\.(status|value|pending_reason)",
    rf"({SECTIONS})\.(status|frozen_at|contract_sha256)",
    r"planned_freeze_sections\.(partitions|calibration_metrics_lod|workflow_matrix"
    r"|failure_and_deviation_reporting|independent_review)\.status",
    r"freeze\.(freeze_manifest|frozen_at)",
    r"pre_freeze_runs\.(inventory|inventory_sha256)([.\[].*)?",
    r"frozen_inputs\.(contract_sha256|recorded_deviations)",
    rf"execution_readiness\.training_blockers\[({RESOLVED_BLOCKERS})\].*",
]
ALLOWED_FILES = {
    PROTOCOL,
    "analysis/v3_validation/deviations.yaml",  # append-only: check_ledger_append_only.py
    "analysis/v3_validation/freeze_sidecar.yaml",
    "analysis/v3_artifact_inventory.tsv",
    "claims/registry.json",
    "PLAN.md",
    ".living/decisions.md",
}


def flatten(node: Any, prefix: str = "") -> dict[str, Any]:
    """Leaf paths of a document; lists of ``id``-keyed mappings are addressed ``[id]``, others by index."""
    if isinstance(node, dict):
        return {
            k: v
            for key, val in node.items()
            for k, v in flatten(val, f"{prefix}.{key}" if prefix else str(key)).items()
        }
    if isinstance(node, list) and node and all(isinstance(i, dict) for i in node):
        keyed = all("id" in i for i in node)
        return {
            k: v
            for n, item in enumerate(node)
            for k, v in flatten(item, f"{prefix}[{item['id'] if keyed else n}]").items()
        }
    return {prefix: node}


def protocol_violations(base: dict[str, Any], head: dict[str, Any]) -> list[str]:
    """Paths that differ between two protocol documents and are not an allowed G3 edit."""
    before, after = flatten(base), flatten(head)
    allowed = [re.compile(p) for p in ALLOWED_PROTOCOL_PATHS]
    bad = [
        path
        for path in sorted(set(before) | set(after))
        if before.get(path, object()) != after.get(path, object())
        and not any(rx.fullmatch(path) for rx in allowed)
    ]
    old = (base.get("frozen_inputs") or {}).get("recorded_deviations") or []
    new = (head.get("frozen_inputs") or {}).get("recorded_deviations") or []
    if list(new[: len(old)]) != list(old):
        bad.append("frozen_inputs.recorded_deviations (not append-only)")
    return bad


def file_violations(names: Sequence[str]) -> list[str]:
    return [n for n in names if n not in ALLOWED_FILES]


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout


def main(argv: Sequence[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    base, head = argv
    names = _git("diff", "--name-only", f"{base}..{head}").split()
    bad = [f"file not allowed at G3: {n}" for n in file_violations(names)]
    if PROTOCOL in names:
        docs = [yaml.safe_load(_git("show", f"{rev}:{PROTOCOL}")) for rev in (base, head)]
        bad += [
            f"protocol path changed outside the G3 allowlist: {p}"
            for p in protocol_violations(*docs)
        ]
    print("\n".join(bad) if bad else f"ok: only the G3 freeze edits differ ({len(names)} files)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
