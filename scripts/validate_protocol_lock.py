#!/usr/bin/env python3
"""LOCK-02 software: validate a discovery-study ``protocol-v1.lock.json`` against its schema.

    python scripts/validate_protocol_lock.py protocol-v1.lock.json [--root DIR] [--allow-test-only]
    python scripts/validate_protocol_lock.py --template > protocol-v1.lock.json

The package ships the container (``schemas/v3/protocol_lock.schema.json``), never the values: the template is
all ``null``. LOCK-01 calibration writes the values. Checks, in order:

1. schema (Draft 2020-12; the lock is its own schema, not ``validation_protocol.schema.json``);
2. ``test_only: true`` (top level, or inside ``adjudication`` where ``adjudicate_candidates.py`` reads it) is
   refused without ``--allow-test-only``, and the two spellings must agree;
3. every ``X_sha256`` with a sibling ``X_path`` (``sha256`` with ``path``) is recomputed from the file; paths
   are optional, relative ones resolve against ``--root`` (default: the lock's directory);
4. any ``null`` leaf (or empty ``calibration_inputs``) makes the lock ``incomplete``, not ``valid``.

Exit: 0 valid, 1 invalid, 2 unreadable input or usage, 3 incomplete (schema-correct, values not yet written).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

try:
    from scripts.governance_utils import json_schema_errors, load_json, sha256_file
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from governance_utils import (  # type: ignore[no-redef]
        json_schema_errors,
        load_json,
        sha256_file,
    )

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEMA = REPO_ROOT / "schemas" / "v3" / "protocol_lock.schema.json"

EXIT_VALID, EXIT_INVALID, EXIT_USAGE, EXIT_INCOMPLETE = 0, 1, 2, 3


def template() -> dict[str, Any]:
    """Skeleton lock: fixed policy labels from the schema, every value ``null``."""
    return {
        "schema_version": "1.0.0",
        "protocol_version": None,
        "created_at": None,
        "status": "exploratory",
        "code": dict.fromkeys(
            ("repo_sha", "package_version", "python", "kb_python", "snakemake", "env_digest")
        ),
        "panel": {
            "gtf_corpus_manifest_sha256": None,
            "panel_manifest_sha256": None,
            "index_sha256": None,
            "t2g_sha256": None,
            "dlist": {"identity": None, "sha256": None},
            "human_reference": {"identity": None, "sha256": None},
            "ncbi_snapshot_date": None,
            "taxonomy_snapshot_date": None,
        },
        "package_validation_protocol_sha256": None,
        "arms": {
            "primary_pair": ["combined_off", "combined_artefact"],
            "two_step": {
                "policy": "only_for_libraries_with_l1_plus_candidate",
                "skipped_arm_policy": "log_reason_per_library",
            },
            "secondary": {
                "name": "max_tier_anellovirus",
                "label": "secondary_sensitivity",
                "changes_primary_panel": False,
            },
        },
        "strand_policy": {
            "mode": "fixed_unstranded",
            "calibration_rule": {
                "adopt_reverse_only_if": [
                    "lower_expected_false_l2_per_library",
                    "keeps_min_positive_viral_molecule_fraction",
                ],
                "min_positive_viral_molecule_fraction": None,
            },
            "reverse_chemistries": None,
        },
        "cell_calling": {
            "method": "emptydrops",
            "emptydrops_seed": None,
            "emptydrops_niters": None,
        },
        "adjudication": {
            "distributed_route": dict.fromkeys(
                ("min_unique_molecules", "min_positive_cells", "min_independent_loci")
            ),
            "max_sibling_allocated_fraction": None,
            "min_virus_best_fraction": None,
            "min_complex_body_fraction": None,
            "l3_min_programme_cells": None,
        },
        "false_discovery": {
            "expected_false_l2_per_library": {"target": None, "interval": None},
            "decoy": {"enabled": None, "seeds": None, "strata": None},
        },
        "sensitivity": {
            "minimum_detectable_burden_reference": None,
            "minimum_detectable_burden_sha256": None,
        },
        "calibration_inputs": [],
        "holdout_manifest_sha256": None,
    }


def _walk(node: Any, where: str = "") -> Iterator[tuple[str, Any]]:
    """Yield ``(dotted path, value)`` for every leaf (an empty list is a leaf)."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, f"{where}.{key}" if where else key)
    elif isinstance(node, list) and node:
        for i, value in enumerate(node):
            yield from _walk(value, f"{where}[{i}]")
    else:
        yield where, node


def _dicts(node: Any, where: str = "") -> Iterator[tuple[str, dict[str, Any]]]:
    if isinstance(node, dict):
        yield where, node
        for key, value in node.items():
            yield from _dicts(value, f"{where}.{key}" if where else key)
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _dicts(value, f"{where}[{i}]")


def sha_errors(doc: dict[str, Any], root: Path) -> list[str]:
    """Recompute each ``X_sha256`` / ``sha256`` that has a sibling ``X_path`` / ``path``."""
    errors: list[str] = []
    for where, d in _dicts(doc):
        for key, claimed in d.items():
            if key != "sha256" and not key.endswith("_sha256"):
                continue
            path_key = "path" if key == "sha256" else key[: -len("sha256")] + "path"
            if path_key not in d or not isinstance(claimed, str):
                continue
            label = f"{where}.{key}" if where else key
            target = Path(d[path_key])
            target = target if target.is_absolute() else root / target
            if not target.is_file():
                errors.append(f"{label}: {path_key} {target} is not a file")
            elif (observed := sha256_file(target)) != claimed:
                errors.append(
                    f"{label}: sha256 mismatch for {target}: lock {claimed}, observed {observed}"
                )
    return errors


def incomplete_fields(doc: dict[str, Any]) -> list[str]:
    """Dotted paths of null leaves and of an empty ``calibration_inputs``; decoy fields skipped when disabled."""
    decoy_off = (doc.get("false_discovery") or {}).get("decoy", {}).get("enabled") is False
    unset = [
        where
        for where, value in _walk(doc)
        if value is None
        and not (
            decoy_off and where in ("false_discovery.decoy.seeds", "false_discovery.decoy.strata")
        )
    ]
    return unset + ([] if doc.get("calibration_inputs") else ["calibration_inputs"])


def test_only_refusals(doc: dict[str, Any], allow: bool) -> list[str]:
    top = doc.get("test_only") is True
    adj = (doc.get("adjudication") or {}).get("test_only") is True
    errors = []
    if (top or adj) and not allow:
        errors.append("test_only lock refused (pass --allow-test-only for fixtures)")
    if top != adj and (top or adj):
        errors.append(
            "test_only must be set both at the top level and in adjudication (the adjudicator reads the latter)"
        )
    return errors


def check(
    doc: Any, root: Path, *, allow_test_only: bool = False, schema: Path = DEFAULT_SCHEMA
) -> tuple[list[str], list[str]]:
    """``(errors, incomplete_fields)``; the lock is valid only when both are empty."""
    errors = json_schema_errors(doc, schema)
    if errors:  # later checks assume the schema shape
        return errors, []
    return test_only_refusals(doc, allow_test_only) + sha_errors(doc, root), incomplete_fields(doc)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("lock", nargs="?", type=Path, help="protocol-v1.lock.json")
    parser.add_argument(
        "--template", action="store_true", help="print a skeleton with null values and exit"
    )
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument(
        "--root", type=Path, help="base for relative *_path entries (default: the lock's directory)"
    )
    parser.add_argument(
        "--allow-test-only", action="store_true", help="accept fixtures marked test_only"
    )
    args = parser.parse_args(argv)

    if args.template:
        print(json.dumps(template(), indent=2))
        return EXIT_VALID
    if args.lock is None:
        parser.error("a lock file or --template is required")
    try:
        doc = load_json(args.lock)
    except (OSError, ValueError) as exc:
        print(f"validate_protocol_lock: cannot read {args.lock}: {exc}", file=sys.stderr)
        return EXIT_USAGE

    errors, incomplete = check(
        doc,
        args.root or args.lock.resolve().parent,
        allow_test_only=args.allow_test_only,
        schema=args.schema,
    )
    if errors:
        print("invalid", file=sys.stderr)
        for line in errors:
            print(f"  {line}", file=sys.stderr)
        return EXIT_INVALID
    if incomplete:
        print(f"incomplete: {len(incomplete)} unset field(s)")
        for line in incomplete:
            print(f"  {line}")
        return EXIT_INCOMPLETE
    print("valid")
    return EXIT_VALID


if __name__ == "__main__":
    sys.exit(main())
