#!/usr/bin/env python3
"""Compare the scientific outputs of two successful legacy-v2/v3 attempts."""

from __future__ import annotations

import argparse
import csv
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any


class AttemptComparisonError(ValueError):
    """Raised when an attempt does not satisfy the comparison input contract."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AttemptComparisonError(f"cannot read {path.name}: {error}") from error
    if not isinstance(payload, dict):
        raise AttemptComparisonError(f"{path.name} is not a JSON object")
    return payload


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.with_suffix(path.suffix + ".tmp")
    staging.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staging.replace(path)


def _read_tsv(path: Path) -> dict[str, Any]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            rows = list(reader)
    except OSError as error:
        raise AttemptComparisonError(f"cannot read {path.name}: {error}") from error
    if not fields:
        raise AttemptComparisonError(f"{path.name} has no header")
    return {"columns": fields, "rows": rows}


def _v3_dir(attempt_dir: Path) -> Path:
    attempt_dir = attempt_dir.expanduser().resolve()
    wrapper_status = _read_json(attempt_dir / "status.json")
    if wrapper_status.get("status") != "success":
        raise AttemptComparisonError(
            f"wrapper status is {wrapper_status.get('status')!r}, not success"
        )
    v3_dir = attempt_dir / "v3"
    v3_status = _read_json(v3_dir / "status.json")
    if v3_status.get("status") != "success":
        raise AttemptComparisonError(f"v3 status is {v3_status.get('status')!r}, not success")
    return v3_dir


def _scientific_provenance(v3_dir: Path) -> dict[str, Any]:
    hashes = _read_json(v3_dir / "hashes.json")
    result = _read_json(v3_dir / "result.json")
    status = _read_json(v3_dir / "status.json")
    before = _read_json(v3_dir / "input_fingerprints_before.json")
    after = _read_json(v3_dir / "input_fingerprints_after.json")
    parameter_hashes = {
        hashes.get("scientific_parameter_hash"),
        result.get("scientific_parameter_hash"),
        status.get("scientific_parameter_hash"),
        before.get("scientific_parameter_hash"),
        after.get("scientific_parameter_hash"),
    }
    if None in parameter_hashes or len(parameter_hashes) != 1:
        raise AttemptComparisonError(
            "scientific_parameter_hash is missing or internally inconsistent"
        )
    inputs = hashes.get("inputs")
    if (
        not isinstance(inputs, dict)
        or inputs != hashes.get("inputs_before")
        or inputs != hashes.get("inputs_after")
        or inputs != before.get("inputs")
        or inputs != after.get("inputs")
    ):
        raise AttemptComparisonError("input fingerprints are missing or internally inconsistent")
    code = hashes.get("code")
    grouping = hashes.get("virus_grouping_sha256")
    if not isinstance(code, dict) or not code or not isinstance(grouping, str):
        raise AttemptComparisonError("code or virus-grouping fingerprints are missing")
    return {
        "scientific_parameter_hash": next(iter(parameter_hashes)),
        "input_fingerprints": inputs,
        "code_fingerprints": code,
        "virus_grouping_sha256": grouping,
    }


def _matrix_mismatch(first: Any, second: Any) -> str | None:
    import numpy as np
    from scipy import sparse

    if first.shape != second.shape:
        return f"shape differs: {first.shape} != {second.shape}"
    first_values = first.data if sparse.issparse(first) else np.asarray(first)
    second_values = second.data if sparse.issparse(second) else np.asarray(second)
    if (
        not np.isfinite(first_values).all()
        or not np.isfinite(second_values).all()
        or (first_values < 0).any()
        or (second_values < 0).any()
    ):
        return "matrix contains non-finite or negative values"
    first_csr = sparse.csr_matrix(first)
    second_csr = sparse.csr_matrix(second)
    difference = first_csr - second_csr
    difference.eliminate_zeros()
    if difference.nnz:
        return f"{difference.nnz} matrix entries differ"
    return None


def _h5ad_mismatches(first_path: Path, second_path: Path) -> list[dict[str, str]]:
    import anndata as ad
    import pandas as pd

    first = ad.read_h5ad(first_path)
    second = ad.read_h5ad(second_path)
    mismatches: list[dict[str, str]] = []
    for name, first_frame, second_frame in (
        ("h5ad.obs", first.obs, second.obs),
        ("h5ad.var", first.var, second.var),
    ):
        try:
            pd.testing.assert_frame_equal(
                first_frame,
                second_frame,
                check_dtype=True,
                check_categorical=True,
                check_like=False,
            )
        except AssertionError as error:
            mismatches.append({"check": name, "message": str(error).splitlines()[0]})
    matrix_error = _matrix_mismatch(first.X, second.X)
    if matrix_error:
        mismatches.append({"check": "h5ad.X", "message": matrix_error})
    first_layers = set(first.layers)
    second_layers = set(second.layers)
    if first_layers != second_layers:
        mismatches.append(
            {
                "check": "h5ad.layers",
                "message": (
                    f"layer names differ: {sorted(first_layers)} != {sorted(second_layers)}"
                ),
            }
        )
    for layer in sorted(first_layers & second_layers):
        matrix_error = _matrix_mismatch(first.layers[layer], second.layers[layer])
        if matrix_error:
            mismatches.append({"check": f"h5ad.layers.{layer}", "message": matrix_error})
    return mismatches


def compare_attempts(
    first_attempt: Path,
    second_attempt: Path,
    report_path: Path,
) -> dict[str, Any]:
    """Compare prespecified deterministic scientific artifacts and write a report."""

    checks: list[dict[str, str]] = []
    mismatches: list[dict[str, str]] = []

    def record(check: str, first: Any, second: Any) -> None:
        if first == second:
            checks.append({"check": check, "status": "match"})
        else:
            row = {
                "check": check,
                "status": "mismatch",
                "message": "canonical scientific values differ",
            }
            checks.append(row)
            mismatches.append({"check": check, "message": row["message"]})

    try:
        first_v3 = _v3_dir(first_attempt)
        second_v3 = _v3_dir(second_attempt)
        first_provenance = _scientific_provenance(first_v3)
        second_provenance = _scientific_provenance(second_v3)
        for key in (
            "scientific_parameter_hash",
            "input_fingerprints",
            "code_fingerprints",
            "virus_grouping_sha256",
        ):
            record(key, first_provenance[key], second_provenance[key])
        record(
            "count_audit",
            _read_json(first_v3 / "count_audit.json"),
            _read_json(second_v3 / "count_audit.json"),
        )
        record(
            "per_virus",
            _read_tsv(first_v3 / "per_virus.tsv"),
            _read_tsv(second_v3 / "per_virus.tsv"),
        )
        record(
            "per_cell",
            _read_tsv(first_v3 / "per_cell.tsv"),
            _read_tsv(second_v3 / "per_cell.tsv"),
        )
        h5ad_mismatches = _h5ad_mismatches(first_v3 / "adata_v3.h5ad", second_v3 / "adata_v3.h5ad")
        h5ad_mismatch_names = {row["check"] for row in h5ad_mismatches}
        for name in ("h5ad.obs", "h5ad.var", "h5ad.X", "h5ad.layers"):
            layer_mismatch = name == "h5ad.layers" and any(
                mismatch_name.startswith("h5ad.layers") for mismatch_name in h5ad_mismatch_names
            )
            if name not in h5ad_mismatch_names and not layer_mismatch:
                checks.append({"check": name, "status": "match"})
        for row in h5ad_mismatches:
            checks.append({**row, "status": "mismatch"})
            mismatches.append(row)
    except Exception as error:
        mismatch = {"check": "attempt_contract", "message": str(error)}
        checks.append({**mismatch, "status": "mismatch"})
        mismatches.append(mismatch)

    report = {
        "schema_version": "1.0.0",
        "status": "mismatch" if mismatches else "match",
        "first_attempt": f"{first_attempt.parent.name}/{first_attempt.name}",
        "second_attempt": f"{second_attempt.parent.name}/{second_attempt.name}",
        "checks": checks,
        "mismatches": mismatches,
        "ignored_nondeterministic_fields": [
            "attempt_id",
            "timings",
            "resource_metrics",
            "scratch_paths",
            "raw_h5ad_file_bytes",
        ],
    }
    _write_json(report_path, report)
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first_attempt", type=Path)
    parser.add_argument("second_attempt", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = compare_attempts(args.first_attempt, args.second_attempt, args.report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "match" else 1


if __name__ == "__main__":
    raise SystemExit(main())
