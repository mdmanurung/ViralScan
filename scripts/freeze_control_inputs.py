#!/usr/bin/env python3
"""Freeze one fail-closed five-control FASTQ manifest for both package stacks."""

from __future__ import annotations

import argparse
import csv
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any


class ControlInputError(RuntimeError):
    """Raised when the five-control input packet is incomplete or inconsistent."""


CONTROL_METADATA = {
    "SRR12682296": ("ebv_control", "EBV", "ena"),
    "SRR12682297": ("ebv_control", "EBV", "ena"),
    "SRR12682298": ("ebv_control", "EBV", "ena"),
    "SRR6825024": ("hiv_control", "HIV-1", "retained_archive"),
    "SRR6825025": ("hiv_control", "HIV-1", "retained_archive"),
}
EXPECTED_IDS = frozenset(CONTROL_METADATA)
TRACKED_FIELDS = (
    "sample_id",
    "sample_class",
    "expected_target",
    "chemistry",
    "source_kind",
    "pair_records",
    "read1_name",
    "read1_storage_bytes",
    "read1_storage_md5",
    "read1_storage_sha256",
    "read1_content_bytes",
    "read1_content_md5",
    "read1_content_sha256",
    "read1_min_sequence_bases",
    "read1_max_sequence_bases",
    "read2_name",
    "read2_storage_bytes",
    "read2_storage_md5",
    "read2_storage_sha256",
    "read2_content_bytes",
    "read2_content_md5",
    "read2_content_sha256",
    "read2_min_sequence_bases",
    "read2_max_sequence_bases",
)
RAW_FIELDS = (*TRACKED_FIELDS, "read1_path", "read2_path")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ControlInputError(f"missing TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_one_tsv(path: Path) -> dict[str, str]:
    rows = _read_tsv(path)
    if len(rows) != 1:
        raise ControlInputError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def _read_audit(path: Path, sample_id: str) -> dict[str, Any]:
    if not path.is_file():
        raise ControlInputError(f"missing pair audit: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "valid":
        raise ControlInputError(f"pair audit is not valid: {path}")
    if payload.get("sample_id") != sample_id:
        raise ControlInputError(f"pair-audit sample mismatch: {path}")
    if payload.get("chemistry") != "10xv2":
        raise ControlInputError(f"pair-audit chemistry mismatch: {path}")
    if payload.get("required_r1_bases") != 26:
        raise ControlInputError(f"pair-audit geometry mismatch: {path}")
    if not isinstance(payload.get("records"), int) or payload["records"] <= 0:
        raise ControlInputError(f"invalid paired-record count: {path}")
    for mate in ("read1", "read2"):
        data = payload.get(mate)
        if not isinstance(data, dict):
            raise ControlInputError(f"missing {mate} audit: {path}")
        for field in (
            "bytes",
            "stored_bytes",
            "md5",
            "sha256",
            "min_sequence_bases",
            "max_sequence_bases",
        ):
            if field not in data:
                raise ControlInputError(f"missing {mate}.{field}: {path}")
        if data["min_sequence_bases"] <= 0:
            raise ControlInputError(f"invalid {mate} sequence geometry: {path}")
        if len(data["md5"]) != 32 or len(data["sha256"]) != 64:
            raise ControlInputError(f"invalid {mate} content digest: {path}")
    if payload["read1"]["min_sequence_bases"] < 26:
        raise ControlInputError(f"R1 is incompatible with 10xv2 geometry: {path}")
    return payload


def _mate_fields(
    *,
    mate: int,
    path: Path,
    audit: dict[str, Any],
    storage_bytes: int,
    storage_md5: str,
    storage_sha256: str,
) -> dict[str, object]:
    key = f"read{mate}"
    data = audit[key]
    if not path.is_file():
        raise ControlInputError(f"missing control FASTQ: {path}")
    if path.stat().st_size != storage_bytes:
        raise ControlInputError(f"stored byte count drifted: {path}")
    if data["stored_bytes"] != storage_bytes:
        raise ControlInputError(f"pair-audit stored byte count mismatch: {path}")
    return {
        f"{key}_name": path.name,
        f"{key}_storage_bytes": storage_bytes,
        f"{key}_storage_md5": storage_md5,
        f"{key}_storage_sha256": storage_sha256,
        f"{key}_content_bytes": data["bytes"],
        f"{key}_content_md5": data["md5"],
        f"{key}_content_sha256": data["sha256"],
        f"{key}_min_sequence_bases": data["min_sequence_bases"],
        f"{key}_max_sequence_bases": data["max_sequence_bases"],
        f"{key}_path": str(path.resolve()),
    }


def _write_tsv(path: Path, rows: list[dict[str, object]], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.with_name(f"{path.name}.tmp")
    with staging.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    staging.replace(path)


def freeze_control_inputs(
    *,
    ena_root: Path,
    ebv_audit_root: Path,
    retained_manifest: Path,
    retained_audit_root: Path,
    raw_output: Path,
    tracked_output: Path,
) -> list[dict[str, object]]:
    """Validate and freeze identical inputs for the v2.2.0 and v3 control runs."""

    ena_rows = _read_tsv(ena_root / "manifest.tsv")
    ena_by_key = {(row["sample_id"], int(row["mate"])): row for row in ena_rows}
    expected_ena_keys = {
        (sample_id, mate)
        for sample_id in CONTROL_METADATA
        if CONTROL_METADATA[sample_id][2] == "ena"
        for mate in (1, 2)
    }
    if len(ena_rows) != len(expected_ena_keys) or set(ena_by_key) != expected_ena_keys:
        raise ControlInputError("ENA manifest does not contain the frozen six mates")

    retained_rows = _read_tsv(retained_manifest)
    retained_by_id = {row["sample_id"]: row for row in retained_rows}
    expected_retained_ids = {
        sample_id
        for sample_id in CONTROL_METADATA
        if CONTROL_METADATA[sample_id][2] == "retained_archive"
    }
    if (
        len(retained_rows) != len(expected_retained_ids)
        or set(retained_by_id) != expected_retained_ids
    ):
        raise ControlInputError("retained manifest does not contain the frozen two pairs")

    frozen: list[dict[str, object]] = []
    for sample_id in sorted(EXPECTED_IDS):
        sample_class, expected_target, source_kind = CONTROL_METADATA[sample_id]
        audit_root = ebv_audit_root if source_kind == "ena" else retained_audit_root
        audit = _read_audit(audit_root / f"{sample_id}.json", sample_id)
        row: dict[str, object] = {
            "sample_id": sample_id,
            "sample_class": sample_class,
            "expected_target": expected_target,
            "chemistry": "10xv2",
            "source_kind": source_kind,
            "pair_records": audit["records"],
        }
        for mate in (1, 2):
            if source_kind == "ena":
                declared = ena_by_key[(sample_id, mate)]
                status = _read_one_tsv(ena_root / "status" / f"{sample_id}_{mate}.tsv")
                if (
                    status.get("sample_id") != sample_id
                    or status.get("mate") != str(mate)
                    or status.get("status") != "valid"
                    or status.get("bytes") != declared["bytes"]
                    or status.get("md5") != declared["md5"]
                ):
                    raise ControlInputError(
                        f"ENA download status mismatch: {sample_id} mate {mate}"
                    )
                path = ena_root / "files" / f"{sample_id}_{mate}.fastq.gz"
                storage_bytes = int(status["bytes"])
                storage_md5 = status["md5"]
                storage_sha256 = status["sha256"]
            else:
                declared = retained_by_id[sample_id]
                path = Path(declared[f"read{mate}"])
                mate_audit = audit[f"read{mate}"]
                storage_bytes = int(mate_audit["stored_bytes"])
                storage_md5 = str(mate_audit["md5"])
                storage_sha256 = str(mate_audit["sha256"])
            row.update(
                _mate_fields(
                    mate=mate,
                    path=path,
                    audit=audit,
                    storage_bytes=storage_bytes,
                    storage_md5=storage_md5,
                    storage_sha256=storage_sha256,
                )
            )
        frozen.append(row)

    _write_tsv(raw_output, frozen, RAW_FIELDS)
    _write_tsv(tracked_output, frozen, TRACKED_FIELDS)
    tracked_text = tracked_output.read_text(encoding="utf-8")
    if any(
        row[path_key] in tracked_text for row in frozen for path_key in ("read1_path", "read2_path")
    ):
        raise ControlInputError("tracked manifest leaked an input path")
    return frozen


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ena-root", type=Path, required=True)
    parser.add_argument("--ebv-audit-root", type=Path, required=True)
    parser.add_argument("--retained-manifest", type=Path, required=True)
    parser.add_argument("--retained-audit-root", type=Path, required=True)
    parser.add_argument("--raw-output", type=Path, required=True)
    parser.add_argument("--tracked-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    freeze_control_inputs(
        ena_root=args.ena_root,
        ebv_audit_root=args.ebv_audit_root,
        retained_manifest=args.retained_manifest,
        retained_audit_root=args.retained_audit_root,
        raw_output=args.raw_output,
        tracked_output=args.tracked_output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
