#!/usr/bin/env python3
"""DSR-16: append-only ledger of exploratory (round 2) runs.

    explore_ledger.py append --dataset D --sample S --arm A --code-sha SHA --reference-id R
    explore_ledger.py check  [--deny-list FILE] [--expect-head HASH]

Each row carries the SHA-256 of the previous row's hash plus its own fields, so editing or
removing a row anywhere but the tail breaks the chain. Dropping the *last* rows is only caught
when the caller remembers the head hash (`--expect-head`; `append` prints it). `check` also fails
when any field names an accession in the Track B deny-list; `append` refuses such a row.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

FIELDS = ["timestamp", "dataset", "sample", "arm", "code_sha", "reference_id", "label"]
HEADER = "\t".join([*FIELDS, "prev_hash", "row_hash"])
GENESIS = "0" * 64
REPO = Path(__file__).resolve().parents[1]
DEFAULT_DENY = REPO / "config" / "track_b_denylist.txt"
DEFAULT_ROUND = "/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round2"


def default_ledger() -> Path:
    return Path(os.environ.get("ROUND_DIR", DEFAULT_ROUND)) / "explore_ledger.tsv"


def row_hash(prev: str, fields: list[str]) -> str:
    return hashlib.sha256("\t".join([prev, *fields]).encode()).hexdigest()


def read_deny(path: Path) -> list[str]:
    lines = (ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines())
    return [ln.upper() for ln in lines if ln]


def denied(fields: list[str], deny: list[str]) -> list[str]:
    text = "\t".join(fields).upper()
    return [d for d in deny if d in text]


def read_rows(ledger: Path) -> list[list[str]]:
    with ledger.open() as fh:
        fcntl.flock(fh, fcntl.LOCK_SH)
        text = fh.read()
    lines = text.splitlines()
    if not text.endswith("\n"):
        raise ValueError(f"{ledger}: missing final newline")
    if not lines or lines[0] != HEADER:
        raise ValueError(f"{ledger}: unexpected header")
    return [ln.split("\t") for ln in lines[1:]]


def verify(rows: list[list[str]]) -> list[str]:
    """Chain problems, one message per broken row."""
    errors, prev = [], GENESIS
    for i, row in enumerate(rows, start=1):
        if len(row) != len(FIELDS) + 2:
            errors.append(f"row {i}: {len(row)} columns")
        elif row[-2] != prev or row[-1] != row_hash(prev, row[: len(FIELDS)]):
            errors.append(f"row {i}: hash chain broken (edited or removed row)")
        prev = row[-1]
    return errors


def append(ledger: Path, record: dict[str, str], deny: list[str]) -> str:
    fields = [record[k] for k in FIELDS]
    if any(c in value for value in fields for c in "\t\r\n"):
        raise ValueError("ledger fields must not contain tabs or newlines")
    if bad := denied(fields, deny):
        raise ValueError(f"refusing to log a deny-listed accession: {', '.join(bad)}")
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger, "a+") as fh:  # a+ creates; the lock serialises concurrent array tasks
        fcntl.flock(fh, fcntl.LOCK_EX)
        fh.seek(0)
        text = fh.read()
        if text and not text.endswith("\n"):
            raise ValueError(f"{ledger}: missing final newline")
        lines = text.splitlines()
        if lines and lines[0] != HEADER:
            raise ValueError(f"{ledger}: unexpected header")
        rows = [ln.split("\t") for ln in lines[1:]]
        if errors := verify(rows):
            raise ValueError(f"{ledger} already broken: {errors[0]}")
        prev = rows[-1][-1] if rows else GENESIS
        digest = row_hash(prev, fields)
        fh.write(("" if lines else HEADER + "\n") + "\t".join([*fields, prev, digest]) + "\n")
    return digest


def check(ledger: Path, deny: list[str], expect_head: str | None = None) -> list[str]:
    rows = read_rows(ledger)
    errors = verify(rows)
    for i, row in enumerate(rows, start=1):
        if bad := denied(row[: len(FIELDS)], deny):
            errors.append(f"row {i}: deny-listed accession {', '.join(bad)}")
    head = rows[-1][-1] if rows else GENESIS
    if expect_head and head != expect_head:
        errors.append(
            f"head {head[:12]} != expected {expect_head[:12]} (rows removed from the tail?)"
        )
    return errors


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ledger", type=Path, default=default_ledger())
    ap.add_argument("--deny-list", type=Path, default=DEFAULT_DENY)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("append")
    for k in ("dataset", "sample", "arm", "code_sha", "reference_id"):
        p.add_argument(f"--{k.replace('_', '-')}", required=True)
    p.add_argument("--label", default="explore")
    c = sub.add_parser("check")
    c.add_argument("--expect-head")
    args = ap.parse_args(argv)
    try:
        deny = read_deny(args.deny_list)
        if args.cmd == "append":
            stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            rec = {k: getattr(args, k) for k in FIELDS if k != "timestamp"} | {"timestamp": stamp}
            print(append(args.ledger, rec, deny))
            return 0
        errors = check(args.ledger, deny, args.expect_head)
    except (OSError, ValueError) as exc:
        print(f"explore_ledger: {exc}", file=sys.stderr)
        return 1
    for e in errors:
        print(f"explore_ledger: {e}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
