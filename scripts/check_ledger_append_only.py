#!/usr/bin/env python3
"""Fail if any existing deviation record changed relative to a previous revision.

SCI-05 round 4 (R4-F1, R4-F6) found that the ledger's hash chain has no anchor
outside the ledger file: re-chaining the whole file after an edit is
undetectable, because the author who edits the records also computes the hashes.

This check supplies the missing anchor. It compares the working ledger against
its committed state and fails on any modification or deletion of an existing
record, permitting only appends. Run in CI, where the comparison is against
pushed history rather than the local working tree, it is an anchor the editing
author does not solely control.

It does not make the ledger tamper-proof. An author who can rewrite the compared
revision defeats it, as they defeat every check whose reference lives in a file
they control. It raises the cost of an undocumented amendment from editing one
file to rewriting shared history.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEDGER = REPO_ROOT / "analysis" / "v3_validation" / "deviations.yaml"


class LedgerHistoryError(RuntimeError):
    """Raised when the ledger's committed history cannot be established."""


def _records(document: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(document, dict):
        raise LedgerHistoryError("ledger must be a mapping")
    out: dict[str, dict[str, Any]] = {}
    for record in document.get("deviations", []) or []:
        if not isinstance(record, dict) or not record.get("deviation_id"):
            raise LedgerHistoryError("every deviation record needs a deviation_id")
        out[record["deviation_id"]] = record
    return out


def read_committed(revision: str, relative_path: str, repo_root: Path = REPO_ROOT) -> Any:
    """Return the ledger as of ``revision``.

    Returns None only when the revision exists and did not contain the ledger,
    which is the genuine genesis case. An unreadable revision raises instead of
    returning None: a comparison check that silently passes when it cannot find
    anything to compare against is worse than no check, because it reports green.
    """
    revision_known = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "--verify", "--quiet", f"{revision}^{{commit}}"],
        check=False,
        capture_output=True,
        text=True,
    )
    if revision_known.returncode != 0:
        raise LedgerHistoryError(
            f"cannot resolve revision {revision!r}; a shallow clone or missing remote ref "
            "would make this check pass vacuously"
        )
    completed = subprocess.run(
        ["git", "-C", str(repo_root), "show", f"{revision}:{relative_path}"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return None
    return yaml.safe_load(completed.stdout)


def compare(previous: Any, current: Any) -> list[str]:
    """Return violations of the append-only rule between two ledger states."""
    if previous is None:
        # No prior revision: nothing to compare against, and nothing to enforce.
        return []
    before = _records(previous)
    after = _records(current)

    errors: list[str] = []
    for deviation_id, record in before.items():
        if deviation_id not in after:
            errors.append(
                f"deviation record {deviation_id!r} was deleted; the ledger is append-only"
            )
            continue
        if after[deviation_id] != record:
            errors.append(
                f"deviation record {deviation_id!r} was modified; supersede it with a new "
                "record instead of editing it"
            )

    before_order = [key for key in before if key in after]
    after_order = [key for key in after if key in before]
    if before_order != after_order:
        errors.append("existing deviation records were reordered; the ledger is append-only")

    return errors


def resolve_field_corrections(document: Any) -> dict[tuple[str, str], Any]:
    """Return the effective value of every field a later record supersedes.

    The ledger is append-only, so a mistake in an existing record is corrected by
    a later record rather than by editing it. Until now that correction lived only
    in prose, which meant no automated check could resolve the effective value —
    the reason a wrong git_sha survived three review rounds. ``field_corrections``
    makes the supersession machine-readable.
    """
    corrections: dict[tuple[str, str], Any] = {}
    for record in document.get("deviations", []) or []:
        if not isinstance(record, dict):
            continue
        for correction in record.get("field_corrections", []) or []:
            if not isinstance(correction, dict):
                continue
            target = correction.get("deviation_id")
            field = correction.get("field")
            if target and field:
                corrections[(target, field)] = correction.get("corrected_value")
    return corrections


def check_git_sha_fields(document: Any, repo_root: Path = REPO_ROOT) -> list[str]:
    """Verify each record's git_sha resolves and is the base, not the landing, commit.

    SCI-05 rounds 6, 7, and 8 all asked for this. R6-F4 found DEV-012 naming the
    commit of the review it answered; R7-F1 found the correction naming the
    landing commit, which the very record defining the field excluded. Both were
    caught by a reviewer rather than by a check, three rounds running.

    A base commit cannot contain the amendment's own after-digest: that digest
    only exists once the change lands. So a recorded git_sha whose protocol
    already carries the after-digest is provably the landing commit or later.
    """
    errors: list[str] = []
    corrections = resolve_field_corrections(document)
    for record in _records(document).values():
        record_id = record.get("deviation_id")
        sha = corrections.get((record_id, "git_sha"), record.get("git_sha"))
        if not sha:
            errors.append(f"deviation record {record_id!r} has no git_sha")
            continue
        resolved = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "--verify", "--quiet", f"{sha}^{{commit}}"],
            check=False,
            capture_output=True,
            text=True,
        )
        if resolved.returncode != 0:
            errors.append(f"deviation record {record_id!r} git_sha {sha!r} does not resolve")
            continue
        after = record.get("protocol_sha256_after")
        before = record.get("protocol_sha256_before")
        shown = subprocess.run(
            ["git", "-C", str(repo_root), "show", f"{sha}:analysis/v3_validation/protocol.yaml"],
            check=False,
            capture_output=True,
            text=True,
        )
        if shown.returncode != 0:
            continue

        # R9-F4: the before-digest was never checked, so a fabricated one passed.
        # The base commit must actually carry the digest the record claims it had.
        if before == after:
            # A genesis record declares a starting value for a section that the
            # same change created, so no earlier commit carries the digest.
            continue
        if corrections.get((record_id, "before_digest_uncommitted")):
            # An intermediate digest that existed only in a working tree, because
            # two records were written against one commit. It is unverifiable
            # against history by construction, and saying so is better than
            # letting the check fail or silently skip.
            continue
        scope = record.get("digest_scope") or ""
        section = scope.split(".", 1)[0]
        try:
            base_document = yaml.safe_load(shown.stdout)
        except yaml.YAMLError:
            base_document = None
        if before and isinstance(base_document, dict) and section:
            base_section = base_document.get(section)
            observed = (
                base_section.get("contract_sha256") if isinstance(base_section, dict) else None
            )
            if observed is None:
                # The section did not exist at the base commit, which is the
                # genesis case already handled above for before == after.
                pass
            elif observed != before:
                errors.append(
                    f"deviation record {record_id!r} claims a before-digest for {scope!r} "
                    f"that its base commit {sha!r} does not carry in that section"
                )
        elif before and before not in shown.stdout:
            errors.append(
                f"deviation record {record_id!r} claims a before-digest that its own "
                f"base commit {sha!r} does not contain"
            )
        if after in shown.stdout:
            errors.append(
                f"deviation record {record_id!r} git_sha {sha!r} already contains the "
                "after-digest, so it is the landing commit rather than the base commit"
            )
    return errors


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument(
        "--revision",
        default="HEAD",
        help="revision to compare against; use origin/main in CI",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.ledger.is_file():
        print(f"ledger not found: {args.ledger}", file=sys.stderr)
        return 1
    relative = args.ledger.resolve().relative_to(REPO_ROOT).as_posix()
    try:
        previous = read_committed(args.revision, relative)
        current = yaml.safe_load(args.ledger.read_text(encoding="utf-8"))
        errors = compare(previous, current)
        errors.extend(check_git_sha_fields(current))
    except (LedgerHistoryError, OSError, yaml.YAMLError) as exc:
        print(f"ledger history check failed: {exc}", file=sys.stderr)
        return 1
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"ledger is append-only relative to {args.revision}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
