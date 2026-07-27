from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from scripts.check_ledger_append_only import check_git_sha_fields

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "analysis" / "v3_validation" / "deviations.yaml"


@pytest.fixture
def ledger() -> dict:
    return yaml.safe_load(LEDGER.read_text(encoding="utf-8"))


def test_every_recorded_git_sha_resolves_and_is_a_base_commit(ledger: dict) -> None:
    assert check_git_sha_fields(ledger) == []


def test_a_landing_commit_recorded_as_the_base_is_detected(ledger: dict) -> None:
    """R8-F5: this exact error took rounds 6, 7, and 8 to catch by eye."""
    tampered = deepcopy(ledger)
    # DEV-012's git_sha is superseded via field_corrections, so tamper the
    # correction itself: that is the value the check actually resolves.
    for record in tampered["deviations"]:
        for correction in record.get("field_corrections", []) or []:
            if correction.get("deviation_id") == "DEV-012" and correction["field"] == "git_sha":
                correction["corrected_value"] = "5cad263"  # the R7-F1 mistake

    errors = check_git_sha_fields(tampered)

    assert any("landing commit rather than the base commit" in error for error in errors)


def test_an_unresolvable_git_sha_is_detected(ledger: dict) -> None:
    tampered = deepcopy(ledger)
    tampered["deviations"][-1]["git_sha"] = "deadbeef"

    errors = check_git_sha_fields(tampered)

    assert any("does not resolve" in error for error in errors)


def test_a_fabricated_before_digest_is_detected(ledger: dict) -> None:
    """R9-F4: the before-digest was never verified, so any value passed."""
    tampered = deepcopy(ledger)
    for record in tampered["deviations"]:
        if record["deviation_id"] == "DEV-005":
            record["protocol_sha256_before"] = "a" * 64

    errors = check_git_sha_fields(tampered)

    # Round 10 made the check scope-aware, so the message names the section.
    assert any(
        "does not carry in that section" in error or "does not contain" in error for error in errors
    )


def test_field_corrections_are_machine_readable(ledger: dict) -> None:
    """A supersession recorded only in prose cannot be resolved by any check."""
    from scripts.check_ledger_append_only import resolve_field_corrections

    corrections = resolve_field_corrections(ledger)

    assert corrections[("DEV-012", "git_sha")] == "62d2bcc"
    assert corrections[("DEV-002", "before_digest_uncommitted")] is True


def test_a_digest_from_the_wrong_section_is_detected(ledger: dict) -> None:
    """R10-F2: substring matching accepted a real digest from another section."""
    import subprocess

    base = subprocess.run(
        ["git", "show", "62d2bcc:analysis/v3_validation/protocol.yaml"],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
    ).stdout
    other_section_digest = yaml.safe_load(base)["calibration"]["contract_sha256"]

    tampered = deepcopy(ledger)
    for record in tampered["deviations"]:
        if record["deviation_id"] == "DEV-012":
            record["protocol_sha256_before"] = other_section_digest

    errors = check_git_sha_fields(tampered)

    assert any("does not carry in that section" in error for error in errors)


def test_a_before_digest_for_an_absent_section_is_not_accepted(ledger: dict) -> None:
    """R11-F2: 'the section is missing' used to pass as if it were genesis.

    Genesis is before == after and returns earlier. Reaching the section lookup
    with a distinct before-digest and finding no section means the claim cannot
    be verified, which is the case this check exists to catch.
    """
    tampered = deepcopy(ledger)
    for record in tampered["deviations"]:
        if record["deviation_id"] == "DEV-012":
            record["digest_scope"] = "a_section_that_never_existed.subkey"

    errors = check_git_sha_fields(tampered)

    assert any("carries no value at that path" in error for error in errors)


def test_an_unreadable_base_protocol_is_not_accepted(ledger: dict, monkeypatch) -> None:
    """R11-F2: a failed `git show` used to skip the record silently."""
    import subprocess as _subprocess

    from scripts import check_ledger_append_only

    real_run = _subprocess.run

    def _fail_protocol_show(cmd, **kwargs):
        if len(cmd) > 3 and cmd[3] == "show" and "protocol.yaml" in cmd[-1]:
            return _subprocess.CompletedProcess(cmd, 128, stdout="", stderr="fatal")
        return real_run(cmd, **kwargs)

    monkeypatch.setattr(check_ledger_append_only.subprocess, "run", _fail_protocol_show)

    errors = check_git_sha_fields(ledger)

    assert errors
    assert any("cannot be verified against history" in error for error in errors)


def test_a_missing_digest_scope_is_not_accepted(ledger: dict) -> None:
    """R11-F2 (1): an absent scope fell through to a whole-file substring search."""
    tampered = deepcopy(ledger)
    for record in tampered["deviations"]:
        if record["deviation_id"] == "DEV-012":
            record.pop("digest_scope", None)

    errors = check_git_sha_fields(tampered)

    assert any("no digest_scope" in error for error in errors)


def test_a_nested_digest_scope_resolves_every_segment() -> None:
    """R11-F2 (3): only the first segment was resolved, so nested paths read None."""
    from scripts.check_ledger_append_only import resolve_digest_scope

    document = {"a": {"b": {"contract_sha256": "beef"}}}

    assert resolve_digest_scope(document, "a.b.contract_sha256") == "beef"
    assert resolve_digest_scope(document, "a.contract_sha256") is None
    assert resolve_digest_scope(document, "missing.contract_sha256") is None
