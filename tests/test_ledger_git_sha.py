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

    assert any("does not contain" in error for error in errors)


def test_field_corrections_are_machine_readable(ledger: dict) -> None:
    """A supersession recorded only in prose cannot be resolved by any check."""
    from scripts.check_ledger_append_only import resolve_field_corrections

    corrections = resolve_field_corrections(ledger)

    assert corrections[("DEV-012", "git_sha")] == "62d2bcc"
    assert corrections[("DEV-002", "before_digest_uncommitted")] is True
