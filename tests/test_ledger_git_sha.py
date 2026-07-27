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
    for record in tampered["deviations"]:
        if record["deviation_id"] == "DEV-012":
            record["git_sha"] = "5cad263"  # the landing commit, the R7-F1 mistake

    errors = check_git_sha_fields(tampered)

    assert any("landing commit rather than the base commit" in error for error in errors)


def test_an_unresolvable_git_sha_is_detected(ledger: dict) -> None:
    tampered = deepcopy(ledger)
    tampered["deviations"][-1]["git_sha"] = "deadbeef"

    errors = check_git_sha_fields(tampered)

    assert any("does not resolve" in error for error in errors)
