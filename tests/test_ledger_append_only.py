from __future__ import annotations

import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from scripts.check_ledger_append_only import (
    LedgerHistoryError,
    compare,
    read_committed,
)

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "analysis" / "v3_validation" / "deviations.yaml"


@pytest.fixture
def ledger() -> dict:
    return yaml.safe_load(LEDGER.read_text(encoding="utf-8"))


def test_an_unchanged_ledger_passes(ledger: dict) -> None:
    assert compare(deepcopy(ledger), ledger) == []


def test_appending_a_record_is_allowed(ledger: dict) -> None:
    appended = deepcopy(ledger)
    appended["deviations"].append({"deviation_id": "DEV-999", "reason": "new"})

    assert compare(ledger, appended) == []


def test_editing_an_existing_record_is_rejected(ledger: dict) -> None:
    """R4-F1: re-chaining after an edit defeats the in-file chain entirely."""
    edited = deepcopy(ledger)
    edited["deviations"][0]["reason"] = "rewritten after the fact"

    errors = compare(ledger, edited)

    assert any("was modified" in error for error in errors)


def test_deleting_a_record_is_rejected(ledger: dict) -> None:
    """R4-F3: genesis records could be deleted and the chain rebuilt."""
    trimmed = deepcopy(ledger)
    removed = trimmed["deviations"].pop(0)

    errors = compare(ledger, trimmed)

    assert any(removed["deviation_id"] in error and "deleted" in error for error in errors)


def test_reordering_records_is_rejected(ledger: dict) -> None:
    reordered = deepcopy(ledger)
    reordered["deviations"][0], reordered["deviations"][1] = (
        reordered["deviations"][1],
        reordered["deviations"][0],
    )

    errors = compare(ledger, reordered)

    assert any("reordered" in error for error in errors)


def test_a_genesis_state_with_no_prior_revision_is_permitted(ledger: dict) -> None:
    assert compare(None, ledger) == []


def test_an_unresolvable_revision_fails_closed_rather_than_passing(tmp_path: Path) -> None:
    """A comparison that silently passes when it cannot compare reports green falsely."""
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)

    with pytest.raises(LedgerHistoryError, match="cannot resolve revision"):
        read_committed("origin/nonexistent", "some/ledger.yaml", repo_root=tmp_path)


def test_the_committed_ledger_is_append_only_against_head() -> None:
    """The real check, against real history."""
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "check_ledger_append_only.py"),
            "--revision",
            "HEAD",
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )

    assert completed.returncode == 0, completed.stderr
