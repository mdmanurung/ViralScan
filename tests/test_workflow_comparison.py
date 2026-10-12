"""DOC-08: the workflow pros/cons page is generated from the protocol and is balanced."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "gen_workflow_comparison", ROOT / "scripts" / "gen_workflow_comparison.py"
)
gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gen)


def test_page_is_current():
    assert gen.main(["--check"]) == 0


def test_five_families_cover_every_non_legacy_workflow_once():
    ids = [w["id"] for w in gen._load()["workflows"] if not w["id"].startswith("W10")]
    covered = [i for _, members in gen.FAMILIES for i in members]
    assert len(gen.FAMILIES) == 5
    assert sorted(covered) == sorted(ids)


def test_every_family_is_balanced_and_prints_no_performance_number():
    rows = [r for r in gen.render().splitlines() if r.startswith("| ") and "`W" in r]
    assert len(rows) == 5
    for row in rows:
        cells = [c.strip() for c in row.strip().strip("|").split(" | ")]
        assert cells[2] and cells[3], row  # strengths and limitations
    assert "%" not in gen.render()
