"""PANEL-01 WP2: scripts/panel_promote.py edits only the overlay columns and refuses unknown accessions."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "panel_promote", Path(__file__).resolve().parents[1] / "scripts" / "panel_promote.py"
)
pp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pp)


def _rows():
    return [
        {"accession": "A1", "species": "x", "panel": "max", "sibling_group": ""},
        {"accession": "B2", "species": "y", "panel": "broad", "sibling_group": "G"},
    ]


def test_apply_changes_only_overlays_and_keeps_blank_sibling_group() -> None:
    rows = _rows()
    changes = pp.apply(
        rows,
        [
            {"accession": "A1", "panel": "shipped", "sibling_group": ""},
            {"accession": "B2", "panel": "shipped", "sibling_group": ""},
        ],
    )
    assert [r["panel"] for r in rows] == ["shipped", "shipped"]
    assert rows[1]["sibling_group"] == "G" and rows[0]["species"] == "x"
    assert changes == [("A1", "panel", "max", "shipped"), ("B2", "panel", "broad", "shipped")]
    assert pp.apply(rows, [{"accession": "A1", "panel": "shipped"}]) == []  # idempotent


def test_unknown_accession_is_an_error() -> None:
    with pytest.raises(SystemExit):
        pp.apply(_rows(), [{"accession": "Z9", "panel": "shipped"}])
