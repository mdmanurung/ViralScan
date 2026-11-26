"""TABLE-SCHEMA-01: a default tiny run writes exactly the columns of ``table_columns.json``."""

from __future__ import annotations

import json
from importlib import resources

import pytest

from .test_tiny_end_to_end import _sample_dir, completed_run  # noqa: F401  (fixture)

pytestmark = pytest.mark.integration


def test_tiny_run_headers_match_the_registry(completed_run) -> None:  # noqa: F811
    registry = json.loads(
        resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    )["tables"]
    sample = _sample_dir(completed_run)
    mismatched = {}
    for name, columns in registry.items():
        path = completed_run / name if name == "combined_summary.tsv" else sample / "results" / name
        header = path.read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
        if header != columns:
            mismatched[name] = header
    assert not mismatched, f"headers differ from table_columns.json: {sorted(mismatched)}"
    optional = json.loads(
        resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    )["optional_tables"]  # TABLE-SCHEMA-02: a default run writes none, a switched-on arm may
    written = {p.name for p in (sample / "results").glob("*.tsv")}
    assert written <= set(registry) | set(optional), (
        f"unregistered tables: {sorted(written - set(registry) - set(optional))}"
    )
    for name in written & set(optional):
        header = (sample / "results" / name).read_text(encoding="utf-8").split("\n", 1)[0]
        assert header.split("\t") == optional[name], name
