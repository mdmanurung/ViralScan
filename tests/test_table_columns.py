"""TABLE-SCHEMA-01: the published tables' columns are pinned and documented.

``schemas/v3/table_columns.json`` lists the columns of each table in a default run. This test
checks the registry against ``docs/output_reference.md``; the integration test
``tests/integration/test_table_headers.py`` checks it against a real tiny run.
"""

import fnmatch
import json
import re
from importlib import resources
from pathlib import Path

import pytest

DOC = (Path(__file__).resolve().parents[1] / "docs" / "output_reference.md").read_text("utf-8")


def _registry() -> dict:
    return json.loads(
        resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    )


def _documented() -> tuple[set[str], list[str]]:
    """Backticked names in the doc, and the ``*``-patterns among them (``*_called``)."""
    tokens = set(re.findall(r"`([A-Za-z0-9_*]+)`", DOC))
    return {t for t in tokens if "*" not in t}, [t for t in tokens if "*" in t]


def test_registry_is_well_formed() -> None:
    registry = _registry()
    assert registry["schema_version"] == "3.0.0"
    for name, columns in registry["tables"].items():
        assert name.endswith(".tsv")
        assert columns and len(columns) == len(set(columns)), f"{name}: empty or duplicate columns"


def test_combined_summary_is_the_viral_summary_with_a_sample_column() -> None:
    tables = _registry()["tables"]
    assert tables["combined_summary.tsv"] == ["sample", *tables["viral_summary.tsv"]]


@pytest.mark.parametrize("table", sorted(_registry()["tables"]))
def test_every_column_is_documented(table: str) -> None:
    exact, patterns = _documented()
    missing = [
        column
        for column in _registry()["tables"][table]
        if column not in exact and not any(fnmatch.fnmatch(column, p) for p in patterns)
    ]
    assert not missing, f"{table}: columns not in docs/output_reference.md: {missing}"
