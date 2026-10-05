"""MECH-F: off-list drops at barcode correction must reach the count audit."""

import json
from pathlib import Path

import pytest

from viralscan.scripts import multimap
from viralscan.validation import require_schema_valid

OLD_AUDIT = {
    "input_molecules": 1,
    "resolved_molecules": 1,
    "unique_molecules": 1,
    "ambiguous_molecules": 0,
    "unresolved_molecules": 0,
    "ignored_read_multiplicity": 0,
    "allocated_ambiguous_mass": 0.0,
}


@pytest.fixture(autouse=True)
def _no_kb_lookup(monkeypatch):
    """Hermetic: with kb on PATH the real resolver runs `kb info` through the fake."""
    monkeypatch.setattr(multimap, "tool_path", lambda name: name)


def _fake_inspect(totals: dict[str, tuple[int, int]]):
    def run(command, check):
        assert command[1] == "inspect" and check is True
        records, reads = totals[Path(command[-1]).name]
        Path(command[command.index("-o") + 1]).write_text(
            json.dumps({"numRecords": records, "numReads": reads})
        )

    return run


def test_bus_totals_parses_inspect_json(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(multimap.subprocess, "run", _fake_inspect({"output.bus": (10, 12)}))
    assert multimap.bus_totals(tmp_path / "output.bus") == (10, 12)


def test_correction_audit_records_raw_corrected_and_dropped(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(
        multimap.subprocess,
        "run",
        _fake_inspect({"raw.bus": (100, 200), "corr.bus": (85, 170)}),
    )
    audit = multimap.correction_audit(tmp_path / "raw.bus", tmp_path / "corr.bus", "kb")
    assert audit == {
        "barcode_correction": "kb",
        "bus_records_raw": 100,
        "bus_records_after_correction": 85,
        "offlist_dropped_records": 15,
        "bus_reads_raw": 200,
        "bus_reads_after_correction": 170,
        "offlist_dropped_reads": 30,
        "bus_totals_source": "bustools_inspect",
    }
    require_schema_valid({**OLD_AUDIT, **audit}, "count_audit.schema.json")


def test_no_correction_is_zero_drop_without_running_bustools(tmp_path: Path) -> None:
    audit = multimap.correction_audit(tmp_path / "raw.bus", tmp_path / "raw.bus", "none")
    assert audit["offlist_dropped_records"] == 0
    assert audit["bus_totals_source"] == "not_corrected"
    require_schema_valid({**OLD_AUDIT, **audit}, "count_audit.schema.json")


def test_old_audits_still_validate() -> None:
    require_schema_valid(OLD_AUDIT, "count_audit.schema.json")


def test_schema_copies_identical() -> None:
    root = Path(__file__).resolve().parents[1]
    name = "v3/count_audit.schema.json"
    assert (root / "schemas" / name).read_bytes() == (
        root / "src" / "viralscan" / "schemas" / name
    ).read_bytes()
