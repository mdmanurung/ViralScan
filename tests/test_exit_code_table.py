"""CLI-EXITS-01: the documented exit codes are the ones the code uses."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = (ROOT / "docs" / "cli_reference.md").read_text(encoding="utf-8")
SRC = ROOT / "src" / "viralscan"


def _documented_codes() -> set[int]:
    section = DOC.split("## Exit codes", 1)[1].split("\n## ", 1)[0]
    return {int(code) for code in re.findall(r"^\| `(\d+)` \|", section, re.M)}


def test_exit_code_table_lists_exactly_the_supported_codes() -> None:
    assert _documented_codes() == {0, 1, 2, 130, 143}


def test_the_signal_codes_are_set_where_the_table_says() -> None:
    assert "sys.exit(130)" in (SRC / "menu.py").read_text(encoding="utf-8")
    assert "sys.exit(143)" in (SRC / "run_safety.py").read_text(encoding="utf-8")
