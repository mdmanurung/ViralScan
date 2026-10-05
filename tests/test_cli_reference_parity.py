"""docs/cli_reference.md must list every parser flag (generated regions)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "scripts" / "gen_cli_reference.py"


def test_generated_cli_tables_are_current() -> None:
    result = subprocess.run([sys.executable, str(GEN), "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_every_parser_flag_is_documented() -> None:
    from viralscan.menu import build_parser

    doc = (ROOT / "docs" / "cli_reference.md").read_text(encoding="utf-8")
    missing: list[str] = []

    def walk(parser) -> None:
        for action in parser._actions:
            if hasattr(action, "choices") and isinstance(action.choices, dict):
                for sub in action.choices.values():
                    walk(sub)
            for opt in action.option_strings:
                if opt not in ("--help", "--version") and opt not in doc:
                    missing.append(opt)

    walk(build_parser())
    assert not missing, f"undocumented flags: {sorted(set(missing))}"
