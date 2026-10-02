#!/usr/bin/env python3
"""Render flag/default/help tables from the argparse parser into docs/cli_reference.md.

Only text between ``<!-- BEGIN GENERATED: <name> -->`` and ``<!-- END GENERATED -->``
is rewritten; hand-written prose stays. ``--check`` exits 1 if the doc is stale.

    PYTHONPATH=src python scripts/gen_cli_reference.py [--check]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DOC = Path(__file__).resolve().parents[1] / "docs" / "cli_reference.md"
REGION = re.compile(
    r"(<!-- BEGIN GENERATED: (?P<name>[\w -]+) -->\n)(?P<body>.*?)(<!-- END GENERATED -->)",
    re.DOTALL,
)


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    return {}


def _cell(text: str) -> str:
    return re.sub(r"\s+", " ", text).replace("|", "\\|").strip()


def _table(parser: argparse.ArgumentParser) -> str:
    fmt = parser._get_formatter()
    rows = []
    for a in parser._actions:
        if isinstance(
            a, (argparse._SubParsersAction, argparse._HelpAction, argparse._VersionAction)
        ):
            continue
        flag = ", ".join(a.option_strings) if a.option_strings else a.dest
        if a.option_strings and a.nargs != 0:
            flag += f" {a.metavar or a.dest.upper()}"
        if a.required:
            default = "*required*"
        elif a.default is None:
            default = "*(none)*"
        elif a.default is argparse.SUPPRESS:
            default = ""
        else:
            default = f"`{a.default}`"
        helptext = fmt._expand_help(a) if a.help and a.help != argparse.SUPPRESS else ""
        rows.append(f"| `{_cell(flag)}` | {_cell(default)} | {_cell(helptext)} |")
    return "| Flag | Default | Help |\n|------|---------|------|\n" + "\n".join(rows) + "\n"


def tables() -> dict[str, str]:
    from viralscan.menu import build_parser

    root = build_parser()
    out = {"main": _table(root)}

    def walk(prefix: str, parser: argparse.ArgumentParser) -> None:
        for name, sub in _subparsers(parser).items():
            key = f"{prefix} {name}".strip()
            out[key] = _table(sub)
            walk(key, sub)

    walk("", root)
    return out


def render(text: str, generated: dict[str, str]) -> str:
    def sub(m: re.Match) -> str:
        name = m.group("name")
        if name not in generated:
            raise SystemExit(f"unknown generated region: {name}")
        return m.group(1) + "\n" + generated[name] + "\n" + m.group(4)

    return REGION.sub(sub, text)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="exit 1 if the doc is out of date")
    args = ap.parse_args()
    current = DOC.read_text(encoding="utf-8")
    new = render(current, tables())
    if args.check:
        if new != current:
            print(f"{DOC} is stale; run: python scripts/gen_cli_reference.py", file=sys.stderr)
            return 1
        return 0
    DOC.write_text(new, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
