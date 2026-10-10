#!/usr/bin/env python3
"""Report remote GitHub Actions lacking an immutable commit pin (REL-12)."""

import argparse
import re
from pathlib import Path


def unpinned_actions(text: str) -> list[str]:
    return [
        value
        for value in re.findall(r"^\s*(?:-\s*)?uses:\s*([^\s#]+)", text, re.M)
        if not value.startswith("./")
        and not re.fullmatch(r"[^@]+@[a-fA-F0-9]{40}", value)
        and not re.fullmatch(r"docker://.+@sha256:[a-fA-F0-9]{64}", value)
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(".github/workflows"))
    args = parser.parse_args()
    paths = sorted([*args.directory.glob("*.yml"), *args.directory.glob("*.yaml")])
    if not paths:
        parser.error("No workflow files found")
    errors = [f"{p}: {v}" for p in paths for v in unpinned_actions(p.read_text())]
    print("\n".join(errors) if errors else "All remote GitHub Actions are commit-pinned.")
    return int(bool(errors))


if __name__ == "__main__":
    raise SystemExit(main())
