#!/usr/bin/env python3
"""Compare package, container, recipe and citation versions without a release."""

import argparse
import re
from pathlib import Path

import yaml


def version_errors(root: Path) -> list[str]:
    patterns = {
        "src/viralscan/__init__.py": r'__version__\s*=\s*"([^"]+)"',
        "Dockerfile": r'version="([^"]+)"',
        "Singularity.def": r"^\s*Version\s+(\S+)",
        "conda-recipe/meta.yaml": r'\{% set version = "([^"]+)" %\}',
    }
    versions = {}
    for rel, pattern in patterns.items():
        match = re.search(pattern, (root / rel).read_text(), re.M)
        if not match:
            return [f"Version not found in {rel}"]
        versions[rel] = match.group(1)
    citation = yaml.safe_load((root / "CITATION.cff").read_text())
    versions["CITATION.cff"] = str(citation.get("version", ""))
    versions["CITATION.cff preferred-citation"] = str(
        citation.get("preferred-citation", {}).get("version", "")
    )
    expected = versions["src/viralscan/__init__.py"]
    return [
        f"{rel}: {value!r}, expected {expected!r}"
        for rel, value in versions.items()
        if value != expected
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    errors = version_errors(parser.parse_args().root)
    print(
        "\n".join(errors)
        if errors
        else "Software versions agree; release metadata and sdist checksum remain maintainer-owned."
    )
    return int(bool(errors))


if __name__ == "__main__":
    raise SystemExit(main())
