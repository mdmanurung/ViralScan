"""Digests of the deterministic outputs of one tiny run (PLAN GOLDEN-01).

A refactor that must not change a result is checked by comparing these digests with the pinned
``golden/tiny_run.json``. Only files that two identical runs reproduce byte for byte are listed:
the result tables and gene lists, ``combined_summary.tsv`` and ``summary.txt`` without its
``Runtime:`` line. JSON provenance files carry paths and timestamps and are left out.
"""

from __future__ import annotations

import json
from pathlib import Path

from viralscan.demo import output_digests

GOLDEN_FILE = Path(__file__).parent / "golden" / "tiny_run.json"


def golden_digests(run: Path) -> dict[str, str]:
    """Map ``<sample>/<relative path>`` (and ``combined_summary.tsv``) to a sha256.

    The function ``viralscan demo`` uses (PLAN DEMO-01), so the demo and this test cannot drift.
    """
    return output_digests(run)


def load_golden() -> dict[str, str]:
    return json.loads(GOLDEN_FILE.read_text(encoding="utf-8"))


def write_golden(run: Path) -> None:
    GOLDEN_FILE.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN_FILE.write_text(
        json.dumps(golden_digests(run), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
