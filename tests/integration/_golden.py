"""Digests of the deterministic outputs of one tiny run (PLAN GOLDEN-01).

A refactor that must not change a result is checked by comparing these digests with the pinned
``golden/tiny_run.json``. Only files that two identical runs reproduce byte for byte are listed:
the result tables and gene lists, ``combined_summary.tsv`` and ``summary.txt`` without its
``Runtime:`` line. JSON provenance files carry paths and timestamps and are left out.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

GOLDEN_FILE = Path(__file__).parent / "golden" / "tiny_run.json"


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def golden_digests(run: Path) -> dict[str, str]:
    """Map ``<sample>/<relative path>`` (and ``combined_summary.tsv``) to a sha256."""
    samples = sorted(p.parent for p in run.glob("*/config.yaml"))
    assert len(samples) == 1, f"expected one sample directory, found {samples}"
    sample = samples[0]
    digests: dict[str, str] = {}
    combined = run / "combined_summary.tsv"
    if combined.is_file():
        digests["combined_summary.tsv"] = _digest(combined.read_bytes())
    for path in sorted((sample / "results").iterdir()):
        if path.is_file() and path.suffix in {".tsv", ".txt"}:
            digests[f"results/{path.name}"] = _digest(path.read_bytes())
    summary = (sample / "summary.txt").read_text(encoding="utf-8").splitlines(keepends=True)
    kept = "".join(line for line in summary if not line.startswith("Runtime:"))
    digests["summary.txt"] = _digest(kept.encode("utf-8"))
    return digests


def load_golden() -> dict[str, str]:
    return json.loads(GOLDEN_FILE.read_text(encoding="utf-8"))


def write_golden(run: Path) -> None:
    GOLDEN_FILE.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN_FILE.write_text(
        json.dumps(golden_digests(run), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
