"""GOLDEN-01: the tiny run's deterministic outputs stay byte-identical.

Every refactor in the API-consolidation plan must leave ``golden/tiny_run.json`` unchanged. A unit
that changes an output on purpose says so in its PLAN row and regenerates the file in the same PR:

    VIRALSCAN_UPDATE_GOLDEN=1 PYTHONPATH=src:. python -m pytest -m integration \
        tests/integration/test_golden_outputs.py

The run is the ``completed_run`` fixture of ``test_tiny_end_to_end`` (needs kallisto, bustools, kb and
snakemake on PATH; skipped otherwise).
"""

from __future__ import annotations

import os

import pytest

from ._golden import golden_digests, load_golden, write_golden
from .test_tiny_end_to_end import completed_run  # noqa: F401  (fixture)

pytestmark = pytest.mark.integration


def test_tiny_run_outputs_match_the_golden_digests(completed_run) -> None:  # noqa: F811
    if os.environ.get("VIRALSCAN_UPDATE_GOLDEN") == "1":
        write_golden(completed_run)
        pytest.skip("golden digests rewritten")
    got = golden_digests(completed_run)
    want = load_golden()
    changed = sorted(name for name in want.keys() & got.keys() if want[name] != got[name])
    assert not changed, f"outputs changed: {changed}"
    assert got.keys() == want.keys(), (
        f"files added {sorted(got.keys() - want.keys())}, removed {sorted(want.keys() - got.keys())}"
    )
