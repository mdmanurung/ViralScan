"""QA-02: a missing ``bustools`` must fail clearly, not put ``None`` in the argv."""

from __future__ import annotations

from pathlib import Path

import pytest

from viralscan.scripts import multimap


def test_missing_bustools_fails_clearly(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(multimap, "tool_path", lambda name: None)
    with pytest.raises(RuntimeError, match="bustools was not found"):
        multimap.bus_totals(tmp_path / "output.bus")
