"""EVID-ONLIST-01: the replay finds the on-list kb count left in kb-python/."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from viralscan.kb_outputs import KbCountOutputs
from viralscan.scripts.evidence_run import _replay_whitelist


def _kb(tmp_path, *names: str) -> KbCountOutputs:
    kb = KbCountOutputs(tmp_path)
    kb.root.mkdir(parents=True)
    for name in names:
        (kb.root / name).write_text("AAAA\n")
    return kb


def test_named_10x_onlist_is_used(tmp_path) -> None:
    kb = _kb(tmp_path, "10x_version3_whitelist.txt")
    got = _replay_whitelist(SimpleNamespace(whitelist=None), kb)
    assert got == str(kb.root / "10x_version3_whitelist.txt")


def test_bustools_whitelist_wins_over_a_named_list(tmp_path) -> None:
    kb = _kb(tmp_path, "whitelist.txt", "10x_version2_whitelist.txt")
    assert _replay_whitelist(SimpleNamespace(whitelist=None), kb) == str(kb.root / "whitelist.txt")


def test_explicit_whitelist_wins(tmp_path) -> None:
    kb = _kb(tmp_path, "10x_version3_whitelist.txt")
    assert _replay_whitelist(SimpleNamespace(whitelist="/x/w.txt"), kb) == "/x/w.txt"


def test_no_onlist_means_no_correction(tmp_path) -> None:
    assert _replay_whitelist(SimpleNamespace(whitelist=None), _kb(tmp_path)) is None


def test_two_named_lists_are_refused(tmp_path) -> None:
    kb = _kb(tmp_path, "10x_version2_whitelist.txt", "10x_version3_whitelist.txt")
    with pytest.raises(ValueError, match="2 on-lists"):
        _replay_whitelist(SimpleNamespace(whitelist=None), kb)
