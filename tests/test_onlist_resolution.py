"""ONLIST-01: one on-list rule; a missing ngs_tools never looks like "no list".

Before, ``anello_align`` and ``discovery`` caught ``ValueError``. ``ChemistryError`` subclasses it, so
without kb-python's ``ngs_tools`` those arms ran without barcode correction instead of stopping.
"""

import pytest

from viralscan import chemistry, discovery
from viralscan.chemistry import ChemistryError, onlist_path_or_none
from viralscan.scripts import anello_align


def _no_ngs_tools(name: str) -> str:
    raise ChemistryError("kb-python's on-lists (ngs_tools) are not importable")


@pytest.mark.parametrize("technology", ["dropseq", "0,0,12:0,12,20:1,0,0"])
def test_a_chemistry_without_a_list_gives_none(technology: str) -> None:
    assert onlist_path_or_none(technology) is None


def test_a_missing_ngs_tools_propagates(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(chemistry, "onlist_path", _no_ngs_tools)
    with pytest.raises(ChemistryError, match="ngs_tools"):
        onlist_path_or_none("10xv3")


def test_the_discovery_arm_stops_instead_of_skipping_correction(monkeypatch) -> None:
    monkeypatch.setattr(chemistry, "onlist_path", _no_ngs_tools)
    with pytest.raises(ChemistryError):
        discovery.resolve_onlist("10xv3", None)


def test_the_anello_arm_uses_the_same_rule() -> None:
    assert anello_align.onlist_path_or_none is onlist_path_or_none
    assert not hasattr(anello_align, "_onlist")


def test_the_discovery_arm_still_gives_none_for_a_list_free_chemistry() -> None:
    assert discovery.resolve_onlist("dropseq", None) is None
