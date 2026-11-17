"""One gem-group rule for every module (PLAN BARCODE-01)."""

import pytest

from viralscan.translated import normalize_barcode
from viralscan.utils import strip_gem_group


@pytest.mark.parametrize(
    ("barcode", "expected"),
    [
        ("ACGT-1", "ACGT"),
        ("ACGT-2", "ACGT"),
        ("ACGT-10", "ACGT"),
        ("ACGT", "ACGT"),
        ("ACGT-1-2", "ACGT-1"),
        ("ACGT-1GCTA-1", "ACGT-1GCTA"),
        ("ACGT-", "ACGT-"),
        ("", ""),
    ],
)
def test_strip_gem_group(barcode, expected):
    assert strip_gem_group(barcode) == expected


def test_translated_uses_the_same_rule():
    assert normalize_barcode(" acgt-2 ") == strip_gem_group("ACGT-2") == "ACGT"
