"""Tests for shared constants."""

from viralscan.constants import VIRUS_NAME_MAP


def test_virus_names_do_not_contain_replacement_character() -> None:
    assert all("\ufffd" not in name for name in VIRUS_NAME_MAP.values())


def test_onyong_name_is_readable() -> None:
    assert VIRUS_NAME_MAP["ONYONG"] == "O'nyong-nyong virus"
