"""Tests for shared constants."""

from viralscan.constants import VIRUS_NAME_MAP


def test_virus_names_do_not_contain_replacement_character() -> None:
    assert all("\ufffd" not in name for name in VIRUS_NAME_MAP.values())


def test_onyong_name_is_readable() -> None:
    assert VIRUS_NAME_MAP["ONYONG"] == "O'nyong-nyong virus"


def test_andet05_aliases_name_the_right_virus() -> None:
    from viralscan.anellovirus import merged_name_map
    from viralscan.constants import VIRUS_GENE_ID_ALIASES
    from viralscan.virus_grouping import virus_name_for_gene

    assert virus_name_for_gene("VARVgp1") == "Variola virus"
    assert VIRUS_GENE_ID_ALIASES["VARV"] == "Variola virus"
    assert VIRUS_NAME_MAP["UUKU"] == "Uukuniemi virus"
    nm = merged_name_map()
    assert virus_name_for_gene("TTVgp1", nm) == virus_name_for_gene("NC_002076.2_TTVgp1", nm)


def test_legacy_name_map_gives_one_name_per_virus() -> None:
    """AUD-20 / D-12: IMMUNO/HIV1 and TLYMPHO/HTLV1 were two display names for one virus."""
    from viralscan.constants import VIRUS_GENE_ID_ALIASES

    for maps in (VIRUS_NAME_MAP, VIRUS_GENE_ID_ALIASES):
        assert maps["IMMUNO"] == VIRUS_NAME_MAP["HIV1"] == "Human immunodeficiency virus 1"
        assert maps["TLYMPHO"] == VIRUS_NAME_MAP["HTLV1"]
