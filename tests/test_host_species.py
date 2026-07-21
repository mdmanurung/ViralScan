"""G2: non-human host support in build-ref (mouse + others via Ensembl)."""

from __future__ import annotations

import pytest

from viralscan.constants import ENSEMBL_SPECIES
from viralscan.scripts.build_reference import _ensembl_species_key


def test_mouse_and_common_hosts_are_supported():
    assert ENSEMBL_SPECIES["mouse"] == ("mus_musculus", "GRCm39")
    for host in ("human", "mouse", "rat", "macaque"):
        assert _ensembl_species_key(host) == host


def test_unknown_host_raises_with_supported_list():
    with pytest.raises(ValueError, match="Supported values"):
        _ensembl_species_key("unicorn")
