"""API-BOUNDARY-01: the documented stable names are exactly the modules' ``__all__``."""

import importlib
import re
from pathlib import Path

import pytest

DOC = (Path(__file__).resolve().parents[1] / "docs" / "api.md").read_text(encoding="utf-8")

STABLE = {
    "viralscan.enrichment": {"cell_type_enrichment", "write_cell_type_enrichment"},
    "viralscan.data_fetch": {"cache_valid", "ensure_viral_data", "fetch_viral_data"},
    "viralscan.errors": {"ViralScanError"},
    "viralscan.run": {"Run", "open_run"},
    "viralscan.compute": {
        "CellCallingError",
        "check_sibling_crossmapping",
        "compute_stats",
        "external_cells",
    },
    "viralscan.utils": {"configure_logging", "setup_script_logging"},
    "viralscan.constants": {"ENSEMBL_SPECIES", "VIRUS_NAME_MAP"},
    "viralscan.scripts.ncbi_fetch": {"fetch_reference"},
    "viralscan.scripts.build_reference": {"build_combined_reference", "fetch_host_cdna"},
}


@pytest.mark.parametrize("module_name, names", sorted(STABLE.items()))
def test_all_matches_the_stable_list(module_name: str, names: set[str]) -> None:
    module = importlib.import_module(module_name)
    assert set(module.__all__) == names
    for name in names:
        assert hasattr(module, name)


def test_every_documented_stable_name_has_a_heading() -> None:
    headings = set(re.findall(r"^### `(\w+)`", DOC, re.M))
    stable = set().union(*STABLE.values())
    assert stable <= headings | {"cache_valid"}, stable - headings


def test_no_private_or_missing_names_are_documented() -> None:
    headings = set(re.findall(r"^### `(\w+)`", DOC, re.M))
    assert not {h for h in headings if h.startswith("_")}
    assert "load_config" not in headings
    assert "## Stability" in DOC


def test_compute_names_are_the_worker_objects() -> None:
    """The facade re-exports; it must not copy."""
    from viralscan import compute
    from viralscan.scripts import cellcalling, detection

    assert compute.compute_stats is detection.compute_stats
    assert compute.check_sibling_crossmapping is detection.check_sibling_crossmapping
    assert compute.external_cells is cellcalling.external_cells
    assert compute.CellCallingError is cellcalling.CellCallingError
