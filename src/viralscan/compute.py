"""Stable import path for the compute helpers that notebooks and tests call.

Each name is re-exported from its worker module (no copy, no behaviour change), so
``viralscan.compute.compute_stats is viralscan.scripts.detection.compute_stats``. The
worker modules stay workflow internals; import these names instead of ``scripts.*``.
"""

from viralscan.scripts.cellcalling import CellCallingError, external_cells
from viralscan.scripts.detection import check_sibling_crossmapping, compute_stats

__all__ = [
    "CellCallingError",
    "check_sibling_crossmapping",
    "compute_stats",
    "external_cells",
]
