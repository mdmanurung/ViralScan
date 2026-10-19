"""CPU budget policy: ``--cores`` is the total concurrent CPU budget (total_budget_v1)."""

KB_WORKERS = 8  # legacy ``kb count`` / ``kb ref`` request
HOMOLOGY_WORKERS = 3  # legacy minimap2 host-homology request
CDHIT_MEMORY_MB = 8000
RESOURCE_POLICY = "total_budget_v1"


def bounded_workers(legacy_request: int, cores: int | None) -> int:
    """Native worker count: the legacy request, never more than the budget.

    ``cores=None`` means no budget was given (build commands without ``--cores``),
    which keeps the legacy request unchanged.
    """
    if cores is None:
        return legacy_request
    return max(1, min(legacy_request, int(cores)))


def cdhit_threads(cores: int | None) -> int:
    """cd-hit-est ``-T``: 0 (all visible CPUs, legacy) without a budget, else the budget."""
    return 0 if cores is None else max(1, int(cores))
