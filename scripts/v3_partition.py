#!/usr/bin/env python3
"""Training/holdout split shared by the VAL-01 generator and the VAL-06 scorer (PLAN VAL-06).

The generator assigns samples to partitions and the scorer reads them back; one implementation
means two correct-looking copies cannot pick different holdout members
(docs/plans/2026-10-05-val01-generator-design.md, step 3 and section 8).

Global largest-remainder apportionment: the holdout seat total is ``round(fraction * N)``, each
stratum's quota is ``size * seats / N`` (the quotas sum to the total), strata take their floors, and
the leftover seats go to the largest fractional remainders, ties broken by stratum key in ASCII
order. Within a stratum the first ``k`` sample IDs in ASCII order are holdout. The split consumes
no randomness: the split seed enters only through :func:`sample_id`.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Iterable, Mapping

TRAINING, HOLDOUT = "training", "holdout"


def sample_id(split_seed: int, dataset: str, stratum_key: str, replicate: int) -> str:
    """Opaque sample ID: ``s`` plus 12 hex characters, carrying no factor or replicate order."""
    text = "|".join([str(split_seed), dataset, stratum_key, str(replicate)])
    return "s" + hashlib.sha256(text.encode()).hexdigest()[:12]


def apportion(sizes: Mapping[str, int], fraction: float) -> dict[str, int]:
    """Holdout seats per stratum by global largest remainder."""
    if not 0 < fraction < 1:
        raise ValueError(f"holdout fraction must be in (0, 1), got {fraction}")
    if any(size < 0 for size in sizes.values()):
        raise ValueError("stratum sizes must be non-negative")
    total = sum(sizes.values())
    if total == 0:
        raise ValueError("no samples to split")
    seats = math.floor(fraction * total + 0.5)  # round half up, not banker's rounding
    quota = {key: size * seats / total for key, size in sizes.items()}
    awarded = {key: math.floor(q) for key, q in quota.items()}
    leftover = seats - sum(awarded.values())
    by_remainder = sorted(sizes, key=lambda key: (-(quota[key] - awarded[key]), key))
    for key in by_remainder[:leftover]:
        awarded[key] += 1
    return awarded


def split(members: Mapping[str, Iterable[str]], fraction: float) -> dict[str, str]:
    """Partition label of every sample, given sample IDs grouped by stratum key."""
    ordered = {key: sorted(ids) for key, ids in members.items()}
    seats = apportion({key: len(ids) for key, ids in ordered.items()}, fraction)
    labels: dict[str, str] = {}
    for key, ids in ordered.items():
        for position, sample in enumerate(ids):
            if sample in labels:
                raise ValueError(f"sample {sample} appears in more than one stratum")
            labels[sample] = HOLDOUT if position < seats[key] else TRAINING
    return labels
