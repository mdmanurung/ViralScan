#!/usr/bin/env python3
"""Training/holdout split shared by the VAL-01 generator and the VAL-06 scorer (PLAN VAL-06).

The generator assigns samples to partitions and the scorer reads them back; one implementation
means two correct-looking copies cannot pick different holdout members
(docs/plans/2026-10-05-val01-generator-design.md, step 3 and section 8).

Global largest-remainder apportionment, as ``partitions.algorithm`` freezes it: each stratum's quota
is ``size * fraction``, strata take their floors, and the remaining seats go to the largest
fractional remainders (ties by stratum key in ASCII order) until the global count reaches the cap,
``ceil(N * fraction)`` (SCI-05 r13 M3; the code had used ``round(fraction * N)`` over rescaled quotas).
The arithmetic is exact (``Fraction``), so ``10 * 0.3`` is 3 and not 3.0000000000000004. Within a stratum the first ``k`` sample IDs in ASCII order are holdout. The split consumes
no randomness: the split seed enters only through :func:`sample_id`.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Iterable, Mapping
from fractions import Fraction

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
    f = Fraction(str(fraction))  # the decimal as written, not its binary approximation
    quota = {key: size * f for key, size in sizes.items()}
    awarded = {key: math.floor(q) for key, q in quota.items()}
    leftover = math.ceil(total * f) - sum(awarded.values())
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


def unit_id(split_seed: int, library_id: str) -> str:
    """Opaque allocation-unit ID of a background donor library (protocol ``partitions.algorithm``).

    A real-background dataset splits donor libraries, not samples: the library id is hashed with the
    split seed so the sort order that picks holdout members carries no information about the library.
    """
    text = "|".join([str(split_seed), "library", library_id])
    return "u" + hashlib.sha256(text.encode()).hexdigest()[:12]


def split_libraries(
    libraries: Mapping[str, Iterable[str]], split_seed: int, fraction: float
) -> dict[str, str]:
    """Partition label of every donor library, drawn once per chemistry over the union of libraries.

    ``libraries`` maps chemistry -> library ids. Every sample and twin made from a library follows
    its label, so a library cannot be training in one dataset and holdout in another.
    """
    units = {c: {unit_id(split_seed, lib): lib for lib in libs} for c, libs in libraries.items()}
    # a chemistry with one library is an unsplittable stratum: training, as for samples
    splittable = {c: list(u) for c, u in units.items() if len(u) > 1}
    labels = split(splittable, fraction) if splittable else {}
    return {lib: labels.get(uid, TRAINING) for u in units.values() for uid, lib in u.items()}
