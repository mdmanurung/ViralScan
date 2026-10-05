"""SW-24: the molecule-mass audit tolerates float rounding at 10M-molecule scale."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from viralscan.multimapping import MoleculeAudit


def _audit(n: int) -> MoleculeAudit:
    return MoleculeAudit(
        input_molecules=n,
        resolved_molecules=n,
        unique_molecules=0,
        ambiguous_molecules=n,
        unresolved_molecules=0,
        ignored_read_multiplicity=0,
    )


def test_thirds_summed_over_millions_of_molecules_validate():
    n = 3_333_334
    shares = sparse.csr_matrix(np.full((n, 3), 1 / 3))
    _audit(n).validate(float(shares.sum()))


def test_observed_ebv_rounding_is_accepted_but_exceeds_absolute_tolerance():
    # Real EBV run: diff 1.86e-9 (relative 1.8e-16) failed under rtol=0, atol=1e-9.
    n, mass = 10_250_998, 10_250_998.000000002
    assert not np.isclose(mass, n, rtol=0.0, atol=1e-9)
    _audit(n).validate(mass)


def test_real_mass_discrepancy_still_raises():
    n = 10_250_998
    with pytest.raises(ValueError, match="Allocated ambiguous mass"):
        _audit(n).validate(n + 0.5)
