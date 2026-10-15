"""AUD3-07: pin three documented-but-unchanged behaviours (P9-01, P3-05, P1-01).

These tests record CURRENT behaviour so that any future change is deliberate; none of them
asserts that the behaviour is ideal.
"""

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from viralscan.defaults import DEFAULTS
from viralscan.evidence_tiers import (
    STATUS_SELECTED,
    TIER_CANDIDATE,
    TIER_PROBABLE,
    Calibration,
    Cutoffs,
    TierEvidence,
    assign_tier,
)
from viralscan.multimapping import build_multimap_layers, em_gene_abundances
from viralscan.scripts.detection import compute_stats, detect_genes
from viralscan.utils import infected_masks

# --- P9-01: em-global stop rule leaves a residue on a zero-unique-support gene ------------------

HOST_UNIQUE, AMBIGUOUS = (
    1e6,
    4e6,
)  # the audited extreme ratio; EM takes float counts, so this is instant


def _em_residue(tol=1e-6, unique=HOST_UNIQUE, ambiguous=AMBIGUOUS):
    diag: dict = {}
    theta = em_gene_abundances(
        {(0, 1): ambiguous}, np.array([unique, 0.0]), 1.0, 100, tol, diagnostics=diag
    )
    return ambiguous * theta[1] / theta.sum(), diag  # allocated viral molecules


def test_em_global_stops_at_iteration_two_and_leaves_residue():
    mass, diag = _em_residue()
    assert (diag["iterations"], diag["converged"]) == (2, True)
    assert 2.0 < mass < 3.0  # audited 2.56; viral gene has zero unique support
    # the residue is real but bounded: it never approaches a second molecule of signal
    assert mass < 5.0


def test_em_global_residue_shrinks_only_with_tighter_tolerance():
    loose, _ = _em_residue(tol=1e-4)
    default, _ = _em_residue(tol=1e-6)
    tight, diag = _em_residue(tol=1e-9)
    assert loose == default  # tol 1e-4..1e-6 all stop at iteration 2
    assert diag["iterations"] > 2
    assert tight < 0.02 < default  # audited 0.0077 at 1e-9


def test_em_global_residue_absent_at_realistic_crossmap_sizes():
    mass, _ = _em_residue(ambiguous=1e3)
    assert mass < 1e-6


def test_em_global_residue_does_not_reach_detected_true():
    mass, _ = _em_residue()
    n_cells = 5000  # audited end-to-end: max per-cell mass 5.9e-4, no cell reaches 1 molecule
    rng = np.random.default_rng(0)
    per_cell = rng.dirichlet(np.ones(n_cells)) * mass
    assert per_cell.max() < 0.01
    X = sparse.csr_matrix(np.column_stack([np.full(n_cells, 5.0), per_cell]))
    genes = ["host", "viral"]
    found = detect_genes(
        genes, X, ["viral"], DEFAULTS["detection_threshold"], groups={"V": ["viral"]}
    )
    assert found == {"viral": float(per_cell.sum())}  # the sample screen passes ...
    infected, subthreshold = infected_masks(per_cell, DEFAULTS["infected_min_molecules"])
    assert not infected.any() and subthreshold.sum() == n_cells  # ... but no cell is infected
    adata = ad.AnnData(
        X=X.astype(float),
        obs=pd.DataFrame(index=[f"C{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=genes),
    )
    stats, _ = compute_stats(adata, found, {"V": ["viral"]}, list(found))
    assert stats["V"]["detected"] is False  # CANDIDATE ONLY, not a detection
    assert stats["V"]["infected_called"] == 0


# --- P3-05: tier cut-offs carry a 1e-6 relative slack ------------------------------------------


def test_assign_tier_relative_slack_on_molecule_cutoff():
    calib = Calibration(status=STATUS_SELECTED, probable=Cutoffs(2, 1, 0.5))

    def tier(m):
        return assign_tier(TierEvidence(m, 1, 1.0), calib)

    assert tier(1.999999) == TIER_PROBABLE  # inside 2 * (1 - 1e-6): reaches probable at M=2
    assert tier(1.99999) == TIER_CANDIDATE  # outside the slack


# --- P1-01: ECs with no gene are skipped from the molecule intersection ------------------------


def test_geneless_and_absent_ecs_are_skipped_not_empty_sets():
    # EC 1 = viral gene; EC 3 has no gene; EC 9 is absent from ec_map. Skipping (not "empty set")
    # means a molecule seen in EC 1 + a gene-less EC stays unique viral instead of an unresolved collision.
    ec_map = {1: [1], 2: [0], 3: []}
    rows = [("A", "M1", 1, 1), ("A", "M1", 3, 1), ("A", "M2", 1, 1), ("A", "M2", 9, 1)]
    rows += [("A", "M3", 3, 1)]  # only gene-less EC: unresolved, no allocation
    rows += [("A", "M4", 2, 1), ("A", "M4", 3, 1)]  # host + gene-less EC: still unique host
    frame = pd.DataFrame(rows, columns=["barcode", "umi", "ec", "count"])
    layers = build_multimap_layers(frame, {"A": 0}, ec_map, 1, 2, {1}, None)
    unique = layers.unique.toarray()
    assert unique.tolist() == [[1.0, 2.0]]  # M4 -> host, M1 + M2 -> viral
    assert layers.corrected.sum() == 0
    assert layers.audit.unresolved_molecules == 1  # M3 only
