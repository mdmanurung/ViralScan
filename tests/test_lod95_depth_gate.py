"""LOD95-DEPTH-01: pin what the LOD95 depth term is today, per index kind.

``detection.py`` passes ``sum(adata.X)`` over all barcodes as ``depth`` to
``sensitivity_record``. Under a virus-only index that sum is the viral molecule total,
so ``lod95_per_10k`` is per 10k viral molecules, not per 10k host molecules (the
yardstick in ``sensitivity.py``). Which denominator is right is a user decision
(PLAN LOD95-DEPTH-01), so these tests describe current behaviour and change with it.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.runconfig import RunConfig
from viralscan.scripts.detection import build_sensitivity_records
from viralscan.sensitivity import LOD95_MOLECULES, lod95_per_10k

HOST_MOLECULES = 3_000_000.0
VIRAL_MOLECULES = 3.0


def _adata(index_kind, viral=VIRAL_MOLECULES):
    """One barcode; virus V holds *viral* molecules; a host gene holds 3M unless virus-only."""
    virus_only = index_kind == "virus_only"
    X = [[viral]] if virus_only else [[HOST_MOLECULES, viral]]
    adata = ad.AnnData(
        X=sp.csr_matrix(np.array(X, dtype=np.float32)),
        obs=pd.DataFrame(index=["c0"]),
        var=pd.DataFrame(index=["v1"] if virus_only else ["h1", "v1"]),
    )
    if index_kind:
        adata.uns["index_kind"] = index_kind
    return adata


def _record(index_kind, viral=VIRAL_MOLECULES):
    stats = {"V": {"viral_molecules_total_est": viral}}
    (rec,) = build_sensitivity_records(
        _adata(index_kind, viral), stats, RunConfig(), index_viruses=["V"]
    )
    return rec


def test_depth_is_the_whole_matrix_sum_so_it_changes_with_index_kind():
    # the record keeps no depth field; recover it from the per-10k LOD (capture 1, threshold 1)
    def depth_of(rec):
        return LOD95_MOLECULES * 1e4 / rec.lod95_per_10k

    assert depth_of(_record(None)) == pytest.approx(HOST_MOLECULES + VIRAL_MOLECULES)
    assert depth_of(_record("virus_only")) == pytest.approx(VIRAL_MOLECULES)


def test_virus_only_lod95_per_10k_is_per_viral_molecule_not_per_host_molecule():
    # Pinned current behaviour (unguarded): the virus-only record divides by the viral total.
    rec = _record("virus_only")
    assert rec.lod95_per_10k == pytest.approx(lod95_per_10k(VIRAL_MOLECULES))
    assert rec.lod95_per_10k > 1e3
    assert not rec.depth_sufficient


def test_combined_lod95_per_10k_uses_host_plus_viral_depth():
    rec = _record(None)
    assert rec.lod95_per_10k == pytest.approx(lod95_per_10k(HOST_MOLECULES + VIRAL_MOLECULES))
    assert rec.depth_sufficient


def test_lod95_molecules_does_not_depend_on_depth():
    # 2.996 * threshold / capture: only the per-10k figure and the depth gate inherit
    # the denominator question.
    assert _record(None).lod95_molecules == _record("virus_only").lod95_molecules
    assert _record(None).lod95_molecules == pytest.approx(LOD95_MOLECULES)


def test_a_deep_virus_only_run_clears_the_depth_gate_on_viral_molecules_alone():
    # The unsafe direction from review R3 F1: enough viral molecules clears the gate with
    # no host depth in view. Pinned as current behaviour, not endorsed.
    rec = _record("virus_only", viral=4_000_000.0)
    assert rec.depth_sufficient and rec.depth_only_sufficient
