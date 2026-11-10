"""LOD95-DEPTH-01: the LOD95 depth term per index kind (user decision 2026-11-09, option a).

``detection.py`` passes ``sum(adata.X)`` as ``depth``. Under a virus-only index, or when
host and viral columns cannot be told apart, that sum is not a host-molecule total, so
``lod95_per_10k``, ``depth_sufficient``, ``depth_only_sufficient`` and
``lod_interpretation`` are NA / not certifiable and no row is an informative negative.
``lod95_molecules`` does not use depth and is unchanged. A combined index is unchanged.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.runconfig import RunConfig
from viralscan.scripts.detection import build_sensitivity_records, build_sensitivity_table
from viralscan.sensitivity import LOD95_MOLECULES, SENSITIVITY_COLUMNS, lod95_per_10k

HOST_MOLECULES = 3_000_000.0
VIRAL_MOLECULES = 3.0
NA_FIELDS = (
    "lod95_per_10k",
    "depth_sufficient",
    "depth_only_sufficient",
    "lod_interpretation",
    "expected_molecules_at_1_per_10k",
    "p_detect_at_1_per_10k",
)


def _adata(index_kind, viral=VIRAL_MOLECULES, identity_source=None):
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
    if identity_source:
        adata.uns["gene_identity_source"] = identity_source
    return adata


def _record(index_kind, viral=VIRAL_MOLECULES, identity_source=None, detected=True):
    stats = {"V": {"viral_molecules_total_est": viral}} if detected else {}
    (rec,) = build_sensitivity_records(
        _adata(index_kind, viral, identity_source), stats, RunConfig(), index_viruses=["V"]
    )
    return rec


@pytest.mark.parametrize(
    "kwargs, reason",
    [
        ({"index_kind": "virus_only"}, "virus_only"),
        ({"index_kind": None, "identity_source": "unavailable"}, "unavailable"),
    ],
)
def test_depth_fields_are_na_when_depth_is_not_host_depth(kwargs, reason):
    rec = _record(**kwargs)
    row = rec.as_row()
    assert all(row[f] == "NA" for f in NA_FIELDS)
    assert rec.depth_sufficient is None and rec.depth_only_sufficient is None
    assert np.isnan(rec.lod95_per_10k)
    assert rec.lod_interpretation == "NA"
    assert reason in rec.depth_na_reason
    assert "NA, not certifiable" in row["notes"] and reason in row["notes"]


def test_virus_only_deep_run_is_still_not_certifiable():
    # The unsafe direction (review R3 F1): enough viral molecules no longer clears a gate.
    rec = _record("virus_only", viral=4_000_000.0)
    assert rec.depth_sufficient is None and rec.depth_only_sufficient is None
    assert rec.lod_interpretation == "NA"


@pytest.mark.parametrize(
    "kwargs",
    [{"index_kind": "virus_only"}, {"index_kind": None, "identity_source": "unavailable"}],
)
def test_informative_negative_is_never_true_when_depth_is_na(kwargs):
    rec = _record(viral=4_000_000.0, detected=False, **kwargs)
    assert not rec.informative and not rec.sensitivity_eligible
    assert rec.as_row()["informative_negative"] is False
    assert "depth-not-certifiable" in rec.negative_blockers


def test_lod95_molecules_does_not_depend_on_depth():
    # 2.996 * threshold / capture: the depth gate does not touch it.
    for rec in (_record(None), _record("virus_only"), _record(None, identity_source="unavailable")):
        assert rec.lod95_molecules == pytest.approx(LOD95_MOLECULES)
        assert rec.as_row()["lod95_molecules"] == pytest.approx(LOD95_MOLECULES)


def test_combined_index_is_unchanged():
    # host columns present: depth is host + viral, the gate and bands work as before
    for source in (None, "virus_identity"):
        rec = _record("combined", identity_source=source)
        row = rec.as_row()
        assert rec.depth_na_reason is None
        assert rec.lod95_per_10k == pytest.approx(lod95_per_10k(HOST_MOLECULES + VIRAL_MOLECULES))
        assert row["lod95_per_10k"] == rec.lod95_per_10k
        assert row["depth_sufficient"] is True and row["depth_only_sufficient"] is True
        assert row["lod_interpretation"] == "informative"
        assert "depth-not-certifiable" not in rec.negative_blockers


def test_combined_undetected_without_capture_is_blocked_by_capture_only():
    rec = _record("combined", detected=False)
    assert rec.depth_sufficient and not rec.informative
    assert rec.negative_blockers == ("capture-not-measured",)


def test_no_index_kind_recorded_is_unchanged():
    # legacy matrices carry no index_kind: treated as before
    rec = _record(None)
    assert rec.depth_na_reason is None and rec.depth_sufficient is True


def test_table_keeps_its_columns_and_shows_na():
    table = build_sensitivity_table(
        _adata("virus_only"),
        {"V": {"viral_molecules_total_est": 3.0}},
        RunConfig(),
        index_viruses=["V"],
    )
    assert tuple(table.columns) == SENSITIVITY_COLUMNS
    assert table.loc[0, list(NA_FIELDS)].tolist() == ["NA"] * len(NA_FIELDS)


def test_yardstick_fields_are_na_for_virus_only_and_numeric_for_combined():
    na = _record("virus_only")
    assert np.isnan(na.expected_molecules_if_present_at_1_per_10k)
    assert np.isnan(na.p_detect_if_present_at_1_per_10k)
    ok = _record("combined")
    row = ok.as_row()
    assert row["expected_molecules_at_1_per_10k"] == pytest.approx(
        (HOST_MOLECULES + VIRAL_MOLECULES) * 1e-4
    )
    assert 0 < row["p_detect_at_1_per_10k"] <= 1
