"""PLAN AUD-17 / M-02: per-marker polyA capture fields in the catalogue (default unknown)."""

from __future__ import annotations

import pandas as pd
import pytest

from viralscan.gene_programs import (
    POLYA_CAPTURE_VALUES,
    catalogue_virus_facts,
    load_catalogue,
    summarise_programs,
    validate_catalogue,
)

HSV1 = "Human herpesvirus 1"
EBV = "Epstein-Barr virus"
_EMPTY = pd.DataFrame(columns=["barcode", "virus_name", "state"])


def test_every_row_has_a_valid_value_and_a_basis_unless_unknown() -> None:
    rows = load_catalogue()
    for row in rows:
        assert row["polya_capture"] in POLYA_CAPTURE_VALUES
        if row["polya_capture"] != "unknown":
            assert row["capture_basis"].strip()
    assert validate_catalogue(rows) == []


def test_only_verified_markers_leave_unknown() -> None:
    by_key = {(r["virus"], r["refseq_gene"]): r for r in load_catalogue()}
    assert by_key[(HSV1, "LAT")]["polya_capture"] == "unknown"
    assert "38445887" in by_key[(HSV1, "LAT")]["capture_basis"]
    # ID-only HHV-6B records and the HSV-2 ncRNA have nothing to derive from
    assert {
        r["polya_capture"] for r in load_catalogue() if r["virus"] == "Human herpesvirus 6b"
    } == {"unknown"}
    assert by_key[("Human herpesvirus 2", "LAT")]["polya_capture"] == "unknown"
    assert by_key[(EBV, "BZLF1")]["polya_capture"] == "yes"


def test_hsv1_lat_is_exons_only_and_the_sole_latency_anchor() -> None:
    """PLAN CAT-29 (2026-10-12): the LAT record is the two exons, never the intron, and adds no second anchor."""
    rows = [r for r in load_catalogue() if r["virus"] == HSV1]
    latent = [r for r in rows if r["programme"] == "latent"]
    assert [r["refseq_gene"] for r in latent] == ["LAT"]
    assert "excludes the intron" in latent[0]["capture_basis"]
    assert {r["panel_completeness"] for r in rows} == {"partial"}
    assert {str(r["latency_observable_in_rna"]).lower() for r in rows} == {"false"}


def test_validation_rejects_bad_values_and_unbased_claims() -> None:
    rows = [dict(r) for r in load_catalogue()]
    rows[0]["polya_capture"] = "maybe"
    assert any("polya_capture" in e for e in validate_catalogue(rows))
    rows[0]["polya_capture"] = "no"
    rows[0]["capture_basis"] = ""
    assert any("needs a capture_basis" in e for e in validate_catalogue(rows))


def test_old_catalogue_without_the_columns_still_validates() -> None:
    rows = [
        {k: v for k, v in r.items() if k not in {"polya_capture", "capture_basis"}}
        for r in load_catalogue()
    ]
    assert validate_catalogue(rows) == []
    assert all("capture_uncertain" not in f for f in catalogue_virus_facts(rows).values())


def test_summary_caveat_names_uncaptured_markers_without_changing_state_counts() -> None:
    rows = load_catalogue()
    out = summarise_programs(_EMPTY, rows, viruses=[HSV1, EBV]).set_index("virus_name")
    assert "polyA capture is not established" in out.loc[HSV1, "caveat"]
    assert "LAT" in out.loc[HSV1, "caveat"]
    assert "polyA" not in out.loc[EBV, "caveat"]
    old = [{k: v for k, v in r.items() if k != "polya_capture"} for r in rows]
    before = summarise_programs(_EMPTY, old, viruses=[HSV1, EBV]).set_index("virus_name")
    state_cols = [c for c in out.columns if c.startswith("n_cells_")]
    assert out[state_cols].equals(before[state_cols])


@pytest.mark.parametrize("capture", ["unknown", "no"])
@pytest.mark.parametrize(
    "resolution",
    [
        None,
        {EBV: {"resolved": 17, "unresolved": 0}},
        {EBV: {"resolved": 16, "unresolved": 1, "absence_assessable": False}},
    ],
)
def test_summary_unreachable_state_uses_capture_even_without_resolution(
    capture, resolution
) -> None:
    rows = [dict(r) for r in load_catalogue()]
    for row in rows:
        if row["virus"] == EBV and row["refseq_gene"] == "BZLF1":
            row["polya_capture"] = capture
            row["capture_basis"] = "capture is not established"
    out = summarise_programs(_EMPTY, rows, viruses=[EBV], marker_resolution=resolution).iloc[0]
    assert out["n_cells_reactivating"] == "not_applicable"
    assert "absence-based states" in out["caveat"]
    assert "productive catalogue markers are unresolved" not in out["caveat"]


def test_productive_marker_without_established_capture_blocks_absence_states() -> None:
    """Same rule as an unresolved productive marker (AUD-19): reactivating/latent -> indeterminate."""
    import anndata as ad
    import numpy as np
    from scipy import sparse

    from viralscan.scripts.gene_programs import run_one

    def run(rows):
        ids = [r["gene_id_bundled"] for r in rows if r["virus"] == EBV]
        ie = next(
            r["gene_id_bundled"]
            for r in rows
            if r["virus"] == EBV and r["programme"] == "productive"
            and r["kinetic_class"] == "immediate_early"
        )  # fmt: skip
        dense = np.zeros((1, len(ids)), dtype=np.float32)
        dense[0, ids.index(ie)] = 3.0
        x = sparse.csr_matrix(dense)
        adata = ad.AnnData(x)
        adata.var_names, adata.obs_names = ids, ["c1"]
        adata.layers["counts_unique_viral"] = x.copy()
        return run_one(adata, [EBV], rows, 2, "bundled")["state"].tolist()

    rows = [dict(r) for r in load_catalogue()]
    assert run(rows) == ["reactivating"]
    for r in rows:
        if r["virus"] == EBV and r["refseq_gene"] == "BZLF1":
            r["polya_capture"], r["capture_basis"] = "unknown", ""
    assert run(rows) == ["indeterminate"]
