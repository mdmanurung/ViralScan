"""Audit 4 programme fixes (PLAN AUD-19): D-01, D-02, D-03, D-18.

D-01  ``reactivating`` is gated by the catalogue ``latency_observable_in_rna`` flag.
D-02  a reduced marker set (PROG-14) cannot yield an absence-based call.
D-03  ``n_cells_reactivating`` is ``not_applicable`` where the state is unreachable.
D-18  one catalogue reader for per-virus facts; ``validate_catalogue`` checks them.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from viralscan.gene_programs import (
    Marker,
    call_cell_programme,
    catalogue_virus_facts,
    load_catalogue,
    summarise_programs,
    validate_catalogue,
)
from viralscan.scripts.gene_programs import run_one

CAT = load_catalogue()
EBV = "Epstein-Barr virus"
VIRUSES = sorted({r["virus"] for r in CAT})


def _ids(virus: str) -> list[str]:
    return [r["gene_id_bundled"] for r in CAT if r["virus"] == virus]


def _ie(virus: str) -> list[str]:
    return [
        r["gene_id_bundled"]
        for r in CAT
        if r["virus"] == virus and r["programme"] == "productive"
        and r["kinetic_class"] == "immediate_early"
    ]  # fmt: skip


def _run(virus: str, hot: list[str], index: list[str] | None = None):
    ids = index if index is not None else _ids(virus)
    dense = np.zeros((1, len(ids)), dtype=np.float32)
    for gene in hot:
        dense[0, ids.index(gene)] = 3.0
    x = sparse.csr_matrix(dense)
    adata = ad.AnnData(x)
    adata.var_names = ids
    adata.obs_names = ["c1"]
    adata.layers["counts_unique_viral"] = x.copy()
    return run_one(adata, [virus], CAT, 2, "bundled")


def _summary(df, virus: str):
    return summarise_programs(
        df, CAT, viruses=[virus], marker_resolution=df.attrs["marker_resolution"]
    ).iloc[0]


NOT_OBSERVABLE_WITH_IE = [
    v for v, f in catalogue_virus_facts(CAT).items() if not f["latency_observable"] and _ie(v)
]


class TestD01FlagGatesReactivating:
    def test_the_viruses_named_by_the_audit_are_covered(self) -> None:
        # VZV, HHV-6A/6B, HHV-7 and HCMV: the flag, not a name list, selects them.
        assert len(NOT_OBSERVABLE_WITH_IE) >= 5

    @pytest.mark.parametrize("virus", NOT_OBSERVABLE_WITH_IE)
    def test_ie_only_cell_is_never_reactivating_when_latency_is_unobservable(
        self, virus: str
    ) -> None:
        df = _run(virus, [_ie(virus)[0]])
        assert df["state"].tolist() == ["indeterminate"]

    def test_ebv_ie_only_cell_is_still_reactivating(self) -> None:
        assert _run(EBV, [_ie(EBV)[0]])["state"].tolist() == ["reactivating"]

    def test_unit_level_flag_overrides_everything(self) -> None:
        ie = Marker("ie", "productive", "g1", True, True, "immediate_early")
        late = Marker("late", "productive", "g2", True, True, "late")
        matrix = sparse.csr_matrix(np.array([[1.0, 0.0]], dtype=np.float32))
        assert call_cell_programme(matrix, [ie, late], latency_observable=False)[0]["state"] == (
            "indeterminate"
        )
        assert call_cell_programme(matrix, [ie, late])[0]["state"] == "reactivating"


class TestD02ReducedMarkerSet:
    def test_one_gene_index_cannot_call_reactivating(self) -> None:
        gene = _ie(EBV)[0]
        df = _run(EBV, [gene], index=[gene])
        assert df["state"].tolist() == ["indeterminate"]
        row = _summary(df, EBV)
        assert row["n_cells_reactivating"] == "not_applicable"
        assert "unresolved" in row["caveat"] and "reactivating" in row["caveat"]

    def test_dropping_one_later_productive_marker_degrades_the_call(self) -> None:
        later = next(
            r["gene_id_bundled"]
            for r in CAT
            if r["virus"] == EBV
            and r["programme"] == "productive"
            and r["kinetic_class"] != "immediate_early"
        )
        index = [g for g in _ids(EBV) if g != later]
        df = _run(EBV, [_ie(EBV)[0]], index=index)
        assert df["state"].tolist() == ["indeterminate"]
        assert _summary(df, EBV)["n_cells_reactivating"] == "not_applicable"

    def test_dropping_a_latent_marker_keeps_the_call(self) -> None:
        """Only missing *productive* markers make the absence unobservable."""
        latent = next(
            r["gene_id_bundled"] for r in CAT if r["virus"] == EBV and r["programme"] == "latent"
        )
        index = [g for g in _ids(EBV) if g != latent]
        df = _run(EBV, [_ie(EBV)[0]], index=index)
        assert df["state"].tolist() == ["reactivating"]
        assert _summary(df, EBV)["n_cells_reactivating"] == 1

    def test_latent_is_not_called_without_the_productive_markers(self) -> None:
        latent = next(
            r["gene_id_bundled"] for r in CAT if r["virus"] == EBV and r["programme"] == "latent"
        )
        productive = {
            r["gene_id_bundled"]
            for r in CAT
            if r["virus"] == EBV and r["programme"] == "productive"
        }
        assert _run(EBV, [latent])["state"].tolist() == ["latent"]
        reduced = [g for g in _ids(EBV) if g not in productive]
        assert _run(EBV, [latent], index=reduced)["state"].tolist() == ["indeterminate"]

    def test_positive_productive_call_survives_a_reduced_set(self) -> None:
        productive = [
            r["gene_id_bundled"]
            for r in CAT
            if r["virus"] == EBV and r["programme"] == "productive"
        ]
        df = _run(EBV, productive[:3], index=productive[:3])
        assert df["state"].tolist() == ["productive"]


class TestD03NotApplicableColumn:
    @pytest.mark.parametrize("virus", VIRUSES)
    def test_value_matches_reachability(self, virus: str) -> None:
        facts = catalogue_virus_facts(CAT)[virus]
        row = _summary(_run(virus, []), virus)
        if virus == EBV:
            assert facts["reactivation_reachable"]
            assert row["n_cells_reactivating"] == 0
        else:
            assert not facts["reactivation_reachable"]
            assert row["n_cells_reactivating"] == "not_applicable"

    @pytest.mark.parametrize(
        "virus", ["Human herpesvirus 1", "Human herpesvirus 2", "Human herpesvirus 8"]
    )
    def test_no_productive_ie_marker_viruses(self, virus: str) -> None:
        assert not _ie(virus)
        assert _summary(_run(virus, []), virus)["n_cells_reactivating"] == "not_applicable"

    def test_other_count_columns_stay_integer(self) -> None:
        row = _summary(_run("Human herpesvirus 1", []), "Human herpesvirus 1")
        for col in (
            "n_cells_productive",
            "n_cells_latent",
            "n_cells_mixed",
            "n_cells_indeterminate",
        ):
            assert isinstance(row[col], (int, np.integer))


class TestD18OneCatalogueReader:
    @staticmethod
    def _flipped(column: str, value: str) -> list[dict]:
        cat = [dict(r) for r in CAT if r["virus"] == EBV]
        cat[-1][column] = value
        return cat

    @pytest.mark.parametrize(
        ("column", "value"),
        [("latency_observable_in_rna", "false"), ("panel_completeness", "partial")],
    )
    def test_validate_catalogue_rejects_disagreeing_rows(self, column: str, value: str) -> None:
        problems = validate_catalogue(self._flipped(column, value))
        assert any("disagree" in p and column in p and EBV in p for p in problems), problems

    def test_packaged_catalogue_agrees(self) -> None:
        assert not [p for p in validate_catalogue(CAT) if "disagree" in p]

    def test_run_one_and_summary_both_refuse_an_inconsistent_catalogue(self) -> None:
        cat = self._flipped("latency_observable_in_rna", "false")
        with pytest.raises(ValueError, match="disagree"):
            catalogue_virus_facts(cat)
        with pytest.raises(ValueError, match="disagree"):
            summarise_programs(None, cat, viruses=[EBV])
        x = sparse.csr_matrix(np.ones((1, 1), dtype=np.float32))
        adata = ad.AnnData(x)
        adata.var_names = ["x"]
        adata.layers["counts_unique_viral"] = x.copy()
        with pytest.raises(ValueError, match="disagree"):
            run_one(adata, [EBV], cat, 2, "bundled")
