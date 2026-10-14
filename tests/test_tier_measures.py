"""TIER-02: the m, F and S measures behind the frozen probable grid (calibration.freeze_scope, Q135, Q141)."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan import tier_measures as tm
from viralscan.evidence_tiers import (
    STATUS_SELECTED,
    Calibration,
    Cutoffs,
    TierEvidence,
    assign_tier,
)
from viralscan.molecule_verdict import write_molecule_verdicts
from viralscan.run_safety import sha256_file

GENES = ["a1", "a2", "a3", "b1", "b2", "host1"]  # virus A: a1-a3, sibling B: b1-b2


def _build(rows, x):
    return ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=list(rows)),
        var=pd.DataFrame(index=GENES),
    )


def _evidence(
    out: Path,
    molecules: list[tuple[str, str, int, bool, int | None, int | None]],
    *,
    target: str | None = "A",
    competitor: bool = False,
) -> Path:
    """Evidence dir through the real ``write_molecule_verdicts``: (cb, ub, n_reads, counted, virus AS, host AS)."""
    out.mkdir(parents=True)
    sam: list[str] = []
    lineage = [["read_number", "read_id", "cb", "ub", "assigned_weight"]]
    n = 0
    for cb, ub, reads, counted, virus_as, host_as in molecules:
        for _ in range(reads):
            n += 1
            qname = f"{cb}_{ub}_{n}|r{n}"
            lineage.append([n, f"r{n}", cb, ub, "1.0" if counted else "0.0"])
            for flag, ref, score in ((0, "VIRUS|x", virus_as), (256, "HOST|chr1", host_as)):
                if score is not None:
                    sam.append(
                        "\t".join(
                            [qname, str(flag), ref, "1", "60", "50M", "*", "0", "0", "A", "I"]
                        )
                        + f"\tAS:i:{score}"
                    )
    with gzip.open(out / "read_lineage.tsv.gz", "wt", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(lineage)
    write_molecule_verdicts(sam, out / "read_lineage.tsv.gz", out)
    if target is not None:
        manifest = {"target": {"label": target, "genes": ["a1"]}, "outputs": {}}
        manifest["outputs"]["molecule_verdicts.tsv.gz"] = sha256_file(
            out / "molecule_verdicts.tsv.gz"
        )
        if competitor:
            manifest["competition"] = {"molecule_verdict": "not_produced_in_competitor_mode"}
        (out / "evidence_manifest.json").write_text(json.dumps(manifest))
    return out


# --- m -------------------------------------------------------------------------------------------------------


def test_m_keeps_fractional_weights_and_float32_thirds_meet_one_molecule() -> None:
    adata = _build(
        {"c1": [0.9, 0, 0, 0, 0, 0], "c2": [1 / 3, 1 / 3, 1 / 3, 0, 0, 5], "c3": [0] * 6},
        sp.csr_matrix(
            np.array(
                [[0.9, 0, 0, 0, 0, 0], [1 / 3, 1 / 3, 1 / 3, 0, 0, 5], [0] * 6], dtype=np.float32
            )
        ),
    )
    m = tm.molecules_per_cell(adata, ["a1", "a2", "a3"])
    assert m[0] == pytest.approx(0.9) and m[2] == 0.0  # measured zero, not missing
    one = Cutoffs(1, 1, 0.5)
    m_ok = tm.probable_cells(_measures(m, F=1, S=1.0), one)
    assert m_ok.tolist() == [
        False,
        True,
        False,
    ]  # 0.9 < 1; 1/3*3 in float32 is 0.99999994 and still counts


def _measures(m, F: int, S: float | None, n_genes: int = 3, called=None) -> tm.VirusMeasures:
    m = np.asarray(m, dtype=float)
    reason = None if S is not None else tm.REASON_NO_COUNTED
    return tm.VirusMeasures(
        m,
        tm.FeatureSupport(F, n_genes),
        tm.SupportFraction(S, reason),
        np.ones(m.shape, bool) if called is None else np.asarray(called, bool),
    )


# --- F -------------------------------------------------------------------------------------------------------


def test_f_needs_one_molecule_over_called_cells_not_merely_positive() -> None:
    #            a1   a2   a3   b1   b2  host1
    rows = {
        "c1": [0.5, 0.5, 0.0, 0, 0, 0],  # a1 0.5 + c2 0.5 = 1.0 -> counts; a2 0.5 alone -> does not
        "c2": [0.5, 0.0, 0.0, 0, 0, 0],
        "c3": [0.0, 0.0, 9.0, 0, 0, 0],  # a3 has molecules only in an empty droplet
    }
    adata = _build(rows, sp.csr_matrix(np.array(list(rows.values()), dtype=np.float32)))
    called = [True, True, False]
    fs = tm.features_supported(adata, ["a1", "a2", "a3", "a4_not_in_matrix"], called)
    assert fs == tm.FeatureSupport(
        1, 4
    )  # > 0 would give 2 (a1, a2); a3 is background; a4 still counts as an index gene
    assert tm.cutoff_unattainable(fs, 5) and not tm.cutoff_unattainable(fs, 4)


def test_f_half_half_split_over_two_genes_is_zero_features() -> None:
    adata = _build(
        {"c1": [0.5, 0.5, 0, 0, 0, 0]},
        sp.csr_matrix(np.array([[0.5, 0.5, 0, 0, 0, 0]], np.float32)),
    )
    assert tm.features_supported(adata, ["a1", "a2"], [True]).n_supported == 0


def test_f_and_m_are_per_virus_and_a_sibling_split_gives_each_its_weight() -> None:
    row = [0.5, 0, 0, 0.5, 0, 0]  # one molecule split a1 / b1 between sibling viruses
    adata = _build({"c1": row, "c2": row}, sp.csr_matrix(np.array([row, row], np.float32)))
    assert tm.molecules_per_cell(adata, ["a1", "a2", "a3"]).tolist() == [0.5, 0.5]
    assert tm.molecules_per_cell(adata, ["b1", "b2"]).tolist() == [0.5, 0.5]
    assert tm.features_supported(adata, ["a1", "a2", "a3"], [True, True]).n_supported == 1
    assert tm.features_supported(adata, ["b1", "b2"], [True, True]).n_supported == 1


def test_called_mask_shape_is_checked() -> None:
    adata = _build(
        {"c1": [1, 0, 0, 0, 0, 0]}, sp.csr_matrix(np.array([[1, 0, 0, 0, 0, 0]], np.float32))
    )
    with pytest.raises(ValueError, match="called_mask"):
        tm.features_supported(adata, ["a1"], [True, True])


# --- S -------------------------------------------------------------------------------------------------------


def test_s_counts_a_molecule_once_however_many_reads_carry_it(tmp_path: Path) -> None:
    ev = _evidence(
        tmp_path / "ev",
        [
            ("AAA", "U1", 5, True, 90, 60),  # five PCR copies of ONE virus_best molecule
            ("AAA", "U2", 1, True, 50, 80),  # one host_best molecule
        ],
    )
    s = tm.support_fraction(ev, "A", {"AAA"})
    assert (s.n_virus_best, s.n_host_best, s.n_tie) == (1, 1, 0)
    assert s.value == 0.5  # 5 reads as 1 molecule; read-counting would give 5/6


def test_s_population_is_counted_called_and_viral_aligned(tmp_path: Path) -> None:
    ev = _evidence(
        tmp_path / "ev",
        [
            ("AAA", "U1", 2, True, 90, 60),  # virus_best
            ("AAA", "U2", 1, True, 70, 70),  # tie counts against the virus
            ("AAA", "U3", 1, True, None, 80),  # counted, host-only: outside the population
            ("AAA", "U4", 1, False, 90, 60),  # not counted by the call
            ("EMPTY", "U5", 1, True, 90, 60),  # empty droplet: not a called cell
            (
                "CCC",
                "U6",
                1,
                True,
                90,
                None,
            ),  # virus_best, second called cell (pooled, not averaged)
        ],
    )
    s = tm.support_fraction(ev, "A", {"AAA", "CCC"})
    assert (s.n_virus_best, s.n_host_best, s.n_tie, s.n_excluded_host_only) == (2, 0, 1, 1)
    assert s.value == pytest.approx(2 / 3) and s.reason is None


@pytest.mark.parametrize(
    "molecules",
    [
        [("AAA", "U1", 1, True, None, 80)],  # only host-only molecules
        [("EMPTY", "U1", 1, True, 90, 60)],  # only empty droplets
        [("AAA", "U1", 1, False, 90, 60)],  # nothing counted
    ],
)
def test_s_zero_denominator_is_not_estimable_never_zero(tmp_path: Path, molecules) -> None:
    s = tm.support_fraction(
        _evidence(tmp_path / "ev", molecules), "A", {"AAA"}, obs_names={"AAA", "EMPTY"}
    )
    assert s.value is None and s.reason == tm.REASON_NO_COUNTED


def test_s_missing_evidence_is_reason_coded(tmp_path: Path) -> None:
    assert tm.support_fraction(None, "A", {"AAA"}).reason == tm.REASON_NO_EVIDENCE
    assert tm.support_fraction(tmp_path / "nope", "A", {"AAA"}).reason == tm.REASON_NO_EVIDENCE
    empty = tmp_path / "empty"
    empty.mkdir()
    s = tm.support_fraction(empty, "A", {"AAA"})
    assert s.value is None and s.reason == tm.REASON_NO_VERDICTS
    comp = _evidence(tmp_path / "comp", [("AAA", "U1", 1, True, 90, 60)], competitor=True)
    assert tm.support_fraction(comp, "A", {"AAA"}).reason == tm.REASON_COMPETITOR_MODE


def test_s_refuses_evidence_for_a_sibling_or_another_sample(tmp_path: Path) -> None:
    ev = _evidence(tmp_path / "ev", [("AAA", "U1", 1, True, 90, 60)], target="B")
    with pytest.raises(tm.EvidenceMismatch, match="not 'A'"):
        tm.support_fraction(ev, "A", {"AAA"})
    ev2 = _evidence(tmp_path / "ev2", [("AAA", "U1", 1, True, 90, 60)])
    with pytest.raises(tm.EvidenceMismatch, match="barcode"):
        tm.support_fraction(ev2, "A", {"ZZZ-1"}, obs_names={"ZZZ-1"})
    (ev2 / "molecule_verdicts.tsv.gz").write_bytes(b"tampered")
    with pytest.raises(tm.EvidenceMismatch, match="hash"):
        tm.support_fraction(ev2, "A", {"AAA"})


def test_s_collapses_a_repeated_cb_ub_key(tmp_path: Path) -> None:
    ev = _evidence(tmp_path / "ev", [("AAA", "U1", 1, True, 90, 60)], target=None)
    path = ev / "molecule_verdicts.tsv.gz"
    with gzip.open(path, "rt") as h:
        lines = h.read().splitlines()
    with gzip.open(
        path, "wt"
    ) as h:  # the same corrected CB-UMI twice, once with a better host score
        h.write("\n".join([*lines, lines[1].replace("\t60\t", "\t95\t")]) + "\n")
    s = tm.support_fraction(ev, "A", {"AAA"})
    assert (s.n_virus_best, s.n_host_best) == (0, 1)  # one molecule, host's best score wins


# --- the call ------------------------------------------------------------------------------------------------


def test_probable_cells_grid_and_equivalence_with_assign_tier(tmp_path: Path) -> None:
    adata = _build(
        {"c1": [2, 1, 0, 0, 0, 0], "c2": [1, 0, 0, 0, 0, 0], "c3": [9, 9, 9, 0, 0, 0]},
        sp.csr_matrix(
            np.array([[2, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [9, 9, 9, 0, 0, 0]], np.float32)
        ),
    )
    ev = _evidence(
        tmp_path / "ev",
        [
            ("c1", "U1", 1, True, 90, 60),
            ("c2", "U2", 1, True, 90, 60),
            ("c3", "U3", 1, True, 50, 80),
        ],
    )
    called = [True, True, False]  # c3 (huge, but empty droplet) feeds neither m, F nor S
    meas = tm.measure_virus(adata, "A", ["a1", "a2", "a3"], called, ev)
    assert meas.features == tm.FeatureSupport(2, 3) and meas.support.value == 1.0
    assert tm.probable_cells(meas, Cutoffs(2, 2, 0.9)).tolist() == [True, False, False]
    assert tm.probable_cells(meas, Cutoffs(1, 3, 0.5)).tolist() == [False] * 3  # F = 2 < 3
    for cut in (Cutoffs(1, 2, 0.9), Cutoffs(3, 1, 0.5)):
        cal = Calibration(STATUS_SELECTED, cut)
        sample_called = tm.probable_cells(meas, cut).any()  # D8
        ev_in = TierEvidence(
            float(meas.m[meas.called].max()), meas.features.n_supported, meas.support.value
        )
        assert (assign_tier(ev_in, cal) == "probable") == bool(sample_called)


def test_not_estimable_s_never_reaches_probable() -> None:
    meas = _measures([5.0, 5.0], F=3, S=None)
    assert not tm.probable_cells(meas, Cutoffs(1, 1, 0.5)).any()
    assert tm.probable_cells(
        _measures([5.0], F=3, S=0.5), Cutoffs(5, 3, 0.5)
    ).all()  # inclusive on all three


def test_detection_exposes_f_and_index_gene_count() -> None:
    from viralscan.scripts.detection import compute_stats

    adata = _build(
        {"c1": [1, 0, 0, 0, 0, 0], "c2": [0.5, 0, 0, 0, 0, 0], "c3": [0, 4, 0, 0, 0, 0]},
        sp.csr_matrix(
            np.array([[1, 0, 0, 0, 0, 0], [0.5, 0, 0, 0, 0, 0], [0, 4, 0, 0, 0, 0]], np.float32)
        ),
    )
    stats, _ = compute_stats(
        adata,
        {},
        {"A": ["a1"]},
        [],
        called_mask=np.array([True, True, False]),
        index_genes_by_virus={"A": ["a1", "a2", "a3"]},
    )
    assert stats["A"]["n_features_supported"] == 1 and stats["A"]["n_index_genes"] == 3
