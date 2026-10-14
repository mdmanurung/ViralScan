"""CARD-01: scripts/dsr03_evidence_cards.py on a tiny synthetic round (offline, no real data)."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path

import anndata
import numpy as np
import pytest
import scipy.sparse as sp

from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

_spec = importlib.util.spec_from_file_location(
    "dsr03_evidence_cards",
    Path(__file__).resolve().parents[1] / "scripts" / "dsr03_evidence_cards.py",
)
card = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(card)

EBV = "Human gammaherpesvirus 4"  # the dataset's expected virus
HHV6 = "Human betaherpesvirus 6B"  # not expected in an EBV dataset
CELLS = ["A", "B", "C", "D"]


def _evidence(root, virus, virus_best, host_best, off_list_host=0):
    d = root / f"ebv__S1__combined_off__{card.dossier.slug(virus)}"
    d.mkdir(parents=True)
    (d / "alignment_qc.tsv").write_text(
        "count_layer\treference\treference_class\treads\tcomplex_body_fraction\t"
        "breadth_1x\tmax_50bp_window_fraction\ndeduplicated\tVIRUS|x\tvirus\t20\t0.9\t0.42\t0.12\n"
    )
    (d / "interpretation_flags.tsv").write_text(
        "flag\tstatus\nhost_homology\tnot_flagged\nlow_complexity\tnot_flagged\n"
        "coverage_hotspot\tnot_flagged\neve_or_integration_like_hotspot\tnot_flagged\n"
    )
    (d / "viral_reads.fasta").write_text(">r\nACGT\n")
    rows = [("A", "virus_best")] * virus_best + [("B", "host_best")] * host_best
    rows += [("EMPTY", "host_best")] * off_list_host  # barcodes outside the frozen cells
    with gzip.open(d / "molecule_verdicts.tsv.gz", "wt") as fh:
        fh.write("cb\tub\treads\tcounted\tvirus_as\thost_as\tverdict\n")
        for i, (cb, v) in enumerate(rows):
            fh.write(f"{cb}\tU{i}\t1\tTrue\t\t\t{v}\n")
    (d / "molecule_verdict_summary.tsv").write_text(
        "scope\tverdict\tmolecules\n"
        f"counted\tvirus_best\t{virus_best}\ncounted\thost_best\t{host_best + off_list_host}\n"
        "counted\ttie\t0\ncounted\tunaligned\t0\n"
    )
    (d / "read_verdict_summary.tsv").write_text(
        "verdict\trepeat_class\treads\nhost_best\ttelomeric\t3\nvirus_best\tnone\t7\n"
    )
    return d


def _round(tmp_path, summary_rows):
    r = tmp_path
    inner = r / "runs" / "ebv" / "combined_off" / "S1" / "S1"
    (inner / "results").mkdir(parents=True)
    (inner / "results" / "viral_summary.tsv").write_text(
        "virus_name\tviral_molecules_total_est\tdetected\n"
        + "".join(f"{v}\t{n}\t{det}\n" for v, n, det in summary_rows)
    )
    (inner / "results" / "called_cells.tsv").write_text(
        "barcode\n" + "".join(f"{c}\n" for c in CELLS)
    )
    genes = [("EBV_g1", EBV), ("HHV6_g1", HHV6), ("host_g", "")]
    VirusIdentityTable(
        tuple(
            GeneIdentity(g, "acc", "catalogued", bool(n), f"k:{n}" if n else "", n)
            for g, n in genes
        )
    ).write_tsv(inner / "results" / "virus_identity.tsv")
    (r / "reference_cells").mkdir()
    (r / "reference_cells" / "ebv__S1.tsv").write_text(
        "barcode\n" + "".join(f"{c}\n" for c in CELLS)
    )
    # cells: 4 x ~1000 UMIs, EBV 5/4/3/0; pool: 6 empty barcodes x 100 host UMIs, one EBV UMI
    barcodes = CELLS + [f"E{i}" for i in range(6)]
    x = np.zeros((10, 3))
    x[:4, 2] = 1000
    x[:4, 0] = [5, 4, 3, 0]
    x[4:, 2] = 100
    x[4, 0] = 1
    ad = anndata.AnnData(sp.csr_matrix(x))
    ad.obs_names = barcodes
    ad.var_names = [g for g, _ in genes]
    counts = inner / "kb-python" / "counts_unfiltered"
    counts.mkdir(parents=True)
    ad.write_h5ad(counts / "adata.h5ad")
    return r


def _by_virus(cards):
    return {c["virus"]: c for c in cards}


def test_expected_virus_reaches_read_level_pass_with_all_tiers_reported(tmp_path):
    r = _round(tmp_path, [(EBV, 12.0, "True")])
    _evidence(r / "evidence", EBV, virus_best=20, host_best=1, off_list_host=50)
    got = _by_virus(card.build_cards(r))[EBV]
    assert got["card_status"] == "read_level_pass"
    # 50 off-list host molecules would sink the call in an all-cells tally
    assert got["molecule_scope"] == "restricted" and got["host_best_molecules"] == 1
    assert got["reads_telomeric"] == 3 and got["reads_total"] == 10
    assert got["genome_structure"] == "reported" and got["best_breadth_1x"] == 0.42
    assert got["ambient"] == "reported" and got["observed_in_cells"] == 12.0
    assert got["expected_ambient_in_cells"] == pytest.approx(1 / 601 * 4012, abs=1e-3)
    assert float(got["ambient_poisson_p"]) == pytest.approx(0.04, abs=0.005)  # P(X >= 12 | 6.68)
    assert got["recurrence"] == "n_a" and got["permutation"] == "n_a"


def test_unexpected_virus_is_unverified_even_when_its_read_level_passes(tmp_path):
    r = _round(tmp_path, [(HHV6, 8.0, "True")])
    _evidence(r / "evidence", HHV6, virus_best=20, host_best=0)
    got = card.build_cards(r)[0]
    assert got["read_level"] == "pass"  # the tier itself is shown...
    assert got["card_status"] == "unverified" and "no calibration record" in got["status_reason"]
    # ...and a passed gate under the current rule values lifts the interim rule
    record = tmp_path / "cal.json"
    record.write_text(json.dumps({"rule": card.calibration.current_rule(), "gate_passed": True}))
    assert card.build_cards(r, calibration_record=record)[0]["card_status"] == "read_level_pass"
    record.write_text(json.dumps({"rule": {"min_virus_best_fraction": 0.9}, "gate_passed": True}))
    assert card.build_cards(r, calibration_record=record)[0]["card_status"] == "unverified"


def test_read_level_failure_and_missing_evidence_and_unrestricted(tmp_path):
    r = _round(tmp_path, [(EBV, 12.0, "True"), (HHV6, 8.0, "True")])
    _evidence(r / "evidence", EBV, virus_best=1, host_best=9)
    got = _by_virus(card.build_cards(r))
    assert got[EBV]["card_status"] == "read_level_fail"
    assert got[HHV6]["card_status"] == "unverified" and got[HHV6]["read_level"] == "no_evidence"
    (r / "reference_cells" / "ebv__S1.tsv").unlink()  # nothing frozen: no cell-level claim
    again = _by_virus(card.build_cards(r))[EBV]
    assert again["card_status"] == "not_restricted" and again["ambient"] == "n_a"


def test_candidate_only_rows_get_no_card(tmp_path):
    """AUD-21: a viral_summary row with detected = false is not a call."""
    r = _round(tmp_path, [(EBV, 12.0, "True"), (HHV6, 9.0, "False")])
    assert [c["virus"] for c in card.build_cards(r)] == [EBV]


def test_recurrence_comes_from_the_dossier_table(tmp_path):
    r = _round(tmp_path, [(EBV, 12.0, "True")])
    (r / "dossier").mkdir()
    (r / "dossier" / "recurrence.tsv").write_text(
        f"virus\tn_samples\tn_datasets\tmax_cells_in_one_sample\tsamples\n{EBV}\t3\t2\t10\tx\n"
    )
    got = card.build_cards(r)[0]
    assert (got["recurrence_samples"], got["recurrence_datasets"]) == ("3", "2")


def test_cell_type_permutation_with_labels(tmp_path):
    r = _round(tmp_path, [(EBV, 12.0, "True")])
    (r / "cell_level").mkdir()
    with gzip.open(r / "cell_level" / "ebv__S1.tsv.gz", "wt") as fh:
        fh.write(
            "arm\tbarcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction\n"
        )
        for b in ("A", "B", "C"):
            fh.write(f"combined_off\t{b}\t{EBV}\t3\t1000\t0.003\n")
    labels = tmp_path / "labels.tsv"
    labels.write_text("barcode\tlabel\nA\tB_cell\nB\tB_cell\nC\tB_cell\nD\tT_cell\n")
    got = card.build_cards(r, labels={"ebv/S1": labels}, n_perm=200, seed=1)[0]
    assert got["permutation"] == "reported" and got["perm_top_label"] == "B_cell"
    assert float(got["perm_top_label_p_bonferroni"]) <= 1.0
    again = card.build_cards(r, labels={"ebv/S1": labels}, n_perm=200, seed=1)[0]
    assert again["perm_top_label_p_bonferroni"] == got["perm_top_label_p_bonferroni"]


def test_cli_writes_tsv_and_markdown(tmp_path):
    r = _round(tmp_path, [(EBV, 12.0, "True")])
    _evidence(r / "evidence", EBV, virus_best=20, host_best=1)
    assert card.main([str(r)]) == 0
    with (r / "cards" / "evidence_cards.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [x["card_status"] for x in rows] == ["read_level_pass"]
    assert EBV in (r / "cards" / "evidence_cards.md").read_text()
