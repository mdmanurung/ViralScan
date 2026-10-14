"""ADJ-01: scripts/adjudicate_candidates.py on a synthetic round (offline, no binaries, no real data).

Thresholds come from tests/data/adjudication_lock_TEST_ONLY.json: arbitrary test numbers, not a protocol lock.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path

import pytest

from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "adjudicate_candidates", ROOT / "scripts" / "adjudicate_candidates.py"
)
adj = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(adj)

LOCK = ROOT / "tests" / "data" / "adjudication_lock_TEST_ONLY.json"
EBV, HPV, TTV, H6B = (
    "Human gammaherpesvirus 4",
    "Alphapapillomavirus 9",
    "Torque teno virus 1",
    "HHV-6B",
)
CELLS = ["A", "B", "C", "D"]
LIB = "ds__S1"

# virus -> dict(unique, allocated, loci, cells, vb, hb, family, sibling, programme)
GOOD = dict(
    unique=40, allocated=2, loci=4, cells=6, vb=30, hb=1, family="", sibling="", programme=None
)


def _evidence(root: Path, virus: str, vb: int, hb: int, complex_body: float = 0.9) -> None:
    d = root / f"ds__S1__combined_off__{adj.cards.dossier.slug(virus)}"
    d.mkdir(parents=True)
    (d / "alignment_qc.tsv").write_text(
        "count_layer\treference\treference_class\treads\tcomplex_body_fraction\tbreadth_1x\t"
        f"max_50bp_window_fraction\nraw\tVIRUS|x\tvirus\t20\t{complex_body}\t0.4\t0.1\n"
    )
    (d / "interpretation_flags.tsv").write_text(
        "flag\tstatus\nhost_homology\tnot_flagged\nlow_complexity\tnot_flagged\n"
    )
    (d / "viral_reads.fasta").write_text(">r\nACGT\n")
    rows = [("A", "virus_best")] * vb + [("B", "host_best")] * hb
    with gzip.open(d / "molecule_verdicts.tsv.gz", "wt") as fh:
        fh.write("cb\tub\treads\tcounted\tvirus_as\thost_as\tverdict\n")
        for i, (cb, v) in enumerate(rows):
            fh.write(f"{cb}\tU{i}\t1\tTrue\t\t\t{v}\n")
    (d / "molecule_verdict_summary.tsv").write_text(
        f"scope\tverdict\tmolecules\ncounted\tvirus_best\t{vb}\ncounted\thost_best\t{hb}\n"
        "counted\ttie\t0\ncounted\tunaligned\t0\n"
    )
    (d / "read_verdict_summary.tsv").write_text(
        "verdict\trepeat_class\treads\nhost_best\ttelomeric\t3\nvirus_best\tnone\t7\n"
    )


def _round(
    tmp_path: Path, viruses: dict[str, dict], frozen: bool = True, with_new_cols: bool = True
) -> Path:
    inner = tmp_path / "runs" / "ds" / "combined_off" / "S1" / "S1"
    res = inner / "results"
    res.mkdir(parents=True)
    names = list(viruses)
    (res / "viral_summary.tsv").write_text(
        "virus_name\tviral_molecules_total_est\tdetected\n"
        + "".join(f"{v}\t{viruses[v]['unique'] + viruses[v]['allocated']}\tTrue\n" for v in names)
    )
    cols = "virus_name\tgene_id\tviral_molecules_unique\tviral_molecules_ambiguous_allocated\tn_unique_viral_cells"
    if with_new_cols:
        cols += "\tlocus_id\tn_unique_positive_cells_called"
    lines = []
    for v in names:
        n = max(viruses[v]["loci"], 1)
        for i in range(n):
            line = (
                f"{v}_g{i}\t{viruses[v]['unique'] / viruses[v]['loci'] if viruses[v]['loci'] else 0}\t"
                f"{viruses[v]['allocated'] / n}\t{viruses[v]['cells']}"
            )
            if (
                with_new_cols
            ):  # `locus_of` maps gene index -> locus (default: one locus per gene); `pos` = union
                line += f"\t{v}_L{viruses[v].get('locus_of', list(range(n)))[i]}\t{viruses[v].get('pos', viruses[v]['cells'])}"
            lines.append(f"{v}\t{line}\n")
    (res / "multimap_evidence.tsv").write_text(cols + "\n" + "".join(lines))
    VirusIdentityTable(
        tuple(
            GeneIdentity(
                f"{v}_g0",
                "acc",
                "catalogued",
                True,
                viruses[v].get("key", f"k:{v}"),
                v,
                family=viruses[v]["family"],
                sibling_group=viruses[v]["sibling"],
            )  # fmt: skip
            for v in names
        )
    ).write_tsv(res / "virus_identity.tsv")
    prog = [(v, viruses[v]["programme"]) for v in names if viruses[v]["programme"] is not None]
    (res / "gene_program_summary.tsv").write_text(
        "virus_name\tn_cells_productive\tn_cells_latent\tn_cells_mixed\tn_cells_reactivating\t"
        "panel_completeness\n"
        + "".join(
            f"{viruses[v].get('prog_name', v)}\t{n}\t0\t0\tnot_applicable\tcomplete\n"
            for v, n in prog
        )
        + "".join(
            f"{v}\t0\t0\t0\tnot_applicable\tnot_applicable\n"
            for v in names
            if viruses[v]["programme"] is None
        )
    )
    if frozen:
        (tmp_path / "reference_cells").mkdir()
        (tmp_path / "reference_cells" / "ds__S1.tsv").write_text(
            "barcode\n" + "".join(f"{c}\n" for c in CELLS)
        )
    for v in names:
        if viruses[v]["vb"] is not None:
            _evidence(tmp_path / "evidence", v, viruses[v]["vb"], viruses[v]["hb"])
    return tmp_path


def _run(
    tmp_path: Path, viruses: dict[str, dict], lock: Path | None = LOCK, **kw
) -> dict[str, dict]:
    r = _round(tmp_path, viruses, kw.pop("frozen", True), kw.pop("with_new_cols", True))
    rows = adj.adjudicate_round(r, adj.load_lock(lock), **kw)
    return {row["virus_name"]: row for row in rows}


def test_l3_for_a_gene_programme_herpesvirus(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 5}})[EBV]
    assert (got["level"], got["level_reason"]) == ("L3", "meets_l3")
    assert (
        got["flags"] == "telomeric;test_only_thresholds"
    )  # repeat class is a label, gates nothing
    assert got["unique_molecules"] == 40 and got["n_loci"] == 4
    assert got["read_virus_best_fraction"] == pytest.approx(30 / 31, abs=1e-6)
    assert got["thresholds_sha256"] and len(got["thresholds_sha256"]) == 64


def test_l3_when_summary_uses_the_legacy_catalogue_name_q58(tmp_path):
    """Identity-table name differs from gene_programs.tsv ("Epstein-Barr virus"): matched through the taxid key."""
    ebv = {**GOOD, "programme": 5, "key": "taxid:10376", "prog_name": "Epstein-Barr virus"}
    got = _run(tmp_path, {EBV: ebv})[EBV]
    assert (got["level"], got["level_reason"], got["gene_programme_cells"]) == ("L3", "meets_l3", 5)


def test_name_mismatch_off_the_taxid_map_still_fails_closed(tmp_path):
    off = {**GOOD, "programme": 5, "key": "taxid:99999", "prog_name": "Epstein-Barr virus"}
    got = _run(tmp_path, {EBV: off})[EBV]
    assert (got["level"], got["level_reason"]) == ("L2", "cap:no_gene_programme")


def test_taxid_map_agrees_with_the_catalogues():
    """The map covers exactly the viruses of gene_programs.tsv, with taxids that exist in virus_catalog.tsv."""
    data = ROOT / "src" / "viralscan" / "data"
    with (data / "gene_programs.tsv").open(newline="") as fh:
        names = {r["virus"] for r in csv.DictReader(fh, delimiter="\t")}
    with (data / "virus_catalog.tsv").open(newline="") as fh:
        taxids = {r["taxid"] for r in csv.DictReader(fh, delimiter="\t")}
    assert set(adj.PROGRAMME_NAME_BY_KEY.values()) == names
    assert {k.split(":")[1] for k in adj.PROGRAMME_NAME_BY_KEY} <= taxids


def test_flags_mark_test_only_thresholds(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 5}})[EBV]
    assert "test_only_thresholds" in got["flags"].split(";")


def test_herpes_with_too_few_programme_cells_stops_at_l2(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 1}})[EBV]
    assert (got["level"], got["level_reason"]) == ("L2", "below_l3:l3_min_programme_cells")


def test_non_herpesvirus_caps_at_l2_with_a_visible_reason(tmp_path):
    got = _run(tmp_path, {HPV: GOOD})[HPV]
    assert (got["level"], got["level_reason"]) == ("L2", "cap:no_gene_programme")


def test_anellovirus_caps_at_l1_and_is_screening_only(tmp_path):
    got = _run(tmp_path, {TTV: {**GOOD, "family": "Anelloviridae"}})[TTV]
    assert (got["level"], got["level_reason"]) == ("L1", "cap:anellovirus_screening_only")
    assert "screening_only" in got["flags"].split(";")


def test_anellovirus_without_identity_table_still_screening_only(tmp_path):
    """No identity table: the genus-name fallback keeps the cap (fails closed)."""
    r = _round(tmp_path, {TTV: GOOD})
    (r / "runs/ds/combined_off/S1/S1/results/virus_identity.tsv").unlink()
    row = adj.adjudicate_round(r, adj.load_lock(LOCK))[0]
    assert (row["level"], row["virus_key"]) == ("L1", f"name:{TTV}")
    assert "screening_only" in row["flags"]


def test_allocated_mass_alone_never_leaves_l0(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "unique": 0, "allocated": 500, "programme": 5}})[EBV]
    assert (got["level"], got["level_reason"]) == ("L0", "allocated_only")


def test_sibling_bleed_blocks_l1_and_is_flagged(tmp_path):
    bleed = {**GOOD, "unique": 10, "allocated": 30, "sibling": "HHV-6"}  # 75 % allocated
    got = _run(tmp_path, {H6B: bleed})[H6B]
    assert (got["level"], got["level_reason"]) == ("L0", "below_l1:max_sibling_allocated_fraction")
    assert "sibling_bleed" in got["flags"].split(";")
    assert got["max_sibling_allocated_fraction"] == 0.75
    # a virus with no sibling group is not judged on its allocated share
    assert (
        _run(tmp_path / "b", {HPV: {**GOOD, "unique": 10, "allocated": 30}})[HPV]["level"] != "L0"
    )


@pytest.mark.parametrize(
    "vb, hb, reason",
    [
        (5, 5, "below_l1:min_virus_best_fraction"),
        (None, 0, "no_read_level_evidence"),
    ],
)
def test_read_level_gates(tmp_path, vb, hb, reason):
    got = _run(tmp_path, {HPV: {**GOOD, "vb": vb, "hb": hb}})[HPV]
    assert (got["level"], got["level_reason"]) == ("L0", reason)


def test_unrestricted_molecule_scope_cannot_promote(tmp_path):
    got = _run(tmp_path, {HPV: GOOD}, frozen=False)[HPV]
    assert (got["level"], got["level_reason"]) == ("L0", "read_level_not_restricted")


def test_distributed_route_minima_decide_l1_vs_l2(tmp_path):
    for key in ("unique", "cells", "loci"):
        weak = {**GOOD, key: {"unique": 5, "cells": 1, "loci": 1}[key]}
        got = _run(tmp_path / key, {HPV: weak})[HPV]
        assert got["level"] == "L1" and got["level_reason"].startswith("below_l2:"), key


def test_no_lock_is_l0_for_everything_but_metrics_are_filled(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 5}}, lock=None)[EBV]
    assert (got["level"], got["level_reason"], got["thresholds_sha256"]) == (
        "L0", "no_protocol_lock", "",
    )  # fmt: skip
    assert got["unique_molecules"] == 40 and got["read_virus_best_fraction"] != ""


def test_incomplete_lock_fails_closed_and_names_the_key(tmp_path):
    doc = json.loads(LOCK.read_text())
    del doc["distributed_route"]["min_independent_loci"]
    doc["max_sibling_allocated_fraction"] = "0.3"  # a string is not a number
    lock = tmp_path / "lock.json"
    lock.write_text(json.dumps({"adjudication": doc}))  # also accepted nested under `adjudication`
    got = _run(tmp_path / "r", {EBV: {**GOOD, "programme": 5}}, lock=lock)[EBV]
    assert got["level"] == "L0"
    assert got["level_reason"] == "incomplete_protocol_lock:distributed_route.min_independent_loci"
    assert got["thresholds_sha256"]  # still the hash of the file that was supplied


def test_malformed_lock_is_an_operator_error(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    with pytest.raises(SystemExit) as exc:
        adj.load_lock(bad)
    assert exc.value.code == 2


def test_host_expression_never_changes_a_call(tmp_path):
    base = _run(tmp_path / "a", {HPV: GOOD})[HPV]
    r = _round(tmp_path / "b", {HPV: GOOD})
    inner = r / "runs/ds/combined_off/S1/S1/results"
    (inner / "hostresponse").mkdir()
    (inner / "hostresponse" / "isg_score.tsv").write_text("virus_name\tisg\n" + f"{HPV}\t99\n")
    text = (inner / "viral_summary.tsv").read_text().splitlines()
    (inner / "viral_summary.tsv").write_text(
        "\n".join([text[0] + "\thost_expression_score"] + [t + "\t99" for t in text[1:]]) + "\n"
    )
    again = adj.adjudicate_round(r, adj.load_lock(LOCK))[0]
    assert (again["level"], again["level_reason"]) == (base["level"], base["level_reason"])


def test_germline_flag_is_passed_through_not_detected(tmp_path):
    flags = tmp_path / "g.tsv"
    flags.write_text(
        f"library_id\tvirus_key\tgermline_integration_compatible\n{LIB}\tk:{HPV}\ttrue\n"
    )
    got = _run(
        tmp_path / "r", {HPV: GOOD, EBV: {**GOOD, "programme": 5}}, germline=adj.load_flags(flags)
    )
    assert "germline_integration_compatible" in got[HPV]["flags"].split(";")
    assert "germline_integration_compatible" not in got[EBV]["flags"]


def test_cli_contract_metadata_sorting_and_no_l4(tmp_path):
    r = _round(
        tmp_path,
        {HPV: GOOD, EBV: {**GOOD, "programme": 5}, TTV: {**GOOD, "family": "Anelloviridae"}},
    )
    meta = tmp_path / "meta.tsv"
    meta.write_text(
        f"library_id\tsubject_id\ttimepoint\trun_id\tlane\textra\n{LIB}\tP1\tT0\tR9\tL2\tx\n"
    )
    out = tmp_path / "out" / "levels.tsv"
    assert adj.main([str(r), "--lock", str(LOCK), "--metadata", str(meta), "--out", str(out)]) == 0
    with out.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames[:13] == [
        "library_id", "virus_key", "level", "level_reason", "flags", "unique_molecules", "n_loci",
        "read_virus_best_fraction", "thresholds_sha256", "subject_id", "timepoint", "run_id", "lane",
    ]  # fmt: skip
    assert [x["virus_key"] for x in rows] == sorted(x["virus_key"] for x in rows)
    assert {x["subject_id"] for x in rows} == {"P1"} and {x["lane"] for x in rows} == {"L2"}
    assert {x["level"] for x in rows} <= {"L0", "L1", "L2", "L3"}
    assert sorted(x["level"] for x in rows) == ["L1", "L2", "L3"]
    # TIER-01 / Q104: one appended column, the explicit level -> package-tier analogue
    assert reader.fieldnames[-1] == "package_tier_analogue"
    assert {x["level"]: x["package_tier_analogue"] for x in rows} == {
        "L1": "probable",
        "L2": "strong",
        "L3": "strong",
    }


def test_every_level_the_ladder_emits_has_a_package_tier_analogue():
    import inspect
    import re

    from viralscan.evidence_tiers import STUDY_LEVEL_TO_PACKAGE_TIER

    source = inspect.getsource(adj.adjudicate) + inspect.getsource(adj._cap)
    literal = {int(n) for n in re.findall(r"return (\d),|_cap\((\d),", source) for n in n if n}
    assert literal and literal <= {0, 1, 2, 3}
    assert {f"L{n}" for n in range(4)} == set(
        STUDY_LEVEL_TO_PACKAGE_TIER
    )  # L4/L5 are never emitted


def test_metadata_duplicate_library_is_refused(tmp_path):
    meta = tmp_path / "m.tsv"
    meta.write_text("library_id\tsubject_id\nx\ta\nx\tb\n")
    with pytest.raises(SystemExit):
        adj.load_metadata(meta)


def test_candidate_only_rows_are_not_adjudicated(tmp_path):
    r = _round(tmp_path, {HPV: GOOD})
    summary = r / "runs/ds/combined_off/S1/S1/results/viral_summary.tsv"
    summary.write_text(summary.read_text().replace("True", "False"))
    assert adj.adjudicate_round(r, adj.load_lock(LOCK)) == []


def test_overlapping_genes_are_one_locus_and_cells_are_a_union(tmp_path):
    # LOCK-01 D1/D2: 4 genes in 2 loci (not 4 loci); 9 positive cells although no gene has more than 6.
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 5, "locus_of": [0, 0, 1, 1], "pos": 9}})[EBV]
    assert (got["n_loci"], got["positive_cells"]) == (2, 9)
    assert (got["level"], got["level_reason"]) == (
        "L3",
        "meets_l3",
    )  # TEST_ONLY lock: 2 loci, 3 cells suffice


def test_loci_below_the_minimum_after_merging_overlaps(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 5, "locus_of": [0, 0, 0, 0]}})[EBV]
    assert got["n_loci"] == 1
    assert (got["level"], got["level_reason"]) == ("L1", "below_l2:min_independent_loci")


def test_run_without_locus_columns_is_not_guessed(tmp_path):
    got = _run(tmp_path, {EBV: {**GOOD, "programme": 5}}, with_new_cols=False)[EBV]
    assert (got["n_loci"], got["positive_cells"]) == ("", "")
    assert (got["level"], got["level_reason"]) == ("L1", "missing:min_positive_cells")
