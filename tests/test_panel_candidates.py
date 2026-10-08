"""PANEL-01: scripts/panel_candidates.py classifies and excludes candidates on a synthetic catalogue."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "panel_candidates", Path(__file__).resolve().parents[1] / "scripts" / "panel_candidates.py"
)
pc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pc)


def _row(accession, species, panel="max", host="", family="Adenoviridae", risk="", segment=""):
    return {
        "accession": accession,
        "species": species,
        "family": family,
        "genus": "",
        "host": host,
        "segment": segment,
        "panel": panel,
        "risk_class": risk,
    }


def _evonk(accession, name, group="evonk-refseq", twin=""):
    return {"accession": accession, "name": name, "group": group, "genbank_twin": twin}


def _by_acc(rows):
    return {r["accession"]: r for r in rows}


def test_relevance_classes_and_never_automatic_for_unreviewed(tmp_path: Path) -> None:
    catalogue = [
        _row("NC_000001", "Human mastadenovirus A", host="Homo sapiens"),
        _row("NC_000002", "Human mastadenovirus A", host=""),  # same species as an H1 record
        _row("NC_000003", "Human adenovirus 54", host=""),  # name-based only
        _row("NC_000004", "Bat adenovirus 1", host=""),  # nothing says human
        _row("NC_000005", "Already shipped virus", panel="shipped", host="Homo sapiens"),
    ]
    curated_path = tmp_path / "curated.tsv"
    curated_path.write_text(
        "key\tkey_type\tclass\tstatus\tbasis\treviewer\n"
        "NC_900001\taccession\tH3\taccepted\tvector\tme\n"
        "NC_900002\taccession\tH2\tproposed\tmaybe\t\n"
    )
    evonk = [
        _evonk("NC_900001", "AAV"),
        _evonk("NC_900002", "Maybe virus"),
        _evonk("NC_900003", "Animal"),
    ]

    rows = _by_acc(pc.build_candidates(catalogue, evonk, pc.load_curated(curated_path), {}, {}))

    assert "NC_000005" not in rows  # shipped rows are not candidates
    assert rows["NC_000001"]["relevance"] == "H1"
    assert rows["NC_000002"]["relevance"] == "H2_species"
    assert rows["NC_000003"]["relevance"] == "H2_name"
    assert rows["NC_000004"]["relevance"] == "unreviewed"
    assert rows["NC_900001"]["relevance"] == "H3_curated"
    assert rows["NC_900002"]["relevance"] == "proposed"  # proposed is not accepted
    assert rows["NC_900003"]["relevance"] == "unreviewed"


def test_exclusions_and_twin_check_do_not_read_unchecked_as_no_twin() -> None:
    catalogue = [
        _row("NC_000010", "Torque teno virus", family="Anelloviridae", host="Homo sapiens"),
        _row("NC_000011", "Some EVE virus", host="Homo sapiens", risk="eve"),
        _row("NC_000012", "Twin of panel", host="Homo sapiens"),
        _row("NC_000013", "Unchecked", host="Homo sapiens"),
    ]
    evonk = [_evonk("NC_900010", "HPV alias", group="evonk-hpv-alias", twin="GQ000001")]
    panel_md5 = {"AB000001": "aaa"}
    candidate_md5 = {"NC_000012": "aaa"}

    rows = _by_acc(pc.build_candidates(catalogue, evonk, {}, panel_md5, candidate_md5))

    assert rows["NC_000010"]["exclusion_reason"] == "anellovirus_max_not_promoted"
    assert rows["NC_000011"]["exclusion_reason"] == "eve_risk"
    assert rows["NC_000012"]["exclusion_reason"] == "sequence_twin_of:AB000001"
    assert rows["NC_000012"]["twin_checked"] == "yes"
    assert rows["NC_000013"]["twin_checked"] == "no"
    assert rows["NC_000013"]["exclusion_reason"] == ""
    assert rows["NC_900010"]["exclusion_reason"] == "refseq_alias_of_genbank_hpv:GQ000001"


def test_output_is_deterministic_and_fasta_md5_is_case_insensitive(tmp_path: Path) -> None:
    fasta = tmp_path / "x.fa"
    fasta.write_text(">NC_1.1 desc\nacgt\nACGT\n>NC_2.1\nACGTACGT\n")
    md5 = pc.fasta_md5(fasta)
    assert md5["NC_1"] == md5["NC_2"]

    out = []
    for name in ("a.tsv", "b.tsv"):
        path = tmp_path / name
        pc.main(
            [
                "--evonk",
                str(tmp_path / "none.tsv"),
                "--curated",
                str(tmp_path / "none2.tsv"),
                "--out",
                str(path),
            ]
        )
        out.append(path.read_bytes())
    assert out[0] == out[1]
    header = next(csv.reader(out[0].decode().splitlines(), delimiter="\t"))
    assert header == pc.COLUMNS
