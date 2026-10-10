"""CAT-28: the panviral marker validation report (parsing and matching)."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "extras"))

import cat28_panviral_marker_check as cat28  # noqa: E402

# A tiny panviral table with every hazard: a repeated EC gene (BWRF1), the stray bracket,
# a fused name, a case-only near miss, and an accession our panel does not use (HHV7).
FIXTURE = (
    "EC\tGene\tNuccore\tVirus\n"
    "0\tBZLF1\tNC_007605.1\tHHV4\n"
    "1\tBWRF1\tNC_007605.1\tHHV4\n"
    "2\tBWRF1\tNC_007605.1\tHHV4\n"
    "3\tBBLF2/BBLF3\tNC_007605.1\tHHV4\n"
    "4\tBaRF1\tNC_007605.1\tHHV4\n"
    "5\tJvgp6\t[NC_001699.1\tJcpolyomavirus\n"
    "6\tU90\tU43400.1\tHHV7\n"
    "7\tU3\tAF157706.1\tHHV6B\n"
)


def _marker(virus: str, gene: str) -> dict[str, str]:
    return {"virus": virus, "refseq_gene": gene, "gene_id_bundled": f"X_{gene}"}


@pytest.fixture
def table(tmp_path: Path):
    path = tmp_path / "pan.tsv"
    path.write_text(FIXTURE)
    return cat28.read_panviral(path)


def test_bracket_is_repaired_and_counted(table) -> None:
    parsed, repaired = table
    assert repaired == 1
    assert ("Jcpolyomavirus", "NC_001699.1") in parsed
    assert not any("[" in acc for _, acc in parsed)


def test_ec_repeats_collapse_to_one_gene_with_a_count(table) -> None:
    parsed, _ = table
    assert parsed[("HHV4", "NC_007605.1")]["BWRF1"] == 2
    assert len(parsed[("HHV4", "NC_007605.1")]) == 4


def test_malformed_row_is_an_error_not_skipped(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("EC\tGene\tNuccore\tVirus\n0\tBZLF1\tNC_007605.1\n")
    with pytest.raises(ValueError, match="expected 4 fields"):
        cat28.read_panviral(path)


def test_statuses(table) -> None:
    parsed, _ = table
    rows = cat28.check_markers(
        [
            _marker("Epstein-Barr virus", "BZLF1"),  # matched
            _marker("Epstein-Barr virus", "BBLF2"),  # fused name: near miss only
            _marker("Epstein-Barr virus", "BARF1"),  # BARF1 is not BaRF1
            _marker("Epstein-Barr virus", "BaRF1"),  # matched
            _marker("Human herpesvirus 7", "U90"),  # gene there, accession differs
            _marker("Human herpesvirus 2", "UL9"),  # virus not in the table
            _marker("Human herpesvirus 6b", "U3"),  # matched on the pseudocontig
        ],
        parsed,
    )
    by = {(r["virus"], r["marker"]): r for r in rows}
    assert [r["status"] for r in rows] == [
        "matched",
        "unmatched",
        "unmatched",
        "matched",
        "accession_mismatch",
        "unmatched",
        "matched",
    ]
    assert "part of a fused name): BBLF2/BBLF3" in by[("Epstein-Barr virus", "BBLF2")]["detail"]
    assert "U43400.1" in by[("Human herpesvirus 7", "U90")]["detail"]
    assert "not in the panviral table" in by[("Human herpesvirus 2", "UL9")]["detail"]
    assert "pseudocontig" in by[("Human herpesvirus 6b", "U3")]["detail"]


def test_matching_is_case_sensitive(table) -> None:
    parsed, _ = table
    (row,) = cat28.check_markers([_marker("Epstein-Barr virus", "Bzlf1")], parsed)
    assert row["status"] == "unmatched"
    assert "case differs): BZLF1" in row["detail"]


def test_cli_writes_the_report(tmp_path: Path, capsys) -> None:
    pan = tmp_path / "pan.tsv"
    pan.write_text(FIXTURE)
    cat = tmp_path / "cat.tsv"
    cat.write_text("virus\trefseq_gene\tgene_id_bundled\nEpstein-Barr virus\tBZLF1\tX_BZLF1\n")
    out = tmp_path / "out" / "r.tsv"
    assert cat28.main(["--panviral", str(pan), "--catalogue", str(cat), "--out", str(out)]) == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert [r["status"] for r in rows] == ["matched"]
    assert "1 matched, 0 unmatched, 0 accession_mismatch" in capsys.readouterr().out


def test_committed_report_covers_every_catalogue_marker() -> None:
    """Staleness guard: the tracked report has one valid row per catalogue row."""
    report = list(
        csv.DictReader(
            (REPO_ROOT / "analysis" / "panel_expansion" / "cat28_panviral_marker_check.tsv").open(),
            delimiter="\t",
        )
    )
    catalogue = list(
        csv.DictReader(
            (REPO_ROOT / "src" / "viralscan" / "data" / "gene_programs.tsv").open(),
            delimiter="\t",
        )
    )
    assert {r["status"] for r in report} <= {"matched", "unmatched", "accession_mismatch"}
    assert [(r["virus"], r["marker"]) for r in report] == [
        (r["virus"], r["refseq_gene"]) for r in catalogue
    ]
