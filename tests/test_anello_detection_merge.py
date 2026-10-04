"""ANDET-09c: detection reads the branch's output and merges it into viral_summary.tsv.

The pure merge is covered in ``test_anello_align.py``; this pins the seam the
pipeline actually runs — ``anello_evidence`` locating and parsing the rule's
TSV, and ``write_tsv_outputs`` writing the merged table to disk.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from viralscan.anello_align import STATUS_NO_HOST_FILTER, STATUS_OK, SUMMARY_COLUMNS
from viralscan.scripts import detection
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

ANELLO = "Betatorquevirus"
OTHER = "Human gammaherpesvirus 4"


def _gene(gene_id, virus_name, family, accession):
    return GeneIdentity(
        gene_id=gene_id,
        genome_accession=accession,
        status="catalogued",
        viral=True,
        virus_key=f"genus:{virus_name}",
        virus_name=virus_name,
        family=family,
    )


@pytest.fixture
def identity():
    return VirusIdentityTable(
        genes=(
            _gene("AB1_gene1", ANELLO, "Anelloviridae", "AB1.1"),
            _gene("ALPHA_gene1", "Alphatorquevirus", "Anelloviridae", "AB2.1"),
            _gene("EBV_gp1", OTHER, "Orthoherpesviridae", "NC_007605.1"),
        )
    )


class _Config:
    def __init__(self, **kw):
        self.anello_align = kw.get("anello_align", True)
        self.host_index = kw.get("host_index", "/host")
        self.anello_index = kw.get("anello_index", "/anello_star")


def _write_branch_output(run: Path, rows):
    results = run / "results"
    results.mkdir(parents=True, exist_ok=True)
    with open(results / "anello_alignment_by_virus.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["virus_name", *SUMMARY_COLUMNS[2:]])
        for name, values in rows.items():
            writer.writerow([name, *values])


def _stats(name, molecules):
    return {
        "viral_molecules_total_est": molecules,
        "infected_cells": 2,
        "total_cells": 100,
        "pct_infected": 2.0,
        "viral_molecules_per_10k_est": 1.0,
        "n_called_cells": 50,
        "infected_called": 1,
        "pct_infected_called": 2.0,
        "infected_comparable": 1,
        "n_comparable_cells": 40,
        "pct_infected_comparable": 2.5,
        "accession_breadth": 1.0,
        "host_viral_ambig_fraction": 0.0,
    }


def _summary_rows(run: Path):
    with open(run / "results" / "viral_summary.tsv", encoding="utf-8") as fh:
        return {r["virus_name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_merge_writes_alignment_columns_and_an_alignment_only_row(tmp_path, identity):
    # Betatorquevirus: kallisto + alignment. Alphatorquevirus: alignment only.
    _write_branch_output(
        tmp_path,
        {
            # reads, unique, molecules, cells, identity, accessions, starts, homopolymer, splice
            ANELLO: [120, 90, 40, 12, 0.97, 3, 55, 0.01, 4],
            "Alphatorquevirus": [30, 25, 7, 5, 0.91, 1, 20, 0.0, 0],
        },
    )
    stats = {ANELLO: _stats(ANELLO, 11.0), OTHER: _stats(OTHER, 900.0)}
    evidence = detection.anello_evidence(_Config(), str(tmp_path), identity)
    assert evidence is not None
    detection.write_tsv_outputs(
        stats, __import__("pandas").DataFrame(), str(tmp_path), facts={}, anello=evidence
    )

    rows = _summary_rows(tmp_path)
    assert rows[OTHER]["detection_source"] == "kallisto"
    assert rows[OTHER]["alignment_reads"] == ""  # not an anellovirus row

    assert rows[ANELLO]["detection_source"] == "kallisto+alignment"
    assert rows[ANELLO]["alignment_status"] == STATUS_OK
    assert rows[ANELLO]["alignment_molecules_unique"] == "40"
    assert rows[ANELLO]["alignment_median_identity"] == "0.97"

    # The case the branch exists for: kallisto never called it.
    alpha = rows["Alphatorquevirus"]
    assert alpha["detection_source"] == "alignment_only"
    assert alpha["alignment_molecules_unique"] == "7"
    assert float(alpha["viral_molecules_total_est"]) == 0.0
    # Run-level denominators are carried over, not blanked.
    assert alpha["n_called_cells"] == "50"


def test_without_a_host_filter_the_status_is_recorded_and_no_row_is_added(tmp_path, identity):
    _write_branch_output(tmp_path, {"Alphatorquevirus": [30, 25, 7, 5, 0.9, 1, 2, 0.0, 0]})
    stats = {ANELLO: _stats(ANELLO, 11.0)}
    evidence = detection.anello_evidence(_Config(host_index=None), str(tmp_path), identity)
    detection.write_tsv_outputs(
        stats, __import__("pandas").DataFrame(), str(tmp_path), facts={}, anello=evidence
    )
    rows = _summary_rows(tmp_path)
    assert rows[ANELLO]["alignment_status"] == STATUS_NO_HOST_FILTER
    assert rows[ANELLO]["alignment_reads"] == ""
    assert "Alphatorquevirus" not in rows


def test_a_legacy_run_without_an_identity_table_keeps_the_old_schema(tmp_path):
    assert detection.anello_evidence(_Config(), str(tmp_path), None) is None
    detection.write_tsv_outputs(
        {ANELLO: _stats(ANELLO, 5.0)},
        __import__("pandas").DataFrame(),
        str(tmp_path),
        facts={},
        anello=None,
    )
    with open(tmp_path / "results" / "viral_summary.tsv", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    assert "detection_source" not in header
    assert header[-1] == "artifact_risk"


def test_an_enabled_branch_with_no_output_is_an_error_not_an_empty_result(tmp_path, identity):
    """A resumed pre-ANDET-09 run directory must not read as 'no anellovirus'."""
    (tmp_path / "results").mkdir()
    with pytest.raises(FileNotFoundError, match="--no-anello-align"):
        detection.anello_evidence(_Config(), str(tmp_path), identity)
