"""ANDET-02 retains raw host loci without changing detection or EVE-risk semantics."""

from __future__ import annotations

import csv
import json
import subprocess

import pandas as pd
import pytest

from viralscan.scripts import build_reference
from viralscan.scripts.detection import write_tsv_outputs
from viralscan.virus_grouping import (
    HOST_HOMOLOGY_LOCUS_COLUMNS,
    host_homology_loci,
    virus_facts,
)
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

PAF = (
    "chr1\t5000\t10\t60\t+\tV1.1\t100\t0\t50\t45\t50\t60\n"
    "chr2\t9000\t100\t180\t-\tV1.1\t100\t20\t100\t60\t80\t40\n"
)


def test_parser_keeps_every_locus_in_both_coordinate_systems() -> None:
    out = build_reference._parse_host_homology_paf(PAF, {"V1.1": 100, "V2.1": 50})
    assert out["V1.1"]["host_homology_loci"] == [
        {
            "viral_start": 0,
            "viral_end": 50,
            "host_contig": "chr1",
            "host_start": 10,
            "host_end": 60,
            "strand": "+",
            "identity": 0.9,
            "aligned_bases": 50,
            "mapq": 60,
        },
        {
            "viral_start": 20,
            "viral_end": 100,
            "host_contig": "chr2",
            "host_start": 100,
            "host_end": 180,
            "strand": "-",
            "identity": 0.75,
            "aligned_bases": 80,
            "mapq": 40,
        },
    ]
    assert out["V1.1"]["host_homology_best_target"] == "chr2"
    assert out["V1.1"]["host_homology_max_identity"] == 0.75
    assert out["V1.1"]["host_homology_max_query_coverage"] == 0.8
    assert out["V2.1"]["host_homology_loci"] == []


@pytest.mark.parametrize(
    "column,value",
    [(2, "-1"), (3, "5001"), (4, "?"), (7, "100"), (8, "101"), (9, "51"), (10, "0"), (11, "256")],
)
def test_parser_rejects_invalid_locus_coordinates_and_metrics(column, value) -> None:
    fields = PAF.splitlines()[0].split("\t")
    fields[column] = value
    with pytest.raises(ValueError, match="Invalid host-homology PAF"):
        build_reference._parse_host_homology_paf("\t".join(fields), {"V1.1": 100})


def test_loci_survive_annotation_manifest_and_report_without_demoting_calls(tmp_path, monkeypatch):
    viral = tmp_path / "viral.fa"
    viral.write_text(">V1.1\n" + "ACGT" * 25 + "\n>V2.1\n" + "AC" * 25 + "\n")
    host = tmp_path / "host.fa"
    host.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(build_reference.shutil, "which", lambda name: "/mock/minimap2")
    monkeypatch.setattr(
        build_reference.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, stdout=PAF),
    )
    annotations = build_reference.measure_host_homology(viral, host, tmp_path / "annotations.tsv")
    with (tmp_path / "annotations.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert json.loads(rows[0]["host_homology_loci"]) == annotations["V1.1"]["host_homology_loci"]
    assert json.loads(rows[1]["host_homology_loci"]) == []
    monkeypatch.setattr(build_reference, "software_identity", lambda: {})
    manifest = build_reference.write_reference_manifest(
        viral,
        tmp_path / "reference_manifest.json",
        profile="curated",
        host_species="human",
        viral_identifiers={"V1.1", "V2.1"},
        annotations=annotations,
    )
    assert (
        json.loads(manifest.read_text())["sequences"][0]["host_homology_loci"]
        == (annotations["V1.1"]["host_homology_loci"])
    )
    table = VirusIdentityTable(
        (
            GeneIdentity("g1", "V1", "catalogued", True, "taxid:1", "Risky", risk_class="eve"),
            GeneIdentity("g2", "V1", "catalogued", True, "taxid:1", "Risky", risk_class="eve"),
            GeneIdentity("g3", "V2.1", "catalogued", True, "taxid:2", "Safe"),
        )
    )
    loci = host_homology_loci(table, tmp_path / "index.idx")
    assert len(loci) == 2
    assert loci[1]["accession_version"] == "V1.1"
    stats = {
        name: {
            "viral_molecules_total_est": 5,
            "infected_cells": 1,
            "total_cells": 10,
            "pct_infected": 10.0,
            "viral_molecules_per_10k_est": 1.0,
            "detected": True,
        }
        for name in ("Risky", "Safe")
    }
    write_tsv_outputs(
        stats, pd.DataFrame(), str(tmp_path), facts=virus_facts(table), homology_loci=loci
    )
    with (tmp_path / "results" / "host_homology_loci.tsv").open() as handle:
        reported = list(csv.DictReader(handle, delimiter="\t"))
    assert len(reported) == 2
    assert reported[0]["eve_risk"] == "True"
    assert reported[1]["host_start"] == "100" and reported[1]["strand"] == "-"
    summary = pd.read_csv(tmp_path / "results" / "viral_summary.tsv", sep="\t").set_index(
        "virus_name"
    )
    assert summary["detected"].all()
    assert summary.loc["Risky", "eve_risk"]
    assert not summary.loc["Safe", "eve_risk"]
    assert summary["viral_molecules_total_est"].tolist() == [5, 5]
    write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path), facts=virus_facts(table))
    assert (tmp_path / "results" / "host_homology_loci.tsv").read_text() == (
        "\t".join(HOST_HOMOLOGY_LOCUS_COLUMNS) + "\n"
    )


def test_older_or_missing_manifests_have_no_loci(tmp_path):
    table = VirusIdentityTable((GeneIdentity("g1", "V1", "catalogued", True, "taxid:1", "V"),))
    assert host_homology_loci(table, None) == []
    assert host_homology_loci(table, tmp_path / "index.idx") == []
    manifest = tmp_path / "reference_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "sequences": [
                    {"accession_version": "V1.1", "host_homology_status": "measured"},
                ]
            }
        )
    )
    assert host_homology_loci(table, tmp_path / "index.idx") == []
    manifest.write_text("incomplete json")
    assert host_homology_loci(table, tmp_path / "index.idx") == []
