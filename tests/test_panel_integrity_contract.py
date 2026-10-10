"""Independent software gates for panel integrity; no index build or calibration."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from viralscan import anellovirus
from viralscan.scripts import detection
from viralscan.virus_grouping import virus_facts
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "panel_integrity", ROOT / "scripts" / "panel_integrity.py"
)
panel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(panel)
BUILDER_SPEC = importlib.util.spec_from_file_location(
    "bundled_builder", ROOT / "scripts" / "build_bundled_panel_ref.py"
)
builder = importlib.util.module_from_spec(BUILDER_SPEC)
BUILDER_SPEC.loader.exec_module(builder)


def _gtf(seq="NC_000001.1", gene="G1", tx="T1", kind="protein_coding", feature="exon"):
    return (
        f"{seq}\tfixture\t{feature}\t1\t9\t.\t+\t.\t"
        f'gene_id "{gene}"; transcript_id "{tx}"; gene_biotype "{kind}";\n'
    )


def test_gtf_real_structure_and_explicit_whole_genome_tag_only(tmp_path):
    path = tmp_path / "a.gtf"
    path.write_text(_gtf() + _gtf(gene="G2", tx="T2", kind="whole_genome"))
    seq_genes, whole, errors = panel.inspect_gtfs([path])
    assert errors == []
    assert seq_genes == {"NC_000001.1": {"G1", "G2"}}
    assert whole == {"G2"}  # Identical spans alone do not reclassify real CDS.


@pytest.mark.parametrize(
    "change,message",
    [
        (lambda s: s.replace("\t1\t9\t", "\t0\t9\t"), "invalid coordinates"),
        (lambda s: s.replace('gene_id "G1";', ""), "missing gene_id"),
        (lambda s: s.replace("\texon\t", "\tCDS\t"), "has no exon"),
        (lambda s: s + s, "duplicate GTF record"),
        (lambda s: s.replace("\texon\t", "\texon\tEXTRA\t"), "nine GTF columns"),
    ],
)
def test_malformed_duplicate_and_unindexed_transcripts_fail(tmp_path, change, message):
    path = tmp_path / "a.gtf"
    path.write_text(change(_gtf()))
    assert any(message in error for error in panel.inspect_gtfs([path])[2])


def test_duplicate_genes_and_transcripts_across_files_fail(tmp_path):
    first, second = tmp_path / "a.gtf", tmp_path / "b.gtf"
    first.write_text(_gtf())
    second.write_text(_gtf(seq="NC_000002.1"))
    errors = panel.inspect_gtfs([first, second])[2]
    assert any("duplicate gene_id" in error for error in errors)
    assert any("inconsistent transcript" in error for error in errors)


def test_shipped_coverage_role_sibling_and_reverse_catalogue_check():
    rows = [
        {
            "accession": "A1",
            "species": "Human herpesvirus 6A",
            "panel": "shipped",
            "sibling_group": "HHV-6",
            "role": "target",
        },
        {
            "accession": "A2",
            "species": "Human herpesvirus 6b",
            "panel": "shipped",
            "sibling_group": "WRONG",
            "role": "typo",
        },
        {"accession": "A3", "species": "other", "panel": "shipped"},
    ]
    matrix = panel.catalogue_matrix(rows, {"A1.1": {"G1"}, "UNKNOWN.1": {"G2"}}, {"A3"})
    by_acc = {row["accession"]: row for row in matrix}
    assert "inconsistent_sibling_group" in by_acc["A1"]["errors"]
    assert "invalid_role" in by_acc["A2"]["errors"]
    assert "missing_gtf_seqname" in by_acc["A2"]["errors"]
    assert by_acc["A3"]["status"] == "excluded"
    assert by_acc["UNKNOWN"]["errors"] == "gtf_accession_not_catalogued"


@pytest.mark.parametrize(
    "summary,size,discarded",
    [
        ({"max_ec_size": 3, "discarded_ec_count": 2}, 2, 2),
        ({"max_ec_size": 2, "discarded_ec_count": 3}, 2, 2),
        ({"max_ec_size": True, "discarded_ec_count": 0}, 2, 2),
        ({"max_ec_size": 1}, 2, 2),
        ({"max_ec_size": 1, "discarded_ec_count": -1}, 2, 2),
    ],
)
def test_ec_summary_exceeding_budget_or_missing_metrics_fails(summary, size, discarded):
    with pytest.raises(ValueError):
        panel.check_ec_budget(summary, size, discarded)


def test_ec_budget_accepts_exact_boundary():
    panel.check_ec_budget({"max_ec_size": 2, "discarded_ec_count": 0}, 2, 0)


def test_builder_manifest_check_runs_before_fetches_and_fails_closed(tmp_path, monkeypatch):
    calls = []

    def run(command, check):
        calls.append(command)
        assert check is True
        raise RuntimeError("stale manifest")

    monkeypatch.setattr(builder.subprocess, "run", run)
    monkeypatch.setattr(builder, "_find_repo_root", lambda: tmp_path)
    monkeypatch.setattr(builder.sys, "path", list(builder.sys.path))
    monkeypatch.setattr(
        builder.sys,
        "argv",
        [
            "builder",
            "--out",
            str(tmp_path / "out"),
            "--ncbi-email",
            "test@example.test",
            "--gtf-manifest",
            str(tmp_path / "manifest.tsv"),
        ],
    )
    with pytest.raises(RuntimeError, match="stale"):
        builder.main()
    assert not (tmp_path / "out").exists()
    assert calls[0][1:] == [
        str(tmp_path / "scripts" / "write_gtf_manifest.py"),
        "--data-dir",
        str(tmp_path / "src" / "viralscan" / "data"),
        "--check",
        str(tmp_path / "manifest.tsv"),
    ]


def test_excluded_ttmdv12_copy_has_no_stale_gene_rows():
    rows = anellovirus.load_gene_table()
    assert not any(row["accession"] == "AB303562.1" for row in rows)
    assert any(row["accession"] == "NC_038359.1" for row in rows)


def test_reference_flags_are_shared_by_tsv_and_html_without_filtering(tmp_path):
    table = VirusIdentityTable(
        (
            GeneIdentity(
                "G1",
                "A1",
                "catalogued",
                True,
                virus_key="k1",
                virus_name="Vector",
                role="contaminant",
            ),
            GeneIdentity(
                "G2",
                "A2",
                "catalogued",
                True,
                virus_key="k2",
                virus_name="Endogenous",
                role="endogenous",
            ),
            GeneIdentity(
                "G3",
                "A3",
                "catalogued",
                True,
                virus_key="k3",
                virus_name="HHV6",
                sibling_group="HHV-6",
            ),
        )
    )
    facts = virus_facts(table)
    assert facts["Vector"].reference_risk_flags == ("vector_reagent_reference",)
    assert facts["Endogenous"].reference_risk_flags == ("endogenous_reference",)
    assert facts["HHV6"].reference_risk_flags == ("iciHHV6_possible",)
    stats = {
        name: {
            "viral_molecules_total_est": 4,
            "infected_cells": 1,
            "total_cells": 2,
            "pct_infected": 50,
            "viral_molecules_per_10k_est": 3,
        }
        for name in facts
    }
    detection.write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path), facts=facts)
    tsv = pd.read_csv(tmp_path / "results" / "viral_summary.tsv", sep="\t")
    assert set(tsv["viral_molecules_total_est"]) == {4}
    assert set(tsv["reference_risk_flags"]) == {
        "vector_reagent_reference",
        "endogenous_reference",
        "iciHHV6_possible",
    }
    from jinja2 import Environment, FileSystemLoader

    env = Environment(loader=FileSystemLoader(str(ROOT / "src" / "viralscan" / "templates")))
    html = env.get_template("report.html.j2").render(
        virus_stats=stats,
        virus_facts=facts,
        params={},
        figures={},
        multimap_primary_call="selected-method",
    )
    for flags in tsv["reference_risk_flags"]:
        assert flags in html


def test_catalogue_with_no_shipped_or_unknown_panel_rows_is_rejected():
    with pytest.raises(ValueError, match="no shipped rows"):
        panel.catalogue_matrix([{"accession": "A1", "panel": "max"}], {"A1.1": {"G"}})
    with pytest.raises(ValueError, match="Unknown catalogue panel"):
        panel.catalogue_matrix([{"accession": "A1", "panel": "Shipped"}], {"A1.1": {"G"}})


def _fixture_data_dir(tmp_path, gtf_text):
    data = tmp_path / "data"
    data.mkdir()
    (data / "a.gtf").write_text(gtf_text)
    (data / "index_exclusions.tsv").write_text(
        "accession\treason\tdecided_by\nX\tfixture\ttester\n"
    )
    catalogue = tmp_path / "catalogue.tsv"
    catalogue.write_text("accession\tspecies\tpanel\nNC_000001\tVirus\tshipped\n")
    return data, catalogue


def test_cli_exit_code_is_nonzero_for_bad_gtf_and_zero_for_clean_panel(tmp_path):
    data, catalogue = _fixture_data_dir(tmp_path, _gtf())
    out = tmp_path / "m.tsv"
    args = ["--data-dir", str(data), "--catalogue", str(catalogue), "--output", str(out)]
    assert panel.main(args) == 0
    (data / "a.gtf").write_text(_gtf(feature="CDS"))  # CDS transcript without an exon
    assert panel.main(args) == 1
    (data / "a.gtf").write_text(_gtf(seq="NC_000009.1"))  # GTF seqname missing from catalogue
    assert panel.main(args) == 1


def test_cli_budget_without_ec_summary_is_a_usage_error_not_a_silent_pass(tmp_path):
    data, catalogue = _fixture_data_dir(tmp_path, _gtf())
    args = ["--data-dir", str(data), "--catalogue", str(catalogue), "--output", str(tmp_path / "m")]
    with pytest.raises(SystemExit):
        panel.main([*args, "--max-ec-size", "1"])


def test_cli_ec_summary_over_budget_fails_the_run(tmp_path):
    data, catalogue = _fixture_data_dir(tmp_path, _gtf())
    summary = tmp_path / "ec.json"
    summary.write_text('{"max_ec_size": 9, "discarded_ec_count": 0}')
    args = ["--data-dir", str(data), "--catalogue", str(catalogue), "--output", str(tmp_path / "m")]
    args += ["--ec-summary", str(summary), "--max-discarded-ecs", "0"]
    assert panel.main([*args, "--max-ec-size", "9"]) == 0
    assert panel.main([*args, "--max-ec-size", "8"]) == 1
