"""TABLE-SCHEMA-02: tables written by optional arms and post-hoc verbs have pinned columns.

``schemas/v3/table_columns.json`` lists them under ``optional_tables``. Each header is read from
the file the real writer produces on a small synthetic fixture that an existing test already
uses; this test pins what exists and renames nothing. It also checks every column against
``docs/output_reference.md``.
"""

from __future__ import annotations

import csv
import fnmatch
import json
import shutil
from importlib import resources
from pathlib import Path

import numpy as np
import pytest

from tests.test_discovery_kb_evidence import MOLECULES
from tests.test_discovery_kb_evidence import _setup as _kb_setup
from tests.test_discovery_read_support import _run as _read_support_run
from tests.test_discovery_star_evidence import RECORDS
from tests.test_discovery_star_evidence import _run as _star_run
from tests.test_strand_qc import TestWrite as _StrandWrite
from tests.test_table_columns import _documented
from viralscan import discovery, read_support
from viralscan.gene_programs import load_catalogue
from viralscan.runconfig import RunConfig
from viralscan.scripts import auto_evidence
from viralscan.scripts import gene_programs as gp_script
from viralscan.strand_qc import FIELDS as STRAND_FIELDS
from viralscan.strand_qc import write_strand_qc
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

EBV = "Epstein-Barr virus"
HAS_SAMTOOLS = shutil.which("samtools") is not None


def _registry() -> dict[str, list[str]]:
    text = resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    return json.loads(text)["optional_tables"]


def _header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def _gene_programs(tmp_path: Path) -> dict[str, Path]:
    """``main()`` of the programme script on one cell that carries an EBV immediate-early marker."""
    import anndata as ad
    from scipy import sparse

    catalogue = load_catalogue()
    ids = [r["gene_id_bundled"] for r in catalogue if r["virus"] == EBV]
    hot = next(
        r["gene_id_bundled"]
        for r in catalogue
        if r["virus"] == EBV
        and r["programme"] == "productive"
        and r["kinetic_class"] == "immediate_early"
    )
    counts = np.zeros((1, len(ids)), dtype=np.float32)
    counts[0, ids.index(hot)] = 3.0
    matrix = sparse.csr_matrix(counts)
    adata = ad.AnnData(matrix)
    adata.var_names = ids
    adata.obs_names = ["c1"]
    adata.layers["counts_unique_viral"] = matrix.copy()
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "called_cells.tsv").write_text("barcode\nc1\n")
    summary = tmp_path / "results" / "viral_summary.tsv"
    summary.write_text(f"virus_name\tviral_molecules_total_est\n{EBV}\t5\n")
    adata.write_h5ad(tmp_path / "adata.h5ad")
    gp_script.main(RunConfig(), str(tmp_path), str(tmp_path / "adata.h5ad"), str(summary), "")
    return {
        "gene_program_summary.tsv": tmp_path / "results" / "gene_program_summary.tsv",
        "gene_program_cells.tsv": tmp_path / "results" / "gene_program_cells.tsv",
    }


def _discovery_evidence(tmp_path: Path) -> dict[str, Path]:
    """The three route tables, their concatenation and the candidate-review table."""
    kb_dir = tmp_path / "kb"
    kb_dir.mkdir()
    disc, sample = _kb_setup(kb_dir, MOLECULES)
    kb = discovery.build_kb_evidence(disc, sample)
    star_dir = tmp_path / "star"
    star_dir.mkdir()
    disc_s, sample_s, reader = _star_run(star_dir, RECORDS)
    star = discovery.build_star_evidence(disc_s, sample_s, record_reader=reader)
    rs_dir = tmp_path / "rs"
    rs_dir.mkdir()
    _, read_support_tsv, _, _ = _read_support_run(rs_dir)
    combined = discovery._concat_evidence([kb, star, read_support_tsv], tmp_path / "evidence.tsv")
    assert combined is not None

    run_dir = tmp_path / "review_run"
    (run_dir / "results").mkdir(parents=True)
    VirusIdentityTable(
        (
            GeneIdentity(
                gene_id="g1",
                genome_accession="A1",
                status="x",
                viral=True,
                taxid="1",
                virus_name="Torque teno virus",
                family="Anelloviridae",
            ),
        )
    ).write_tsv(run_dir / "results" / "virus_identity.tsv")
    (run_dir / "results" / "viral_summary.tsv").write_text(
        "virus_name\tdetected\tdetection_source\nTorque teno virus\tTrue\tunsupported_for_test\n"
    )
    auto_evidence.review_candidates(run_dir, auto_evidence.RunConfig(), tmp_path / "review")
    return {
        "evidence_kb.tsv": kb,
        "evidence_star.tsv": star,
        "evidence_read_support.tsv": read_support_tsv,
        "evidence.tsv": combined,
        "evidence_review.tsv": tmp_path / "review" / "evidence_review.tsv",
    }


@pytest.fixture(scope="module")
def produced(tmp_path_factory) -> dict[str, Path]:
    tables: dict[str, Path] = {}
    tables |= _gene_programs(tmp_path_factory.mktemp("programmes"))
    tables |= _discovery_evidence(tmp_path_factory.mktemp("evidence"))
    return tables


@pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not on PATH")
def test_strand_and_read_support_headers(tmp_path: Path) -> None:
    out, bam = _StrandWrite()._evidence_dir(tmp_path)
    write_strand_qc(
        bam,
        out,
        technology="10xv3",
        host_gtf=tmp_path / "host.gtf",
        viral_gtfs=[tmp_path / "viral.gtf"],
    )
    assert _header(out / "strand_qc.tsv") == _registry()["strand_qc.tsv"]
    read_support.write_molecule_support(bam, out)
    assert _header(out / "molecule_support.tsv") == _registry()["molecule_support.tsv"]
    assert _header(out / "read_support_summary.tsv") == _registry()["read_support_summary.tsv"]
    # The cross-evidence-directory roll-up uses the summary's writer fields, so the two agree.
    (out / "read_lineage.tsv.gz").unlink()  # the fixture's lineage has no assigned_weight column
    read_support.main([str(out), "-o", str(tmp_path / "calls")])
    assert (
        _header(tmp_path / "calls" / "read_support_calls.tsv")
        == _registry()["read_support_calls.tsv"]
    )
    assert _registry()["read_support_calls.tsv"] == _registry()["read_support_summary.tsv"]


@pytest.mark.parametrize(
    "table",
    [
        "gene_program_summary.tsv",
        "gene_program_cells.tsv",
        "evidence_kb.tsv",
        "evidence_star.tsv",
        "evidence_read_support.tsv",
        "evidence.tsv",
        "evidence_review.tsv",
    ],
)
def test_header_matches_registry(table: str, produced: dict[str, Path]) -> None:
    assert _header(produced[table]) == _registry()[table]


def test_registry_matches_the_writers_constants() -> None:
    tables = _registry()
    assert tables["strand_qc.tsv"] == list(STRAND_FIELDS)
    assert tables["molecule_support.tsv"] == list(read_support.ROW_FIELDS)
    assert tables["read_support_summary.tsv"] == list(read_support.SUMMARY_FIELDS)
    for name in (
        "evidence_kb.tsv",
        "evidence_star.tsv",
        "evidence_read_support.tsv",
        "evidence.tsv",
    ):
        assert tables[name] == list(discovery.EVIDENCE_COLUMNS)
    assert tables["evidence_review.tsv"] == list(auto_evidence.REVIEW_FIELDS)


def test_optional_tables_are_disjoint_from_the_default_tables() -> None:
    default = json.loads(
        resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    )["tables"]
    assert not set(_registry()) & set(default)
    for name, columns in _registry().items():
        assert name.endswith(".tsv")
        assert columns and len(columns) == len(set(columns)), f"{name}: empty or duplicate columns"


@pytest.mark.parametrize("table", sorted(_registry()))
def test_every_column_is_documented(table: str) -> None:
    exact, patterns = _documented()
    missing = [
        column
        for column in _registry()[table]
        if column not in exact and not any(fnmatch.fnmatch(column, p) for p in patterns)
    ]
    assert not missing, f"{table}: columns not in docs/output_reference.md: {missing}"
