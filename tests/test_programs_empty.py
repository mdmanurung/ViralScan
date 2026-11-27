"""PROGRAMS-EMPTY-01: gene-programme inference finishes when layer 1 detected no virus."""

from __future__ import annotations

import csv
import json
from importlib import resources
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from viralscan.gene_programs import SUMMARY_COLUMNS, load_catalogue, summarise_programs
from viralscan.runconfig import RunConfig
from viralscan.scripts import gene_programs as gp_script

EBV = "Epstein-Barr virus"
SUMMARY_TEXT = {
    "all_candidate_only": f"virus_name\tdetected\tviral_molecules_total_est\n{EBV}\tFalse\t3\n",
    "header_only": "virus_name\tdetected\tviral_molecules_total_est\n",
    "zero_bytes": "",
}


def _registry() -> dict[str, list[str]]:
    text = resources.files("viralscan.schemas.v3").joinpath("table_columns.json").read_text("utf-8")
    return json.loads(text)["optional_tables"]


def _header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def _zero_virus_run(tmp_path: Path, summary_text: str) -> tuple[Path, Path]:
    catalogue = load_catalogue()
    ids = [r["gene_id_bundled"] for r in catalogue if r["virus"] == EBV][:3]
    matrix = sparse.csr_matrix(np.ones((2, len(ids)), dtype=np.float32))
    adata = ad.AnnData(matrix)
    adata.var_names = ids
    adata.obs_names = ["c1", "c2"]
    adata.layers["counts_unique_viral"] = matrix.copy()
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "called_cells.tsv").write_text("barcode\nc1\nc2\n")
    summary = tmp_path / "results" / "viral_summary.tsv"
    summary.write_text(summary_text)
    adata.write_h5ad(tmp_path / "adata.h5ad")
    return summary, tmp_path / "adata.h5ad"


@pytest.mark.parametrize("case", sorted(SUMMARY_TEXT))
def test_zero_detected_viruses_write_empty_tables_with_their_schema(
    tmp_path: Path, case: str, caplog: pytest.LogCaptureFixture
) -> None:
    run = tmp_path / "run"
    run.mkdir()
    summary, h5ad = _zero_virus_run(run, SUMMARY_TEXT[case])
    out = tmp_path / "out"
    with caplog.at_level("INFO"):
        gp_script.main(RunConfig(), str(run), str(h5ad), str(summary), "", results_root=str(out))
    registry = _registry()
    summary_tsv = out / "results" / "gene_program_summary.tsv"
    cells_tsv = out / "results" / "gene_program_cells.tsv"
    assert _header(summary_tsv) == registry["gene_program_summary.tsv"]
    assert _header(cells_tsv) == registry["gene_program_cells.tsv"]
    assert len(summary_tsv.read_text().splitlines()) == 1
    assert len(cells_tsv.read_text().splitlines()) == 1
    assert "no virus detected by layer 1" in caplog.text


def test_summary_columns_match_a_populated_summary() -> None:
    catalogue = load_catalogue()
    populated = summarise_programs(None, catalogue, viruses=[EBV])
    empty = summarise_programs(None, catalogue, viruses=[])
    assert tuple(populated.columns) == SUMMARY_COLUMNS == tuple(empty.columns)
    assert empty.empty
