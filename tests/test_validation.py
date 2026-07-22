from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from viralscan.run_safety import build_run_manifest, prepare_output_directory
from viralscan.validation import (
    REQUIRED_V3_SCHEMAS,
    doctor_report,
    packaged_schema_resource,
    validate_json_schema,
    validate_run,
)


def _valid_run(tmp_path: Path) -> Path:
    r1, r2 = tmp_path / "R1.fastq", tmp_path / "R2.fastq"
    r1.write_text("r1")
    r2.write_text("r2")
    args = argparse.Namespace(
        sample1=str(r1),
        sample2=str(r2),
        output=str(tmp_path / "run"),
        multimap_method="equal",
        cell_calling="emptydrops",
        called_cells_file=None,
        resume=False,
        overwrite=False,
        yes=False,
        verbose=False,
        quiet=False,
    )
    run = Path(args.output)
    prepare_output_directory(
        run, build_run_manifest(args), resume=False, overwrite=False, yes=False
    )
    target = run / "sample" / "kb-python" / "counts_unfiltered"
    target.mkdir(parents=True)
    unique = sparse.csr_matrix([[1.0, 0.0]])
    ambiguous = sparse.csr_matrix([[0.5, 0.5]])
    adata = ad.AnnData(
        X=unique + ambiguous,
        obs=pd.DataFrame(index=["BC1"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    adata.layers["counts_unique"] = unique
    adata.layers["counts_ambiguous_allocated"] = ambiguous
    adata.uns["count_schema_version"] = "3.0.0"
    adata.uns["molecule_audit"] = {
        "input_molecules": 2,
        "unique_molecules": 1,
        "ambiguous_molecules": 1,
        "unresolved_molecules": 0,
        "allocated_ambiguous_mass": 1.0,
    }
    adata.write_h5ad(target / "adata_multimap.h5ad")
    return run


def test_validate_run_accepts_conserved_v3_output(tmp_path: Path) -> None:
    report = validate_run(_valid_run(tmp_path))
    assert report["ok"] is True
    assert report["issues"] == []


def test_validate_run_detects_layer_drift(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    path = next(run.rglob("adata_multimap.h5ad"))
    adata = ad.read_h5ad(path)
    adata.X = sparse.csr_matrix(np.zeros((1, 2)))
    adata.write_h5ad(path)
    report = validate_run(run)
    assert report["ok"] is False
    assert "x_layer_mismatch" in {issue["code"] for issue in report["issues"]}


def test_doctor_pip_profile_has_no_external_tool_gate() -> None:
    report = doctor_report("pip")
    assert report["tools"] == {}
    assert set(report["python"]) >= {"anndata", "snakemake"}
    assert set(report["schemas"]) == set(REQUIRED_V3_SCHEMAS)
    assert all(report["schemas"].values())
    assert report["schema_errors"] == {}


def test_all_required_v3_schemas_are_packaged() -> None:
    for name in REQUIRED_V3_SCHEMAS:
        schema = json.loads(packaged_schema_resource(name).read_text(encoding="utf-8"))
        assert isinstance(schema, dict), name


def test_packaged_v3_schemas_match_public_contracts() -> None:
    public_schema_dir = Path(__file__).resolve().parents[1] / "schemas" / "v3"
    for name in REQUIRED_V3_SCHEMAS:
        assert (
            packaged_schema_resource(name).read_bytes() == (public_schema_dir / name).read_bytes()
        )


def test_doctor_fails_closed_when_packaged_schemas_are_missing() -> None:
    with patch(
        "viralscan.validation.packaged_schema_resource",
        side_effect=FileNotFoundError("missing packaged schema"),
    ):
        report = doctor_report("pip")
    assert report["ok"] is False
    assert not any(report["schemas"].values())
    assert set(report["schema_errors"]) == set(REQUIRED_V3_SCHEMAS)


def test_validate_run_fails_closed_when_packaged_schema_is_missing(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    with patch(
        "viralscan.validation.packaged_schema_resource",
        side_effect=FileNotFoundError("missing packaged schema"),
    ):
        report = validate_run(run)
    assert report["ok"] is False
    assert "schema_unavailable" in {issue["code"] for issue in report["issues"]}


def test_validate_json_schema_reports_contract_errors(tmp_path: Path) -> None:
    schema = tmp_path / "schema.json"
    schema.write_text('{"type":"object","required":["value"]}')
    issues = validate_json_schema({}, schema)
    assert issues and issues[0].code == "schema_validation"
