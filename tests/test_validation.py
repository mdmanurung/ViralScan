from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from viralscan.run_safety import (
    RUN_MANIFEST,
    build_run_manifest,
    prepare_output_directory,
    write_run_complete,
)
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
    # required_uns from h5ad_contract.json. Before SW-02 this fixture omitted
    # quantification_unit and multimap_method and still passed, because nothing
    # read the contract.
    adata.uns["quantification_unit"] = "bustools-resolved-cb-umi-molecule"
    adata.uns["multimap_method"] = "equal"
    # All seven fields count_audit.schema.json requires.
    adata.uns["molecule_audit"] = {
        "input_molecules": 2,
        "resolved_molecules": 2,
        "unique_molecules": 1,
        "ambiguous_molecules": 1,
        "unresolved_molecules": 0,
        "ignored_read_multiplicity": 0,
        "allocated_ambiguous_mass": 1.0,
    }
    adata.write_h5ad(target / "adata_multimap.h5ad")
    (run / "sample" / "config.yaml").write_text("output: x\n")
    write_run_complete(run)
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


# ── SW-02: every shipped schema is enforced at some boundary ─────────────────


def test_h5ad_contract_is_not_mistaken_for_a_json_schema() -> None:
    """h5ad_contract.json declares no JSON Schema keywords.

    Handing it to a validator would accept every document while looking like
    enforcement, so validate_json_schema must refuse it outright.
    """
    issues = validate_json_schema(
        {"anything": True}, packaged_schema_resource("h5ad_contract.json")
    )

    assert [issue.code for issue in issues] == ["not_a_json_schema"]


def test_matrix_checks_are_driven_by_the_packaged_contract() -> None:
    """The contract file and the enforced checks must not be two sources of truth."""
    from viralscan.validation import h5ad_contract

    contract = h5ad_contract()

    assert set(contract["required_layers"]) == {"counts_unique", "counts_ambiguous_allocated"}
    assert set(contract["required_uns"]) == {
        "count_schema_version",
        "quantification_unit",
        "multimap_method",
        "molecule_audit",
    }


def test_validate_run_detects_a_missing_required_uns_key(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    path = next(run.rglob("adata_multimap.h5ad"))
    adata = ad.read_h5ad(path)
    del adata.uns["quantification_unit"]
    adata.write_h5ad(path)

    report = validate_run(run)

    assert report["ok"] is False
    assert "missing_uns" in {issue["code"] for issue in report["issues"]}


def test_validate_run_detects_an_incomplete_count_audit(tmp_path: Path) -> None:
    """count_audit.schema.json had no reader before SW-02."""
    run = _valid_run(tmp_path)
    path = next(run.rglob("adata_multimap.h5ad"))
    adata = ad.read_h5ad(path)
    audit = dict(adata.uns["molecule_audit"])
    del audit["ignored_read_multiplicity"]
    adata.uns["molecule_audit"] = audit
    adata.write_h5ad(path)

    report = validate_run(run)

    assert report["ok"] is False
    assert "schema_validation" in {issue["code"] for issue in report["issues"]}


def test_validate_run_detects_unique_mass_that_contradicts_the_audit(tmp_path: Path) -> None:
    """The contract's non-overlapping-partition invariant, enforced from disk."""
    run = _valid_run(tmp_path)
    path = next(run.rglob("adata_multimap.h5ad"))
    adata = ad.read_h5ad(path)
    audit = dict(adata.uns["molecule_audit"])
    audit["unique_molecules"] = 5
    audit["input_molecules"] = 6
    audit["resolved_molecules"] = 6
    adata.uns["molecule_audit"] = audit
    adata.write_h5ad(path)

    report = validate_run(run)

    assert report["ok"] is False
    assert "audit_unique_mismatch" in {issue["code"] for issue in report["issues"]}


def test_validate_run_detects_a_resolved_molecule_miscount(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    path = next(run.rglob("adata_multimap.h5ad"))
    adata = ad.read_h5ad(path)
    audit = dict(adata.uns["molecule_audit"])
    audit["resolved_molecules"] = 99
    adata.uns["molecule_audit"] = audit
    adata.write_h5ad(path)

    report = validate_run(run)

    assert report["ok"] is False
    assert "audit_resolved_mismatch" in {issue["code"] for issue in report["issues"]}


def test_validate_run_validates_a_reference_manifest_it_finds(tmp_path: Path) -> None:
    """reference_manifest.schema.json had no reader at any boundary before SW-02."""
    run = _valid_run(tmp_path)
    (run / "reference_manifest.json").write_text(
        json.dumps({"schema_version": "3.0.0"}), encoding="utf-8"
    )

    report = validate_run(run)

    assert report["ok"] is False
    assert "schema_validation" in {issue["code"] for issue in report["issues"]}


def test_validate_run_validates_an_evidence_manifest_it_finds(tmp_path: Path) -> None:
    run = _valid_run(tmp_path)
    (run / "evidence_manifest.json").write_text(json.dumps({"target": {}}), encoding="utf-8")

    report = validate_run(run)

    assert report["ok"] is False
    assert "schema_validation" in {issue["code"] for issue in report["issues"]}


def test_an_empty_manifest_no_longer_skips_schema_validation(tmp_path: Path) -> None:
    """The `if manifest:` guard reported a clean run for a manifest with no fields."""
    run = _valid_run(tmp_path)
    (run / RUN_MANIFEST).write_text("{}", encoding="utf-8")

    report = validate_run(run, verify_inputs=False)

    assert report["ok"] is False
    assert "schema_validation" in {issue["code"] for issue in report["issues"]}


def test_writing_an_artifact_that_violates_its_schema_raises(tmp_path: Path) -> None:
    """Write boundary raises; ViralScan authors these, so a violation is a defect."""
    from viralscan.validation import SchemaContractError, require_schema_valid

    with pytest.raises(SchemaContractError, match="refusing to write"):
        require_schema_valid({"schema_version": "3.0.0"}, "reference_manifest.schema.json")


def test_writing_a_conformant_artifact_is_allowed() -> None:
    from viralscan.validation import require_schema_valid

    require_schema_valid(
        {
            "input_molecules": 2,
            "resolved_molecules": 2,
            "unique_molecules": 1,
            "ambiguous_molecules": 1,
            "unresolved_molecules": 0,
            "ignored_read_multiplicity": 0,
            "allocated_ambiguous_mass": 1.0,
        },
        "count_audit.schema.json",
    )


def test_numpy_scalars_do_not_trip_the_write_boundary() -> None:
    """Counters arrive as np.int64 in production but as int in every fixture."""
    from viralscan.validation import require_schema_valid

    require_schema_valid(
        {
            "input_molecules": np.int64(2),
            "resolved_molecules": np.int64(2),
            "unique_molecules": np.int64(1),
            "ambiguous_molecules": np.int64(1),
            "unresolved_molecules": np.int64(0),
            "ignored_read_multiplicity": np.int64(0),
            "allocated_ambiguous_mass": np.float64(1.0),
        },
        "count_audit.schema.json",
    )
