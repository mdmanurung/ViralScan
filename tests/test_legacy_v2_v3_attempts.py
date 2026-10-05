from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.sparse as sp

from scripts.compare_legacy_v2_v3_attempts import compare_attempts, main


def _write_attempt(root: Path, *, attempt_id: str) -> None:
    v3 = root / "v3"
    v3.mkdir(parents=True)
    (root / "status.json").write_text(
        json.dumps(
            {
                "status": "success",
                "attempt_id": attempt_id,
                "total_wall_seconds": 12.3,
            }
        ),
        encoding="utf-8",
    )
    (v3 / "status.json").write_text(
        json.dumps(
            {
                "status": "success",
                "attempt_id": attempt_id,
                "scientific_parameter_hash": "p" * 64,
                "scratch_path": f"/tmp/{attempt_id}",
                "total_wall_seconds": 10.0,
            }
        ),
        encoding="utf-8",
    )
    (v3 / "result.json").write_text(
        json.dumps(
            {
                "scientific_parameter_hash": "p" * 64,
                "attempt_id": attempt_id,
                "setup_seconds": 1.0,
                "allocation_seconds": 2.0,
            }
        ),
        encoding="utf-8",
    )
    fingerprints = {
        "bus": {"size_bytes": 3, "sha256": "a" * 64},
        "t2g": {"size_bytes": 4, "sha256": "b" * 64},
    }
    (v3 / "hashes.json").write_text(
        json.dumps(
            {
                "scientific_parameter_hash": "p" * 64,
                "inputs": fingerprints,
                "inputs_before": fingerprints,
                "inputs_after": fingerprints,
                "code": {
                    "benchmark": {"size_bytes": 5, "sha256": "c" * 64},
                    "multimapping": {"size_bytes": 6, "sha256": "d" * 64},
                },
                "virus_grouping_sha256": "e" * 64,
            }
        ),
        encoding="utf-8",
    )
    (v3 / "input_fingerprints_before.json").write_text(
        json.dumps(
            {
                "attempt_id": attempt_id,
                "scientific_parameter_hash": "p" * 64,
                "inputs": fingerprints,
            }
        ),
        encoding="utf-8",
    )
    (v3 / "input_fingerprints_after.json").write_text(
        json.dumps(
            {
                "attempt_id": attempt_id,
                "scientific_parameter_hash": "p" * 64,
                "inputs": fingerprints,
            }
        ),
        encoding="utf-8",
    )
    (v3 / "count_audit.json").write_text(
        json.dumps(
            {
                "input_molecules": 3,
                "unique_molecules": 2,
                "ambiguous_molecules": 1,
                "unresolved_molecules": 0,
                "allocated_ambiguous_mass": 1.0,
            }
        ),
        encoding="utf-8",
    )
    (v3 / "per_virus.tsv").write_text(
        "virus_name\tn_features\tv3_unique\tv3_equal\tv3_host_conservative\n"
        "Epstein-Barr virus\t1\t2\t3\t2.5\n",
        encoding="utf-8",
    )
    (v3 / "per_cell.tsv").write_text(
        "barcode\tv3_unique_viral\tv3_equal_viral\tv3_host_conservative_viral\nBC1\t2\t3\t2.5\n",
        encoding="utf-8",
    )
    adata = ad.AnnData(
        X=sp.csr_matrix([[2.5, 0.0]]),
        obs=pd.DataFrame({"batch": ["A"]}, index=["BC1"]),
        var=pd.DataFrame(
            {"gene_name": ["EBNA", "HOST"], "is_viral": [True, False]},
            index=["EPSTEIN_A", "HOST_A"],
        ),
    )
    adata.layers["counts_unique"] = sp.csr_matrix([[2.0, 0.0]])
    adata.layers["counts_ambiguous_allocated"] = sp.csr_matrix([[0.5, 0.0]])
    adata.layers["counts_multimap_equal"] = sp.csr_matrix([[1.0, 0.0]])
    adata.uns["attempt_id"] = attempt_id
    adata.uns["timing_seconds"] = 10.0 if attempt_id == "attempt-one" else 20.0
    adata.write_h5ad(v3 / "adata_v3.h5ad")


def test_identical_science_matches_despite_attempt_metadata_and_h5ad_bytes(
    tmp_path: Path,
) -> None:
    first = tmp_path / "pilot-attempt-1"
    second = tmp_path / "pilot-attempt-2"
    _write_attempt(first, attempt_id="attempt-one")
    _write_attempt(second, attempt_id="attempt-two")
    report_path = tmp_path / "comparison.json"

    report = compare_attempts(first, second, report_path)

    assert report["status"] == "match"
    assert report["mismatches"] == []
    assert report_path.is_file()
    assert all(check["status"] == "match" for check in report["checks"])


def test_scientific_table_or_matrix_change_fails_closed(tmp_path: Path) -> None:
    first = tmp_path / "pilot-attempt-1"
    second = tmp_path / "pilot-attempt-2"
    _write_attempt(first, attempt_id="attempt-one")
    _write_attempt(second, attempt_id="attempt-two")
    (second / "v3" / "per_cell.tsv").write_text(
        "barcode\tv3_unique_viral\tv3_equal_viral\tv3_host_conservative_viral\nBC1\t2\t99\t2.5\n",
        encoding="utf-8",
    )
    changed = ad.read_h5ad(second / "v3" / "adata_v3.h5ad")
    changed.layers["counts_unique"] = sp.csr_matrix([[7.0, 0.0]])
    changed.write_h5ad(second / "v3" / "adata_v3.h5ad")
    report_path = tmp_path / "comparison.json"

    report = compare_attempts(first, second, report_path)

    assert report["status"] == "mismatch"
    mismatch_checks = {row["check"] for row in report["mismatches"]}
    assert "per_cell" in mismatch_checks
    assert "h5ad.layers.counts_unique" in mismatch_checks
    assert main([str(first), str(second), "--report", str(report_path)]) == 1


def test_missing_or_failed_attempt_is_a_reported_mismatch(tmp_path: Path) -> None:
    first = tmp_path / "pilot-attempt-1"
    second = tmp_path / "pilot-attempt-2"
    _write_attempt(first, attempt_id="attempt-one")
    second.mkdir()
    (second / "status.json").write_text('{"status":"failed"}\n', encoding="utf-8")

    report = compare_attempts(first, second, tmp_path / "comparison.json")

    assert report["status"] == "mismatch"
    assert any(row["check"] == "attempt_contract" for row in report["mismatches"])
