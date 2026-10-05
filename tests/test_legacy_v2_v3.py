from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scripts.compare_legacy_v2_v3 import (
    InventoryContract,
    InventoryContractError,
    inventory,
    reconstruct_legacy,
    run_row,
    summarize_results,
    validate_results,
)


def _write_tree(
    root: Path,
    relative_run: str,
    *,
    sample1: str,
    sample2: str,
    technology: str,
) -> Path:
    run_dir = root / relative_run
    kb_dir = run_dir / "kb-python"
    counts_dir = kb_dir / "counts_unfiltered"
    reference_dir = root.parent / "reference"
    counts_dir.mkdir(parents=True)
    reference_dir.mkdir(exist_ok=True)

    index = reference_dir / "index.idx"
    t2g = reference_dir / "t2g.txt"
    transcriptome = reference_dir / "transcriptome.fa"
    index.write_bytes(b"index")
    t2g.write_text("tx\tgene\n", encoding="utf-8")
    transcriptome.write_text(">tx\nACGT\n", encoding="utf-8")

    (run_dir / "config.yaml").write_text(
        "\n".join(
            (
                f"index: {index}",
                f"transcripts: {t2g}",
                f"sample1: {sample1}",
                f"sample2: {sample2}",
                f"technology: {technology}",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    (run_dir / "summary.txt").write_text("summary\n", encoding="utf-8")
    (kb_dir / "output.bus").write_bytes(b"bus")
    (kb_dir / "matrix.ec").write_text("0\t0\n", encoding="utf-8")
    (kb_dir / "transcripts.txt").write_text("tx\n", encoding="utf-8")
    (kb_dir / f"10x_version{technology[-1]}_whitelist.txt").write_text("ACGT\n", encoding="utf-8")
    (kb_dir / "kb_info.json").write_text(
        '{"version":"0.29.5","kallisto":{"version":"0.51.1"},"bustools":{"version":"0.45.1"}}\n',
        encoding="utf-8",
    )
    (counts_dir / "adata_multimap.h5ad").write_bytes(b"h5ad")
    (counts_dir / "adata.h5ad").write_bytes(b"adata")
    (counts_dir / "cells_x_genes.barcodes.txt").write_text("ACGT-1\n", encoding="utf-8")
    (counts_dir / "cells_x_genes.genes.txt").write_text("gene\n", encoding="utf-8")
    (counts_dir / "cells_x_genes.genes.names.txt").write_text("Gene\n", encoding="utf-8")
    (kb_dir / "run_info.json").write_text('{"n_pseudoaligned":1}\n', encoding="utf-8")
    (run_dir / "log").mkdir()
    (run_dir / "log" / "analysis.txt").write_text("gene\n", encoding="utf-8")
    return run_dir


def test_inventory_emits_raw_and_sanitized_manifests(tmp_path: Path) -> None:
    source = tmp_path / "archive"
    _write_tree(
        source,
        "SRR12682296/SRR12682296",
        sample1="/private/input/SRR12682296_1.fastq",
        sample2="/private/input/SRR12682296_2.fastq",
        technology="10xv2",
    )
    raw_manifest = tmp_path / "benchmark_runs" / "raw.tsv"
    cohort_manifest = tmp_path / "analysis" / "cohort.tsv"

    rows = inventory(
        source,
        raw_manifest,
        cohort_manifest,
        contract=InventoryContract(
            technical_rows=1,
            logical_inputs=1,
            chemistry_counts={"10xv2": 1, "10xv3": 0},
        ),
    )

    assert len(rows) == 1
    with cohort_manifest.open(newline="", encoding="utf-8") as handle:
        sanitized = list(csv.DictReader(handle, delimiter="\t"))
    with raw_manifest.open(newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle, delimiter="\t"))

    assert sanitized[0]["run_id"] == "SRR12682296__SRR12682296"
    assert sanitized[0]["logical_id"] == "SRR12682296"
    assert sanitized[0]["chemistry"] == "10xv2"
    assert sanitized[0]["sample_class"] == "ebv_control"
    assert sanitized[0]["expected_target"] == "EBV"
    assert sanitized[0]["output_bus_path"] == ("SRR12682296/SRR12682296/kb-python/output.bus")
    assert sanitized[0]["adata_path"].endswith("kb-python/counts_unfiltered/adata.h5ad")
    assert sanitized[0]["analysis_path"].endswith("log/analysis.txt")
    assert sanitized[0]["run_info_sha256"]
    assert sanitized[0]["output_bus_sha256"]
    assert str(tmp_path) not in cohort_manifest.read_text(encoding="utf-8")
    assert raw[0]["output_bus_path"] == str(source / "SRR12682296/SRR12682296/kb-python/output.bus")


def test_inventory_enforces_frozen_44_tree_42_input_contract(tmp_path: Path) -> None:
    source = tmp_path / "archive"
    for srr in ("SRR12682296", "SRR12682297", "SRR12682298", "SRR6825024"):
        _write_tree(
            source,
            f"{srr}/{srr}",
            sample1=f"/inputs/{srr}_1.fastq",
            sample2=f"/inputs/{srr}_2.fastq",
            technology="10xv2",
        )
    for relative_run in ("SRR6825025/SRR6825025", "old/SRR6825025"):
        _write_tree(
            source,
            relative_run,
            sample1="/inputs/SRR6825025_1.fastq",
            sample2="/inputs/SRR6825025_2.fastq",
            technology="10xv2",
        )

    repeated_skin = "WS_SKN_KCL10525740_S3_L002"
    _write_tree(
        source,
        "WS_SKN__KCL10525740/S3/L001/WS",
        sample1=f"{repeated_skin}_R1_001.fastq.gz",
        sample2=f"{repeated_skin}_R2_001.fastq.gz",
        technology="10xv3",
    )
    _write_tree(
        source,
        "WS_SKN_KCL10525740/S3/L002/WS",
        sample1=f"{repeated_skin}_R1_001.fastq.gz",
        sample2=f"{repeated_skin}_R2_001.fastq.gz",
        technology="10xv3",
    )
    for number in range(36):
        input_id = f"WS_SKN_KCL_TEST_S{number:02d}_L001"
        _write_tree(
            source,
            f"skin/sample-{number:02d}/WS",
            sample1=f"{input_id}_R1_001.fastq.gz",
            sample2=f"{input_id}_R2_001.fastq.gz",
            technology="10xv3",
        )

    rows = inventory(
        source,
        tmp_path / "benchmark_runs/raw.tsv",
        tmp_path / "analysis/cohort.tsv",
    )

    assert len(rows) == 44
    assert len({row["logical_id"] for row in rows}) == 42
    assert Counter(row["chemistry"] for row in rows) == {
        "10xv2": 6,
        "10xv3": 38,
    }
    repeat_groups = Counter(
        row["technical_repeat_group"] for row in rows if row["technical_repeat_group"]
    )
    assert repeat_groups == {
        "repeat__SRR6825025": 2,
        "repeat__WS_SKN_KCL10525740_S3_L002": 2,
    }


def test_inventory_fails_before_writing_when_cohort_counts_drift(
    tmp_path: Path,
) -> None:
    source = tmp_path / "archive"
    _write_tree(
        source,
        "SRR12682296/SRR12682296",
        sample1="SRR12682296_1.fastq",
        sample2="SRR12682296_2.fastq",
        technology="10xv2",
    )
    raw_manifest = tmp_path / "benchmark_runs/raw.tsv"
    cohort_manifest = tmp_path / "analysis/cohort.tsv"

    with pytest.raises(
        InventoryContractError,
        match=r"technical rows: expected 44, found 1.*logical inputs: expected 42, found 1",
    ):
        inventory(source, raw_manifest, cohort_manifest)

    assert not raw_manifest.exists()
    assert not cohort_manifest.exists()


def test_legacy_reconstruction_adds_original_and_corrected_before_strict_call(
    tmp_path: Path,
) -> None:
    h5ad_path = tmp_path / "adata_multimap.h5ad"
    summary_path = tmp_path / "summary.txt"
    accessions_path = tmp_path / "analysis.txt"
    adata = ad.AnnData(
        X=sp.csr_matrix(np.zeros((2, 2))),
        obs=pd.DataFrame(index=["BC1", "BC2"]),
        var=pd.DataFrame(index=["EPSTEIN_HHV4_EBNA-2", "EPSTEIN_HHV4_LOW"]),
    )
    adata.layers["counts_original"] = sp.csr_matrix([[1.0, 0.5], [0.0, 0.0]])
    adata.layers["counts_corrected"] = sp.csr_matrix([[1.0, 0.5], [0.0, 0.0]])
    adata.write_h5ad(h5ad_path)
    accessions_path.write_text("EPSTEIN_HHV4_EBNA-2\nEPSTEIN_HHV4_LOW\n", encoding="utf-8")
    summary_path.write_text(
        "Found viral Gene IDs including the count:\n"
        "Gene ID; Gene Count\n"
        "EPSTEIN_HHV4_EBNA-2;2.0\n\n"
        "Epstein-Barr virus has a viral load of: 2.0 UMIs.\n\n"
        "Total amount of viral load found: 2.0\n",
        encoding="utf-8",
    )

    result = reconstruct_legacy(
        h5ad_path,
        summary_path,
        accessions_path,
        virus_name_map={"EPSTEIN": "Epstein-Barr virus"},
    )

    assert result.called_gene_totals == {"EPSTEIN_HHV4_EBNA-2": 2.0}
    assert result.reconstructed_virus_totals == {"Epstein-Barr virus": 2.0}
    assert result.summary.gene_totals == {"EPSTEIN_HHV4_EBNA-2": 2.0}
    assert result.summary.virus_totals == {"Epstein-Barr virus": 2.0}
    assert all(row["status"] == "match" for row in result.comparisons)


def test_legacy_reconstruction_retains_matches_mismatches_and_missing_values(
    tmp_path: Path,
) -> None:
    h5ad_path = tmp_path / "adata_multimap.h5ad"
    summary_path = tmp_path / "summary.txt"
    accessions_path = tmp_path / "analysis.txt"
    genes = ["EPSTEIN_A", "EPSTEIN_B"]
    adata = ad.AnnData(
        X=sp.csr_matrix((1, 2)),
        obs=pd.DataFrame(index=["BC1"]),
        var=pd.DataFrame(index=genes),
    )
    adata.layers["counts_original"] = sp.csr_matrix([[1.0, 1.0]])
    adata.layers["counts_corrected"] = sp.csr_matrix([[1.0, 2.0]])
    adata.write_h5ad(h5ad_path)
    accessions_path.write_text("\n".join(genes) + "\n", encoding="utf-8")
    summary_path.write_text(
        "Found viral Gene IDs including the count:\n"
        "Gene ID; Gene Count\n"
        "EPSTEIN_A;2.0000005\n"
        "EPSTEIN_MISSING_FROM_H5AD;7.0\n\n"
        "Epstein-Barr virus has a viral load of: 10.0 UMIs.\n\n"
        "Total amount of viral load found: 5.0000005\n",
        encoding="utf-8",
    )

    result = reconstruct_legacy(
        h5ad_path,
        summary_path,
        accessions_path,
        virus_name_map={"EPSTEIN": "Epstein-Barr virus"},
    )
    audit = {(row["record_type"], row["identifier"]): row for row in result.comparisons}

    assert audit[("gene", "EPSTEIN_A")]["status"] == "match"
    assert audit[("gene", "EPSTEIN_B")]["status"] == "missing_in_summary"
    assert audit[("gene", "EPSTEIN_MISSING_FROM_H5AD")]["status"] == "missing_in_reconstruction"
    assert audit[("virus", "Epstein-Barr virus")]["status"] == "mismatch"
    assert audit[("total", "total_viral_load")]["status"] == "match"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_run_row_fixture(tmp_path: Path) -> tuple[Path, Path, str]:
    sample_dir = tmp_path / "source" / "sample"
    counts_dir = sample_dir / "kb-python" / "counts_unfiltered"
    log_dir = sample_dir / "log"
    reference_dir = tmp_path / "reference"
    counts_dir.mkdir(parents=True)
    log_dir.mkdir()
    reference_dir.mkdir()
    paths = {
        "config": sample_dir / "config.yaml",
        "summary": sample_dir / "summary.txt",
        "output_bus": sample_dir / "kb-python" / "output.bus",
        "matrix_ec": sample_dir / "kb-python" / "matrix.ec",
        "transcripts": sample_dir / "kb-python" / "transcripts.txt",
        "whitelist": sample_dir / "kb-python" / "10x_version3_whitelist.txt",
        "adata_multimap": counts_dir / "adata_multimap.h5ad",
        "barcodes": counts_dir / "cells_x_genes.barcodes.txt",
        "genes": counts_dir / "cells_x_genes.genes.txt",
        "gene_names": counts_dir / "cells_x_genes.genes.names.txt",
        "adata": counts_dir / "adata.h5ad",
        "run_info": sample_dir / "kb-python" / "run_info.json",
        "analysis": log_dir / "analysis.txt",
        "index": reference_dir / "index.idx",
        "t2g": reference_dir / "t2g.txt",
        "transcriptome": reference_dir / "transcriptome.fa",
    }
    for name, path in paths.items():
        if name not in {"adata_multimap", "adata"}:
            path.write_text(name + "\n", encoding="utf-8")
    adata = ad.AnnData(
        X=sp.csr_matrix([[0.0]]),
        obs=pd.DataFrame(index=["BC1"]),
        var=pd.DataFrame(index=["EPSTEIN_GENE"]),
    )
    adata.layers["counts_original"] = sp.csr_matrix([[1.0]])
    adata.layers["counts_corrected"] = sp.csr_matrix([[1.0]])
    adata.write_h5ad(paths["adata_multimap"])
    ad.AnnData(
        X=sp.csr_matrix([[1.0]]),
        obs=pd.DataFrame(index=["BC1-1"]),
        var=pd.DataFrame(index=["EPSTEIN_GENE"]),
    ).write_h5ad(paths["adata"])
    paths["config"].write_text(
        f"index: {paths['index']}\ntranscripts: {paths['t2g']}\ntechnology: 10xv3\n",
        encoding="utf-8",
    )
    paths["summary"].write_text(
        "Found viral Gene IDs including the count:\n"
        "Gene ID; Gene Count\n"
        "EPSTEIN_GENE;2.0\n\n"
        "Epstein-Barr virus has a viral load of: 2.0 UMIs.\n\n"
        "Total amount of viral load found: 2.0\n",
        encoding="utf-8",
    )
    paths["analysis"].write_text("EPSTEIN_GENE\n", encoding="utf-8")

    run_id = "skin__sample"
    row = {
        "run_id": run_id,
        "logical_id": "skin_sample",
        "chemistry": "10xv3",
    }
    for name, path in paths.items():
        row[f"{name}_path"] = str(path)
        row[f"{name}_sha256"] = _sha256(path)
        row[f"{name}_bytes"] = str(path.stat().st_size)
    manifest = tmp_path / "raw_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    return manifest, sample_dir, run_id


def test_run_row_verifies_manifest_writes_legacy_and_invokes_v3(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, sample_dir, run_id = _write_run_row_fixture(tmp_path)
    output_root = tmp_path / "outputs"
    captured: dict[str, object] = {}

    from scripts import benchmark_v3_multimap

    def fake_v3_main(argv):
        captured["argv"] = argv
        output = Path(argv[argv.index("--output") + 1])
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text('{"status":"success"}\n', encoding="utf-8")
        (output.parent / "status.json").write_text('{"status":"success"}\n', encoding="utf-8")
        return 0

    monkeypatch.setattr(benchmark_v3_multimap, "main", fake_v3_main)
    source_hashes = {path: _sha256(path) for path in sample_dir.rglob("*") if path.is_file()}

    status = run_row(
        manifest,
        run_id,
        output_root,
        source_root=tmp_path / "source",
        threads=3,
    )

    assert status["status"] == "success"
    row_output = output_root / run_id
    legacy_rows = list(
        csv.DictReader(
            (row_output / "legacy" / "legacy_audit.tsv").open(newline="", encoding="utf-8"),
            delimiter="\t",
        )
    )
    assert {row["status"] for row in legacy_rows} == {"match"}
    argv = captured["argv"]
    assert argv[argv.index("--sample-id") + 1] == run_id
    assert Path(argv[argv.index("--sample-dir") + 1]) == sample_dir
    assert Path(argv[argv.index("--whitelist") + 1]).name == ("10x_version3_whitelist.txt")
    assert Path(argv[argv.index("--t2g") + 1]).name == "t2g.txt"
    assert argv[argv.index("--threads") + 1] == "3"
    assert json.loads((row_output / "hash_check_before.json").read_text())["status"] == "match"
    assert json.loads((row_output / "hash_check_after.json").read_text())["status"] == "match"
    assert {
        path: _sha256(path) for path in sample_dir.rglob("*") if path.is_file()
    } == source_hashes


def test_run_row_retains_v3_failure_and_detects_post_run_source_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, sample_dir, run_id = _write_run_row_fixture(tmp_path)
    output_root = tmp_path / "outputs"

    from scripts import benchmark_v3_multimap

    def mutating_failure(_argv):
        (sample_dir / "kb-python" / "output.bus").write_text("changed\n", encoding="utf-8")
        raise RuntimeError("allocator failed")

    monkeypatch.setattr(benchmark_v3_multimap, "main", mutating_failure)

    status = run_row(
        manifest,
        run_id,
        output_root,
        source_root=tmp_path / "source",
    )

    assert status["status"] == "failed"
    assert [error["stage"] for error in status["errors"]] == [
        "v3_benchmark",
        "verify_after",
    ]
    row_output = output_root / run_id
    assert (row_output / "legacy" / "legacy_audit.tsv").is_file()
    assert (row_output / "failure.json").is_file()
    after = json.loads((row_output / "hash_check_after.json").read_text())
    changed = {record["artifact"]: record["status"] for record in after["artifacts"]}
    assert changed["output_bus"] in {"size_mismatch", "sha256_mismatch"}


def test_run_row_fails_closed_before_processing_a_hash_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest, sample_dir, run_id = _write_run_row_fixture(tmp_path)
    (sample_dir / "summary.txt").write_text("tampered\n", encoding="utf-8")

    from scripts import benchmark_v3_multimap

    def forbidden(_argv):
        raise AssertionError("v3 must not run after a preflight mismatch")

    monkeypatch.setattr(benchmark_v3_multimap, "main", forbidden)
    output_root = tmp_path / "outputs"

    status = run_row(
        manifest,
        run_id,
        output_root,
        source_root=tmp_path / "source",
    )

    assert status["status"] == "failed"
    assert status["errors"][0]["stage"] == "verify_before"
    assert not (output_root / run_id / "legacy").exists()
    assert not (output_root / run_id / "v3").exists()


def test_run_row_requires_exactly_one_manifest_match_before_output(
    tmp_path: Path,
) -> None:
    manifest, _sample_dir, run_id = _write_run_row_fixture(tmp_path)
    lines = manifest.read_text(encoding="utf-8").splitlines()
    manifest.write_text("\n".join([*lines, lines[-1]]) + "\n", encoding="utf-8")

    with pytest.raises(
        ValueError,
        match="must select exactly one raw-manifest row; found 2",
    ):
        run_row(
            manifest,
            run_id,
            tmp_path / "outputs",
            source_root=tmp_path / "source",
        )
    assert not (tmp_path / "outputs").exists()


def test_run_row_rejects_source_containment_before_creating_output(
    tmp_path: Path,
) -> None:
    manifest, sample_dir, run_id = _write_run_row_fixture(tmp_path)
    unsafe_root = sample_dir / "comparison-output"

    with pytest.raises(ValueError, match="outside the read-only archive source root"):
        run_row(
            manifest,
            run_id,
            unsafe_root,
            source_root=tmp_path / "source",
        )

    assert not unsafe_root.exists()


def test_run_row_rejects_manifest_path_inconsistency_before_output(
    tmp_path: Path,
) -> None:
    manifest, _sample_dir, run_id = _write_run_row_fixture(tmp_path)
    with manifest.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    alternate = tmp_path / "alternate.bus"
    alternate.write_text("output_bus\n", encoding="utf-8")
    row["output_bus_path"] = str(alternate)
    row["output_bus_sha256"] = _sha256(alternate)
    row["output_bus_bytes"] = str(alternate.stat().st_size)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    with pytest.raises(ValueError, match="output_bus_path"):
        run_row(
            manifest,
            run_id,
            tmp_path / "outputs",
            source_root=tmp_path / "source",
        )

    assert not (tmp_path / "outputs").exists()


def test_run_row_refuses_to_reuse_a_named_pilot_output(
    tmp_path: Path,
) -> None:
    manifest, _sample_dir, run_id = _write_run_row_fixture(tmp_path)
    row_output = tmp_path / "pilot-run-1" / run_id
    row_output.mkdir(parents=True)
    sentinel = row_output / "existing-status.json"
    sentinel.write_text('{"status":"old"}\n', encoding="utf-8")

    with pytest.raises(ValueError, match="choose a clean named output root"):
        run_row(
            manifest,
            run_id,
            tmp_path / "pilot-run-1",
            source_root=tmp_path / "source",
        )

    assert sentinel.read_text(encoding="utf-8") == '{"status":"old"}\n'


def _write_planned_manifest(path: Path, run_ids: list[str]) -> None:
    fields = ["run_id", "logical_id", "sample_class", "expected_target"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for run_id in run_ids:
            writer.writerow(
                {
                    "run_id": run_id,
                    "logical_id": f"logical-{run_id}",
                    "sample_class": "skin",
                    "expected_target": "",
                }
            )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_summarize_retains_success_failed_and_unsupported_planned_rows(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "raw.tsv"
    run_ids = ["run-success", "run-failed", "run-unsupported", "run-missing"]
    _write_planned_manifest(manifest, run_ids)
    run_root = tmp_path / "runs"
    success = run_root / "run-success"
    failed = run_root / "run-failed"
    unsupported = run_root / "run-unsupported"
    (success / "legacy").mkdir(parents=True)
    (success / "v3").mkdir()
    failed.mkdir(parents=True)
    unsupported.mkdir(parents=True)
    (success / "status.json").write_text(
        '{"status":"success","legacy_audit_status_counts":{"match":3}}\n',
        encoding="utf-8",
    )
    (success / "legacy" / "legacy_reconstruction.json").write_text(
        '{"reconstructed_total_viral_load":2.0}\n', encoding="utf-8"
    )
    (success / "legacy" / "legacy_audit.tsv").write_text(
        "run_id\tlogical_id\trecord_type\tidentifier\treconstructed_value\t"
        "summary_value\tdelta\tstatus\n"
        "run-success\tlogical-run-success\tgene\tEPSTEIN_A\t2\t2\t0\tmatch\n",
        encoding="utf-8",
    )
    (success / "v3" / "result.json").write_text(
        '{"selected_matrix_mass":5.0,"unique_molecule_mass":3.0,'
        '"allocated_ambiguous_mass":2.0,"audit":{"input_molecules":7,'
        '"unique_molecules":3,"ambiguous_molecules":2,"unresolved_molecules":2}}\n',
        encoding="utf-8",
    )
    (success / "v3" / "per_virus.tsv").write_text(
        "virus_name\tn_features\tv3_unique\tv3_equal\tv3_host_conservative\n"
        "Epstein-Barr virus\t1\t3\t5\t4\n",
        encoding="utf-8",
    )
    failure_payload = {
        "status": "failed",
        "errors": [
            {
                "stage": "v3_benchmark",
                "error_type": "RuntimeError",
                "message": "boom",
            }
        ],
    }
    (failed / "status.json").write_text(json.dumps(failure_payload), encoding="utf-8")
    (failed / "failure.json").write_text(json.dumps(failure_payload), encoding="utf-8")
    (unsupported / "status.json").write_text(
        '{"status":"unsupported","reason":"chemistry"}\n', encoding="utf-8"
    )

    result = summarize_results(manifest, run_root, tmp_path / "summary")

    run_rows = _read_tsv(Path(result["run_metrics"]))
    assert [(row["run_id"], row["status"]) for row in run_rows] == [
        ("run-failed", "failed"),
        ("run-missing", "missing"),
        ("run-success", "success"),
        ("run-unsupported", "unsupported"),
    ]
    virus_rows = _read_tsv(Path(result["virus_metrics"]))
    assert {row["run_id"] for row in virus_rows} == set(run_ids)
    legacy_rows = _read_tsv(Path(result["legacy_reproduction"]))
    assert {row["run_id"] for row in legacy_rows} == set(run_ids)
    failures = _read_tsv(Path(result["failures"]))
    assert {(row["run_id"], row["status"]) for row in failures} == {
        ("run-failed", "failed"),
        ("run-missing", "missing"),
        ("run-unsupported", "unsupported"),
    }


def _write_validated_success_run(run_root: Path, run_id: str) -> Path:
    row_dir = run_root / run_id
    legacy = row_dir / "legacy"
    v3 = row_dir / "v3"
    legacy.mkdir(parents=True)
    v3.mkdir()
    (row_dir / "status.json").write_text(
        json.dumps({"run_id": run_id, "status": "success"}), encoding="utf-8"
    )
    for name in ("hash_check_before.json", "hash_check_after.json"):
        (row_dir / name).write_text('{"status":"match"}\n', encoding="utf-8")
    (legacy / "legacy_reconstruction.json").write_text(
        '{"reconstructed_total_viral_load":2.0}\n', encoding="utf-8"
    )
    (legacy / "legacy_audit.tsv").write_text(
        "run_id\tlogical_id\trecord_type\tidentifier\treconstructed_value\t"
        "summary_value\tdelta\tstatus\n"
        f"{run_id}\tlogical\tgene\tEPSTEIN_A\t2\t2\t0\tmatch\n",
        encoding="utf-8",
    )
    result = {
        "sample_id": run_id,
        "n_cells": 1,
        "n_genes": 1,
        "unique_molecule_mass": 3.0,
        "allocated_ambiguous_mass": 2.0,
        "selected_matrix_mass": 5.0,
        "comparison_endpoints": {
            "v3-unique": {},
            "v3-equal": {},
            "v3-host-conservative": {},
        },
        "audit": {
            "input_molecules": 7,
            "unique_molecules": 3,
            "ambiguous_molecules": 2,
            "unresolved_molecules": 2,
        },
    }
    (v3 / "result.json").write_text(json.dumps(result), encoding="utf-8")
    (v3 / "status.json").write_text(
        json.dumps({"sample_id": run_id, "status": "success"}), encoding="utf-8"
    )
    (v3 / "count_audit.json").write_text(json.dumps(result["audit"]), encoding="utf-8")
    exported = ad.AnnData(
        X=sp.csr_matrix([[5.0]]),
        obs=pd.DataFrame(index=["BC1"]),
        var=pd.DataFrame(index=["EPSTEIN_A"]),
    )
    exported.layers["counts_unique"] = sp.csr_matrix([[3.0]])
    exported.layers["counts_ambiguous_allocated"] = sp.csr_matrix([[2.0]])
    exported.layers["counts_multimap_equal"] = sp.csr_matrix([[2.0]])
    exported.layers["counts_multimap_host_conservative"] = sp.csr_matrix([[1.0]])
    exported.write_h5ad(v3 / "adata_v3.h5ad")
    (v3 / "per_cell.tsv").write_text(
        "barcode\tv3_unique_viral\tv3_equal_viral\tv3_host_conservative_viral\nBC1\t3\t5\t4\n",
        encoding="utf-8",
    )
    (v3 / "per_virus.tsv").write_text(
        "virus_name\tn_features\tv3_unique\tv3_equal\tv3_host_conservative\n"
        "Epstein-Barr virus\t1\t3\t5\t4\n",
        encoding="utf-8",
    )
    output_names = [
        "result.json",
        "status.json",
        "count_audit.json",
        "adata_v3.h5ad",
        "per_cell.tsv",
        "per_virus.tsv",
    ]
    hashes = {
        "outputs": {
            name: {
                "size_bytes": (v3 / name).stat().st_size,
                "sha256": _sha256(v3 / name),
            }
            for name in output_names
        }
    }
    (v3 / "hashes.json").write_text(json.dumps(hashes), encoding="utf-8")
    return row_dir


def test_validate_accepts_complete_outputs_and_explicit_failures(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "raw.tsv"
    _write_planned_manifest(manifest, ["success", "failed"])
    run_root = tmp_path / "runs"
    _write_validated_success_run(run_root, "success")
    failed = run_root / "failed"
    failed.mkdir()
    failure = {
        "run_id": "failed",
        "status": "failed",
        "errors": [{"stage": "v3_benchmark", "error_type": "RuntimeError", "message": "boom"}],
    }
    (failed / "status.json").write_text(json.dumps(failure), encoding="utf-8")
    (failed / "failure.json").write_text(json.dumps(failure), encoding="utf-8")

    report = validate_results(manifest, run_root)

    assert report["status"] == "valid"
    assert report["planned_rows"] == 2
    assert report["status_counts"] == {"failed": 1, "success": 1}
    assert report["errors"] == []


def test_validate_fails_closed_on_conservation_hash_and_missing_row(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "raw.tsv"
    _write_planned_manifest(manifest, ["bad", "missing"])
    row_dir = _write_validated_success_run(tmp_path / "runs", "bad")
    result_path = row_dir / "v3" / "result.json"
    payload = json.loads(result_path.read_text())
    payload["selected_matrix_mass"] = 99
    result_path.write_text(json.dumps(payload), encoding="utf-8")

    report = validate_results(manifest, tmp_path / "runs")

    assert report["status"] == "invalid"
    checks = {(error["run_id"], error["check"]) for error in report["errors"]}
    assert ("bad", "matrix_conservation") in checks
    assert ("bad", "output_hashes") in checks
    assert ("missing", "planned_row") in checks
