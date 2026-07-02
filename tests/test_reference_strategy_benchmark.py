from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from viralscan.reference_strategy import (
    BenchmarkContractError,
    benchmark_rows,
    parse_starsolo_metrics,
    parse_viralscan_metrics,
    parsed_benchmark_rows,
    validate_results,
    write_slurm_array,
)

FIELDS = [
    "dataset",
    "srr",
    "technology",
    "method",
    "reference_strategy",
    "target_virus",
    "reference_hash_id",
    "command_id",
    "slurm_job_id",
    "count_layer",
    "barcode_universe",
    "denominator_fixed_barcodes",
    "denominator_method_called_cells",
    "denominator_shared_anchor_cells",
    "target_positive_cells_fixed",
    "target_umi_fixed",
    "related_off_target_positive_cells_fixed",
    "related_off_target_umi_fixed",
    "combined_vs_two_step_delta_target_umi",
    "status",
    "failure_reason",
]


def _write_results(path: Path, mutate=None) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for base in benchmark_rows():
            row = {
                "dataset": base["dataset"],
                "srr": base["srr"],
                "technology": base["technology"],
                "method": base["method"],
                "reference_strategy": base["reference_strategy"],
                "target_virus": base["target_virus"],
                "reference_hash_id": "sha256:abc",
                "command_id": base["row_id"],
                "slurm_job_id": "12345_0",
                "count_layer": "GeneFull" if base["method"] == "starsolo" else "counts_corrected",
                "barcode_universe": "combined_called_cells",
                "denominator_fixed_barcodes": "1000",
                "denominator_method_called_cells": "900",
                "denominator_shared_anchor_cells": "850",
                "target_positive_cells_fixed": "10",
                "target_umi_fixed": "25",
                "related_off_target_positive_cells_fixed": "1",
                "related_off_target_umi_fixed": "2",
                "combined_vs_two_step_delta_target_umi": "0",
                "status": "complete",
                "failure_reason": "",
            }
            if mutate:
                row = mutate(row)
            writer.writerow(row)


def test_validate_results_accepts_12_row_contract(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"
    _write_results(path)
    rows = validate_results(path)
    assert len(rows) == 12


def test_validate_results_requires_count_layer(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"

    def mutate(row):
        if row["dataset"] == "ebv" and row["method"] == "starsolo":
            row["count_layer"] = ""
        return row

    _write_results(path, mutate=mutate)
    with pytest.raises(BenchmarkContractError, match="count_layer"):
        validate_results(path)


def test_validate_results_rejects_viral_only_two_step_denominator(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"

    def mutate(row):
        if row["reference_strategy"] == "two_step":
            row["barcode_universe"] = "viral_only_called_cells"
        return row

    _write_results(path, mutate=mutate)
    with pytest.raises(BenchmarkContractError, match="two-step viral-only denominator"):
        validate_results(path)


def test_validate_results_rejects_inconsistent_fixed_denominator(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"

    def mutate(row):
        if row["dataset"] == "ebv" and row["method"] == "starsolo":
            row["denominator_fixed_barcodes"] = "999"
        return row

    _write_results(path, mutate=mutate)
    with pytest.raises(BenchmarkContractError, match="inconsistent fixed barcode denominator"):
        validate_results(path)


def test_validate_results_requires_complete_row_provenance(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"

    def mutate(row):
        if row["dataset"] == "hsv1" and row["method"] == "viralscan":
            row["reference_hash_id"] = ""
        return row

    _write_results(path, mutate=mutate)
    with pytest.raises(BenchmarkContractError, match="reference_hash_id"):
        validate_results(path)


def test_validate_results_rejects_historical_single_virus_tokens(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"

    def mutate(row):
        if row["dataset"] == "ebv":
            row["failure_reason"] = "used GRCh38_EBV by mistake"
        return row

    _write_results(path, mutate=mutate)
    with pytest.raises(BenchmarkContractError, match="forbidden"):
        validate_results(path)


def test_parse_starsolo_metrics_counts_target_and_off_target(tmp_path: Path) -> None:
    raw = tmp_path / "Solo.out" / "GeneFull" / "raw"
    filtered = tmp_path / "Solo.out" / "GeneFull" / "filtered"
    raw.mkdir(parents=True)
    filtered.mkdir(parents=True)
    (raw / "features.tsv").write_text(
        "gene1\tEpstein Barr virus gene\tGene Expression\n"
        "gene2\tHuman herpesvirus 1 gene\tGene Expression\n"
        "gene3\tACTB\tGene Expression\n",
        encoding="utf-8",
    )
    (raw / "barcodes.tsv").write_text("BC1\nBC2\n", encoding="utf-8")
    (raw / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 5\n2 1 2\n1 2 1\n",
        encoding="utf-8",
    )
    (filtered / "barcodes.tsv").write_text("BC1\n", encoding="utf-8")

    metrics = parse_starsolo_metrics(
        raw,
        filtered / "barcodes.tsv",
        target_regex=r"(?i)epstein|ebv",
        off_target_regex=r"(?i)herpesvirus 1",
    )

    assert metrics.status == "complete"
    assert metrics.counts_by_barcode["BC1"] == (5.0, 2.0)
    assert metrics.counts_by_barcode["BC2"] == (1.0, 0.0)
    assert metrics.anchor_barcodes == {"BC1"}


def test_parse_viralscan_metrics_counts_target_and_off_target(tmp_path: Path) -> None:
    results = tmp_path / "SRR" / "results"
    results.mkdir(parents=True)
    (results / "per_cell_viral.tsv").write_text(
        "barcode\tvirus_name\tviral_umi\ttotal_umi\tviral_fraction\n"
        "BC1\tEpstein Barr virus\t3\t3\t1\n"
        "BC1\tHuman herpesvirus 1\t2\t5\t0.4\n"
        "BC2\tEpstein Barr virus\t1\t1\t1\n",
        encoding="utf-8",
    )
    (results / "viral_summary.tsv").write_text(
        "virus_name\ttotal_umi\tinfected_cells\ttotal_cells\tpct_infected\tumi_per_10k\n"
        "Epstein Barr virus\t4\t2\t100\t2\t400\n",
        encoding="utf-8",
    )

    metrics = parse_viralscan_metrics(
        tmp_path,
        target_regex=r"(?i)epstein|ebv",
        off_target_regex=r"(?i)herpesvirus 1",
    )

    assert metrics.status == "complete"
    assert metrics.counts_by_barcode["BC1"] == (3.0, 2.0)
    assert metrics.counts_by_barcode["BC2"] == (1.0, 0.0)
    assert metrics.anchor_barcodes == {"BC1", "BC2"}


def test_parsed_benchmark_rows_preserves_run_status_for_missing_outputs(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "reference_audit.tsv").write_text("key\tpath\nexample\tref\n", encoding="utf-8")
    with (run_dir / "commands.jsonl").open("w", encoding="utf-8") as handle:
        for row in benchmark_rows():
            handle.write(json.dumps(row) + "\n")
    (run_dir / "run_status.tsv").write_text(
        "dataset\tsrr\tmethod\treference_strategy\tstatus\tfailure_reason\n"
        "hhv6b\tSRR20710641\tviralscan\ttwo_step\tblocked\tkallisto missing\n",
        encoding="utf-8",
    )

    rows = parsed_benchmark_rows(run_dir)
    row = next(
        row
        for row in rows
        if row["dataset"] == "hhv6b"
        and row["method"] == "viralscan"
        and row["reference_strategy"] == "two_step"
    )

    assert row["status"] == "blocked"
    assert row["failure_reason"] == "kallisto missing"


def test_write_slurm_array_preflights_required_runtime_tools(tmp_path: Path) -> None:
    script = write_slurm_array(tmp_path)
    text = script.read_text(encoding="utf-8")

    assert "KB_PYTHON_BIN_DIR=$(python - <<'PY'" in text
    assert 'export PATH="$KB_PYTHON_BIN_DIR:$PATH"' in text
    assert "for tool in python kb snakemake kallisto bustools" in text
    assert "missing required benchmark runtime tools" in text
    assert "exit 127" in text
