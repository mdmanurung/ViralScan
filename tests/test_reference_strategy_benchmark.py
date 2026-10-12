from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pytest

from viralscan.reference_strategy import (
    DATASETS,
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


# ---------------------------------------------------------------------------
# target_regex vs STARsolo gene naming (regression guard)
# ---------------------------------------------------------------------------
# STARsolo's Serratus gene_ids use the `HHV<n>` form for human herpesviruses
# (HHV1 == HSV-1, HHV4 == EBV). A benchmark bug counted STARsolo HSV-1 as 0
# because the HSV-1 target_regex omitted the `hhv-?1` alias that the EBV
# (`hhv-?4`) and HHV-6B (`hhv-?6b`) regexes already had. These tests pin the fix.
_DATASETS_BY_VIRUS = {ds["target_virus"]: ds for ds in DATASETS}


def test_hsv1_target_regex_matches_starsolo_hhv1_gene_ids():
    regex = _DATASETS_BY_VIRUS["HSV-1"]["target_regex"]
    for gene_id in ("HHV1gp00s01", "HHV1gp00p76", "NC_001806"):
        assert re.search(regex, gene_id), (
            f"HSV-1 target_regex {regex!r} must match STARsolo gene id {gene_id!r}"
        )


def test_hsv1_target_regex_does_not_match_other_herpesviruses():
    regex = _DATASETS_BY_VIRUS["HSV-1"]["target_regex"]
    # `hhv-?1(?![0-9])` must not leak onto EBV (HHV4), HHV-6B, or a hypothetical HHV-1x.
    for gene_id in ("HHV4_EBNA-1.1", "HHV6Bgp01", "HHV12gp01"):
        assert not re.search(regex, gene_id), (
            f"HSV-1 target_regex {regex!r} should not match {gene_id!r}"
        )


def test_ebv_target_regex_matches_starsolo_hhv4_gene_ids():
    regex = _DATASETS_BY_VIRUS["EBV"]["target_regex"]
    assert re.search(regex, "HHV4_EBNA-1.1")


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


_V3_PER_CELL = (
    "barcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction\t"
    "is_called_cell\n"
)


def _write_viralscan_results(root: Path, per_cell: str, called: str | None = None) -> None:
    results = root / "results"
    results.mkdir(parents=True)
    (results / "per_cell_viral.tsv").write_text(per_cell, encoding="utf-8")
    (results / "viral_summary.tsv").write_text(
        "virus_name\ttotal_cells\nEpstein Barr virus\t100\n", encoding="utf-8"
    )
    if called is not None:
        (results / "called_cells.tsv").write_text(called, encoding="utf-8")


def test_parse_viralscan_metrics_reads_the_v3_schema(tmp_path: Path) -> None:
    """AUD-20 / D-06: real v3 rows carry viral_molecules_total_est, not viral_umi; the anchor
    is the called-cell universe (called_cells.tsv), not the infected cells."""
    _write_viralscan_results(
        tmp_path,
        _V3_PER_CELL
        + "BC1\tEpstein Barr virus\t3\t1000\t0.003\ttrue\n"
        + "BC1\tHuman herpesvirus 1\t2.5\t1000\t0.0025\ttrue\n"
        + "BC2\tEpstein Barr virus\t1\t500\t0.002\tfalse\n",
        called="barcode\nBC1\nBC7\nBC8\n",
    )
    metrics = parse_viralscan_metrics(
        tmp_path, target_regex=r"(?i)epstein|ebv", off_target_regex=r"(?i)herpesvirus 1"
    )
    assert metrics.status == "complete", metrics.failure_reason
    assert metrics.counts_by_barcode["BC1"] == (3.0, 2.5)
    assert metrics.counts_by_barcode["BC2"] == (1.0, 0.0)
    assert metrics.anchor_barcodes == {"BC1", "BC7", "BC8"}


def test_parse_viralscan_metrics_anchor_falls_back_to_is_called_cell(tmp_path: Path) -> None:
    _write_viralscan_results(
        tmp_path,
        _V3_PER_CELL
        + "BC1\tEpstein Barr virus\t3\t1000\t0.003\ttrue\n"
        + "BC2\tEpstein Barr virus\t1\t500\t0.002\tfalse\n",
    )
    metrics = parse_viralscan_metrics(tmp_path, target_regex="(?i)epstein", off_target_regex="x^")
    assert metrics.anchor_barcodes == {"BC1"}


def test_parse_viralscan_metrics_tolerates_extra_and_reordered_columns(tmp_path: Path) -> None:
    """Named-column reads: new columns (e.g. allocated molecules) and a different order."""
    _write_viralscan_results(
        tmp_path,
        "extra\tvirus_name\tbarcode\tviral_molecules_total_est\tallocated_molecules\n"
        "z\tEpstein Barr virus\tBC1\t3\t3.0\n",
    )
    metrics = parse_viralscan_metrics(tmp_path, target_regex="(?i)epstein", off_target_regex="x^")
    assert metrics.status == "complete"
    assert metrics.counts_by_barcode == {"BC1": (3.0, 0.0)}


def test_parse_viralscan_metrics_accepts_real_detection_output(tmp_path: Path) -> None:
    """The parser must accept what detection.compute_stats actually writes."""
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    from viralscan.scripts.detection import compute_stats

    adata = ad.AnnData(sp.csr_matrix(np.array([[10, 0, 100], [0, 5, 200], [1, 0, 3]], float)))
    adata.obs_names = ["bc0", "bc1", "bc2"]
    adata.var_names = ["v1", "v2", "h1"]
    _, per_cell = compute_stats(
        adata, {}, {"Epstein Barr virus": ["v1", "v2"]}, [], called_mask=np.array([1, 1, 0], bool)
    )
    _write_viralscan_results(tmp_path, per_cell.to_csv(sep="\t", index=False))
    metrics = parse_viralscan_metrics(tmp_path, target_regex="(?i)epstein", off_target_regex="x^")
    assert metrics.status == "complete", metrics.failure_reason
    assert metrics.counts_by_barcode["bc0"] == (10.0, 0.0)
    assert metrics.anchor_barcodes == {"bc0", "bc1"}


def _metrics(counts, status="complete"):
    from viralscan.reference_strategy import ParsedMetrics

    return ParsedMetrics(counts, set(counts), status, "" if status == "complete" else "boom")


def _parse_with(tmp_path: Path, table: dict) -> dict:
    from unittest.mock import patch

    from viralscan import reference_strategy as rs

    run = tmp_path / "run"
    run.mkdir()
    (run / "reference_audit.tsv").write_text("key\n", encoding="utf-8")
    with (run / "commands.jsonl").open("w", encoding="utf-8") as handle:
        for row in benchmark_rows():
            handle.write(json.dumps(row) + "\n")
    empty = _metrics({}, "incomplete")
    with patch.object(rs, "_row_metrics", lambda _d, row: table.get(row["row_id"], empty)):
        rows = parsed_benchmark_rows(run)
    return {r["command_id"]: r for r in rows}


def test_incomplete_combined_row_does_not_shrink_the_barcode_universe(tmp_path: Path) -> None:
    """AUD-20 / D-07: with STARsolo-combined failed the fixed universe was built from the
    ViralScan infected cells only (2 of 3 barcodes) and the rows still read 'complete'."""
    rows = _parse_with(
        tmp_path,
        {
            "ebv__starsolo__combined": _metrics({}, "incomplete"),
            "ebv__starsolo__two_step": _metrics(
                {"A": (5.0, 0.0), "B": (5.0, 0.0), "C": (5.0, 0.0)}
            ),
            "ebv__viralscan__combined": _metrics({"A": (4.0, 0.0), "B": (4.0, 0.0)}),
            "ebv__viralscan__two_step": _metrics({"A": (4.0, 0.0), "B": (4.0, 0.0)}),
        },
    )
    for cid in ("ebv__starsolo__two_step", "ebv__viralscan__combined", "ebv__viralscan__two_step"):
        assert rows[cid]["status"] == "incomplete", cid
        assert "fixed barcode universe unavailable" in rows[cid]["failure_reason"]
        assert "ebv__starsolo__combined" in rows[cid]["failure_reason"]
        assert rows[cid]["combined_vs_two_step_delta_target_umi"] == "NA"


def test_combined_vs_two_step_delta_is_combined_minus_two_step_on_both_rows(
    tmp_path: Path,
) -> None:
    table = {
        "ebv__starsolo__combined": _metrics({"A": (6.0, 0.0), "B": (4.0, 0.0)}),
        "ebv__starsolo__two_step": _metrics({"A": (5.0, 0.0), "B": (1.0, 0.0)}),
        "ebv__viralscan__combined": _metrics({"A": (6.0, 0.0), "B": (4.0, 0.0)}),
        "ebv__viralscan__two_step": _metrics({"A": (6.0, 0.0), "B": (4.0, 0.0)}),
    }
    rows = _parse_with(tmp_path, table)
    assert rows["ebv__starsolo__combined"]["status"] == "complete"
    for cid in ("ebv__starsolo__combined", "ebv__starsolo__two_step"):
        assert rows[cid]["combined_vs_two_step_delta_target_umi"] == "4.0"
    for cid in ("ebv__viralscan__combined", "ebv__viralscan__two_step"):
        assert rows[cid]["combined_vs_two_step_delta_target_umi"] == "0.0"


def test_validate_results_accepts_na_delta(tmp_path: Path) -> None:
    path = tmp_path / "reference_strategy_benchmark.tsv"

    def mutate(row):
        row["combined_vs_two_step_delta_target_umi"] = "NA"
        return row

    _write_results(path, mutate=mutate)
    assert len(validate_results(path)) == 12


def test_hsv1_and_ebv_regexes_do_not_match_non_human_herpesviruses() -> None:
    """AUD-20 / D-08: 'Macacine alphaherpesvirus 1' (B virus) is not HSV-1."""
    by_id = {d["dataset"]: d for d in DATASETS}
    hsv1, ebv = by_id["hsv1"], by_id["ebv"]
    assert re.search(hsv1["target_regex"], "Human herpesvirus 1")
    assert re.search(hsv1["target_regex"], "Human alphaherpesvirus 1")
    assert re.search(hsv1["target_regex"], "HHV1gp001")
    assert not re.search(hsv1["target_regex"], "Macacine alphaherpesvirus 1")
    assert not re.search(hsv1["target_regex"], "Human herpesvirus 12")
    assert re.search(hsv1["off_target_regex"], "Human herpesvirus 2")
    assert not re.search(hsv1["off_target_regex"], "Macacine alphaherpesvirus 2")
    assert re.search(ebv["off_target_regex"], "Human herpesvirus 8")
    assert not re.search(ebv["off_target_regex"], "Macacine alphaherpesvirus 1")
    assert not re.search(ebv["off_target_regex"], "Papiine betaherpesvirus 3")


def test_slurm_array_script_never_chdirs_to_a_placeholder(tmp_path: Path) -> None:
    """AUD-20 / D-10: the array script chdir'd to /path/to/ViralScan."""
    text = write_slurm_array(tmp_path).read_text(encoding="utf-8")
    assert "/path/to" not in text
    assert 'os.chdir(os.environ["VS_BENCH_REPO"])' in text
    assert 'VS_BENCH_REPO="${VS_BENCH_REPO:-' in text


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


def test_viralscan_two_step_row_parses_with_the_real_cli(tmp_path: Path) -> None:
    """PLAN SW-25: the harness argv must be accepted by menu.py, not just built."""
    from viralscan.menu import build_parser
    from viralscan.reference_strategy import commands_for_row, default_manifest

    row = next(
        r
        for r in benchmark_rows()
        if r["method"] == "viralscan" and r["reference_strategy"] == "two_step"
    )
    (cmd,) = commands_for_row(row, default_manifest(), tmp_path)
    ns = build_parser().parse_args([a.replace("__VS_THREADS__", "1") for a in cmd[2:]])

    assert ns.host_filter == "starsolo"
    assert ns.host_index.endswith("human_GRCh38_2024A")
