"""Output layout: v3 paths, layout_version, and the run-level combined summary."""

import argparse

from viralscan.outputs import LAYOUT_VERSION, V3Outputs
from viralscan.run_safety import build_run_manifest, write_combined_summary


def test_v3_outputs_paths(tmp_path):
    out = V3Outputs.from_config_output(f"{tmp_path}/S/")
    assert out.results_dir == tmp_path / "S" / "results"
    assert out.read_filtered_dir == tmp_path / "S" / "read_filtered"
    assert out.count_audit == tmp_path / "S" / "results" / "count_audit.tsv"


def test_manifest_carries_layout_version():
    args = argparse.Namespace(multimap_method="equal", cell_calling="none", called_cells_file=None)
    assert build_run_manifest(args)["layout_version"] == LAYOUT_VERSION


def _sample(root, name, rows):
    results = root / name / "results"
    results.mkdir(parents=True)
    (root / name / "config.yaml").write_text("output: x\n")
    (results / "viral_summary.tsv").write_text(
        "virus\tn\n" + "".join(f"{v}\t{n}\n" for v, n in rows)
    )


def test_combined_summary_stacks_samples(tmp_path):
    _sample(tmp_path, "GSM1", [("HHV-6B", 10)])
    _sample(tmp_path, "GSM2", [("HHV-6B", 0)])
    assert write_combined_summary(tmp_path)
    lines = (tmp_path / "combined_summary.tsv").read_text().splitlines()
    assert lines == ["sample\tvirus\tn", "GSM1\tHHV-6B\t10", "GSM2\tHHV-6B\t0"]


def test_combined_summary_header_only_for_zero_virus_runs(tmp_path):
    _sample(tmp_path, "GSM1", [])
    (tmp_path / "combined_summary.tsv").write_text("sample\tvirus\tn\nOLD\tX\t1\n")
    assert write_combined_summary(tmp_path)
    assert (tmp_path / "combined_summary.tsv").read_text().splitlines() == ["sample\tvirus\tn"]


def test_combined_summary_skipped_without_summaries(tmp_path):
    (tmp_path / "S").mkdir()
    (tmp_path / "S" / "config.yaml").write_text("output: x\n")
    assert not write_combined_summary(tmp_path)
    assert not (tmp_path / "combined_summary.tsv").exists()
