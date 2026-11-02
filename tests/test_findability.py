"""FINDABILITY-01: start-of-run README, report file table, summary footer, overwrite guard, docs."""

import argparse
import json
import re
from pathlib import Path

import pytest

from viralscan import menu
from viralscan.artifact_inventory import write_layout_docs, write_start_docs
from viralscan.outputs import RESULT_FILES, result_files, safe_to_delete
from viralscan.run_safety import RunSafetyError, require_fresh_output

ROOT = Path(__file__).resolve().parents[1]


def _run(tmp_path: Path, options: dict) -> Path:
    (tmp_path / "run_manifest.json").write_text(
        json.dumps({"options": options, "run_fingerprint": "abc123", "allocation_method": "em"})
    )
    return tmp_path


def test_start_readme_says_not_finished_and_lists_only_enabled_options(tmp_path: Path) -> None:
    run = _run(tmp_path, {"umap": True, "gene_programs": False})
    write_start_docs(run, ["S1"])
    sample = (run / "S1" / "README.txt").read_text()
    assert "NOT FINISHED" in sample
    assert "plots/umap_binary.html" in sample
    assert "gene_program_summary.tsv" not in sample
    assert "not detected" in sample and "results/sensitivity.tsv" in sample
    for path, _text in safe_to_delete():
        assert path in sample
    top = (run / "README.txt").read_text()
    assert "NOT FINISHED" in top and "abc123" in top and "multimap method   em" in top


def test_completion_rewrites_readme_as_finished(tmp_path: Path) -> None:
    run = _run(tmp_path, {})
    write_start_docs(run, ["S1"])
    (run / "S1" / "config.yaml").write_text("x")
    write_layout_docs(run)
    assert "finished" in (run / "S1" / "README.txt").read_text()
    assert "NOT FINISHED" not in (run / "README.txt").read_text()
    assert "plots/umap_binary.html" not in (run / "S1" / "README.txt").read_text()


def test_report_table_comes_from_registry_not_template() -> None:
    assert {p for p, _ in result_files(lambda _o: False)}.isdisjoint(
        {"plots/umap_binary.html", "results/gene_program_cells.tsv"}
    )
    assert len(result_files(lambda _o: True)) == len(RESULT_FILES)
    template = (ROOT / "src/viralscan/templates/report.html.j2").read_text()
    assert "umap_binary" not in template and "output_files" in template


def test_summary_footer_is_not_duplicated_on_resume(tmp_path: Path) -> None:
    (tmp_path / "summary.txt").write_text("headline\n")
    menu._write_sample_summary(str(tmp_path), 1.0, 10, 3)
    menu._write_sample_summary(str(tmp_path), 2.0, 10, 3)
    text = (tmp_path / "summary.txt").read_text()
    assert text.count("Runtime:") == 1 and "3.0000" in text and text.startswith("headline\n")


def test_require_fresh_output(tmp_path: Path) -> None:
    require_fresh_output(tmp_path / "absent")
    require_fresh_output(tmp_path)  # empty
    (tmp_path / "f").write_text("x")
    with pytest.raises(RunSafetyError, match="--force"):
        require_fresh_output(tmp_path)
    with pytest.raises(RunSafetyError) as exc:
        require_fresh_output(tmp_path, forceable=False)
    assert "--force" not in str(exc.value)
    require_fresh_output(tmp_path, force=True)


def test_post_hoc_verb_refuses_existing_output_without_force(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    (out / "old.json").write_text("{}")
    with pytest.raises(SystemExit):
        menu._require_fresh(out, argparse.Namespace(force=False))
    menu._require_fresh(out, argparse.Namespace(force=True))


def test_docs_describe_every_multimap_layer() -> None:
    code = (ROOT / "src/viralscan/scripts/multimap.py").read_text()
    layers = set(re.findall(r'adata\.layers\["(counts_[a-z_]+)"\]', code))
    doc = (ROOT / "docs/output_reference.md").read_text()
    missing = {name for name in layers if f"`{name}`" not in doc}
    # the four per-method host_viral_selected variants share one table row
    assert {m for m in missing if not m.startswith("counts_host_viral_selected_")} == set()
