"""LAYOUT-V3-01: v3 path fallbacks, artifact inventory and docs drift guards."""

import json
import re
from importlib import resources
from pathlib import Path

import jsonschema

from viralscan.artifact_inventory import (
    INVENTORY_REL,
    inventory_paths,
    is_ignored,
    write_layout_docs,
)
from viralscan.outputs import (
    QUESTIONS,
    SAMPLE_LAYOUT,
    analysis_genes_path,
    found_genes_path,
    is_posthoc,
)

ROOT = Path(__file__).resolve().parents[1]


def _touch(base: Path, rel: str, text: str = "x") -> None:
    path = base / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _sample(tmp_path: Path) -> Path:
    sample = tmp_path / "S1"
    for rel in (
        "config.yaml",
        "summary.txt",
        "report.html",
        "results/viral_summary.tsv",
        "results/analysis_genes.txt",
        "results/superexpressor/SuperExpressor_EBV.tsv",
        "plots/EBV_histogram.png",
        "log/kb.done",
        "log/kb_count.log",
        "kb-python/output.bus",
        "kb-python/counts_unfiltered/adata_multimap.h5ad",
        "host_filtered/R1.fastq.gz",
        "posthoc/discovery/candidates.tsv",
        ".snakemake/log/x.log",
    ):
        _touch(sample, rel)
    return sample


def test_gene_lists_fall_back_to_v2_paths(tmp_path: Path) -> None:
    assert analysis_genes_path(tmp_path) == tmp_path / "results" / "analysis_genes.txt"
    _touch(tmp_path, "log/analysis.txt")
    _touch(tmp_path, "log/found_genes.txt")
    assert analysis_genes_path(tmp_path) == tmp_path / "log" / "analysis.txt"
    assert found_genes_path(tmp_path) == tmp_path / "log" / "found_genes.txt"
    _touch(tmp_path, "results/analysis_genes.txt")
    assert analysis_genes_path(tmp_path) == tmp_path / "results" / "analysis_genes.txt"


def test_is_posthoc_sees_both_levels(tmp_path: Path) -> None:
    assert is_posthoc(tmp_path / "posthoc" / "x", tmp_path)  # v2 / rerun-multimap
    assert is_posthoc(tmp_path / "S1" / "posthoc" / "x", tmp_path)  # v3, from the run root
    assert not is_posthoc(tmp_path / "S1" / "results" / "x", tmp_path)


def test_inventory_lists_everything_not_explicitly_ignored(tmp_path: Path) -> None:
    sample = _sample(tmp_path)
    write_layout_docs(tmp_path)
    inventory = json.loads((sample / INVENTORY_REL).read_text())
    listed = {a["path"] for a in inventory["artifacts"]}
    on_disk = {p.relative_to(sample).as_posix() for p in sample.rglob("*") if p.is_file()}
    assert on_disk - listed == {p for p in on_disk if is_ignored(p)}
    assert not listed & {p for p in on_disk if is_ignored(p)}
    assert listed == set(inventory_paths(sample))
    assert "kb-python/counts_unfiltered/adata_multimap.h5ad" in listed
    assert "kb-python/output.bus" not in listed and "log/kb.done" not in listed
    assert "README.txt" in listed and INVENTORY_REL not in listed
    schema = json.loads(
        resources.files("viralscan.schemas.v3")
        .joinpath("artifact_inventory.schema.json")
        .read_text()
    )
    jsonschema.Draft202012Validator(schema).validate(inventory)
    run_readme = (tmp_path / "README.txt").read_text()
    assert "S1/report.html" in run_readme


def test_docs_tree_names_every_sample_folder_the_code_creates() -> None:
    doc = (ROOT / "docs" / "output_reference.md").read_text()
    tree = doc.split("## Directory layout", 1)[1].split("```", 2)[1]
    for path, _role, _text in SAMPLE_LAYOUT:
        assert path in tree or path.rstrip("/") + "/" in tree, f"{path} missing from the doc tree"
    for _question, path in QUESTIONS:
        assert Path(path).name in doc, f"{path} missing from docs/output_reference.md"
    # Every top-level folder a Snakefile rule writes is in the registry.
    snakefile = (ROOT / "src" / "viralscan" / "Snakefile").read_text()
    rule_dirs = set(re.findall(r"\{config\['output'\]\}([A-Za-z_.-]+)/", snakefile))
    known = {p.rstrip("/") for p, _r, _t in SAMPLE_LAYOUT}
    assert rule_dirs <= known, f"Snakefile writes unregistered folders: {rule_dirs - known}"
