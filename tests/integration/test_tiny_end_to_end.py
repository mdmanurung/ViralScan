"""Integration test — the documented CLI, end to end, on a tiny fixture (SW-10).

This is the ``G1`` exit clause "the full tiny workflow is green". It drives the
real ``viralscan`` entry point over ``tests/data/evidence_tiny`` — two read pairs
against a two-sequence reference — through preflight, quantification, molecule
allocation, cell calling, summaries, and ``validate-run``.

Run with::

    PYTHONPATH=src python -m pytest -m integration tests/integration/test_tiny_end_to_end.py

Requires ``kallisto``, ``bustools``, ``kb`` and ``snakemake`` on PATH; skipped
otherwise. ``kb`` is slow to start, so expect a few minutes.

Two things this pinned down that no unit test could
---------------------------------------------------
1. ``run_manifest.json`` sits at the **root** of the result tree, beside the
   per-sample directories rather than inside one. A review finding claimed the
   opposite and a fix was written to match it; running the pipeline showed
   ``out/run_manifest.json`` next to ``out/<sample>/config.yaml`` and the change
   was reverted. See ``test_run_manifest_is_at_the_tree_root``.
2. A default run cannot complete at all without ``-gtf``, because the bundled
   viral panel is fetched from an unregistered Zenodo DOI (``REF-11``). The
   fixture therefore ships its own minimal GTF. When ``REF-11`` is resolved, a
   variant of this test should drop ``-gtf`` and exercise the default path.

What this does NOT cover
------------------------
The evidence/BAM/BLAST/IGV leg, which needs ``blastn``, ``makeblastdb`` and
``minimap2``. Those are exercised separately by ``test_exact_lineage.py`` and
``test_evidence_chain.py``.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from viralscan.evidence import have_tools

pytestmark = pytest.mark.integration

FIXTURE = Path(__file__).parents[1] / "data" / "evidence_tiny"
REQUIRED_TOOLS = ["kallisto", "bustools", "kb", "snakemake"]

#: Artifacts a completed run must publish, relative to the sample directory.
EXPECTED_ARTIFACTS = (
    "config.yaml",
    "count_audit.tsv",
    "summary.txt",
    "report.html",
    "results/virus_identity.tsv",
    "results/viral_summary.tsv",
    "results/per_cell_viral.tsv",
    "results/reference_provenance.json",
    "kb-python/counts_unfiltered/adata_multimap.h5ad",
    "kb-python/output.resolved.sorted.bus",
    "log/create_config.done",
    "log/kb.done",
    "log/multimap.done",
    "log/detection.done",
)


def _viralscan(*args: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "viralscan.menu", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
        timeout=1800,
    )


@pytest.fixture(scope="module")
def completed_run(tmp_path_factory) -> Path:
    """Build an index and take one tiny sample through the documented commands."""
    missing = have_tools(REQUIRED_TOOLS)
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

    work = tmp_path_factory.mktemp("sw10")
    repo = Path(__file__).parents[2]
    index = work / "index.idx"
    subprocess.run(
        ["kallisto", "index", "-i", str(index), str(FIXTURE / "transcripts.fasta")],
        check=True,
        capture_output=True,
        timeout=600,
    )

    out = work / "out"
    _viralscan(
        "-o",
        str(out),
        "-i",
        str(index),
        "-t",
        str(FIXTURE / "t2g.txt"),
        "-s1",
        str(FIXTURE / "R1.fastq"),
        "-s2",
        str(FIXTURE / "R2.fastq"),
        # The bundled panel is unfetchable until REF-11 is resolved.
        "-gtf",
        str(FIXTURE / "viral.gtf"),
        "-x",
        "10xv3",
        "-c",
        "2",
        "--no-visual",
        "--cell-calling",
        "none",
        "--yes",
        cwd=repo,
    )
    return out


def _sample_dir(run: Path) -> Path:
    samples = sorted(p.parent for p in run.glob("*/config.yaml"))
    assert len(samples) == 1, f"expected one sample directory, found {samples}"
    return samples[0]


class TestTinyWorkflowCompletes:
    @pytest.mark.parametrize("relative", EXPECTED_ARTIFACTS)
    def test_every_documented_artifact_is_published(self, completed_run, relative) -> None:
        assert (_sample_dir(completed_run) / relative).is_file()

    def test_run_manifest_is_at_the_tree_root(self, completed_run) -> None:
        """Pins the layout `rerun-multimap` depends on. See the module docstring."""
        assert (completed_run / "run_manifest.json").is_file()
        assert not (_sample_dir(completed_run) / "run_manifest.json").exists()

    def test_the_manifest_records_the_method_and_unit_actually_used(self, completed_run) -> None:
        manifest = json.loads((completed_run / "run_manifest.json").read_text(encoding="utf-8"))

        assert manifest["schema_version"] == "3.0.0"
        assert manifest["quantification_unit"] == "bustools-resolved-cb-umi-molecule"
        assert manifest["allocation_method"] == manifest["options"]["multimap_method"]

    def test_summary_leads_with_the_molecule_totals(self, completed_run) -> None:
        """The block multimap used to write and detection used to truncate away."""
        summary = (_sample_dir(completed_run) / "summary.txt").read_text(encoding="utf-8")

        assert "Viral molecules in unique-count matrix:" in summary
        assert "Total viral molecules (selected method):" in summary
        assert "Cells with viral reads:" in summary

    def test_summary_headline_matches_viral_summary_total(self, completed_run) -> None:
        """SW-16: the headline said 0 beside 906,202 EBV molecules in viral_summary.tsv."""
        sample = _sample_dir(completed_run)
        summary = (sample / "summary.txt").read_text(encoding="utf-8")
        headline = next(
            line
            for line in summary.splitlines()
            if line.startswith("Total viral molecules (selected method):")
        )
        headline_total = float(headline.split(":", 1)[1].strip())
        with open(sample / "results" / "viral_summary.tsv", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        table_total = sum(float(r["viral_molecules_total_est"]) for r in rows)

        assert table_total > 0
        assert headline_total == pytest.approx(table_total)

    def test_virus_identity_table_resolves_the_fixture_index(self, completed_run) -> None:
        """MECH-A step 3: the fixture's 3-column t2g takes the legacy fallback."""
        path = _sample_dir(completed_run) / "results" / "virus_identity.tsv"
        with open(path, encoding="utf-8") as fh:
            rows = {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

        assert rows["VIRUS_TARGET"]["status"] == "legacy_prefix"
        assert rows["VIRUS_TARGET"]["viral"] == "true"
        assert rows["HOST_GENE"]["status"] == "host"


class TestValidateRunAcceptsTheResult:
    def test_validate_run_reports_no_issues(self, completed_run, tmp_path) -> None:
        """Exercises every SW-02 schema check against a real artifact."""
        report_path = tmp_path / "validate.json"
        _viralscan(
            "validate-run",
            str(completed_run),
            "--json-output",
            str(report_path),
            cwd=Path(__file__).parents[2],
        )
        report = json.loads(report_path.read_text(encoding="utf-8"))

        assert report["ok"] is True, report["issues"]
        assert report["issues"] == []
        assert len(report["h5ad_files"]) == 1
