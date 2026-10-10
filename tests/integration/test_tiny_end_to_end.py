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
2. A run without ``-gtf`` and without an index build manifest cannot complete,
   because the bundled viral panel is fetched from an unregistered Zenodo DOI
   (``REF-11``). The fixture therefore ships its own minimal GTF. Once the
   index has a build manifest (``CAT-06``), ``analysis`` reads the viral gene
   IDs from it and ``-gtf`` is unnecessary: see ``test_manifest_run_*``.

The whole documented sequence, one fixture (SW-10)
--------------------------------------------------
``TestDocumentedSequence`` chains every documented CLI step on the same tiny
pair, in order, each step consuming the previous step's real output:

1. ``viralscan doctor --profile full`` (needs STAR, blastn, makeblastdb,
   minimap2, samtools, Rscript, cd-hit-est, kb, snakemake on PATH).
2. Reference build through the documented ``-ref -fasta -gtf`` path
   (``kb ref``), not a hand-run ``kallisto index``. ``viralscan build-ref`` is
   not usable here because it downloads Ensembl/NCBI sequence; the offline
   equivalent is ``-ref`` on ``reference.fasta`` + ``combined.gtf``. Fixture
   rule: the viral genome record must NOT share its name with a transcript ID,
   or ``virus_identity`` reads the gene as host cDNA (``structural_host``); the
   host record is named after its transcript on purpose (the host-cDNA shape).
3. A run **without** ``--no-visual``: plots and ``report.html`` must exist.
4. ``--cell-calling external --called-cells-file called_cells.txt``.
5. ``viralscan evidence --run-dir <real sample dir> ... --blast`` against the
   real run (not the hand-faked layout of ``test_exact_lineage.py``), asserting
   the 15 artifacts that test lists.
6. ``viralscan validate-run`` again after evidence.

The pre-existing ``completed_run`` tests above keep the cheaper
``kallisto index`` + ``--no-visual`` path.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from viralscan.evidence import have_tools
from viralscan.validation import tool_path

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


def _run_tiny(work: Path, *extra: str, manifest: bool = False) -> Path:
    """Build an index and take one tiny sample through the documented commands.

    ``manifest=True`` writes the index build manifest and drops ``-gtf`` (CAT-06).
    """
    missing = have_tools(REQUIRED_TOOLS)
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

    repo = Path(__file__).parents[2]
    index = work / "index.idx"
    subprocess.run(
        # The binary kb count will read the index with (REL-16): a conda kallisto
        # of the same version builds an index kb's bundled one spins on forever.
        [tool_path("kallisto"), "index", "-i", str(index), str(FIXTURE / "transcripts.fasta")],
        check=True,
        capture_output=True,
        timeout=600,
    )

    if manifest:
        from viralscan.virus_identity import write_build_manifest

        write_build_manifest(index, ["HOST_GENE"], ["VIRUS_TARGET"], {"builder": "test"})
        gtf_args: tuple[str, ...] = ()
    else:
        # Without a manifest the bundled panel is unfetchable until REF-11 is resolved.
        gtf_args = ("-gtf", str(FIXTURE / "viral.gtf"))

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
        *gtf_args,
        "-x",
        "10xv3",
        "-c",
        "2",
        "--no-visual",
        "--cell-calling",
        "none",
        "--yes",
        *extra,
        cwd=repo,
    )
    return out


@pytest.fixture(scope="module")
def completed_run(tmp_path_factory) -> Path:
    return _run_tiny(tmp_path_factory.mktemp("sw10"))


@pytest.fixture(scope="module")
def strand_run(tmp_path_factory) -> Path:
    return _run_tiny(tmp_path_factory.mktemp("sw10_strand"), "--strand", "unstranded")


@pytest.fixture(scope="module")
def read_filter_run(tmp_path_factory) -> Path:
    return _run_tiny(tmp_path_factory.mktemp("def01"), "--read-filter", "artefact")


@pytest.fixture(scope="module")
def manifest_run(tmp_path_factory) -> Path:
    return _run_tiny(tmp_path_factory.mktemp("cat06"), manifest=True)


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


class TestManifestRunWithoutGtf:
    """CAT-06 / SW-10: no ``-gtf``, no Zenodo panel; the build manifest names the viruses."""

    def test_manifest_run_completes_and_detects_the_virus(self, manifest_run) -> None:
        sample = _sample_dir(manifest_run)
        for relative in EXPECTED_ARTIFACTS:
            assert (sample / relative).is_file(), relative
        assert "VIRUS_TARGET" in (sample / "log" / "analysis.txt").read_text().split()
        with open(sample / "results" / "viral_summary.tsv", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert sum(float(r["viral_molecules_total_est"]) for r in rows) > 0


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


class TestCompletionMarker:
    """SW-06: a finished run carries run_complete.json; deleting it fails validate-run."""

    def test_marker_exists_and_matches_manifest(self, completed_run) -> None:
        marker = json.loads((completed_run / "run_complete.json").read_text(encoding="utf-8"))
        manifest = json.loads((completed_run / "run_manifest.json").read_text(encoding="utf-8"))
        assert manifest["completion_marker"] is True
        assert marker["run_fingerprint"] == manifest["run_fingerprint"]
        assert any(k.endswith("results/viral_summary.tsv") for k in marker["artifacts"])

    def test_validate_run_fails_without_the_marker(self, completed_run, tmp_path) -> None:
        import shutil

        copy = tmp_path / "copy"
        shutil.copytree(completed_run, copy)
        (copy / "run_complete.json").unlink()
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "viralscan.menu",
                "validate-run",
                str(copy),
                "--no-verify-inputs",
            ],
            cwd=Path(__file__).parents[2],
            capture_output=True,
            text=True,
            timeout=600,
        )
        assert result.returncode != 0
        assert "missing_completion_marker" in result.stdout


class TestStrandIsPassedToKb:
    """DEF-02: ``--strand`` reaches ``kb count`` and is recorded for resume safety."""

    def test_run_completes_and_kb_log_shows_the_flag(self, strand_run) -> None:
        sample = _sample_dir(strand_run)
        assert (sample / "log" / "kb.done").is_file()
        log = (sample / "kb-python" / "kb_info.json").read_text()  # records the kb call
        assert "--strand unstranded" in log, log[:2000]
        assert "kallisto bus" in log and "--unstranded" in log

    def test_manifest_records_strand(self, strand_run) -> None:
        manifest = json.loads((strand_run / "run_manifest.json").read_text())
        assert manifest["options"]["strand"] == "unstranded"


class TestReadFilterFeedsKb:
    """DEF-01: ``--read-filter artefact`` runs before ``kb count`` and is audited."""

    def test_kb_counts_the_filtered_pair(self, read_filter_run) -> None:
        sample = _sample_dir(read_filter_run)
        assert (sample / "log" / "kb.done").is_file()
        kb_call = (sample / "kb-python" / "kb_info.json").read_text()
        assert "read_filtered/R2.fastq.gz" in kb_call, kb_call[:2000]

    def test_audit_balances_and_manifest_records_the_mode(self, read_filter_run) -> None:
        sample = _sample_dir(read_filter_run)
        with (sample / "read_filtered" / "read_filter_audit.tsv").open() as handle:
            audit = {
                r["category"]: int(float(r["fragments"]))
                for r in csv.DictReader(handle, delimiter="\t")
                if not r["category"].startswith("param:")
            }
        removed = sum(v for k, v in audit.items() if k.startswith("removed_"))
        assert audit["input"] == audit["retained"] + removed > 0
        manifest = json.loads((read_filter_run / "run_manifest.json").read_text())
        assert manifest["options"]["read_filter"] == "artefact"

    def test_validate_run_accepts_a_filtered_run(self, read_filter_run, tmp_path) -> None:
        report_path = tmp_path / "validate.json"
        _viralscan(
            "validate-run",
            str(read_filter_run),
            "--json-output",
            str(report_path),
            cwd=Path(__file__).parents[2],
        )
        report = json.loads(report_path.read_text(encoding="utf-8"))
        assert report["ok"] is True, report["issues"]


# --------------------------------------------------------------------------- #
# SW-10: the whole documented sequence on one fixture
# --------------------------------------------------------------------------- #
FULL_TOOLS = REQUIRED_TOOLS + [
    "STAR",
    "blastn",
    "makeblastdb",
    "minimap2",
    "samtools",
    "Rscript",
    "cd-hit-est",
]
#: The artifacts ``test_exact_lineage.py`` asserts for the faked run layout.
EVIDENCE_ARTIFACTS = (
    "read_lineage.tsv.gz",
    "evidence_manifest.json",
    "competitive_reads.raw.bam",
    "competitive_reads.raw.bam.bai",
    "competitive_reads.umi_dedup.bam",
    "coverage.raw.tsv",
    "coverage.deduplicated.tsv",
    "alignment_qc.tsv",
    "per_cell_alignment_qc.tsv",
    "coverage.raw_vs_deduplicated.png",
    "interpretation_flags.tsv",
    "blast_identity.tsv",
    "blast_sampling.json",
    "read_start_profile.tsv",
    "viralscan_evidence.igv.xml",
)


def _validate(run: Path, report_path: Path, repo: Path) -> dict:
    _viralscan("validate-run", str(run), "--json-output", str(report_path), cwd=repo)
    return json.loads(report_path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def documented_sequence(tmp_path_factory) -> dict:
    missing = have_tools(FULL_TOOLS)
    if missing:
        pytest.skip(f"Full-profile binaries not on PATH: {', '.join(missing)}")
    repo = Path(__file__).parents[2]
    work = tmp_path_factory.mktemp("sw10_full")
    steps: dict = {}

    # 1. preflight
    doctor = _viralscan("doctor", "--profile", "full", "--json", cwd=repo)
    steps["doctor"] = json.loads(doctor.stdout)

    # 2 + 3 + 4. reference build, visuals on, external cell calling: one command.
    out = work / "out"
    _viralscan(
        "-o", str(out),
        "-ref",
        "-fasta", str(FIXTURE / "reference.fasta"),
        "-gtf", str(FIXTURE / "combined.gtf"),
        "-s1", str(FIXTURE / "R1.fastq"),
        "-s2", str(FIXTURE / "R2.fastq"),
        "-x", "10xv3",
        "-c", "2",
        "--cell-calling", "external",
        "--called-cells-file", str(FIXTURE / "called_cells.txt"),
        "--yes",
        cwd=repo,
    )  # fmt: skip
    steps["out"] = out
    steps["sample"] = _sample_dir(out)

    # 6a. validate before evidence
    steps["validate_before"] = _validate(out, work / "validate_before.json", repo)

    # 5. evidence against the real run
    evidence = work / "evidence"
    _viralscan(
        "evidence",
        "--run-dir", str(steps["sample"]),
        "-o", str(evidence),
        "--virus", "VIRUS_TARGET",
        "--viral-fasta", str(FIXTURE / "viral.fasta"),
        "--host-fasta", str(FIXTURE / "host.fasta"),
        "--blast",
        "--read-start-profile",
        "--cell-tags",
        "--dedup", "umi",
        "--bin-size", "1",
        "--sampling-seed", "11",
        "-c", "1",
        cwd=repo,
    )  # fmt: skip
    steps["evidence"] = evidence

    # 6b. validate again after evidence
    steps["validate_after"] = _validate(out, work / "validate_after.json", repo)
    return steps


class TestDocumentedSequence:
    def test_doctor_full_profile_reports_every_tool(self, documented_sequence) -> None:
        tools = documented_sequence["doctor"]["tools"]
        for tool in FULL_TOOLS:
            assert tools.get(tool), f"doctor did not resolve {tool}"

    def test_reference_was_built_by_kb_ref(self, documented_sequence) -> None:
        index = documented_sequence["out"] / "index"
        for name in ("index.idx", "t2g.txt", "index.idx.build_manifest.json"):
            assert (index / name).is_file(), name
        manifest = json.loads((index / "index.idx.build_manifest.json").read_text())
        assert manifest["provenance"]["builder"] == "viralscan --reference"
        assert manifest["viral_gene_ids"] == ["VIRUS_TARGET"]
        assert manifest["host_gene_ids"] == ["HOST_GENE"]

    def test_visual_run_publishes_plots_and_report(self, documented_sequence) -> None:
        sample = documented_sequence["sample"]
        assert sorted((sample / "plots").glob("*.png")), "no plots without --no-visual"
        assert (sample / "report.html").is_file()
        assert (sample / "log" / "detection.done").is_file()

    def test_external_cell_calling_is_recorded_and_used(self, documented_sequence) -> None:
        manifest = json.loads(
            (documented_sequence["out"] / "run_manifest.json").read_text(encoding="utf-8")
        )
        assert manifest["cell_calling"]["method"] == "external"
        path = documented_sequence["sample"] / "results" / "per_cell_viral.tsv"
        with open(path, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        called = [r["barcode"] for r in rows if r["is_called_cell"] == "True"]
        assert called == ["AAACCCAAGAAACACT"]

    @pytest.mark.parametrize("relative", EVIDENCE_ARTIFACTS)
    def test_evidence_on_the_real_run_publishes_every_artifact(
        self, documented_sequence, relative
    ) -> None:
        assert (documented_sequence["evidence"] / relative).is_file()

    def test_evidence_lineage_names_only_the_viral_fragment(self, documented_sequence) -> None:
        import gzip

        with gzip.open(documented_sequence["evidence"] / "read_lineage.tsv.gz", "rt") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [r["read_id"] for r in rows] == ["viral_read"]

    def test_validate_run_is_clean_before_and_after_evidence(self, documented_sequence) -> None:
        for key in ("validate_before", "validate_after"):
            report = documented_sequence[key]
            assert report["ok"] is True, (key, report["issues"])
            assert report["issues"] == []
