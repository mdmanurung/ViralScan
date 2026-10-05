"""DAG-structure tests for the ViralScan Snakefile.

The pure-Python ``TestRuleOrdering`` class runs in the default suite without
needing the snakemake binary.  ``TestHostFilterDag`` runs ``snakemake -n``
and is marked ``@pytest.mark.integration``.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

from viralscan.evidence import have_tools

SNAKEFILE = Path(__file__).parent.parent / "src" / "viralscan" / "Snakefile"

# Minimal config values needed to keep the Snakefile shell-block template
# from raising KeyError during a dry-run (kb_count's shell references these).
_BASE_CONFIG = [
    "output=/tmp/vs_dag_test/",
    "index=/fake/index.idx",
    "transcripts=/fake/t2g.txt",
    "technology=10xv3",
    "whitelist=",
    "host_index=",
    "host_filter_aligner=",
    "cores=2",
]


class TestRuleOrdering:
    """rule all must be the first rule in the Snakefile.

    Snakemake uses the first rule it sees as the default target.  When
    ``host_filter`` was defined before ``rule all`` (conditional on
    ``config.get("host_index")``), running with ``--host-filter`` caused
    Snakemake to use ``host_filter`` as the default target — stopping the
    pipeline after host_filter.done with exit 0, no viral_summary.tsv.
    This test locks the fix in place.
    """

    def test_rule_all_is_first_rule(self) -> None:
        text = SNAKEFILE.read_text()
        matches = list(re.finditer(r"^rule\s+(\w+)\s*:", text, re.MULTILINE))
        assert matches, "No rules found in Snakefile"
        first = matches[0].group(1)
        assert first == "all", (
            f"First rule is '{first}', expected 'all'. "
            "When host_filter is defined first it becomes the default Snakemake "
            "target, halting the pipeline before kb_count (PLAN S2)."
        )

    def test_kb_count_shell_never_templates_raw_whitelist_config(self) -> None:
        """An unset whitelist must not reach kb as the literal ``-w None``.

        Snakemake parses ``--config whitelist=`` as None, and ``{config[whitelist]:q}``
        renders that as ``None``; kb reads ``-w None`` as "bypass barcode error
        correction". The whitelist has to be resolved in Python (a rule param).
        """
        assert "{config[whitelist]" not in SNAKEFILE.read_text()

    def test_kb_count_shell_never_templates_raw_strand_config(self) -> None:
        """An unset strand must not reach kb as the literal ``--strand None``."""
        assert "{config[strand]" not in SNAKEFILE.read_text()

    def test_kb_count_lists_filtered_fastqs_as_inputs_when_host_filter_set(self) -> None:
        """_kb_count_inputs must yield R1/R2 filter-FASTQ paths — not only the .done file."""
        text = SNAKEFILE.read_text()
        # The function body must reference both filtered FASTQ paths explicitly.
        assert "host_filtered/R1.fastq.gz" in text
        assert "host_filtered/R2.fastq.gz" in text


@pytest.mark.integration
class TestHostFilterDag:
    """Snakemake dry-run DAG tests.

    These require the ``snakemake`` binary on PATH and are excluded from the
    default pytest run (use ``-m integration`` to include them).
    """

    def _dryrun(self, extra_config: list[str]) -> str:
        """Return combined stdout+stderr of ``snakemake -n``.

        Skips when the binary is absent. Every other integration test here skips
        on a missing tool; these two failed instead, which made the required WP1
        gate red on any machine without snakemake and so made a green gate
        unavailable rather than merely inconvenient.
        """
        missing = have_tools(["snakemake"])
        if missing:
            pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

        cmd = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--dryrun",
            "--quiet",
            "all",
            "--config",
            *_BASE_CONFIG,
            *extra_config,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.stdout + result.stderr

    def test_no_host_filter_plans_full_pipeline(self) -> None:
        """Without host_index all core rules must appear in the dry-run plan."""
        output = self._dryrun(
            [
                "sample1=/fake/R1.fastq.gz",
                "sample2=/fake/R2.fastq.gz",
                "kb_r1=/fake/R1.fastq.gz",
                "kb_r2=/fake/R2.fastq.gz",
            ]
        )
        for rule in ("kb_count", "analysis", "multimap", "detection", "umap"):
            assert rule in output, (
                f"Rule '{rule}' missing from dry-run plan (no host filter). Full output:\n{output}"
            )

    def _kb_shell(self, extra_config: list[str]) -> str:
        """Return the rendered kb_count shell command from ``snakemake -n -p``."""
        missing = have_tools(["snakemake"])
        if missing:
            pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
        cmd = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "all",
            "--dryrun",
            "--printshellcmds",
            "--forcerun",
            "kb_count",
            "--config",
            *_BASE_CONFIG,
            "sample1=/fake/R1.fastq.gz",
            "sample2=/fake/R2.fastq.gz",
            "kb_r1=/fake/R1.fastq.gz",
            "kb_r2=/fake/R2.fastq.gz",
            *extra_config,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.stdout + result.stderr

    def test_unset_whitelist_renders_empty_not_none(self) -> None:
        """No -w given → WL is empty, so kb uses its packaged on-list and corrects."""
        output = self._kb_shell([])
        # Snakemake's `:q` renders "" as nothing, i.e. a bare `WL=` assignment.
        assert re.search(r"^\s*WL=\s*$", output, re.MULTILINE), output
        assert "WL=None" not in output, output

    def test_dropseq_without_w_bypasses_correction(self) -> None:
        """SW-21: no official on-list, so kb gets -w None instead of a knee allowlist."""
        output = self._kb_shell(["technology=dropseq"])
        assert re.search(r"^\s*WL=None\s*$", output, re.MULTILINE), output

    def test_set_whitelist_is_passed_through(self) -> None:
        output = self._kb_shell(["whitelist=/fake/wl.txt"])
        assert re.search(r"^\s*WL=/fake/wl\.txt\s*$", output, re.MULTILINE), output

    def test_set_strand_is_passed_to_both_kb_branches(self) -> None:
        output = self._kb_shell(["strand=reverse", "whitelist=/fake/wl.txt"])
        assert re.search(r"^\s*ST=reverse\s*$", output, re.MULTILINE), output
        # Both the -w and the no -w kb invocation carry the flag, guarded on ST.
        flagged = [ln for ln in output.splitlines() if ln.lstrip().startswith("kb count")]
        assert len(flagged) == 2, output
        assert all('${ST:+--strand "$ST"}' in ln for ln in flagged), flagged

    def test_unset_strand_renders_empty_not_none(self) -> None:
        output = self._kb_shell([])
        assert re.search(r"^\s*ST=\s*$", output, re.MULTILINE), output
        assert "ST=None" not in output, output

    def test_host_filter_plans_full_pipeline(self) -> None:
        """With host_index ALL rules through umap must appear — PLAN S2 regression guard.

        Before the fix, only create_config + host_filter appeared and the run
        exited 0 without producing viral_summary.tsv.
        """
        output = self._dryrun(
            [
                "sample1=/fake/R1.fastq.gz",
                "sample2=/fake/R2.fastq.gz",
                "host_index=/fake/host.idx",
                "host_filter_aligner=kallisto",
                "kb_r1=/tmp/vs_dag_test/host_filtered/R1.fastq.gz",
                "kb_r2=/tmp/vs_dag_test/host_filtered/R2.fastq.gz",
            ]
        )
        for rule in ("host_filter", "kb_count", "analysis", "multimap", "detection", "umap"):
            assert rule in output, (
                f"Rule '{rule}' missing from host-filter dry-run plan (PLAN S2). "
                f"Full output:\n{output}"
            )
