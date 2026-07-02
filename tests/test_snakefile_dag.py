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
        """Return combined stdout+stderr of ``snakemake -n``."""
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
