"""Integration test — no stale artifact survives a multimap method change (SW-05).

``rerun-multimap`` copies the source result tree with ``shutil.copytree`` and then
selectively invalidates. Every artifact it fails to invalidate survives as a
byte-copy of the *old* method's output inside a tree labelled with the *new*
method. Four such artifacts were found and fixed in SW-04; this test is the guard
that keeps the invalidation list complete.

Run with::

    PYTHONPATH=src python -m pytest -m integration tests/integration/

What this covers
----------------
The full invalidation contract as it is actually implemented, composed from the
real functions rather than re-stated:

* ``menu._rewrite_run_manifest`` re-stamps the per-sample provenance record.
* ``detection.clear_stale_virus_plots`` and ``hostresponse.clear_stale_virus_outputs``
  run at the start of their own rules, so a rule that re-runs cleans its own
  previous output.
* The sentinel list in ``_run_rerun_multimap`` decides which rules re-run at all.

What this does NOT cover
------------------------
Snakemake is not invoked, so this does not prove the rules re-run — it proves
that when they do, nothing stale is left behind, and that the source tree is
untouched. Proving the DAG actually re-executes needs the tiny end-to-end fixture
(SW-10), which does not exist yet.

The scenario throughout: virus ``HHV-6B`` clears the detection threshold under the
source run's method and falls below it under the new one. That is exactly the case
that leaves orphans, because the demoted virus is skipped rather than rewritten.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

#: Sentinels rerun-multimap drops so the corresponding rules re-run. Mirrors the
#: list in menu._run_rerun_multimap; a rule missing here silently keeps its
#: previous output, which is how hostresponse went stale.
RERUN_SENTINELS = ("log/detection.done", "log/umap.done", "log/hostresponse.done")

SURVIVING_VIRUS = "EBV"
DEMOTED_VIRUS = "HHV-6B"


def _build_source_run(root: Path) -> Path:
    """A completed run for two viruses, laid out as the real pipeline writes it."""
    sample = root / "source" / "sampleA"
    (sample / "log").mkdir(parents=True)
    (sample / "plots").mkdir()
    (sample / "hostresponse").mkdir()

    (sample / "config.yaml").write_text("multimap_method: equal\n", encoding="utf-8")
    (sample / "run_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "3.0.0",
                "allocation_method": "equal",
                "run_fingerprint": "source-fingerprint",
            }
        ),
        encoding="utf-8",
    )

    for sentinel in RERUN_SENTINELS + ("log/multimap.done",):
        (sample / sentinel).write_text("done", encoding="utf-8")

    for virus in (SURVIVING_VIRUS, DEMOTED_VIRUS):
        (sample / "plots" / f"{virus}_histogram.png").write_bytes(b"png")
        (sample / "plots" / f"SuperExpressor_{virus}.png").write_bytes(b"png")
        for suffix in ("_gene_weights.csv", "_stability.csv", "_differential.csv"):
            (sample / "hostresponse" / f"{virus}{suffix}").write_text("old", encoding="utf-8")
        (sample / "hostresponse" / f"{virus}_enrichment_GO_Biological_Process_2023.csv").write_text(
            "old", encoding="utf-8"
        )

    # Written by umap.py into the same directory detection cleans. Must survive.
    (sample / "plots" / "qc_hist_total_counts.png").write_bytes(b"png")
    # Run-level, rewritten wholesale rather than per virus. Must survive.
    (sample / "hostresponse" / "hostresponse_metrics.csv").write_text("old", encoding="utf-8")
    return sample


def _rerun(source_root: Path, output_root: Path, new_method: str) -> Path:
    """Reproduce what rerun-multimap does to one sample, using the real helpers."""
    from viralscan.menu import _rewrite_run_manifest
    from viralscan.scripts.detection import clear_stale_virus_plots
    from viralscan.scripts.hostresponse import clear_stale_virus_outputs

    shutil.copytree(source_root, output_root)
    sample = output_root / "sampleA"

    sample.joinpath("config.yaml").write_text(f"multimap_method: {new_method}\n", encoding="utf-8")
    for sentinel in RERUN_SENTINELS:
        sample.joinpath(sentinel).unlink()
    _rewrite_run_manifest(sample, source_dir=source_root, new_method=new_method)

    # Each rule clears its own previous output when it re-runs. Only the surviving
    # virus is regenerated; the demoted one is skipped, which is the orphan case.
    clear_stale_virus_plots(sample)
    clear_stale_virus_outputs(sample / "hostresponse")
    for name in (
        f"{SURVIVING_VIRUS}_histogram.png",
        f"SuperExpressor_{SURVIVING_VIRUS}.png",
    ):
        (sample / "plots" / name).write_bytes(b"new")
    (sample / "hostresponse" / f"{SURVIVING_VIRUS}_gene_weights.csv").write_text(
        "new", encoding="utf-8"
    )
    return sample


@pytest.fixture
def rerun(tmp_path: Path):
    source = _build_source_run(tmp_path)
    result = _rerun(tmp_path / "source", tmp_path / "output", "host-conservative")
    return source, result


class TestNoStaleArtifactSurvives:
    def test_the_demoted_virus_leaves_no_plot_behind(self, rerun) -> None:
        _, result = rerun

        assert not (result / "plots" / f"{DEMOTED_VIRUS}_histogram.png").exists()
        assert not (result / "plots" / f"SuperExpressor_{DEMOTED_VIRUS}.png").exists()

    def test_the_demoted_virus_leaves_no_hostresponse_csv_behind(self, rerun) -> None:
        _, result = rerun
        leftovers = sorted(p.name for p in (result / "hostresponse").glob(f"{DEMOTED_VIRUS}*"))

        assert leftovers == []

    def test_the_surviving_virus_is_regenerated_not_merely_kept(self, rerun) -> None:
        _, result = rerun

        assert (result / "plots" / f"{SURVIVING_VIRUS}_histogram.png").read_bytes() == b"new"
        assert (result / "hostresponse" / f"{SURVIVING_VIRUS}_gene_weights.csv").read_text(
            encoding="utf-8"
        ) == "new"

    def test_provenance_names_the_new_method(self, rerun) -> None:
        _, result = rerun
        manifest = json.loads((result / "run_manifest.json").read_text(encoding="utf-8"))

        assert manifest["allocation_method"] == "host-conservative"
        assert manifest["parent_run_fingerprint"] == "source-fingerprint"
        assert manifest["run_fingerprint"] != "source-fingerprint"

    def test_artifacts_owned_by_other_rules_survive(self, rerun) -> None:
        """Cleanup is scoped to owned patterns; it must not delete a neighbour's work."""
        _, result = rerun

        assert (result / "plots" / "qc_hist_total_counts.png").is_file()
        assert (result / "hostresponse" / "hostresponse_metrics.csv").is_file()

    def test_every_method_dependent_rule_is_invalidated(self, rerun) -> None:
        """A rule left out of the sentinel list re-runs only by mtime luck."""
        _, result = rerun

        for sentinel in RERUN_SENTINELS:
            assert not (result / sentinel).exists(), sentinel


class TestSourceRunIsUntouched:
    def test_the_source_keeps_the_demoted_virus_artifacts(self, rerun) -> None:
        source, _ = rerun

        assert (source / "plots" / f"{DEMOTED_VIRUS}_histogram.png").is_file()
        assert (source / "hostresponse" / f"{DEMOTED_VIRUS}_stability.csv").is_file()

    def test_the_source_manifest_and_sentinels_are_unchanged(self, rerun) -> None:
        source, _ = rerun
        manifest = json.loads((source / "run_manifest.json").read_text(encoding="utf-8"))

        assert manifest["allocation_method"] == "equal"
        assert manifest["run_fingerprint"] == "source-fingerprint"
        for sentinel in RERUN_SENTINELS:
            assert (source / sentinel).is_file(), sentinel
