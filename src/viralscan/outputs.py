"""Where ViralScan-only artifacts live under a sample's output directory.

evonk's ViralScan 2.2.0 wrote ``config.yaml``, ``summary.txt``, ``log/``,
``kb-python/`` and ``plots/``. Everything v3 added lives under ``results/`` (and
``read_filtered/``). This module names those paths once, like
:mod:`viralscan.kb_outputs` does for ``kb-python/``, and performs no I/O.

``LAYOUT_VERSION`` is stamped in ``run_manifest.json``; bump it when a published
path moves. Version 2 adds the reserved ``posthoc/`` subfolder, where post-hoc verbs
write by default. Version 3 (PLAN ``LAYOUT-V3-01``) moves the data files out of ``log/``
and ``plots/`` (``log/analysis.txt`` -> ``results/analysis_genes.txt``, ``log/found_genes.txt``
-> ``results/found_genes.txt``, ``plots/SuperExpressor_<virus>.{tsv,json}`` ->
``results/superexpressor/``) and gives post-hoc verbs one level, ``<sample>/posthoc/<verb>/``
(``rerun-multimap`` copies the whole run, so it stays at ``<run>/posthoc/``).
Readers accept every version in ``SUPPORTED_LAYOUT_VERSIONS``; v1 and v2 trees stay valid and
the ``*_path`` helpers below fall back to the old locations.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

LAYOUT_VERSION = 3
SUPPORTED_LAYOUT_VERSIONS = (1, 2, 3)
POSTHOC_DIRNAME = "posthoc"
ANALYSIS_GENES = "results/analysis_genes.txt"  # v1/v2: log/analysis.txt
FOUND_GENES = "results/found_genes.txt"  # v1/v2: log/found_genes.txt


def _new_or_legacy(output: Union[str, Path], new: str, legacy: str) -> Path:
    """The v3 path, or the v1/v2 path when only that one exists (the v3 path when neither)."""
    new_path, legacy_path = Path(output) / new, Path(output) / legacy
    return legacy_path if not new_path.is_file() and legacy_path.is_file() else new_path


def analysis_genes_path(output: Union[str, Path]) -> Path:
    """Viral gene IDs of the index: ``results/analysis_genes.txt`` (v1/v2: ``log/analysis.txt``)."""
    return _new_or_legacy(output, ANALYSIS_GENES, "log/analysis.txt")


def found_genes_path(output: Union[str, Path]) -> Path:
    """Detected genes with counts: ``results/found_genes.txt`` (v1/v2: ``log/found_genes.txt``)."""
    return _new_or_legacy(output, FOUND_GENES, "log/found_genes.txt")


#: Sample-directory registry: ONE definition of the v3 layout. ``(path, role, what it holds)``.
#: ``role`` is ``root`` (first stop), ``results`` (user-facing tables), ``plots``,
#: ``intermediate`` (bulky machine-only work dirs; the inventory skips them, see
#: ``artifact_inventory.py``), ``log`` or ``posthoc``. The README.txt of every sample, the
#: inventory's ignore list and the docs drift test all read this tuple.
#: ponytail: intermediates stay at the sample root (kb-python/ etc.) instead of
#: ``intermediate/``; moving them would touch ~200 path strings, the frozen oracles and
#: the kb staging shell block. Move them here (and in KbCountOutputs) if that is wanted later.
SAMPLE_LAYOUT = (
    ("report.html", "root", "HTML report: start here"),
    ("summary.txt", "root", "plain-text summary: headline, detected viruses, gate lines"),
    ("README.txt", "root", "this file"),
    ("config.yaml", "root", "the settings this sample was run with"),
    ("results/", "results", "user-facing tables (TSV/JSON); see the questions below"),
    ("plots/", "plots", "images and interactive HTML only"),
    ("posthoc/", "posthoc", "post-hoc verbs write here: posthoc/<verb>/ (never part of the run)"),
    ("log/", "log", "logs and .done sentinels only"),
    (
        "kb-python/",
        "intermediate",
        "kb count work dir; counts_unfiltered/adata_multimap.h5ad is the matrix to load",
    ),
    ("read_filtered/", "intermediate", "FASTQs after the read filter"),
    ("host_filtered/", "intermediate", "STARsolo host-filter work dir and filtered FASTQs"),
    ("star_hostvirus/", "intermediate", "STAR host+virus confirmation BAM"),
    ("host_capture/", "intermediate", "host-capture work dir (table: results/host_capture.tsv)"),
)

#: "Which file answers which question" lines of README.txt and the top of the docs.
QUESTIONS = (
    ("Is virus X present?", "results/viral_summary.tsv"),
    ("Which viral genes were detected, with counts?", "results/found_genes.txt"),
    ("Which genes in the index are viral?", "results/analysis_genes.txt"),
    ("Which cells carry virus, and how much?", "results/per_cell_viral.tsv"),
    ("Which barcodes were called as cells?", "results/called_cells.tsv"),
    ("Per-cell counts matrix (AnnData)?", "kb-python/counts_unfiltered/adata_multimap.h5ad"),
    ("Is the detection robust (unique vs multimapped)?", "results/multimap_evidence.tsv"),
    ("Which files does this run contain, with checksums?", "results/artifact_inventory.json"),
)


def sample_readme(sample: str) -> str:
    """Text of ``<sample>/README.txt``, from ``SAMPLE_LAYOUT`` and ``QUESTIONS``."""
    lines = [f"ViralScan sample {sample}: layout version {LAYOUT_VERSION}", "", "What is where"]
    lines += [f"  {path:<17} {text}" for path, _role, text in SAMPLE_LAYOUT]
    lines += ["", "Where do I find ..."]
    lines += [f"  {q:<48} {path}" for q, path in QUESTIONS]
    lines += [
        "",
        "Bulky intermediates (kb-python/ and the other work dirs) are not in the inventory.",
    ]
    return "\n".join(lines) + "\n"


def run_readme(samples: list[str]) -> str:
    """Text of the run-level ``README.txt``."""
    lines = [
        "ViralScan run",
        "",
        "  run_manifest.json     what was run (settings, versions, inputs)",
    ]
    lines += ["  run_complete.json     marker: the run finished; hashes of the key outputs"]
    lines += ["  combined_summary.tsv  every sample's results/viral_summary.tsv, stacked", ""]
    lines += ["Samples (open report.html first; each folder has its own README.txt)"]
    lines += [f"  {sample}/report.html" for sample in samples]
    return "\n".join(lines) + "\n"


def posthoc_dir(run_dir: Union[str, Path], *parts: str) -> Path:
    """``<run_dir>/posthoc/<parts>``: the only place post-hoc verbs write by default."""
    return Path(run_dir).joinpath(POSTHOC_DIRNAME, *parts)


def is_posthoc(path: Path, run_dir: Path) -> bool:
    """True when ``path`` is inside a ``posthoc/`` folder of ``run_dir`` or of one of its samples.

    Matches ``<run_dir>/posthoc/...`` (v2 and ``rerun-multimap``) and
    ``<run_dir>/<sample>/posthoc/...`` (v3). Pass the run root or a sample directory.
    """
    return POSTHOC_DIRNAME in path.relative_to(run_dir).parts[:2]


@dataclass(frozen=True)
class V3Outputs:
    """Named v3 paths beneath one sample's output directory."""

    output: Path

    @classmethod
    def from_config_output(cls, output: Union[str, Path]) -> V3Outputs:
        """Build from a ``config["output"]`` value (which may carry a trailing sep)."""
        return cls(Path(output))

    @property
    def results_dir(self) -> Path:
        return self.output / "results"

    @property
    def posthoc_dir(self) -> Path:
        return posthoc_dir(self.output)

    @property
    def read_filtered_dir(self) -> Path:
        return self.output / "read_filtered"

    @property
    def count_audit(self) -> Path:
        return self.results_dir / "count_audit.tsv"

    @property
    def analysis_genes(self) -> Path:
        return analysis_genes_path(self.output)

    @property
    def found_genes(self) -> Path:
        return found_genes_path(self.output)

    @property
    def superexpressor_dir(self) -> Path:
        return self.results_dir / "superexpressor"

    @property
    def hostresponse_dir(self) -> Path:
        return self.results_dir / "hostresponse"

    @property
    def artifact_inventory(self) -> Path:
        return self.results_dir / "artifact_inventory.json"
