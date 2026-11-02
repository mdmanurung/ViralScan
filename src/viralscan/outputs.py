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

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

from viralscan import __version__

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


#: Files the report lists. ``(path, what it holds, option)``: ``option`` is the run option that
#: writes the file (``None`` = every run). The report table and the README "written by the
#: options you enabled" block both read this, so a gated file shows up only when it exists.
RESULT_FILES = (
    ("summary.txt", "Plain-text summary (legacy, human-readable)", None),
    ("results/viral_summary.tsv", "Per-virus statistics table (TSV)", None),
    (
        "results/sensitivity.tsv",
        "Per-virus detection limit (LOD95) and negative-certifiability (TSV)",
        None,
    ),
    (
        "results/gene_program_summary.tsv",
        "Per-virus gene-programme state counts and breadth (TSV)",
        "gene_programs",
    ),
    (
        "results/gene_program_cells.tsv",
        "Per-cell gene-programme assignment (TSV)",
        "gene_programs",
    ),
    ("results/per_cell_viral.tsv", "Per-barcode viral UMI annotation (TSV)", None),
    ("results/multimap_evidence.tsv", "Unique and ambiguous viral evidence table (TSV)", None),
    ("plots/", "UMI bar plots and super-expressor scatter plots (PNG)", None),
    ("plots/umap_binary.html", "Interactive UMAP coloured by virus", "umap"),
    ("plots/umap_continuous.html", "Interactive UMAP coloured by viral load", "umap"),
    ("kb-python/", "Raw kallisto|bustools count matrices", None),
)


def result_files(is_on: Callable[[str], bool]) -> list[tuple[str, str]]:
    """``(path, text)`` of the ``RESULT_FILES`` this run writes; ``is_on(option)`` gates the rest."""
    return [(path, text) for path, text, opt in RESULT_FILES if opt is None or is_on(opt)]


def safe_to_delete() -> list[tuple[str, str]]:
    """What a user may remove from a sample after a finished run, from ``SAMPLE_LAYOUT``.

    Every ``intermediate`` work dir except ``kb-python/`` (it holds the matrix), plus the
    kb-python files around that matrix and the staging dirs an interrupted ``kb count`` leaves.
    Deleting any of it costs ``--resume`` reuse and the post-hoc verbs that re-read BUS or FASTQs.
    """
    rows = [
        (path, text)
        for path, role, text in SAMPLE_LAYOUT
        if role == "intermediate" and path != "kb-python/"
    ]
    rows.append(("kb-python/ except counts_unfiltered/", "BUS files and run_info.json"))
    rows.append((".kb-count.* and .kb-previous.*", "staging left by an interrupted kb count"))
    return rows


def _on(options: Mapping[str, Any] | None) -> Callable[[str], bool]:
    return lambda name: bool((options or {}).get(name))


def _status(finished: bool) -> list[str]:
    if finished:
        return ["Status: finished (run_complete.json is at the run root)"]
    return [
        "Status: NOT FINISHED. This file was written when the run started.",
        "  If run_complete.json is missing from the run root, the run failed or is still going;",
        "  check .snakemake/log/ in the directory you ran viralscan from, then re-run with",
        "  --resume. It is rewritten on completion.",
    ]


def sample_readme(
    sample: str,
    *,
    finished: bool = True,
    options: Mapping[str, Any] | None = None,
) -> str:
    """Text of ``<sample>/README.txt``, from ``SAMPLE_LAYOUT``, ``QUESTIONS`` and the options on.

    ``options`` is the run's option map (``run_manifest.json`` ``options``); it decides which
    gated files in ``RESULT_FILES`` are listed.
    """
    lines = [f"ViralScan {__version__}, sample {sample}, layout version {LAYOUT_VERSION}", ""]
    lines += _status(finished) + ["", "What is where"]
    lines += [f"  {path:<17} {text}" for path, _role, text in SAMPLE_LAYOUT]
    lines += ["", "Where do I find ..."]
    lines += [f"  {q:<48} {path}" for q, path in QUESTIONS]
    gated = [(p, t) for p, t, opt in RESULT_FILES if opt and _on(options)(opt)]
    if gated:
        lines += ["", "Also written because you enabled the matching option"]
        lines += [f"  {path:<36} {text}" for path, text in gated]
    lines += [
        "",
        "How to read a negative",
        "  A virus that is absent from results/viral_summary.tsv, or has detected=False, was",
        "  not detected. That does not show it is absent: results/sensitivity.tsv gives each",
        "  virus's detection limit (LOD95), and a zero below that limit says little.",
        "",
        "Which settings made this: config.yaml here; the whole run is in ../run_manifest.json.",
        "",
        "Safe to delete after a finished run (you lose --resume reuse and BUS-based post-hoc verbs)",
    ]
    lines += [f"  {path:<38} {text}" for path, text in safe_to_delete()]
    lines += [
        "  Keep kb-python/counts_unfiltered/, results/, plots/, config.yaml and the run root files.",
        "",
        "Bulky intermediates (kb-python/ and the other work dirs) are not in the inventory.",
    ]
    return "\n".join(lines) + "\n"


def run_readme(
    samples: list[str],
    *,
    finished: bool = True,
    manifest: Mapping[str, Any] | None = None,
) -> str:
    """Text of the run-level ``README.txt``; ``manifest`` adds the fingerprint and method."""
    manifest = manifest or {}
    lines = [f"ViralScan {__version__} run", ""] + _status(finished) + [""]
    if manifest.get("run_fingerprint"):
        lines += [f"  run fingerprint   {manifest['run_fingerprint']}"]
    if manifest.get("allocation_method"):
        lines += [f"  multimap method   {manifest['allocation_method']}"]
    lines += [
        "",
        "  run_manifest.json     what was run (settings, versions, inputs)",
        "  run_complete.json     marker: the run finished; hashes of the key outputs",
        "  combined_summary.tsv  every sample's results/viral_summary.tsv, stacked",
        "",
        "Samples (open report.html first; each folder has its own README.txt)",
    ]
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
