"""
The analysis script creates a text file containing all the (viral) gene IDs.
It also checks whether the user has created the index itself, and if so, it
adds the gene IDs as well.

It then builds the Run's Virus Identity table from the index t2g and writes it
to ``results/virus_identity.tsv`` (PLAN ``MECH-A``). ``log/analysis.txt`` is
kept, unchanged, while consumers move to the table.

The module is importable without Snakemake: the magic-global wiring runs only
under the ``if "snakemake" in globals()`` guard at the bottom, and all logic is
reachable through :func:`run` / :func:`obtain_gtf` for direct testing.
"""

# Importing packages
import json
import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from viralscan.anellovirus import candidate_gene_ids
from viralscan.data_fetch import ViralScanDataError, ensure_viral_data
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging, split_comma_paths
from viralscan.virus_identity import (
    manifest_path_for_index,
    manifest_viral_gene_ids,
    write_identity_table,
)

log = setup_script_logging()

_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def _custom_gtf_paths(value: Any) -> list[str]:
    """Return configured custom GTF paths, ignoring unset Snakemake sentinels."""
    if value is None:
        return []
    return [path for path in split_comma_paths(str(value)) if path.lower() not in {"none", "null"}]


def extract_gene_ids(lines: Iterable[str]) -> set[str]:
    """Extract ``gene_id`` values from GTF lines.

    Skips comment lines and rows with fewer than 9 tab-separated columns. Pure
    and importable — this is the parsing core the tests exercise directly.
    """
    accessions: set[str] = set()
    for line in lines:
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        m = _GENE_ID_RE.search(cols[8])
        if m:
            accessions.add(m.group(1))
    return accessions


def _bundled_panel_gene_ids(config: RunConfig, custom_gtf_paths: list[str]) -> set[str]:
    """Gene IDs of the bundled panel, for an index without a build manifest."""
    gtf_files: list[Path] = []
    try:
        data_dir = ensure_viral_data(config.data_cache_dir)
        gtf_files = list(data_dir.glob("*.gtf"))
    except ViralScanDataError as exc:
        if not custom_gtf_paths:
            raise RuntimeError(str(exc)) from exc
        # A custom GTF was supplied but the packaged panel is unavailable. Warn
        # loudly rather than only logging: continuing means the bundled panel's
        # viruses are absent from `viral_accessions` while the index still
        # contains them, so they become undetectable and the run looks clean for
        # the wrong reason.
        log.warning(
            "Bundled viral reference panel is unavailable (%s); continuing with the "
            "custom GTF(s) only. Any bundled-panel virus present in the index will "
            "NOT be reported as detected.",
            exc,
        )

    # Serratus viruses (bundled panel)
    accessions: set[str] = set()
    for file in gtf_files:
        with open(file) as f:
            accessions |= extract_gene_ids(f)
    return accessions


def obtain_gtf(config: RunConfig) -> set[str]:
    """
    This function collects the viral gene IDs of the index: from its build
    manifest (``<index>.build_manifest.json``, PLAN ``CAT-06``) when it has one,
    otherwise from the bundled panel GTFs. GTFs added by the user (``-gtf``) are
    unioned in either way.
    ---------------------------------------------------------------------
    Returns:
        viral_accessions (set): a set of (viral) gene IDs

    Anellovirus gene IDs are added from the packaged accession table rather than
    from a GTF file. The expanded anellovirus panel (2,041 accessions) is
    materialized into the *built index* by ``build-reference``, not into the
    packaged panel directory, so globbing the panel GTFs alone left 2,021 of
    those genomes — 99 % of the panel, and the whole human anellovirus sequence
    space — countable but invisible to detection.
    ``anellovirus.candidate_gene_ids`` derives the ``{accession}_geneN`` IDs the
    builder emits, so they are recognised however the reference was built. Set
    ``--no-anellovirus-gene-ids`` to restore the pre-v3 GTF-glob-only behaviour.
    """
    viral_accessions: set[str] = set()

    custom_gtf_paths = _custom_gtf_paths(config.gtf)
    manifest = manifest_path_for_index(config.index) if (config.index or "").strip() else None
    if manifest is not None and manifest.is_file():
        # The index records its own viral genes (PLAN CAT-06): no bundled panel
        # (and no Zenodo fetch) is needed, and a custom GTF is added, not dropped.
        viral_accessions |= manifest_viral_gene_ids(manifest)
        log.info("Viral gene IDs from the index build manifest %s.", manifest)
    else:
        viral_accessions |= _bundled_panel_gene_ids(config, custom_gtf_paths)

    # check if GTF has been added by user. If so, add them to the viral list
    for custom_file in custom_gtf_paths:
        gtf_path = Path(custom_file)
        if not gtf_path.exists():
            raise FileNotFoundError(f"Custom GTF path does not exist: {gtf_path}")
        with open(gtf_path) as f:
            viral_accessions |= extract_gene_ids(f)

    # Expanded anellovirus panel, derived from the authoritative accession table.
    n_before = len(viral_accessions)
    if getattr(config, "anellovirus_gene_ids", True):
        try:
            viral_accessions |= candidate_gene_ids()
        except Exception as exc:  # noqa: BLE001 - a missing table must not block a run
            log.warning(
                "Could not derive anellovirus gene IDs from the packaged accession "
                "table (%s). Anellovirus genomes present in the index will not be "
                "reported as detected. Pass --no-anellovirus-gene-ids to silence "
                "this warning once the gap is understood.",
                exc,
            )
    log.info(
        "Viral gene IDs: %d total (%d from manifest/GTF, %d anellovirus candidates added)",
        len(viral_accessions),
        n_before,
        len(viral_accessions) - n_before,
    )

    # write list to file
    with open(f"{config.output}log/analysis.txt", "w") as f:
        for v in sorted(viral_accessions):
            f.write(v + "\n")
    return viral_accessions


def _chemistry_sanity(ctx: RunContext) -> None:
    """Persist post-count diagnostics before the optional fail-closed gate."""
    from viralscan.chemistry_check import sanity_report
    from viralscan.run_safety import RUN_MANIFEST, record_manifest_block

    kb_dir = Path(ctx.config.output) / "kb-python"
    host_in_index = not ctx.config.host_filter_aligner
    report = sanity_report(kb_dir, host_in_index=host_in_index)
    results = Path(ctx.config.output) / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "chemistry_sanity.json").write_text(json.dumps(report, indent=2) + "\n")
    run_root = Path(ctx.config.output).resolve().parent
    if (run_root / RUN_MANIFEST).is_file():
        record_manifest_block(run_root, "chemistry_sanity", Path(ctx.config.output).name, report)
    for finding in report["findings"]:
        log.warning("chemistry sanity (%s): %s", finding["check"], finding["message"])
    if ctx.config.require_chemistry_sanity and (
        report["status"] != "checked" or any(f["level"] == "error" for f in report["findings"])
    ):
        raise RuntimeError("Library sanity gate failed; see results/chemistry_sanity.json.")


def run(ctx: RunContext) -> set[str]:
    """Entry point: obtain viral accessions and the Virus Identity table for one Run."""
    _chemistry_sanity(ctx)
    accessions = obtain_gtf(ctx.config)
    write_identity_table(ctx.config, accessions)
    log.info("Analysis is done!")
    return accessions


if "snakemake" in globals():
    run(RunContext.from_yaml(snakemake.params.configfile))  # noqa: F821 (snakemake magic global)
