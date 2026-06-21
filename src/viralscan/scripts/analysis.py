"""
The analysis script creates a text file containing all the (viral) gene IDs.
It also checks whether the user has created the index itself, and if so, it
adds the gene IDs as well.

The module is importable without Snakemake: the magic-global wiring runs only
under the ``if "snakemake" in globals()`` guard at the bottom, and all logic is
reachable through :func:`run` / :func:`obtain_gtf` for direct testing.
"""

# Importing packages
import re
from pathlib import Path
from typing import Any, Iterable

from viralscan.data_fetch import ViralScanDataError, ensure_viral_data
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging, split_comma_paths

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


def obtain_gtf(config: RunConfig) -> set[str]:
    """
    This function obtains the GTF files out of the data directory of the package
    and checks whether GTFs (by using kb ref) have been added by the user.
    ---------------------------------------------------------------------
    Returns:
        viral_accessions (set): a set of (viral) gene IDs
    """
    viral_accessions: set[str] = set()

    custom_gtf_paths = _custom_gtf_paths(config.gtf)
    gtf_files = []
    try:
        data_dir = ensure_viral_data(config.data_cache_dir)
        gtf_files = list(data_dir.glob("*.gtf"))
    except ViralScanDataError as exc:
        if not custom_gtf_paths:
            raise RuntimeError(str(exc)) from exc
        log.warning("Bundled viral reference panel is unavailable: %s", exc)

    # Serratus viruses (bundled panel)
    for file in gtf_files:
        with open(file, "r") as f:
            viral_accessions |= extract_gene_ids(f)

    # check if GTF has been added by user. If so, add them to the viral list
    for file in custom_gtf_paths:
        gtf_path = Path(file)
        if not gtf_path.exists():
            raise FileNotFoundError(f"Custom GTF path does not exist: {gtf_path}")
        with open(gtf_path, "r") as f:
            viral_accessions |= extract_gene_ids(f)

    # write list to file
    with open(f"{config.output}log/analysis.txt", "w") as f:
        for v in viral_accessions:
            f.write(v + "\n")
    return viral_accessions


def run(ctx: RunContext) -> set[str]:
    """Entry point: obtain viral accessions for one Run."""
    accessions = obtain_gtf(ctx.config)
    log.info("Analysis is done!")
    return accessions


if "snakemake" in globals():
    run(RunContext.from_yaml(snakemake.params.configfile))  # noqa: F821 (snakemake magic global)
