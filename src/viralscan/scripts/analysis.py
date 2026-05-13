"""
The analysis script creates a text file containing all the (viral) gene IDs.
It also checks whether the user has created the index itself, and if so, it
adds the gene IDs as well.
"""

# Importing packages
import os
import re
import logging
from pathlib import Path

from viralscan.data_fetch import ViralScanDataError, ensure_viral_data
from viralscan.utils import load_config, setup_script_logging, split_comma_paths

log = setup_script_logging()

# Loading configfile of Snakefile
configfile = snakemake.params.configfile
config = load_config(configfile)


def _custom_gtf_paths(value):
    """Return configured custom GTF paths, ignoring unset Snakemake sentinels."""
    if value is None:
        return []
    return [path for path in split_comma_paths(str(value)) if path.lower() not in {"none", "null"}]


def obtain_gtf():
    """
    This function obtains the GTF files out of the data directory of the package
    and checks whether GTFs (by using kb ref) have been added by the user.
    ---------------------------------------------------------------------
    Returns:
        viral_accessions (set): a set of (viral) gene IDs
    """
    viral_accessions = set()

    custom_gtf_paths = _custom_gtf_paths(config.get("gtf"))
    gtf_files = []
    try:
        data_dir = ensure_viral_data()
        gtf_files = list(data_dir.glob("*.gtf"))
    except ViralScanDataError as exc:
        if not custom_gtf_paths:
            raise RuntimeError(str(exc)) from exc
        log.warning("Bundled viral reference panel is unavailable: %s", exc)

    # Serratus viruses
    for file in gtf_files:
        with open(file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    cols = line.split("\t")
                    if len(cols) < 9:
                        continue
                    info = cols[8]
                    m = re.search(r'gene_id "([^"]+)"', info)
                    if m:
                        viral_accessions.add(m.group(1))

    # check if GTF has been added by user. If so, add them to the viral list
    for file in custom_gtf_paths:
        gtf_path = Path(file)
        if not gtf_path.exists():
            raise FileNotFoundError(f"Custom GTF path does not exist: {gtf_path}")
        with open(gtf_path, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    cols = line.split("\t")
                    if len(cols) < 9:
                        continue
                    info = cols[8]
                    m = re.search(r'gene_id "([^"]+)"', info)
                    if m:
                        viral_accessions.add(m.group(1))

    # write list to file
    with open(f"{config['output']}log/analysis.txt", "w") as f:
        for v in viral_accessions:
            f.write(v + "\n")
    return viral_accessions


def main():
    obtain_gtf()
    log.info("Analysis is done!")


main()
