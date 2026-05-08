"""
The analysis script creates a text file containing all the (viral) gene IDs.
It also checks whether the user has created the index itself, and if so, it
adds the gene IDs as well.
"""

# Importing packages
from pathlib import Path
import os
import logging

from viralscan.utils import load_config, setup_script_logging

log = setup_script_logging()

# Loading configfile of Snakefile
configfile = snakemake.params.configfile
config = load_config(configfile)


def obtain_gtf():
    """
    This function obtains the GTF files out of the data directory of the package
    and checks whether GTFs (by using kb ref) have been added by the user.
    ---------------------------------------------------------------------
    Returns:
        viral_accessions (set): a set of (viral) gene IDs
    """
    viral_accessions = set()

    # get the standard gtf files from this package
    project_root = Path(__file__).resolve().parent.parent

    # Point to folder and get the gtf files
    data_dir = project_root / "data"
    gtf_files = list(data_dir.glob("*.gtf"))

    # Serratus viruses
    for file in gtf_files:
        with open(file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    info = line.split("\t")[8]
                    gene_id = info.split('"')[1]
                    viral_accessions.add(gene_id)

    # check if GTF has been added by user. If so, add them to the viral list
    if config.get("gtf"):
        if os.path.exists(config["gtf"]):
            gtf_files = config["gtf"].split(",")
            for file in gtf_files:
                with open(file, "r") as f:
                    for line in f:
                        if not line.startswith("#"):
                            info = line.split("\t")[8]
                            gene_id = info.split('"')[1]
                            viral_accessions.add(gene_id)

    # write list to file
    with open(f"{config['output']}log/analysis.txt", "w") as f:
        for v in viral_accessions:
            f.write(v + "\n")
    return viral_accessions


def main():
    obtain_gtf()
    log.info("Analysis is done!")


main()
