"""
This script creates the config (yaml) file which is used throughout the workflow
to obtain information about the users' input
"""

# importing packages
import os
import logging
import yaml

from viralscan.utils import setup_script_logging

log = setup_script_logging()

# Read parameters from Snakefile
log_done, config_yaml = snakemake.output
cfg_in = snakemake.config

# Ensure log directory exists
os.makedirs(f"{cfg_in['output']}log", exist_ok=True)

# Write new config.yaml
cfg = {
    "output": cfg_in["output"],
    "index": cfg_in["index"],
    "transcripts": cfg_in["transcripts"],
    "sample1": cfg_in["sample1"],
    "sample2": cfg_in["sample2"],
    "overwrite": "yes",
    "gtf": cfg_in["gtf"] or None,
    "fasta": cfg_in["fasta"] or None,
    "visual": bool(cfg_in["visual"]),
    "f1": cfg_in["f1"] or None,
    "reference": bool(cfg_in["reference"]),
    "umap": bool(cfg_in["umap"]),
    "technology": cfg_in["technology"],
    "whitelist": cfg_in["whitelist"] or None,
    "multimapping": bool(cfg_in["multimapping"]),
    # Reporting / threshold parameters (PR 11 A4)
    "se_threshold": int(cfg_in.get("se_threshold", 10)),
    "detection_threshold": int(cfg_in.get("detection_threshold", 1)),
    "min_counts": int(cfg_in.get("min_counts", 1000)),
    "min_genes": int(cfg_in.get("min_genes", 200)),
    "cell_types": cfg_in.get("cell_types") or None,
}

with open(config_yaml, "w") as out:
    yaml.dump(cfg, out)

log.info("Creating the config is done!")

# Touch the done file
with open(log_done, "w") as f:
    pass
