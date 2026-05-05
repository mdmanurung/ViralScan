"""
This script creates the config (yaml) file which is used throughout the workflow
to obtain information about the users' input
"""

# importing packages
import os
import yaml

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
    "gtf": cfg_in["gtf"],
    "fasta": cfg_in["fasta"],
    "visual": cfg_in["visual"],
    "f1": cfg_in["f1"],
    "reference": cfg_in["reference"],
    "umap": cfg_in["umap"],
    "technology": cfg_in["technology"],
    "whitelist": cfg_in["whitelist"],
    "multimapping": cfg_in["multimapping"],
}

with open(config_yaml, "w") as out:
    yaml.dump(cfg, out)

print("\033[32mCreating the config is done!\033[0m")

# Touch the done file
with open(log_done, "w") as f:
    pass
