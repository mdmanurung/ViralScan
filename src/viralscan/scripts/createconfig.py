"""Snakemake worker: write the per-Run ``config.yaml``.

This rule is the sole writer of ``config.yaml``. All coercion and validation
lives in :class:`viralscan.runconfig.RunConfig`; this script just wires the
Snakemake globals to it. See ``CONTEXT.md`` ("Run Config").
"""

import os

from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging

log = setup_script_logging()

# Read parameters from Snakefile
log_done, config_yaml = snakemake.output  # noqa: F821 (snakemake magic global)
cfg_in = snakemake.config  # noqa: F821

# Ensure log directory exists
os.makedirs(f"{cfg_in['output']}log", exist_ok=True)

# Build + validate once, then write the trusted YAML downstream rules read.
RunConfig.from_snakemake_config(cfg_in).to_yaml(config_yaml)
log.info("Creating the config is done!")

# Touch the done file
with open(log_done, "w") as f:
    pass
