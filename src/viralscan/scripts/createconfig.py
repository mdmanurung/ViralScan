"""
This script creates the config (yaml) file which is used throughout the workflow
to obtain information about the users' input
"""

# importing packages
import os
import yaml

from viralscan.defaults import DEFAULTS
from viralscan.utils import setup_script_logging

log = setup_script_logging()

# Read parameters from Snakefile
log_done, config_yaml = snakemake.output
cfg_in = snakemake.config

# Ensure log directory exists
os.makedirs(f"{cfg_in['output']}log", exist_ok=True)

# Write new config.yaml
_detection_threshold = int(cfg_in.get("detection_threshold", 1))
if _detection_threshold < 1:
    raise ValueError(
        f"detection_threshold must be >= 1, got {_detection_threshold}. "
        "A threshold of 0 or below would flag every viral accession as detected."
    )
_multimap_pseudocount = float(
    cfg_in.get("multimap_pseudocount", DEFAULTS["multimap_pseudocount"])
)
if _multimap_pseudocount <= 0:
    raise ValueError(
        f"multimap_pseudocount must be > 0, got {_multimap_pseudocount}."
    )

cfg = {
    **DEFAULTS,
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
    "se_threshold": int(cfg_in.get("se_threshold", DEFAULTS["se_threshold"])),
    "detection_threshold": _detection_threshold,
    "min_counts": int(cfg_in.get("min_counts", DEFAULTS["min_counts"])),
    "min_genes": int(cfg_in.get("min_genes", DEFAULTS["min_genes"])),
    "hvg_min_mean": float(cfg_in.get("hvg_min_mean", DEFAULTS["hvg_min_mean"])),
    "hvg_max_mean": float(cfg_in.get("hvg_max_mean", DEFAULTS["hvg_max_mean"])),
    "hvg_min_disp": float(cfg_in.get("hvg_min_disp", DEFAULTS["hvg_min_disp"])),
    "umap_n_neighbors": int(cfg_in.get("umap_n_neighbors", DEFAULTS["umap_n_neighbors"])),
    "multimap_method": cfg_in.get("multimap_method", DEFAULTS["multimap_method"]),
    "multimap_pseudocount": _multimap_pseudocount,
    "multimap_primary_call": cfg_in.get(
        "multimap_primary_call", DEFAULTS["multimap_primary_call"]
    ),
    "cell_types": cfg_in.get("cell_types") or None,
}

# Host pre-subtraction (optional)
host_index = cfg_in.get("host_index") or None
host_filter_aligner = cfg_in.get("host_filter_aligner") or None

# Precompute the FASTQ paths that kb_count will actually consume.
# When host subtraction is active the host_filter rule writes filtered FASTQs
# to a predictable location; we record that path now so the Snakefile shell
# block never needs a conditional expression.
if host_index:
    _out = cfg_in["output"]
    kb_r1 = os.path.join(_out, "host_filtered", "R1.fastq.gz")
    kb_r2 = os.path.join(_out, "host_filtered", "R2.fastq.gz")
else:
    kb_r1 = cfg_in["sample1"]
    kb_r2 = cfg_in["sample2"]

cfg["host_index"] = host_index
cfg["host_filter_aligner"] = host_filter_aligner
cfg["kb_r1"] = kb_r1
cfg["kb_r2"] = kb_r2

with open(config_yaml, "w") as out:
    yaml.dump(cfg, out)

log.info("Creating the config is done!")

# Touch the done file
with open(log_done, "w") as f:
    pass
