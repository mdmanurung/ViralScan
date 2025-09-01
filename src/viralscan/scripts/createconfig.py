import os
import shutil
import yaml

# Outputs are in snakemake.output (list or dict)
output_files = snakemake.output
log_done, config_yaml = output_files

# Config dictionary from Snakemake
cfg_in = snakemake.config

# Clean up existing output
if os.path.exists(cfg_in["output"]):
    for filename in os.listdir(cfg_in["output"]):
        file_path = os.path.join(cfg_in["output"], filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

# Extra cleanup in output/log
directory = f"{cfg_in['output']}log"
file_to_keep = "create_config.done"

# Ensure directory exists, create if it doesn't
os.makedirs(directory, exist_ok=True)

if os.path.isdir(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and filename != file_to_keep:
            try:
                os.remove(file_path)
            except Exception:
                pass

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
    "umap": cfg_in["umap"]
}

with open(config_yaml, "w") as out:
    yaml.dump(cfg, out)

# Make sure log directory exist
os.makedirs(f"{cfg_in['output']}/log/", exist_ok=True)

# Touch the done file
with open(log_done, "w") as f:
    pass
