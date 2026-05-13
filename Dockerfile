# ViralScan Docker image
#
# Build:
#   docker build -t viralscan:2.3.0 .
#
# Run (interactive):
#   docker run --rm -it \
#     -v "$PWD/data:/data" \
#     viralscan:2.3.0 \
#     viralscan -t /data/t2g.txt -i /data/index.idx \
#               -o /data/output/ \
#               -s1 /data/R1.fastq.gz -s2 /data/R2.fastq.gz
#
# The image uses Miniforge (conda-forge) so kb-python and snakemake are
# installed from bioconda without needing a separate base image.

FROM condaforge/miniforge3:24.3.0-0

LABEL maintainer="emma.vonk@hotmail.nl" \
      version="2.3.0" \
      description="ViralScan — viral load quantification from single-cell RNA-seq"

# --------------------------------------------------------------------------
# System dependencies
# --------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------------------------------
# Conda environment
# --------------------------------------------------------------------------
COPY environment.yml /tmp/environment.yml

RUN mamba env create -f /tmp/environment.yml \
    && mamba clean --all -f -y

# Activate the viralscan conda environment for all subsequent RUN / CMD steps.
SHELL ["conda", "run", "-n", "viralscan", "/bin/bash", "-c"]

# Install the local checkout into the environment.
COPY . /opt/ViralScan
WORKDIR /opt/ViralScan
RUN python -m pip install .

# Verify the installation.
RUN viralscan --help

# --------------------------------------------------------------------------
# Default working directory and entry point
# --------------------------------------------------------------------------
WORKDIR /data

# Use conda run so the environment is active in non-interactive containers.
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "viralscan", "viralscan"]
CMD ["--help"]
