# Installation

## Requirements

- Python ≥ 3.9
- [kb-python](https://www.kallistobus.tools/) (`kb` command)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0

## Option 1 — Conda (recommended)

A ready-to-use conda environment file is provided:

```bash
conda env create -f environment.yml
conda activate viralscan
```

This installs Python, `kb-python`, `snakemake`, and all Python dependencies
in one step.

## Option 2 — pip

```bash
pip install ViralScan
```

You still need `kb-python` and `snakemake` on your PATH.  Install them via
conda:

```bash
conda install -c bioconda -c conda-forge kb-python snakemake
```

## Option 3 — Container (Docker / Singularity)

Pre-built containers bundle every dependency including external tools.

### Docker

```bash
docker build -t viralscan:2.2.0 .
docker run --rm -it -v "$PWD:/data" viralscan:2.2.0 --help
```

### Singularity / Apptainer (HPC)

```bash
singularity build viralscan_2.2.0.sif Singularity.def
singularity exec viralscan_2.2.0.sif viralscan --help
```

## Verify the installation

```bash
viralscan --help
```

You should see the ViralScan banner and a list of all available options.

## Download the bundled viral reference panel

ViralScan's default viral GTF annotation panel is distributed separately on
Zenodo to keep the Python package small. Fetch it once after installation:

```bash
viralscan data fetch
```

This downloads the archive for DOI `10.5281/zenodo.20112332`, verifies the
Zenodo checksum, and unpacks the GTF files under
`~/.cache/viralscan/data/`. Custom `-gtf` workflows do not require this cache.
