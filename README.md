# ViralScan

[![CI](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml/badge.svg)](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mdmanurung/ViralScan/branch/main/graph/badge.svg)](https://codecov.io/gh/mdmanurung/ViralScan)
[![PyPI](https://img.shields.io/pypi/v/ViralScan)](https://pypi.org/project/ViralScan/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue)](https://pypi.org/project/ViralScan/)

## Introduction

**ViralScan** is a Snakemake-driven Python CLI that quantifies viral load at
single-cell resolution from paired-end FASTQ data.  It uses
[kallisto | bustools](https://www.kallistobus.tools/) (via `kb-python`) to
pseudo-align reads against a viral or combined host+virus reference, then
produces per-cell and per-virus UMI counts, an interactive HTML report, and
optional UMAP visualisations. For human samples, the recommended workflow is a
combined host+virus reference so host and viral targets compete in one
quantification step. ViralScan is developed and maintained at the
[Leiden University Medical Centre (LUMC)](https://www.lumc.nl/).

---

## Features

- **195 pre-bundled viral reference annotations** (GTF) covering common
  human-infecting viruses, distributed separately through Zenodo
- **Host-aware reference building** with `viralscan build-ref` for combined
  host + virus kallisto indexes
- **Multimapping correction** uses host-conservative allocation by default,
  with legacy equal splitting still available by flag
- **Per-cell and per-virus summary tables** (`viral_summary.tsv`,
  `per_cell_viral.tsv`) plus a self-contained HTML report
- **Cell-type enrichment analysis** (`--cell-types`) — Fisher exact test with
  BH correction for each virus × cell-type combination
- **NCBI auto-fetch + local cache** (`--ncbi-accession`) — downloads FASTA +
  GTF via E-utilities and caches under `~/.cache/viralscan/ncbi/`
- **`viralscan build-ref`** — builds a combined host + virus kallisto index
  in one command
- **Docker and Singularity containers** provided for fully reproducible runs
- **Snakemake backend** — each step is a named rule with logged output

---

## Installation

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install -e .
viralscan data fetch
```

This installs all runtime dependencies including `kb-python` and `snakemake`,
then installs the local ViralScan checkout into the environment.

### pip

```bash
pip install ViralScan
viralscan data fetch
```

> **Note:** `kb-python` and `snakemake` must be installed separately via
> conda or another mechanism; pip does not guarantee the native binaries
> (`kb`, `snakemake`) are on `PATH`.

### Container

```bash
docker run --rm -v "$PWD":/data ghcr.io/mdmanurung/viralscan:latest \
    viralscan --help
```

A `Singularity.def` is also provided for HPC environments.

---

## User Guide

All output — counts, logs, plots, and the HTML report — is written under the
directory given by `--output`. Each FASTQ pair gets its own sample directory
named from the R1 filename before the first underscore; `sample_R1.fastq.gz`
writes to `output/sample/`.

### Mode 1 — Pre-built index

Use this when you already have a kallisto index and transcript-to-gene map
built with `kb ref`:

```bash
viralscan \
  --index index.idx \
  --transcripts t2g.txt \
  --output output/ \
  --sample1 sample_R1.fastq.gz \
  --sample2 sample_R2.fastq.gz
```

Multiple samples can be processed together by passing comma-separated paths:

```bash
viralscan \
  --index index.idx --transcripts t2g.txt \
  --output output/ \
  --sample1 s1_R1.fastq.gz,s2_R1.fastq.gz \
  --sample2 s1_R2.fastq.gz,s2_R2.fastq.gz
```

### Mode 2 — Build index from FASTA + GTF

Pass `--reference` with one or more FASTA and GTF files.  The index is
written to `output/index/`:

```bash
viralscan \
  --reference \
  --fasta viral.fasta \
  --gtf viral.gtf \
  --output output/ \
  --sample1 sample_R1.fastq.gz \
  --sample2 sample_R2.fastq.gz
```

Multiple FASTA/GTF files are accepted as comma-separated values.

### Mode 3 — NCBI accession auto-fetch

Provide one or more RefSeq accession numbers.  ViralScan downloads the FASTA
and GTF from NCBI, builds the index, and runs quantification in one step.
Supply `--ncbi-email` or set the `NCBI_EMAIL` environment variable:

```bash
export NCBI_EMAIL=you@example.org
viralscan \
  --ncbi-accession NC_002021.3 \
  --output output/ \
  --sample1 sample_R1.fastq.gz \
  --sample2 sample_R2.fastq.gz
```

Multiple accessions are comma-separated:

```bash
viralscan \
  --ncbi-accession NC_002021.3,NC_001512.1 \
  --output output/ \
  --sample1 sample_R1.fastq.gz \
  --sample2 sample_R2.fastq.gz
```

Downloads are cached under `~/.cache/viralscan/ncbi/`; re-runs do not
re-download.

### Viral annotation panel cache

The standard 195-virus GTF panel is hosted on Zenodo
(`10.5281/zenodo.20112332`) and must be fetched once before running workflows
that rely on the built-in panel:

```bash
viralscan data fetch
```

The files are cached under `~/.cache/viralscan/data/`. Workflows that pass
custom `-gtf` files still use those annotations directly.

For shared HPC caches, fetch to a shared root and pass the same root during
quantification:

```bash
viralscan data fetch --cache-dir /shared/viralscan-cache
viralscan --data-cache-dir /shared/viralscan-cache ...
```

You can also set `VIRALSCAN_CACHE=/shared/viralscan-cache` for both commands.

### Recommended — Build a combined host + virus reference

`viralscan build-ref` downloads a host transcriptome and one or more viral
sequences from NCBI, then builds a single kallisto index. This is the preferred
host-aware workflow because reads compete against host and viral targets in the
same alignment space:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_002021.3 NC_045512.2 \
  --output viralscan_ref/ \
  --ncbi-email you@example.org
```

Use `--list-species` to print all supported host species names.

Then quantify with the generated files. The default
`--multimap-method host-conservative` keeps host-virus ambiguous
equivalence-class mass out of primary viral counts:

```bash
viralscan \
  --index viralscan_ref/index.idx \
  --transcripts viralscan_ref/t2g.txt \
  --output output/ \
  --sample1 sample_R1.fastq.gz \
  --sample2 sample_R2.fastq.gz
```

Optional host pre-subtraction is still available with `--host-filter starsolo`
or `--host-filter kallisto`, but it is an advanced extra filter rather than a
required first step.

---

## Output

For an input named `sample_R1.fastq.gz`, key results are written under
`output/sample/`:

| File | Description |
|------|-------------|
| `results/viral_summary.tsv` | Per-virus totals: `total_umi`, `infected_cells`, `pct_infected`, `umi_per_10k` |
| `results/per_cell_viral.tsv` | Per-barcode × per-virus: `viral_umi`, `total_umi`, `viral_fraction`; UMI values may be fractional after multimapping correction |
| `results/multimap_evidence.tsv` | Unique, ambiguous, and host-virus ambiguous viral evidence |
| `report.html` | Self-contained interactive HTML report |
| `kb-python/counts_unfiltered/adata_multimap.h5ad` | AnnData with multimapping-corrected counts |
| `plots/` | PNG plots and optional UMAP HTML files |

See [docs/output_reference.md](docs/output_reference.md) for the full output
schema including optional files (`cell_type_enrichment.tsv`, UMAP plots, etc.).

---

## Citation

If you use ViralScan in your work, please cite it using the metadata in
[CITATION.cff](CITATION.cff).  ViralScan is developed at the Leiden University
Medical Centre (LUMC), Leiden, the Netherlands.

---

## License

This project is licensed under the MIT License.  See the [LICENSE](LICENSE)
file for details.
