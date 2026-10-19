# ViralScan

**ViralScan** is a Snakemake-driven Python bioinformatics CLI that quantifies viral
load from paired-end FASTQ samples using
[kb-python](https://www.kallistobus.tools/) (kallisto + bustools). It is developed at
the Leiden University Medical Centre (LUMC).

[![CI](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml/badge.svg)](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mdmanurung/ViralScan/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mdmanurung/ViralScan)
[![PyPI](https://img.shields.io/pypi/v/ViralScan)](https://pypi.org/project/ViralScan/)

---

## Introduction

ViralScan pseudoaligns reads with kallisto | bustools against a viral or combined
host + virus reference, then reports per-cell and per-virus molecule estimates, an
HTML report and optional UMAP plots. It can run `kb ref` and `kb count` for you, so the
same tool covers reference building and scalable detection of viral load.

## Installation

```bash
conda env create -f environment.yml   # recommended: kb-python, snakemake, dependencies
conda activate viralscan
python -m pip install -e .
viralscan --help                       # all commands and options
```

See [installation](installation.md) for the pip tier, containers and the tested toolchain.

---

## Input Data

### Reference index

ViralScan needs an index created with `kb ref`. Either provide your own pre-built
index and transcript-to-gene map, or let ViralScan create one from FASTA + GTF files,
from NCBI accessions, or with `viralscan build-ref` (combined host + virus). The viral
GTF panel is distributed through Zenodo (`viralscan data fetch`).

### Samples

Paired-end FASTQ files, optionally gzipped. Both the forward (`-s1`) and the reverse
(`-s2`) file are required.

---

## User Guide

All output (counts, logs, plots, HTML report) is written under `-o/--output`, one
directory per FASTQ pair. There are two ways to run ViralScan.

**1. You already have a reference index**

```bash
viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
```

**2. ViralScan creates the index** (written to `output/index/`)

```bash
viralscan -o output/ -reference True -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -fasta fasta.fasta -gtf gtf.gtf
```

**Multiple samples or reference files** - comma-separated, no spaces:

```bash
viralscan -t transcripts.txt -i index.idx -o output/ -s1 a_1.fastq.gz,b_1.fastq.gz -s2 a_2.fastq.gz,b_2.fastq.gz
```

Stages (multimapping, UMAP, host filtering, evidence, gene programmes, host response)
are switched on or off with flags; see the [CLI reference](cli_reference.md) and the
[quickstart](quickstart.md). Output files are described in the
[output reference](output_reference.md); executable examples are in the
[vignettes](vignettes/README.md). Upgrading from an earlier release? See
[migration](migration.md). Project policies:
[SUPPORT.md](https://github.com/mdmanurung/ViralScan/blob/main/SUPPORT.md) (version and
support policy) and
[SECURITY.md](https://github.com/mdmanurung/ViralScan/blob/main/SECURITY.md)
(reporting vulnerabilities).

## License

MIT; see the LICENSE file.

```{toctree}
:maxdepth: 2
:caption: Home

installation
quickstart
```

```{toctree}
:maxdepth: 1
:caption: Tutorials

vignettes/README
vignettes/quickstart_fastq_to_viral_load
vignettes/building_a_reference
vignettes/multimapping_correction
vignettes/cell_calling_denominators
vignettes/cell_type_enrichment
vignettes/specificity_true_negative
vignettes/qc_and_read_evidence
vignettes/host_response_depth_control
```

```{toctree}
:maxdepth: 1
:caption: API

api
api_reference
```

```{toctree}
:maxdepth: 1
:caption: CLI and outputs

cli_reference
output_reference
```

```{toctree}
:maxdepth: 1
:caption: More

reference_panel
faq
workflow_comparison
migration
release_policy
reference_maintenance
v3_counting_contract
panel_integrity
changelog
```
