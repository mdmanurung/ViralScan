# ViralScan

> **V3 development notice:** pre-v3 corrected/combined counts are scientifically
> incompatible with the molecule-counting contract under development and are not
> eligible for v3 quantitative claims. Rebuild from retained corrected BUS/reference
> inputs or FASTQ; do not numerically migrate old H5AD values.

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
produces per-cell and per-virus molecule estimates, an interactive HTML report, and
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
- **Molecule-aware multimapping** uses `host-conservative` by default so
  host–virus ambiguity does not create viral evidence; `equal`,
  `unique-weighted`, `em-global`, and hierarchical `em-cell` are explicit alternatives
- **Per-cell and per-virus summary tables** (`viral_summary.tsv`,
  `per_cell_viral.tsv`) plus a self-contained HTML report
- **Cell-type enrichment analysis** (`--cell-types`) — Fisher exact test with
  BH correction for each virus × cell-type combination
- **NCBI auto-fetch + local cache** (`--ncbi-accession`) — downloads FASTA +
  GTF via E-utilities and caches under `~/.cache/viralscan/ncbi/`
- **`viralscan build-ref`** — builds a combined host + virus kallisto index
  in one command
- **Host-response analysis** (`viralscan hostresponse` or `--host-h5ad`) —
  associates per-virus presence with host gene expression via L2 logistic
  regression + Lasso stability selection, with optional gget pathway enrichment
- **Docker and Singularity containers** provided for fully reproducible runs
- **Snakemake backend** — each step is a named rule with logged output

---

## Installation

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install .          # use "-e ." for a development checkout
viralscan data fetch
```

This installs all runtime dependencies including `kb-python` and `snakemake`
from conda first, then installs ViralScan into the environment. Because
`snakemake` is already satisfied by conda, pip does not rebuild it here.

### pip

```bash
pip install ViralScan
viralscan data fetch
```

> **Note:** `kb-python` and `snakemake` must be on `PATH` for ViralScan to run;
> pip does not provide the native binaries (`kb`, `snakemake`). Install them via
> conda/bioconda first. In particular, letting pip resolve `snakemake` from PyPI
> can fail while building the `connection_pool` transitive dependency on older
> `setuptools`; installing `snakemake` from conda (as in the Conda flow above)
> avoids this. The most reproducible option is the container below.

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
The curated default does not include the expanded anellovirus panel. Opt in
with `--anellovirus`; without full-host-genome competition, resulting hits are
screening candidates only. Each build writes a per-sequence
`reference_manifest.json` and fails closed on incomplete panels or duplicate
reference records.
Use `--genome-dlist GRCh38.fa` to mask host-genomic shared k-mers and record raw
per-virus host-homology metrics; this full-genome build is intentionally an
explicit, compute-intensive mode.

Then quantify with the generated files. The default
`--multimap-method host-conservative` keeps host-virus ambiguous
equivalence-class mass out of primary viral estimates (use `equal` for an
explicit equal-allocation comparison):

```bash
viralscan \
  --index viralscan_ref/index.idx \
  --transcripts viralscan_ref/t2g.txt \
  --output output/ \
  --sample1 sample_R1.fastq.gz \
  --sample2 sample_R2.fastq.gz
```

Optional host pre-subtraction is available with `--host-filter starsolo`, but
it is an advanced, irreversible filter rather than a required first step. The
pre-v3 kallisto subtraction mode was removed because BUS output does not retain
exact read IDs and CB–UMI-wide deletion can remove unrelated fragments.

### Associate viral presence with host gene expression

After a completed run, `viralscan hostresponse` trains per-virus L2 logistic
regression models predicting candidate-positive vs. candidate-negative molecule-
support labels from host gene expression, then runs randomized Lasso stability selection to identify
robustly associated genes. Provide a matched host-gene h5ad (same barcodes as
the viralscan run):

```bash
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad
```

The same analysis runs inline during a full `viralscan` run when `--host-h5ad`
is supplied.

> **Depth confounding — read the depth baseline, not the AUC alone.** The default
> `counts >= threshold` label tracks sequencing depth, so the headline model AUC can
> partly reflect library size. Every run reports `depth_alone_auc` and per-gene
> depth-adjusted sensitivity metrics. `--label cpm` changes the label to a
> depth-normalized ratio and `--depth-match` constructs a depth-matched cohort;
> neither option proves that all depth or technical confounding has been removed.
> `%mito` is controlled by default. See the
> [CLI reference](docs/cli_reference.md) for interpretation guidance.

Optional pathway enrichment via gget requires the `[enrichment]` extra:

```bash
pip install "viralscan[enrichment]"
viralscan hostresponse -o output/sample/ --host-h5ad host_genes.h5ad --enrichment
```

**Chemistry preflight.** If a run yields near-zero viral (or host) counts,
verify the barcode chemistry before re-running — a wrong `--technology`/`--whitelist`
makes `bustools` silently discard most reads:

```bash
viralscan check-whitelist -s1 sample_R1.fastq.gz -w whitelist.txt -x 10xv3
```

---

## Output

For an input named `sample_R1.fastq.gz`, key results are written under
`output/sample/`:

| File | Description |
|------|-------------|
| `results/viral_summary.tsv` | Per-virus molecule estimates and called/all-barcode candidate-evidence denominators |
| `results/per_cell_viral.tsv` | Per-barcode × per-virus unique/allocated molecule estimates; allocated values may be fractional |
| `results/multimap_evidence.tsv` | Unique, ambiguous, and host-virus ambiguous viral evidence |
| `report.html` | Self-contained interactive HTML report |
| `kb-python/counts_unfiltered/adata_multimap.h5ad` | AnnData with the complete selected-method molecule matrix and audit layers |
| `plots/` | PNG plots and optional UMAP HTML files |

See [docs/output_reference.md](docs/output_reference.md) for the full output
schema including optional files (`cell_type_enrichment.tsv`, UMAP plots, etc.).

---

## Limitations

- **Formal v3 sensitivity, specificity, and comparative performance are not yet established.**
  The preregistered truth-panel holdout and harmonized comparator program must pass before those
  claims are made. Pre-v3 public-dataset results are historical and cannot validate v3.
- **Multimapper allocations are model-based estimates.** `em-global` fits one sample-wide
  abundance model; `em-cell` fits cell-local models shrunk toward that sample model. Their
  fractional allocations are not raw UMI observations or biological infection calls.
- **Production scope is human paired-end droplet scRNA-seq.** Individual chemistry and workflow
  combinations remain unvalidated until they pass the v3 truth panel. Other hosts, bulk, spatial,
  long-read, and single-nucleus-specific models are experimental or future work.
- **Cross-homology with host sequence** can create candidate viral evidence when a transcriptome
  reference lacks competing intronic or intergenic host sequence. The default
  `host-conservative` method allocates mixed host-virus molecules only among compatible host
  genes and keeps their ambiguity visible. Use full-host-genome competition and
  `viralscan evidence` for read-level diagnostics; neither alone confirms infection.
- **Ambient RNA** from high-burden cells is not corrected. In samples with extreme
  viral-burden heterogeneity, ambient RNA may create molecule support in cells without true
  intracellular viral transcription.
  Run SoupX or CellBender on the host matrix before using `--host-h5ad` if ambient correction
  is needed.

---

## Citation

If you use ViralScan in your work, please cite it using the metadata in
[CITATION.cff](CITATION.cff).  ViralScan is developed at the Leiden University
Medical Centre (LUMC), Leiden, the Netherlands.

---

## License

This project is licensed under the MIT License.  See the [LICENSE](LICENSE)
file for details.
