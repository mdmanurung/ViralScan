# ViralScan

<!-- viralscan-claim:legacy-corrected-counts status=ineligible_legacy -->
> **V3 development notice:** pre-v3 corrected/combined counts are scientifically
> incompatible with the molecule-counting contract under development and are not
> eligible for v3 quantitative claims. Rebuild from retained corrected BUS/reference
> inputs or FASTQ; do not numerically migrate old H5AD values.

[![CI](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml/badge.svg)](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mdmanurung/ViralScan/branch/main/graph/badge.svg)](https://codecov.io/gh/mdmanurung/ViralScan)
[![PyPI](https://img.shields.io/pypi/v/ViralScan)](https://pypi.org/project/ViralScan/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](https://pypi.org/project/ViralScan/)

## Introduction

**ViralScan** is a Snakemake-based computational bioinformatics framework for
quantifying viral load in single-cell RNA-seq samples. It pseudoaligns paired-end
FASTQ files with [kallisto | bustools](https://www.kallistobus.tools/) (via
`kb-python`) and performs the downstream analysis: per-cell and per-virus molecule
estimates, super-expressor detection, bar plots of detected UMIs per virus, an
interactive HTML report and optional UMAP visualisations. For human samples, the
recommended workflow is a combined host+virus reference so host and viral targets
compete in one quantification step.

The tool can use both `kb ref` and `kb count` to perform scalable detection of the
viral load. ViralScan is developed and maintained at the
[Leiden University Medical Centre (LUMC)](https://www.lumc.nl/). The documentation
is built with Sphinx from `docs/` (`pip install ".[docs]"`, then
`python -m sphinx -b html docs docs/_build/html`).

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
- **Cell-type enrichment analysis** (`--cell-type-labels`) — Fisher exact test with
  BH correction for each virus × cell-type combination
- **NCBI auto-fetch + local cache** (`--ncbi-accession`) — downloads FASTA +
  GTF via E-utilities and caches under `~/.cache/viralscan/ncbi/`
- **`viralscan build-ref`** — builds a combined host + virus kallisto index
  in one command
- **Host-response analysis** (`viralscan hostresponse` or `--host-h5ad`) —
  associates per-virus presence with host gene expression via L2 logistic
  regression + randomized-Lasso selection frequency (uncalibrated), with optional gget pathway enrichment
  (default `--cv auto`: group-disjoint CV when a donor/sample column is found or named, otherwise cell-level CV with a donor-leakage warning)
- **Docker and Singularity containers** provided for fully reproducible runs
- **Snakemake backend** — each step is a named rule with logged output

---

## Installation

ViralScan can be installed using `pip` in a Python environment:

```bash
pip install ViralScan            # pip tier: no Snakemake
viralscan doctor --profile pip
```

A conda environment is recommended. It provides `snakemake`, `kb-python` and the
external tools (the **full tier**):

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install .          # use "-e ." for a development checkout
```

The minimal original-style environment is
`conda create -n bioenv -c conda-forge -c bioconda snakemake kb-python`, but only
`environment.yml` is the tested full-workflow toolchain (Python 3.11 on `linux-64`).

ViralScan is command-line based. To view all options, including an example, run:

```bash
viralscan --help
```

Note: the PyPI release can lag this repository. The commands documented here
track the source tree (v3 development); if `pip install ViralScan` installs an
older release that lacks them, install from the checkout instead
(`python -m pip install .` in a clone).

The pip tier covers `build-ref`, `evidence`, `validate-run` and `doctor --profile pip`.
It does not run the quantification workflow. `pip install "viralscan[full]"` adds
the Python workflow dependencies (`snakemake`, `scanpy`, ...), but `kb`, `kallisto`,
`bustools` and `STAR` are native binaries that pip does not provide: use the Conda
flow above (or the container) for the full workflow, then check it with
`viralscan doctor --profile full`. Without `snakemake`/`kb` on `PATH`, a workflow
run stops with an error naming the `full` extra.

### Container

```bash
docker run --rm -v "$PWD":/data ghcr.io/mdmanurung/viralscan:latest \
    --help
```

Release container builds require the tested wheel and committed runtime lock.
For HPC, generate the Apptainer definition from the tested OCI digest; see
[installation](docs/installation.md).

---

## Input Data

### Reference index

ViralScan requires an index created with `kb ref` (from the `kb-python` package).
There are 2 options:

- provide your own pre-built index and transcript-to-gene map, or
- let ViralScan create the index for you from FASTA + GTF files, or from NCBI
  accessions (`viralscan build-ref` builds a combined host + virus index).

The standard 195-virus GTF annotation panel is hosted on Zenodo
(`10.5281/zenodo.20112332`); see "Viral annotation panel cache" below.

### Samples

ViralScan expects paired-end FASTQ files (potentially gunzipped). Both the forward
(`--r1`) and the reverse (`--r2`) file must be provided.

---

## User Guide

All output (counts, logs, plots and the HTML report) is written under the
directory given by `-o/--output`. Each FASTQ pair gets its own sample directory
named from the R1 filename before the first underscore; `sample_R1.fastq.gz`
writes to `output/sample/`. Quantification is `viralscan run`; the original
single-dash spellings (`-s1 -s2 -gtf -fasta -ref ...`) were removed in 3.0 (see
[docs/migration.md](docs/migration.md)).

### Option 1 - You already have a reference index

```bash
viralscan run -t transcripts.txt -i index.idx -o output/ --r1 sample_1.fastq.gz --r2 sample_2.fastq.gz
```

### Option 2 - ViralScan creates the index from FASTA + GTF

Provide the GTF and FASTA file(s). The index is placed in the output directory, in
a subfolder called `index`:

```bash
viralscan run -o output/ --reference --r1 sample_1.fastq.gz --r2 sample_2.fastq.gz --fasta fasta.fasta --gtf gtf.gtf
```

`--reference`, `--visual` and `--multimapping` take no value; turn the last two off with
`--no-visual` / `--no-multimapping`.

### Multiple samples, FASTA and GTF files

Split the names with a comma (without a space in between):

```bash
# multiple samples
viralscan run -t transcripts.txt -i index.idx -o output/ --r1 sample1_1.fastq.gz,sample2_1.fastq.gz --r2 sample1_2.fastq.gz,sample2_2.fastq.gz

# multiple FASTA and GTF files
viralscan run -o output/ --reference --r1 sample_1.fastq.gz --r2 sample_2.fastq.gz --fasta fasta1.fasta,fasta2.fasta --gtf gtf1.gtf,gtf2.gtf
```

### Option 3 - NCBI accession auto-fetch

Provide one or more RefSeq accession numbers (comma-separated). ViralScan downloads
the FASTA and GTF from NCBI, builds the index, and runs quantification in one step.
NCBI asks for a contact email whenever something has to be downloaded: supply
`--ncbi-email` or set the `NCBI_EMAIL` environment variable. A run whose sequences
are all in the local cache (`~/.cache/viralscan/ncbi/`) needs no email:

```bash
export NCBI_EMAIL=you@example.org
viralscan run --ncbi-accession NC_002021.3,NC_001512.1 -o output/ --r1 sample_1.fastq.gz --r2 sample_2.fastq.gz
```

### Viral annotation panel cache

The standard 195-virus GTF panel is hosted on Zenodo
(`10.5281/zenodo.20112332`). It is required only for an index that has no
`<index>.build_manifest.json` (an index built by another tool, or by an older
`build-ref`): such a run reads its viral gene list from this panel. An index built
by `viralscan build-ref` carries its own gene list in that manifest and does not
need the panel. When it is needed, fetch it once:

```bash
viralscan data fetch
```

The files are cached under `~/.cache/viralscan/data/`. Workflows that pass
custom `--gtf` files still use those annotations directly.

For shared HPC caches, fetch to a shared root and pass the same root during
quantification:

```bash
viralscan data fetch --cache-dir /shared/viralscan-cache
viralscan run --cache-dir /shared/viralscan-cache ...
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
<!-- viralscan-claim:v3-reference-homology status=implemented_not_calibrated -->
Use `--genome-dlist GRCh38.fa` to mask host-genomic shared k-mers and record raw
per-virus host-homology metrics; this full-genome build is intentionally an
explicit, compute-intensive mode.

Then quantify with the generated files. The default
`--multimap-method conservative` (recorded as `host-conservative`) keeps host-virus ambiguous
equivalence-class mass out of primary viral estimates (use `equal` for an
explicit equal-allocation comparison):

```bash
viralscan run \
  --index viralscan_ref/index.idx \
  --transcripts viralscan_ref/t2g.txt \
  --output output/ \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz
```

Optional host pre-subtraction is available with `--host-filter starsolo`, but
it is an advanced, irreversible filter rather than a required first step. The
pre-v3 kallisto subtraction mode was removed because BUS output does not retain
exact read IDs and CB–UMI-wide deletion can remove unrelated fragments.

### Associate viral presence with host gene expression

After a completed run, `viralscan hostresponse` trains per-virus L2 logistic
regression models predicting candidate-positive vs. candidate-negative molecule-
support labels from host gene expression, then runs a randomized-Lasso selection-frequency ranking (uncalibrated: it does
not control an error rate, and with hundreds of positives noise genes pass the threshold)
to rank associated genes. Metrics use `--cv auto` by default: group-disjoint CV on `--groups <column>`, or on the one donor/sample-like obs column if exactly one exists; otherwise cell-level CV with a donor-leakage warning (`--cv cell` forces the legacy cell-level scheme, where a donor can span train and test). Provide a matched host-gene h5ad (same barcodes as
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
viralscan hostresponse --run-dir output/sample/ --host-h5ad host_genes.h5ad --enrichment
```

**Chemistry preflight.** If a run yields near-zero viral (or host) counts,
verify the barcode chemistry before re-running — a wrong `--technology`/`--whitelist`
makes `bustools` silently discard most reads:

```bash
viralscan check-whitelist --r1 sample_R1.fastq.gz -w whitelist.txt -x 10xv3
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

## Reproducibility

Validation experiments, benchmarks and manuscript analyses live in [ViralScan-reproducibility](https://github.com/mdmanurung/ViralScan-reproducibility).

---

## Citation

If you use ViralScan in your work, please cite it using the metadata in
[CITATION.cff](CITATION.cff).  ViralScan is developed at the Leiden University
Medical Centre (LUMC), Leiden, the Netherlands.

---

## License

This project is licensed under the MIT License.  See the [LICENSE](LICENSE)
file for details.
