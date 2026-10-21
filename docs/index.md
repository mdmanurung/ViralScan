# ViralScan - viral load quantification from single-cell RNA-seq

[![CI](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml/badge.svg)](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mdmanurung/ViralScan/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mdmanurung/ViralScan)
[![PyPI](https://img.shields.io/pypi/v/ViralScan)](https://pypi.org/project/ViralScan/)

ViralScan is a Snakemake-driven command-line tool for quantifying viral load in
single-cell RNA-seq. It pseudoaligns paired-end FASTQ files with
[kb-python](https://www.kallistobus.tools/) (kallisto | bustools) against a viral or
combined host + virus reference. It reports per-cell and per-virus molecule estimates, an
HTML report and optional UMAP plots. ViralScan can build the reference and run
`kb ref` and `kb count` itself, so one tool covers reference building and detection. It is
developed at the Leiden University Medical Centre (LUMC).

::::{grid} 1 2 2 3
:gutter: 2

:::{grid-item-card} Installation
:link: installation
:link-type: doc

New to ViralScan? Start with the installation guide.
:::

:::{grid-item-card} Quickstart
:link: quickstart
:link-type: doc

Run your first sample from FASTQ to viral load.
:::

:::{grid-item-card} Tutorials
:link: vignettes/README
:link-type: doc

Worked examples, from building a reference to host response.
:::

:::{grid-item-card} CLI and outputs
:link: cli_reference
:link-type: doc

Every command and option, and every file ViralScan writes.
:::

:::{grid-item-card} API
:link: api
:link-type: doc

Use ViralScan from Python.
:::

:::{grid-item-card} Issues and support
:link: https://github.com/mdmanurung/ViralScan/issues

Found a bug? Need help? Open an issue on GitHub.
:::
::::

## Getting started

```bash
conda env create -f environment.yml   # kb-python, snakemake and dependencies
conda activate viralscan
python -m pip install -e .
viralscan --help                       # all commands and options
```

ViralScan needs a reference index made with `kb ref`. Provide a pre-built index and
transcript-to-gene map, or let ViralScan create one from FASTA and GTF files, from NCBI
accessions, or with `viralscan build-ref` for a combined host + virus reference. The viral
GTF panel is distributed through Zenodo (`viralscan data fetch`). Samples are paired-end
FASTQ files, optionally gzipped; both `-s1` (forward) and `-s2` (reverse) are required. See
[installation](installation.md) for the pip tier, containers and the tested toolchain.

All output is written under `-o/--output`, one directory per FASTQ pair.

**With a pre-built index:**

```bash
viralscan -t transcripts.txt -i index.idx -o output/ -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz
```

**With FASTA and GTF files.** ViralScan writes the index to `output/index/`:

```bash
viralscan -o output/ -reference True -s1 sample_1.fastq.gz -s2 sample_2.fastq.gz -fasta fasta.fasta -gtf gtf.gtf
```

**With several samples or reference files.** Separate names with commas, without spaces:

```bash
viralscan -t transcripts.txt -i index.idx -o output/ -s1 a_1.fastq.gz,b_1.fastq.gz -s2 a_2.fastq.gz,b_2.fastq.gz
```

Flags switch the optional stages on or off: multimapping, UMAP, host filtering, evidence,
gene programmes and host response. The [CLI reference](cli_reference.md) lists them, and
the [output reference](output_reference.md) describes the files. Upgrading from an earlier
release? See [migration](migration.md). Project policies are in
[SUPPORT.md](https://github.com/mdmanurung/ViralScan/blob/main/SUPPORT.md) (version and
support policy) and
[SECURITY.md](https://github.com/mdmanurung/ViralScan/blob/main/SECURITY.md)
(reporting vulnerabilities).

## License

ViralScan is released under the MIT license; see the LICENSE file.

```{toctree}
:hidden:
:maxdepth: 1
:caption: General

installation
quickstart
cli_reference
output_reference
discovery
api
api_reference
```

```{toctree}
:hidden:
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
:hidden:
:maxdepth: 1
:caption: About

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
