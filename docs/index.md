# ViralScan Documentation

**ViralScan** is a Snakemake-driven Python bioinformatics CLI that quantifies
viral load from paired-end FASTQ samples using
[kb-python](https://www.kallistobus.tools/) (kallisto + bustools).

[![CI](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml/badge.svg)](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mdmanurung/ViralScan/branch/main/graph/badge.svg)](https://codecov.io/gh/mdmanurung/ViralScan)
[![PyPI](https://img.shields.io/pypi/v/ViralScan)](https://pypi.org/project/ViralScan/)

---

## Where to start

1. Install ViralScan with the [conda workflow](installation.md) unless your
   system already provides `kb`, `snakemake`, and the Python dependencies.
2. Run the [quickstart](quickstart.md) with the reference mode that matches
   the files you have: pre-built `kb ref` outputs, FASTA/GTF files, or NCBI
   accessions.
3. Use the [output reference](output_reference.md) to find the per-sample
   result tables, HTML report, plots, and AnnData file.
4. Open the [vignettes](vignettes/README.md) for executable, task-oriented
   examples — from a full run to multimapping correction, cell-calling
   denominators, enrichment, QC, specificity, and host-response.

```{toctree}
:maxdepth: 2
:caption: User Guide

installation
quickstart
cli_reference
output_reference
reference_panel
faq
```

```{toctree}
:maxdepth: 1
:caption: Vignettes

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
:caption: Developer Reference

api
changelog
```
