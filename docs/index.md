# ViralScan

**ViralScan quantifies viral RNA at single-cell resolution.** It pseudoaligns
paired-end FASTQ files against a viral or combined host + virus kallisto index,
counts corrected (cell barcode, UMI) molecules for each virus and cell, and
writes a report, per-virus and per-cell tables, and provenance receipts.

ViralScan reports *candidate* evidence with explicit denominators. A `detected`
row is not a diagnosis, and passing a threshold does not establish infection.

[![CI](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml/badge.svg)](https://github.com/mdmanurung/ViralScan/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/mdmanurung/ViralScan/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mdmanurung/ViralScan)
[![PyPI](https://img.shields.io/pypi/v/ViralScan)](https://pypi.org/project/ViralScan/)

## Choose your next step

| Task | Start here |
|---|---|
| Install ViralScan | [Installation](installation.md) — the recommended conda route, then the pip and container alternatives |
| Run one sample | [Run ViralScan](quickstart.md) — preflight, reference, quantification, validation |
| Read the results | [Inspect results](inspect_results.md) — report, tables, units, denominators |
| Judge a candidate | [Assess the evidence](assess_evidence.md) — five checks from run integrity to conclusion limits |
| Browse worked examples | [Examples](vignettes/README.md) — notebooks from a full run to host response |

## Inputs at a glance

- **Reads:** paired-end FASTQ, optionally gzipped. R1 carries the cell barcode
  and UMI. 10x 3′ and 5′ chemistries are detected from the reads before counting;
  Drop-seq is supported.
- **Reference:** a kallisto `index.idx` and `t2g.txt`. For host-aware work, build
  a combined host + virus index with `viralscan build-ref`, or reuse a compatible
  index that carries its build manifest.
- **Samples:** one or more FASTQ pairs per run; separate multiple samples with
  commas.
- **Toolchain:** the full workflow needs Linux x86-64, Python 3.11, `kb`,
  `snakemake`, kallisto, bustools and STAR. EmptyDrops cell calling additionally
  needs R with `DropletUtils`.

Reference modes, every flag and every output column are under
[Reference & help](reference_help.md).

## Project information

Project policies: [SUPPORT.md](https://github.com/mdmanurung/ViralScan/blob/main/SUPPORT.md)
(version and support policy) and
[SECURITY.md](https://github.com/mdmanurung/ViralScan/blob/main/SECURITY.md)
(reporting vulnerabilities). Upgrading from an earlier release? See
[migration](migration.md); v3 does not migrate legacy outputs.

```{toctree}
:hidden:
:maxdepth: 1
:caption: General

installation
quickstart
inspect_results
assess_evidence
cli_reference
output_reference
reference_help
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
project_info
changelog
```
