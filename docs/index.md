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

## Three steps

1. [Install](installation.md) the full toolchain with conda or pixi.
2. [Run](quickstart.md) one sample: build a combined host + virus reference,
   quantify, validate.
3. [Read the results](results.md): which file answers which question, what the
   units and denominators mean, and five checks before a candidate becomes a
   claim.

Worked notebooks are under [Examples](vignettes/README.md). Every flag, file and
column is in the [CLI reference](cli_reference.md) and the
[Output reference](output_reference.md); symptom-first troubleshooting is in the
[FAQ](faq.md).

## Inputs at a glance

- **Reads:** paired-end FASTQ, optionally gzipped. R1 carries the cell barcode
  and UMI. 10x 3′ and 5′ chemistries are detected from the reads before counting;
  Drop-seq is supported.
- **Reference:** a kallisto `index.idx` and `t2g.txt`. For host-aware work, build
  a combined host + virus index with `viralscan build-ref`, or reuse a compatible
  index that carries its build manifest. The bundled panel and host-control
  strategies are described under [Reference panel](reference_panel.md).
- **Samples:** one or more FASTQ pairs per run; separate multiple samples with
  commas.
- **Toolchain:** the full workflow needs Linux x86-64, Python 3.11, `kb`,
  `snakemake`, kallisto, bustools and STAR. EmptyDrops cell calling additionally
  needs R with `DropletUtils`.

## Project

Upgrading from an earlier release? See [migration](migration.md); v3 does not
migrate legacy outputs. Supported versions, contact routes and disclosure are in
the [release policy](release_policy.md), with the short statements in
[SUPPORT.md](https://github.com/mdmanurung/ViralScan/blob/main/SUPPORT.md) and
[SECURITY.md](https://github.com/mdmanurung/ViralScan/blob/main/SECURITY.md).

```{toctree}
:hidden:
:maxdepth: 1
:caption: Guide

installation
quickstart
results
vignettes/README
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Reference

reference_panel
cli_reference
output_reference
discovery
faq
api
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Project

migration
release_policy
changelog
```
