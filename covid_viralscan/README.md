# Detecting viral reads

Accession numbers:
- `NC_045512.2` - SARS-CoV-2 reference genome accession for ViralScan reference
  construction.

This is the compact accession to use for a SARS-CoV-2 detection reference. The
`gget virus` documentation uses this accession for the optimized SARS-CoV-2
single-accession workflow:

```bash
gget virus NC_045512.2 --is_accession --is_sars_cov2 --out covid_viralscan/gget_sars_cov2_refseq
```

For a ViralScan host-aware reference:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 \
  --output covid_viralscan/viralscan_ref \
  --ncbi-email mikhael.manurung@gmail.com
```

Path to CellRanger output:
Set `CELLRANGER_OUTS_A` and `CELLRANGER_OUTS_B` to the two local
CellRanger `outs/` directories before running the COVID analysis helpers.

Possible path to FASTQs:
Set `COVID_FASTQ_ARCHIVE` or point the helper scripts at an extracted local
FASTQ directory.

AIM:
- use ViralScan to detect viral reads in single-cell RNA-seq data, specifically for SARS-CoV-2.
