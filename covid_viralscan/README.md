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

Path to cellranger output:
/exports/para-lipg-hpc/Youvika/20250605_scRNAseq_YS/20250814_tino_scRNAseq_batch2_YS/data_raw/202502341a_count_v2

Possible path to fastq:
/exports/para-lipg-hpc/Youvika/Project_202502341_fastq.tar.gz

AIM:
- use ViralScan to detect viral reads in single-cell RNA-seq data, specifically for SARS-CoV-2.