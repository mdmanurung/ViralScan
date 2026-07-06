# Summary Statistics: reference_strategy_fastqs

<!-- Generated: 2026-07-01 -->
<!-- Script: manual (from source_urls.tsv + ls) -->

## Overview

| Property | Value |
|----------|-------|
| Benchmark runs (total) | 3 (SRR12682296, SRR20710641, SRR8315713) |
| Runs staged locally | 1 (SRR12682296) |
| Staged size | ~10.9 GB gzipped (R1 2.27 GB + R2 8.67 GB) |
| Technology | 10xv3 paired-end scRNA-seq |
| Format | FASTQ.gz |

## Staged files (SRR12682296)

| File | Role | Size (gzip) | MD5 |
|------|------|-------------|-----|
| SRR12682296_1.fastq.gz | R1 (barcode + UMI) | 2.27 GB | dd1bfe5861d10c89e0f67b836d3ea262 |
| SRR12682296_2.fastq.gz | R2 (cDNA) | 8.67 GB | 54c5806e695152f7d03a01d8159f8261 |
| source_urls.tsv | provenance (tracked) | 282 B | — |

## Quality flags

- Read counts / per-base quality not computed here (would require decompressing
  10.9 GB). Run `seqkit stats` on the FASTQs if QC metrics are needed for a report.
- Two of three benchmark runs are NOT in this repo — downstream benchmark
  reproducibility depends on the external showcase `fastq_root`.

## Notes

- MD5s here verify the download; they are gzip-file MD5s (ENA-reported), not
  read-content hashes.
