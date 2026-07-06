# reference_strategy_fastqs

Paired-end FASTQ inputs for the reference-strategy benchmark. The FASTQs
(~10.3 GB gzipped) are **gitignored**; only each sample's `source_urls.tsv`
(ENA URL + MD5 + byte count) is tracked, so the data is fully re-fetchable.

## Samples

The benchmark spans three SRA runs (see `srr` column in
`results/reference_strategy_benchmark.tsv`): `SRR12682296`, `SRR20710641`,
`SRR8315713`. **Only `SRR12682296` is staged under `benchmark_inputs/` here**;
the other two are read from the showcase data root
(`/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data`, see
`reference_manifest.json:fastq_root`).

| Run | Location | R1 / R2 |
|-----|----------|---------|
| SRR12682296 | `benchmark_inputs/reference_strategy/SRR12682296/` | `_1.fastq.gz` (2.27 GB), `_2.fastq.gz` (8.67 GB) |
| SRR20710641 | showcase `fastq_root` | not staged here |
| SRR8315713 | showcase `fastq_root` | not staged here |

## How to obtain

Each staged sample directory contains `source_urls.tsv` with the ENA FTP URL,
MD5, and byte count per mate. To re-fetch:

```bash
# Uses the tracked source_urls.tsv provenance
python scripts/fetch_reference_strategy_fastqs.py   # see script for args
# or manually:
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/096/SRR12682296/SRR12682296_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/096/SRR12682296/SRR12682296_2.fastq.gz
md5sum -c   # against source_urls.tsv md5 column
```

Expected MD5s (SRR12682296): R1 `dd1bfe5861d10c89e0f67b836d3ea262`,
R2 `54c5806e695152f7d03a01d8159f8261`.

See `data/metadata/reference_strategy_fastqs/` for schema, provenance, summary.
