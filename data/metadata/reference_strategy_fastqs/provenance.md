# Provenance: reference_strategy_fastqs

## Source

**Type**: database (public SRA/ENA)

**Origin**:
- Database: ENA / SRA. Staged run **SRR12682296** downloaded from ENA FTP
  (`https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/096/SRR12682296/`). Full URLs +
  MD5 + byte counts in `benchmark_inputs/reference_strategy/SRR12682296/source_urls.tsv`.
- Two further benchmark runs (SRR20710641, SRR8315713) are read from the showcase
  `fastq_root` (`/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data`) and are
  not staged in this repo.

**Citation / accession**: SRR12682296, SRR20710641, SRR8315713 (NCBI SRA).

## Acquisition details

**Date acquired**: 2026-06-28 (SRR12682296 FASTQ mtimes)

**Obtained by**: mdmanurung

**Method**: `scripts/fetch_reference_strategy_fastqs.py` (downloads per
source_urls.tsv). MD5s recorded for verification.

**Checksum**: MD5 per mate in `source_urls.tsv` — R1 `dd1bfe5861d10c89e0f67b836d3ea262`,
R2 `54c5806e695152f7d03a01d8159f8261`.

## Access restrictions

**Restriction level**: none (public SRA data)

**Details**: Publicly available via ENA/SRA. Repo stores only URLs/checksums, not bytes.

## Known issues

- Only 1 of 3 benchmark runs is staged locally; the other two depend on the
  external showcase `fastq_root` being present.
- R2 (cDNA) is ~8.7 GB gzipped — large; ensure disk space before re-fetching.

## Contact

**Primary contact**: mdmanurung (mikhael.manurung@gmail.com)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-07-01 | Initial mycelium ingest (SRR12682296 staged; FASTQs gitignored). |
