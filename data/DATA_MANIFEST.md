# Data Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### reference_strategy_fastqs
```yaml
name: reference_strategy_fastqs
type: genomic
source: NCBI SRA / ENA (SRR12682296 staged; SRR20710641, SRR8315713 from showcase fastq_root)
date_acquired: 2026-06-28
format: FASTQ.gz (paired-end 10xv3) + source_urls.tsv provenance
size: "~10.9 GB gzipped (staged sample only)"
raw_path: benchmark_inputs/reference_strategy/   # FASTQs gitignored; source_urls.tsv tracked
metadata_path: data/metadata/reference_strategy_fastqs/
status: raw
known_issues:
  - Only 1 of 3 benchmark runs staged locally; others depend on external showcase fastq_root
  - R2 (cDNA) ~8.7 GB gzipped
access_restrictions: none
tags: [fastq, 10x, scrna-seq, sra, ena, benchmark, input]
```

Paired-end 10x FASTQ inputs for the reference-strategy benchmark. Bytes are
gitignored and re-fetchable from ENA via each sample's tracked `source_urls.tsv`
(URL + MD5 + bytes). Feeds the ViralScan/STARsolo runs summarized in
`results/reference_strategy_benchmark.tsv`.

### reference_strategy_refs
```yaml
name: reference_strategy_refs
type: genomic
source: Derived — 10x GRCh38-2024-A + Serratus+expanded-anellovirus panel; STARsolo & kallisto indices built locally
date_acquired: 2026-06-27
format: STAR genomeDir, kallisto .idx, FASTA, GTF, t2g/TSV (20 artifacts)
size: "~64 GB (in-repo references/); more out-of-repo"
raw_path: data/raw/reference_strategy_refs/   # pointer doc only; bytes are gitignored / out-of-repo
metadata_path: data/metadata/reference_strategy_refs/
status: validated   # audited via reference_audit.tsv (SHA256 + feature counts)
known_issues:
  - HHV-6A (NC_001664.4) vs HHV-6B (HUM_HERP6B contigs) must not be conflated
  - Confirm STARsolo genomeDir completeness before benchmark submission
  - viralscan all_virus and combined share identical viral index/t2g (same SHA256)
access_restrictions: institutional-only
tags: [reference, kallisto, starsolo, benchmark, viral, GRCh38, anellovirus]
```

Reference indices and annotations for the reference-strategy benchmark (three
strategies × two aligners). Too large for git (~64 GB); registered in place and
gitignored, with machine-readable provenance in `reference_manifest.json` and
per-artifact SHA256/feature-count audit in `reference_audit.tsv` (both committed).
Rebuild commands are stored verbatim under `build_commands` in the manifest. This
dataset is the shared input for `results/reference_strategy_benchmark.tsv`.
