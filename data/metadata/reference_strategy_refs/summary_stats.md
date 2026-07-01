# Summary Statistics: reference_strategy_refs

<!-- Generated: 2026-07-01 -->
<!-- Script: derived from reference_audit.tsv (scripts/audit_reference_strategy.py) -->

## Overview

| Property | Value |
|----------|-------|
| Reference artifacts | 20 (across 2 aligners × 3 strategies + shared human/viral inputs) |
| Aligners | ViralScan (kallisto/bustools), STARsolo |
| Strategies | human_only, all_virus, combined |
| Total size (in-repo `references/`) | ~64 GB |
| Largest single file | STARsolo combined.gtf, ~1.69 GB (3,310,456 features) |
| Viral panel | serratus_plus_expanded_anellovirus |
| Format | STAR genomeDir, kallisto .idx, FASTA, GTF, t2g/TSV |

## Artifact summaries (from reference_audit.tsv)

| Artifact | Size | Feature count | Expected viruses |
|----------|------|---------------|------------------|
| starsolo all_virus genome_fasta | 10.2 MB | — | HHV-6B, HSV-1 |
| starsolo all_virus genome_gtf | 4.16 MB | 17,295 | HHV-6B, EBV, HSV-1 |
| starsolo combined genome_fasta | 3.16 GB | — | HHV-6B |
| starsolo combined genome_gtf | 1.69 GB | 3,310,456 | HHV-6B |
| viralscan all_virus kallisto_index | 390 MB | — | — |
| viralscan all_virus t2g | 8.54 MB | 230,901 | HHV-6B, EBV, HSV-1 |
| viralscan all_virus gtf | 1.15 MB | 8,929 | HHV-6B, EBV, HSV-1 |
| viralscan human_only kallisto_index | 409 MB | — | — |
| viralscan human_only t2g | 18.5 MB | 226,005 | — |
| human genome_fasta (GRCh38-2024-A) | 3.15 GB | — | — |
| human genes_gtf | 1.69 GB | 3,293,161 | — |
| anellovirus accession table | 149 KB | 2,043 | not_detected_in_sample |

## Quality flags

- All 20 audited artifacts report `exists=True` in `reference_audit.tsv`.
- Confirm STARsolo `human_only` genomeDir completeness before the benchmark (the
  audit records the directory but not per-file index contents).
- HHV-6B must be scored on HUM_HERP6B contigs — see provenance "Known issues".

## Notes

- Checksums (SHA256) for every artifact are in `reference_audit.tsv`; this summary
  intentionally does not duplicate them.
- `combined` viralscan index/t2g are identical files to `all_virus` (same SHA256):
  the ViralScan combined strategy reuses the all-virus viral index; the human side
  differs by run configuration, not by a separate viral index.
