# reference_strategy_refs

Reference indices and annotations for the **reference-strategy benchmark**
(comparing three reference strategies — `human_only`, `all_virus`, `combined` —
across two aligners, ViralScan/kallisto and STARsolo).

## This directory holds NO data bytes

The reference artifacts total **~64 GB** and are **too large for git**. They live
in place at the paths recorded in the machine-readable manifest and audit at the
repo root:

- `reference_manifest.json` — build commands, source releases, per-artifact paths.
- `reference_audit.tsv` — existence, size, mtime, **SHA256**, and feature counts
  for every artifact (the authoritative checksum record).

Both files are committed; this directory is a pointer only. The bulk data
(`references/`, and the STARsolo genome dirs) is gitignored.

## Where the bytes live

| Component | Location | Notes |
|-----------|----------|-------|
| STARsolo genome indices | `references/starsolo/{all_virus_serratus_plus_anellovirus,combined_GRCh38_2024A_serratus_plus_anellovirus,human_GRCh38_2024A}/` | In-repo, gitignored (`references/`). |
| ViralScan/kallisto index + t2g (viral) | `/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/merged/` | Out-of-repo, shared showcase dir. |
| Human kallisto index + t2g | `/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/human_index/` | Out-of-repo, archive. |
| Human genome FASTA + GTF (GRCh38-2024-A) | `/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/` | 10x Genomics refdata, Ensembl 110. |
| Anellovirus accession table | `src/viralscan/data/anellovirus_accessions.tsv` | In-repo, tracked (2042 accessions). |

## How to rebuild

The exact build commands (STAR `--runMode genomeGenerate`, `kb ref`) are stored
verbatim under `build_commands` in `reference_manifest.json`. To reconstruct:

1. Obtain the human refdata (10x `refdata-gex-GRCh38-2024-A`) and the merged viral
   panel FASTA/GTF (`viral_serratus_plus_anellovirus.{fa,gtf}`) from the showcase dir.
2. Run the `starsolo_*` and `kallisto_*` commands from `build_commands`.
3. Verify against `reference_audit.tsv` SHA256s and feature counts.

## Panel

`serratus_plus_expanded_anellovirus` — Serratus viral panel merged with the
ViralScan expanded anellovirus panel (2042 anellovirus accessions, 0 missing).
Expected detectable viruses in the benchmark samples include HHV-6B, EBV, HSV-1.

See `data/metadata/reference_strategy_refs/` for schema, provenance, and summary.
