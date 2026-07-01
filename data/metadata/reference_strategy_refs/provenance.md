# Provenance: reference_strategy_refs

## Source

**Type**: derived (built locally from public reference sources)

**Origin**:
- Human genome/annotation: 10x Genomics `refdata-gex-GRCh38-2024-A` (genome FASTA
  + GTF), with the matching Ensembl 110 transcriptome used by the local kallisto
  human index. Located at
  `/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/`.
- Viral panel: Serratus viral panel merged with the ViralScan expanded
  anellovirus panel, under
  `/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/merged/`
  (`viral_serratus_plus_anellovirus.{fa,gtf}`).
- STARsolo indices: built locally via `STAR --runMode genomeGenerate` from the
  above FASTA/GTF. The all-virus and combined FASTA/GTF include 97 HUM_HERP6B
  pseudo-contigs synthesized from the kallisto transcriptome.

**Citation / accession**: Serratus (Edgar et al. 2022, Nature); 10x refdata-gex-GRCh38-2024-A; Ensembl release 110. Anellovirus accessions in `src/viralscan/data/anellovirus_accessions.tsv`.

## Acquisition details

**Date acquired**: 2026-06-27 (STARsolo genome dirs built; see mtimes in `reference_audit.tsv`)

**Obtained by**: mdmanurung

**Method**: Built with the exact commands stored under `build_commands` in
`reference_manifest.json` (repo root). Audited with `scripts/audit_reference_strategy.py`,
which produced `reference_audit.tsv`.

**Checksum**: Per-artifact SHA256 recorded in `reference_audit.tsv` (authoritative).
Not duplicated here to avoid drift. Verify with that file after any rebuild.

## Access restrictions

**Restriction level**: institutional-only

**Details**: Data lives on institutional storage (`/exports/...`). Human reference
is public (10x/Ensembl). No DUA on the derived indices, but paths are internal.

## Known issues

- **HHV-6A vs HHV-6B ambiguity**: the Serratus genome source contains NC_001664.4
  labeled *Human betaherpesvirus 6A*; HHV-6B target calls must use the HUM_HERP6B
  pseudo-contigs, not NC_001664.4. Getting this wrong conflates 6A/6B.
- **Preflight status at manifest creation**: `reference_manifest.json` notes STAR
  genomeGenerate index dirs "must be built and audited before benchmark submission."
  Confirm `references/starsolo/*/` are complete indices (not just FASTA/GTF inputs)
  before running the benchmark.
- **Cross-filesystem paths**: artifacts span `/exports/archive/...` and
  `/exports/para-lipg-hpc/...`; a move of either breaks the manifest paths.

## Contact

**Primary contact**: mdmanurung (mikhael.manurung@gmail.com)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-07-01 | Initial mycelium ingest (registered in place; data gitignored). |
