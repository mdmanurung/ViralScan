# F-019 — Most of the published covid Alphatorquevirus signal is poly-G "no-signal" reads, not host homology

**Status:** confirmed (2026-09-29) · **Tags:** anellovirus, artefact, poly-G, homopolymer, CAT-17, CAT-30, F-005 · **Refines:** F-005

## Observation

SLURM job 25670680 (`viralscan_work/f005_readorigin/`, outside the repo) reproduced the published x213 run exactly:

- kallisto 0.51.1 / bustools 0.45.1, the same `index.idx`, `-x 10xv3`, and the CellRanger whitelist;
- `kallisto bus -n`, so every BUS record keeps its read number.

It extracted **every R2 read in an anellovirus EC** whose barcode is in the published matrix: 1,269,176 reads, 1.18M molecules. It then aligned them competitively with minimap2 (`-ax sr`) against GRCh38 2024-A genome + the published viral panel, keeping unmapped reads.

| class (all `anello_only` / `anello_host` ECs) | reads | % |
|---|---:|---:|
| best hit viral, homopolymer ≥15 nt | 1,143,775 | 90.1 |
| best hit host, homopolymer ≥15 nt | 124,918 | 9.8 |
| unmapped, homopolymer | 474 | 0.04 |
| **clean** (no homopolymer, no TSO) | **9** (all host-best) | 0.0007 |

- Of the "viral-best" reads, **1,143,452 hit one reference: NC_001479.1, encephalomyocarditis virus**, on its 115-nt poly(C) tract (position 148). No anellovirus.
- In a 200,000-read sample of the extracted reads, **91 % are ≥80 % poly-G**, e.g. a 90-nt run of G.
- `clean_viral_coverage.tsv` is empty. No clean, full-length read has a viral best hit, so no anellovirus genome has any clean coverage.

## Interpretation

- **Mechanism.** The 1.17M published Alphatorquevirus UMIs are mostly **poly-G reads**. On two-colour Illumina chemistry (NovaSeq/NextSeq), an empty cluster reads as G. These reads pseudoalign to a G/C-rich stretch of anellovirus transcripts in an index with no homopolymer masking; this is the `CAT-17` class.
  - Host-genome homology (F-005's mechanism) explains at most the ~10 % host-best share.
  - The BLAST sample (S100A16 mRNA, TSO/polyA junk) was drawn from the STAR-filtered survivors, so it describes the residue, not the bulk.
- **F-005's conclusion stands; its mechanism is corrected.** The signal is artefact. But the "0 of 4.5M" test sampled the first 5M library reads and counted only primary alignments. This test uses the assigned reads themselves.
- **Absence of true TTV is still not established.** 0 clean viral reads rules out *detectable* TTV among the reads this index assigned. It says nothing about divergent strains the index cannot capture (SENS audit: 1.36 % of anellovirus 31-mer space).

## Action

- CAT-30 retraction wording (grill Q11): "an artefact of poly-G no-signal reads (≈90 %) and host-homologous reads (≈10 %) (F-005, F-019); low-level divergent anellovirus not excluded."
- CAT-17: gate homopolymer/low-complexity reads or k-mers at count time, not only in references. A poly-G read filter before `kb count` would have removed ~90 % of this signal.

## Update 2026-09-30 — tonsil reproduces it; strand and rerun results

- **Tonsil SFL x223** (job 25671631, same kallisto path as covid): 3,020,792 reads were in anellovirus ECs.
  - 88.4 % are artefact-flagged and viral-best. The rest are host-best artefacts.
  - 3 clean full-length reads are host-best. 1 clean read is viral-best, and it is partial.
  - TONSIL-01 found 0 genuine anellovirus reads, so the false TTV signal is a library/pipeline artefact, not covid biology.
- **Current code reproduces the covid numbers** (job 25672275, v3 code, same index and whitelist): 1,083,687 (x213) and 1,499,051 (x216) Alphatorquevirus molecules, against 1,167,103 / 1,605,631 published.
- **The artefact is forward-strand specific** (F-020): under `--strand reverse` it falls 9–40×.
