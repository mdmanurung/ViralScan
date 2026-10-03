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

---

## Update 2026-10-03 — the mechanism has a second, larger component: TSO/poly-A chimeras

F-019 characterised the 1,269,176 reads in anellovirus **equivalence classes**
(91 % ≥80 % poly-G). This update looks at the 19,785 reads that actually
**aligned** to anellovirus references in the competitive minimap2 BAM
(`covid_viralscan/results_hostfilter/LUM-SJ-x213-g/evidence/viral_reads.bam`),
i.e. the strongest-looking subset, the one a reviewer would call real.

They are a different artefact class from the poly-G reads:

- 99.91 % contain a homopolymer run ≥15 nt, but only **0.01 %** are ≥80 % poly-G
  and **0.07 %** ≥80 % poly-A; mean GC 0.361. So they are *not* no-signal reads.
- **56.5 % contain the 10x template-switch oligo verbatim**
  (`AAGCAGTGGTATCAACGCAGAGTAC`). The TSO is a synthetic oligo; it is not viral
  sequence and not human sequence.
- Typical structure: a short low-complexity stretch, then poly-A of 20–45 nt,
  then the TSO.
- Some reads carry the TSO *in the cell-barcode position of R1* — barcodes
  literally `AAGCAGTGGTATCAAC` / `GCAGTGGTATCAACGC`. Those molecules have no
  real barcode at all.

Positional profile per reference (genomes ~3 kb, so ~30 possible 100-nt bins):

| reference | reads | distinct 100-nt bins | top-bin share | ≥15-nt homopolymer | median NM |
|---|---|---|---|---|---|
| MZ286238.1 | 6849 | 2 | 100.0 % | 100.0 % | 0 |
| MW455373.1 | 3684 | 1 | 100.0 % | 100.0 % | 0 |
| MW455378.1 | 3441 | 1 | 100.0 % | 100.0 % | 0 |
| MW455365.1 | 2428 | 1 | 100.0 % | 100.0 % | 0 |
| KP343825.1 | 1390 | 3 | 99.4 % | 99.6 % | 0 |
| MW455439.1 | 1104 | 1 | 100.0 % | 100.0 % | 0 |
| MN771265.1 | 486 | 1 | 100.0 % | 100.0 % | 0 |
| KP343824.1 | 120 | 1 | 100.0 % | 100.0 % | 0 |

Every reference: all reads in **one** 100-nt window, perfect identity (NM = 0,
because a homopolymer matches a homopolymer tract exactly), zero breadth.

**Caveat kept open deliberately.** Single-window pileup is weaker evidence than
it looks: 10x 3′ chemistry is 3′-biased and TTV mRNA is polyadenylated, so a
genuine 3′-end read would also be A-rich and would also pile up near the polyA
site. The decisive discriminators here are the **TSO content** (synthetic, cannot
be viral) and the absence of any read with non-adapter, non-homopolymer viral
sequence — not the pileup on its own. An adversarial review of exactly this point
was commissioned (Biomni task `tsk_010G28jS5K1qC5TzDMva8eZR`, 2026-10-03).

**What this changes.** F-005 said host homology; F-019 said poly-G; this says the
aligned subset is predominantly TSO/poly-A chimera. All three are artefact, but
the actionable filter is adapter/TSO + poly-A trimming plus a complexity gate
before counting (`ANELLO-PRIOR.3`), not a host-homology or N-masking fix.

**Independent corroboration.** The ANDET-09 STARsolo branch, whose filters
require ≥80 % of the read matched at ≤8 % mismatch, returns **0** anellovirus
molecules on 5 M host-unmapped reads of this same library — a chimera can only
align over its poly-A stretch, roughly a third of its length. The cat42b
low-complexity mask independently takes the same library's call from 57,715 UMI
to 0.
