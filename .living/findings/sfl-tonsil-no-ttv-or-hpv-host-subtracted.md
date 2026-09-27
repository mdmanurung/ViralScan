# F-010 — SFL tonsil: no TTV or HPV in host-subtracted reads; the raw "hits" were low-complexity sequence

**Date**: 2026-09-27
**Status**: result recorded (screening; not an informative negative)
**Found by**: `TONSIL-01` (PLAN WP4J)

## Question

Can torquetenovirus (TTV) or HPV be detected in the SFL tonsil CITE-seq pool?
The pool is 24 hashtagged donors, `s1`–`s24`: 18 adults with non-active
tonsillitis and 6 children aged 2–9 with hypertrophy.

## Design

- **Library.** Only `x223` is gene expression; `x225` is ADT. The chemistry is
  10x 5′ v3 R2-only.
- **Host exclusion.** The screen used cellranger's own GRCh38-2024-A alignment.
  Only reads cellranger left unmapped were screened, which excludes F-005-type
  host-homology reads by construction.
- **Pipeline.**
  1. Trim TSO and adapter; filter by complexity.
  2. minimap2 against 2,326 viral genomes: 2,042 anelloviruses and 16 HPV
     types, with the HPV set from WP4E.
  3. Competitive re-check against GRCh38 plus the same viral set.
  4. Remove UMI duplicates; attribute cells to donors through the cellhashr
     singlets.
- **Decision rules.** Written into PLAN.md before any result was seen.

## Result

| Stage | Reads |
|---|---|
| Unmapped GEX reads (valid barcode) | 118,856,604 |
| Pass prefilter | 24,969,278 |
| Any viral alignment | 34,208 |
| ≥ 50 bp aligned, ≥ 85 % identity, DUST ≤ 20 | 1,511 |
| Best alignment viral, not host | 23 |

- **Anelloviridae: 0 reads in every donor.** There were 1,852 raw hits, all
  partial alignments of 20–49 bp on low-complexity sequence:
  - 1,485 are a 28-bp match to a CAG trinucleotide repeat (`AB303556.1:2318`).
  - The rest are poly(A) runs joined to reverse-complement TSO, on A-rich
    anellovirus regions.
  - The one hit that passed the viral filters aligns better to human.
- **HPV (16 types): 0 reads.** There were 2 raw hits, both shorter than 50 bp.
- **The 23 survivors are artifacts or unattributable:**
  - **HCV:** every read falls on the 3′-UTR poly(U/UC) tract.
  - **Macaque *Cercopithecine herpesvirus 2*:** reads fall on two GC-rich
    positions shared across donors.
  - **HSV-1:** 2 reads at one position.
  - **HHV-6B:** 5 host-free reads at 5 positions, in 5 barcodes that are all
    doublets, hashtag-negatives or non-cells. They are below any call and
    cannot be attributed to a donor.
- **Positive-control plant.** 1,000 reads per target, in 5′ R2 geometry, into
  1 M real reads. Recovery was 99.2 % for HPV16 (from p97) and 99.4 % for TTV,
  all attributed to `s1`. No planted read was lost to the host re-check.

## Why this is not an absence claim

1. **Unscreened reads.** The per-sample BAM holds only reads with a valid
   barcode, 1.35 B of 1.61 B GEX reads.
2. **Plant is not an LOD.** The plant used genomes that are in the database, so
   it does not measure sensitivity to divergent anelloviruses. That is
   `ANDET-08`.
3. **Low priors.** The EBV zero was predicted (E[EBV+ cells] ≈ 0.02–0.15). The
   HPV prior in benign tonsil suspensions, which are mostly lymphocytes, is low.
4. **Anellovirus scope.** Any anellovirus result is `screening_only` (`REF-10`).

## What it teaches

The anellovirus false-positive mechanism is visible here at read level:
- **The source.** Low-complexity and repeat sequence gives short, perfect
  partial matches to A-rich and CAG-like stretches of anellovirus genomes.
- **The fix.** Requiring ≥ 50 bp aligned plus a competitive host re-check
  removes all of them.

This is direct support for the `ANDET-01`/`ANDET-02` gates.
