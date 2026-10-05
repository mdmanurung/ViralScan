# F-023 — Viral/host homology in the v3 panel is almost entirely low-complexity

**Status:** preliminary (one aligner, one masker, one host assembly; measured 2026-10-05)
**Task:** `REF-07` (host-homology table), feeding `VAL-01` / `DEF-00`
**Script:** `scripts/ref07_host_homology_table.py`, `scripts/ref07_level_proposal.py`
**Output:** `viral_ref_cat42d/ref07/host_homology_loci.tsv` (5,971 rows) and
`grch38_vs_viral.paf`

## Claim

Across the 2,343-genome v3 panel against GRCh38 (10x `refdata-gex-GRCh38-2024-A`,
sha256 `fb742121…`), nucleotide homology between viral and human sequence is
**99.7 % low-complexity**. There is effectively no genuine shared-gene homology
that could produce a confusable read.

## Evidence

| quantity | value |
|---|---|
| alignments ≥ 50 nt | 5,971, over 51 of 2,343 genomes, in **89** clusters |
| dustmasker coverage of the panel | **4.50 %** of bases (497,395 of 11,048,660; per-genome median 4.4 %, max 25.9 %) |
| dust coverage of **aligned viral bases** | **98.8 %** (1,097,260 of 1,110,425) |
| hits with viral-side dust ≥ 0.2 | 99.7 % (5,954 / 5,971) |
| clusters with mean dust ≥ 0.2 | 88.8 % (79 / 89) |
| hits sharing a 31-mer | 73.9 % (4,413) |
| clusters with a 31-mer **outside** low complexity | **1** of 89 |
| clusters touching a host exon | 23 of 89 |
| hits touching a transcript 3′ window (300 nt) | **12** of 5,971 |
| 31-mer + non-low-complexity + exonic | **0** |

The like-for-like comparison is the base-weighted one: dustmasker covers
**4.5 %** of the panel's bases but **98.8 %** of the bases that align to the
human genome. This is the structure of the sequence, not a masking artefact.
Row counts are dominated by two genomes -- HHV-6A and HHV-6B telomeric hits are
53.4 % of the 5,971 rows -- which is why the cluster and base views are given
beside them.

Top contributors are all repeat-driven: HHV-6B `AF157706.1` (1,751 hits, 100 %
low-complexity), HHV-6A `NC_001664.4` (1,439, 100 %), Cowpox `NC_003663.2` (690,
99 %), HPV19 `X74470.1` (388, 100 %), MPXV `NC_003310.1` (364, 100 %). HHV-6's
hits are its telomeric `TTAGGG` arrays against human telomeres — real biology,
and low-complexity.

## Sensitivity controls

Two, because the first one alone is too weak to carry the claim.

1. **Genome-level presence (weak).** Every genome whose spurious unique calls
   the REF-06 D-list removed appears in the table (8/8: HHV-6, MPXV, MOCV,
   HPV9, HPV77, EBV, CPXV, HHV-1). This only shows the genomes are represented.
   HHV-6 has 1,439 telomeric hits, so it would pass whatever happens at `p23`,
   and every genome in the set produced low-complexity hits in the first place.
2. **Synthetic recovery of non-repeat decoys (the real control).**
   `scripts/ref07_recovery_control.py` takes 20 dust-free human exons of at
   least 300 nt per tier from chr20-22, mutates them to a known identity,
   writes them as a fake viral panel, and runs the identical minimap2 command.

   | planted identity | recovered | median block |
   |---|---|---|
   | 0.95 | 20/20 (100 %) | 516 nt |
   | 0.90 | 20/20 (100 %) | 482 nt |
   | 0.85 | 16/20 (80 %) | 498 nt |
   | 0.80 | 12/20 (60 %) | 803 nt |

   So genuine non-repeat homology at 90-95 % identity would have been found
   without exception, and at 85 % four times in five. The absence of that class
   in the panel is a measurement, not blindness. Below 0.85 the floor falls and
   the finding makes no claim there -- which costs little, since a 90-nt read at
   80 % identity carries no exact 31-mer and pseudoalignment cannot see it
   either.

   Limit: the decoys are >= 300 nt. A short (90-150 nt) diverged match would
   recover worse than the table above, so this is the sensitivity floor for
   *long* homology.

## Why it matters

1. It corroborates [F-005] and [F-019] at the *reference* level rather than the
   read level. The covid Alphatorquevirus artefact was poly-G reads and the raft
   Gammatorquevirus artefact was poly-A capture; this says the panel offers
   almost no other route for a host read to look viral.
2. **`host_virus_homology` and `low_complexity` are not separable factors for
   this panel.** The VAL-01 design treats them as two axes, and the protocol
   strata are abundance × homology × chemistry. On this measurement a
   "homology" challenge read is a low-complexity read. `DEF-00` has to merge
   them, reduce homology to presence/absence, or state that the band structure
   has no material.
3. For a 3′ chemistry the exposure is smaller still: 12 of 5,971 hits touch a
   transcript 3′ window, and the realistic host-read route onto viral sequence
   is the poly-A tail (F-021), not shared genes.

## Limits

- Nucleotide level only, `minimap2 -k15 -w10 -s40 -m20`, blocks ≥ 50 nt.
  Protein-level homology (EBV `BARF1` ~ `CSF1R`, `BaRF1` ~ `RRM2`) is real but
  too diverged in nucleotide space to yield a confusable 90-nt read, which is
  the point rather than a gap.
- Dead end recorded: with the **genome** as the minimap2 reference, `asm10`/
  `asm20` return **zero** alignments for all 2,343 genomes. The asm presets want
  long colinear high-identity blocks. Indexing the 11 MB panel and streaming the
  genome as query is ~70x faster, uses 0.5 GB against 13 GB, and finds the
  signal.
- One host assembly, one masker. A second masker (`sdust`, `tantan`) was not
  run.
- Clusters join on a shared viral segment **or** a shared host locus
  (union-find). Viral-side merging alone put one host locus into both
  partitions: MPXV and Cowpox both hit `chr9:109765432-109765517` from separate
  records, and HHV-6A/6B/7 all hit the same telomeres. After the fix, 0 of
  1,721 1-kb host windows appear in both partitions.
- Realized holdout share over the 89 clusters is 0.213 against a nominal 0.30.
  Hash assignment is unbiased but coarse at this cluster count. These are
  challenge templates rather than samples, so the protocol's largest-remainder
  apportionment does not apply; the deviation is recorded rather than corrected.

## Implications

- `DEF-00`: `factors.host_virus_homology` cannot take identity/length bands as
  specified. See `docs/plans/2026-10-05-def00-protocol-amendment.md` row D1.
- `REF-09` ("prove planted human-homology reads cannot reach probable/strong")
  is now mostly a statement about low-complexity reads, which `DEF-01`'s
  read-artefact filter already targets.
- `VAL-01`'s `synthetic_host_homology` dataset should plant from these 89
  clusters, with their low-complexity status recorded per read.
