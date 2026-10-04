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

---

## Update 2026-10-04 — the adversarial review landed: conclusion upheld, stated reasoning corrected

Biomni task `tsk_010G28jS5K1qC5TzDMva8eZR` (completed 2026-10-03, `model="max"`)
re-derived the read-level analysis independently from the 30 submitted reads and
the 8 implicated reference genomes fetched from NCBI, rather than from our
description. Verdict: **the reads are artefacts, but two of the three headline
arguments are non-diagnostic.**

**The decisive test, which we had not run.** Segmenting each 90-nt read into
`[body][poly-A run][3' tail]`: **0/30 read bodies** — the sequence 5' of the
poly-A run — align to any of the 8 anellovirus genomes (mlen >= 25, NM <= 2). The
bodies are low-complexity tracts, human Alu, or Illumina adapter; never viral. A
genuine 10x R2 read of a TTV mRNA 3' end must be `[complex TTV 3'UTR][untemplated
poly-A]`. It fails.

**Read composition (30-read cross-section, <= 2 mismatches):**

| Class | Fraction | Evidence |
|---|---|---|
| Carries literal reagent sequence | 22/30 (73 %) | 20/30 end in TSO-rc `GTACTCTGCGTTGATACCACTGCTT` at the **fixed** position 65–90, mean 0.1 mismatches; 1/30 body is TruSeq adapter; 2/30 have TSO-derived "cell barcodes" |
| Low-complexity body, no reagent motif | 7/30 (23 %) | C/A-rich, dinucleotide entropy < 3.5 bits |
| Complex body | 1/30 (3 %) | Complex but **non-viral** |

**What does not discriminate** (so we must stop leading with it):

- **Single-window pileup.** TTV mRNAs (3.0/1.2/1.0 kb) are alternatively spliced
  but share a common 3' end; 10x 3' chemistry sees only terminal fragments. A
  genuine infection gives the same single-bin pileup. Worse, the artefact piles
  up at the polyA site *by construction*, because that is the genome's longest
  templated A-tract (MZ286238.1 nt 2835–2865, A31; the other references have
  36–48-nt tracts of their own). Position is evidence for neither hypothesis.
- **Poly-A richness / >= 15-nt homopolymer.** A genuine 3'-end read carries the
  tail by definition.
- **NM = 0, as we reported it.** Re-alignment shows median aligned fraction
  **38 %** (mlen ~ 34 nt), MAPQ <= 4. NM = 0 holds over ~34 nt of A/C-rich
  sequence while ~62 % of the read — including the whole TSO — is soft-clipped.
  **Report aligned fraction alongside NM or the number misleads.**

**What does discriminate:** (1) body match, above; (2) TSO at the read 3' end is
physically impossible in a genuine molecule — in 10x 3' v3 only fragments
carrying the bead-oligo end receive P5, so the TSO end of the cDNA is discarded,
yet 20/30 reads end in verbatim TSO-rc; (3) the poly-A runs align base-for-base
to the genome's own A-tract and stop mid-tract, where a genuine untemplated tail
would extend past it; (4) perfect identity to 8 diverse references at once, where
TTV genotypes differ by > 30 %; (5) splice-junction reads would be decisive
positive evidence and there are none.

**The honest hole, which stays open.** Our homopolymer filters remove exactly the
reads that would prove a low-level genuine component, so the evidence **bounds**
rather than excludes it. Running the body census on all 19,785 aligned reads
would convert the 0/30 into a rule-of-three bound of ~3/19,785 ~ 0.015 % of
aligned reads — against 57,715 claimed UMI. Separately, the 79–100 % prevalence
figure that made us distrust our own negative is *plasma DNA*; PBMC DNA carriage
is ~20 % and transcription is found mainly in activated, not resting,
mononuclear cells, so a negative in a resting-PBMC 10x library is expected. And
our "positive" has the wrong shape: 12 % of cells at ~2.2 UMI each, uniformly
low, where genuine infection should be heavy-tailed.

**Also corrected:** the cat42b low-complexity masking test is weakly circular as
stated — we removed low-complexity k-mers from the reference and the artefact
reads match only low-complexity k-mers. The HPV16 control shows masking does not
break a complex, high-titre genome; it does not show the masked panel retains
sensitivity for a *low-abundance anellovirus*. The fix is the read-side census,
not another reference-side test.

**Decisive orthogonal, if this ever needs closing:** pan-anellovirus TaqMan qPCR
on the **final library** (~$100–300, days). At 57,715 genuine UMI the library
holds ~10^6–10^7 TTV molecules per 30-ng aliquot, Ct ~ 12–18; under the artefact
hypothesis no library molecule contains any TTV ORF sequence, Ct >= 38. Use
pan-anellovirus UTR primers plus several genotype-family ORF1 assays — **not**
primers against MZ286238, since that genotype call is itself a product of the
artefact.

**Literature:** no published report of TSO/poly-A chimeras mistaken for
anellovirus, and no published false-positive rate for anellovirus detection in
scRNA-seq — the gap cuts both ways. Nearest precedent is NIH-CQV/PHV, a
"novel circovirus-like hybrid" traced to contaminated silica spin columns
(Naccache et al., J Virol 2013, `10.1128/JVI.02323-13`). Systematic contamination:
Asplund et al. 2019 (`10.1016/j.cmi.2019.04.028`, > 65 % of viral sequences across
700 virome libraries linkable to lab components); Salter et al. 2014
(`10.1186/s12915-014-0087-z`). TSO/template-switch artefacts: Tang et al. 2012
(`10.1093/nar/gks1128`); Balázs et al. 2019 (`10.1186/s12864-019-6199-7`, template
switching mimicking alternative polyadenylation — directly our poly-A pileup);
Verwilt et al. 2023 (`10.1261/rna.079623.123`). Internal oligo-dT priming at
A-tracts: Nam et al. 2002 (`10.1073/pnas.092140899`); Svoboda et al. 2022
(`10.1093/nargab/lqac035`).

Trace and figure (`ttv_artefact_read_structure.png`):
https://biomni.phylo.bio/projects/prj_011fdGq9tp9AeHdJuQxbSZfD/tasks/tsk_010G28jS5K1qC5TzDMva8eZR

---

## Update 2026-10-04 (later) — full body census on all 19,785 aligned reads

`scripts/anello_body_census.py` runs the branch's own `is_complex_body` and
`has_tso` over every **primary** record (`-F 0x904`) on the 2,041 anellovirus
accessions in the covid x213 BAM. That is exactly 19,785 reads and 2,791
distinct CB+UMI. The 160 extra read names in the BAM have their primary
alignment elsewhere. It then places each read body on the panel itself (both
strands, ungapped, ≤ 2 mismatches, exhaustive pigeonhole seeding) as a label
that does not depend on those measures.

Orientation was checked on real data first. `read_seq()` equals the
pre-alignment `viral_reads.fasta` for 25/25 reads, 20 of them flag 16. No
primary record is hard-clipped.

**The bound.** No read has a body longer than 26 nt that places on any
anellovirus genome over ≥ 90 % of its length.

- **Review's criterion** (a ≥ 25 nt stretch at ≤ 2 mismatches): 81 reads,
  62 CB+UMI. Every one is a low-complexity G/C/A tract that lands on TTV's
  GC-rich region (`GCGGCGGCGG…`, G-tracts). Most place over exactly 25 nt of a
  50–60 nt body.
- **Full-length placement, reagent-free:** 29 reads. 11 have a low-complexity
  body, 5 are TSO variants with mismatches inside the core, and 13 are complex.
  All 13 are 20–26 nt and mostly G/A/C mosaics.
- **Short placements are not evidence.** A 20 nt TSO fragment places 19/20 on
  `MN774952.1`, yet no panel genome carries the TSO (≤ 3 mm). A 20–26 nt query
  against 12 Mb of both-strand sequence matches by chance at this rate.
- **Bound on the genuine component of the aligned subset:** at most 13 reads,
  that is ≤ 13/19,785 = 0.066 % (CP95 upper 0.11 %), or ≤ 13/2,791 CB+UMI =
  0.47 % (CP95 0.80 %). The call claimed 57,715 UMI. The point estimate is
  consistent with 0. The bound covers the **aligned** subset only, not the
  1.27 M equivalence-class reads.

**Validating the measures: they leak, and they leak the reagent.**

| | reads | CB+UMI |
|---|---|---|
| `has_tso` | 11,393 (57.6 %) | 751 |
| `is_complex_body` | 17,564 (88.8 %) | 1,351 |
| kept = complex and no TSO | 6,171 | 786 |
| ↳ carries the TSO 3′ 15-nt core | 4,762 | 345 |
| ↳ carries TruSeq R1 (`CTACACGACGCTCTTCCGATCT`) | 459 | 4 |
| ↳ no reagent found | 950 | 485 |
| kept, but the body places nowhere | 6,102 (98.9 % of kept) | 739 |

- `has_tso` meets its floor (57.6 % ≥ the 56.5 % counted verbatim). It misses
  any TSO **truncated at the read edge**, because `_contains` needs the whole
  25-mer inside the read. The dominant leak is the 5′-truncated read
  `GCAGTGGTATCAACGCAGAGTAC|T{30+}|…`.
- `read_body` returns whatever precedes the first homopolymer run. On a
  `[TSO][poly-T]…` or `[TruSeq R1][…][poly-A]` read that is the reagent itself.
  The TSO scores H = 3.52–3.58 and TruSeq 3.39, so it passes as a "complex body".
  Raising `MIN_BODY_ENTROPY` would not help: those values clear even the
  review's 3.5. The leak is reagent at the read edge, not low complexity.
- **Sensitivity** against the placement label: 69/81 body-mapping reads are kept.
  The 12 dropped all have `complex_body = 0`. They are G/C mosaics placing on
  TTV G-tracts, not genuine bodies.
- The strand split (complex 46 % forward vs 96 % reverse) comes from read
  composition. The reagent-led reads align reverse. It is not an orientation bug.
- The review's 30 reads cannot be identified in the upload (`sample.fa` holds
  146), so no per-read comparison with its classes was possible.

**Implication.** As shipped, `tso_fraction` under-reports and
`complex_body_fraction` over-reports on exactly this artefact class. Both errors
point the same way: the read looks more like a genuine molecule than it is.
Edge-truncated TSO matching is a flag-only fix and fits the locked decision.
Whether TruSeq and other reagent hits belong in `tso_fraction` is a
column-semantics question for the user.
