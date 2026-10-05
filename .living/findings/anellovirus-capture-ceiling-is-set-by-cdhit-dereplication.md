# F-013 — The anellovirus panel's sensitivity ceiling was set by CD-HIT, not by available diversity

**Date**: 2026-09-27
**Status**: mechanism identified and measured; the actionable lever is a reversal of upstream dereplication
**Found by**: tracing the panel's provenance while reviewing whether extra anellovirus genera could raise coverage

## The question that prompted this

"Can we make the TTV genomes more complete? Which extra anellovirus genera
should we add?" The premise — that our low measured capture reflects a sampling
gap we can close by adding genomes from under-represented genera — turns out to
be wrong in an instructive way.

## Where the panel actually comes from

The `source` column of `src/viralscan/data/anellovirus_accessions.tsv` tallies
2,022 `clareaulab` and 20 `viralscan-refseq`. The upstream is
`github.com/clareaulab/human_anellovirus_pangenome`, which built its set as:

1. **3,545** complete anellovirus genomes, NCBI Virus taxon 687329, **human host
   filter**, June 2025 freeze.
2. RepeatMasker against human sequence.
3. **CD-HIT at 95 % ANI / 85 % coverage → 2,023 representatives.**
4. dustmasker hardmasking of low-complexity regions.
5. ORF1 protein phylogeny (MAFFT → trimAl → IQ-TREE, 40 clusters), which
   **resolved genus for ~584 genomes NCBI leaves unclassified**.

## The mechanism

A 31-mer matches only with **zero** mismatches. At the CD-HIT threshold of 95 %
ANI, the expected fraction of 31-mers shared between a discarded genome and the
representative that absorbed it is:

```
0.95 ** 31 = 0.2039
```

Measured median canonical leave-one-out capture across the 2,042-genome panel:

```
0.2162
```

These agree to ~0.012. **CD-HIT's 95 % threshold, not NCBI's holdings, set our
sensitivity ceiling.** The 1,522 genomes dereplication discarded sit in exactly
the divergence band where 31-mer pseudoalignment fails — they were removed
*because* they were ~95 % identical to a kept representative, which is precisely
the identity at which they stop sharing k-mers with it.

**Caveat, not yet resolved:** the binomial assumes substitutions are independent
and uniformly distributed. [[kmer-capture-independence-model-overstates-detection]]
(F-011) showed that clustering of shared k-mers in conserved blocks makes the
empirical fragment capture 0.51 where the analytic model gives 1.00. The two
results must be reconciled before the correspondence above is treated as causal
rather than coincidental. Anellovirus divergence concentrates in hypervariable
ORF1 regions, which should push measured capture *away* from the binomial
expectation, so the closeness of the match deserves suspicion.

## Two things this falsifies

**1. "Our genus labels are unreliable."** I measured 582/2,042 panel genomes
with no genus in their GenBank lineage, and real disagreements between the panel
labels and GenBank (5 "Hetorquevirus" → Betatorquevirus, 2 → Sadetorquevirus,
single records → Yodtorquevirus and Lamedtorquevirus). I read this as our labels
being untrustworthy predictions. It is the opposite: those ~582 are the ~584 the
upstream ORF1 phylogeny **resolved**, and the disagreements are the phylogenetic
correction doing its job. Upstream reports Betatorquevirus at 1,558 by phylogeny
vs 1,360 by NCBI label; measured here independently, 1,542 vs 1,368 — the same
correction reproduced. **Do not "fix" the panel labels against GenBank lineages.**

**2. "~15 Anelloviridae genera are entirely missing from our panel."** True, and
correct. Every absent genus has **zero** complete genomes with
`"Homo sapiens"[Host]` (live NCBI query, 2026-09-27):

| genus | complete genomes | human-host |
|---|---|---|
| Lambdatorquevirus | 408 | 0 |
| Etatorquevirus | 256 | 0 |
| Iotatorquevirus | 139 | 0 |
| Kappatorquevirus | 24 | 0 |
| Rhotorquevirus | 19 | 0 |
| Epsilontorquevirus | 18 | 0 |
| Pitorquevirus | 14 | 0 |
| Thetatorquevirus | 10 | 0 |
| Sigmatorquevirus | 9 | 0 |
| Upsilontorquevirus | 8 | 0 |
| Mutorquevirus | 5 | 0 |
| Xitorquevirus | 4 | 0 |
| Zetatorquevirus | 2 | 0 |
| Delta-, Nu-, Tautorquevirus | 0 | 0 |

They are swine, feline, canine, tupaia and pinniped anelloviruses. The upstream
human-host filter excluded them correctly. Adding them to a human detection
panel buys no sensitivity and adds false-positive surface. The same reasoning
applies to the 6 Gyrovirus genomes already in the panel, which are avian
CAV relatives — they are a *contaminant-explanation* set, not a detection set.

Memtorquevirus and Samektorquevirus are likewise already exhausted: NCBI holds 3
and 6 complete genomes respectively; we hold 4 and 8.

## Ruled out

**Masking is not a confound.** The upstream set is hardmasked and
`measure_kmer_capture.kmers()` skips any window containing `N`, so masking could
in principle depress capture. Measured on the shipped panel FASTA: **728 N bases
in 5,998,625 (0.01 %)**, 16 of 2,042 genomes affected, worst genome 11 %.
Negligible.

## What follows

- The one lever that raises anellovirus capture is **restoring the 1,522
  genomes CD-HIT discarded** (2,023 → 3,545, +75 % genomes, ~10.5 Mb). Quantify
  the gain before committing to it; do not assume it.
- Genomes deposited since the June 2025 freeze are the only genuinely new
  material. 15 months of deposits are unexamined as of 2026-09-27.
- If capture is ceiling-limited by strain divergence rather than panel size,
  then growing the panel is the wrong response to a weak negative, and the tool
  should instead report a measured per-genus detection probability with every
  call. Anellovirus results remain `screening_only` (`REF-10`) regardless.
- `NC_038359.1` and `AB303562.1` are an exact duplicate pair that survived
  CD-HIT into the shipped panel (RefSeq mirror of a GenBank record); both show
  canonical leave-one-out of exactly 1.000000. Logged against `CAT-05`.

## Related

- [[kmer-capture-independence-model-overstates-detection]] — F-011, the
  clustering correction that this finding must be reconciled against.
- [[capture-metric-measured-forward-strand-only]] — F-012, the canonical-k-mer
  bug found in the same session; every number here uses the corrected metric.
