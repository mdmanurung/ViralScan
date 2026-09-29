# F-020 — On 10x 5′ libraries, kallisto's forward-strand default loses ~85–88 % of reads; reverse strand also suppresses the poly-G artefact

**Status:** confirmed (2026-09-30) · **Tags:** 5-prime, strand, chemistry, kallisto, poly-G, HHV-6B, covid, MECH-D, R2.6

## Observation

Job 25672276 (`scripts/slurm_strand_test.sh`) ran `kb count --strand forward|reverse|unstranded` on the first 4M read pairs of each 5′ library:

| library | strand | pseudoaligned | host UMIs | viral UMIs (top) |
|---|---|---:|---:|---|
| covid x213 (5′, 16+12, `-x 10xv3`) | forward (kallisto 10x default) | 6.5 % | 126,504 | 5,662 (MW455439.1 anello 5,652) |
| | reverse | 47.6 % | 1,477,822 | 615 |
| | unstranded | 53.4 % | 1,570,907 | 6,270 |
| covid x216 | forward | 8.9 % | 134,212 | 9,809 (MW455439.1 9,805) |
| | reverse | 45.9 % | 1,339,862 | 243 |
| | unstranded | 54.2 % | 1,445,625 | 10,046 |
| HHV-6B SRR20710641 (5′, 16+10, `-x 10xv2`) | forward | 6.8 % | 103,332 | 137 (HHV-6B 136) |
| | reverse | 53.3 % | 383,057 | 108 |
| | unstranded | 60.1 % | 482,652 | 243 |

## Interpretation

- **The forward default discards most 5′ data.** Host UMIs rise 3.7–11× under reverse or unstranded. The 5′ R2 read is antisense to the transcript, and kallisto sets `--fr-stranded` for 10x technologies when no strand is given.
  - Every 5′ run so far (the published covid runs, the HHV-6B benchmark, tonsil) was quantified on about 7–9 % of its reads.
- **The poly-G anellovirus artefact is forward-strand specific.** It drops 9× (x213) and 40× (x216) under reverse, and comes back under unstranded. Poly-G reads match C-rich forward-strand stretches (F-019). So for 5′ data, reverse is both the correct strand for host signal and a partial artefact suppressor.
- **Viral signal is not strictly strand-bound.** HHV-6B UMIs are highest unstranded (243, 1.8× forward). Viral transcripts, or reads, fall on both strands. A strand choice made for host signal may cost viral sensitivity.
  - This is what the R2.6 pilot and the D-experiments under the R2.3 specificity constraint must weigh.
- **Caveat.** These are 4M-read subsamples, and viral counts are small for HHV-6B.

## Action

- DEF-02 / MECH-D (R2.6): priority. Pass `--strand` to kb, and infer it for 5′ data. 5′ numbers produced before this fix are not comparable.
- Rerun the 5′ positives (HHV-6B, covid, tonsil) with the inferred strand once `--strand` exists.
