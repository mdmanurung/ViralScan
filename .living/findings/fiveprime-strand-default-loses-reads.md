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

## Update 2026-10-03 — full-run reruns (25695057_4-5 HHV-6B; 25695488_0-3 covid, resumed after an Rscript PATH fix)
- Code: d01043c (pinned worktree). Each library keeps its baseline index (covid: `viralscan_ref`, unmasked; HHV-6B: v1 full panel). Cell calling: emptyDrops.
- Baselines: covid = `results_v3_exploratory` (no `--strand`, knee). HHV-6B = `combined_corrected_5p` (no `--strand`, knee).

| library | strand | pseudoaligned | called cells | top viral rows (molecules, % called cells) |
|---|---|---:|---:|---|
| covid x213 | none (baseline) | 6.4 % | 118,061 (knee) | Alphatorquevirus 1,083,687 (62.6 %), all on MW455439.1_gene1 (poly-G, F-019) |
| | reverse | 47.8 % | 30,711 | **Gammatorquevirus 92,216 (64.9 %)**; Alphatorquevirus 32,466 (45.1 %) |
| | unstranded | 53.6 % | 45,220 | Alphatorquevirus 494,564 (91.0 %); Gammatorquevirus 88,639 (54.7 %) |
| covid x216 | none (baseline) | 8.9 % | 99,109 (knee) | Alphatorquevirus 1,499,051 (69.0 %) |
| | reverse | 46.2 % | 15,895 | Alphatorquevirus 22,826 (46.4 %); Gammatorquevirus 14,262 (39.4 %) |
| | unstranded | 54.4 % | 23,780 | Alphatorquevirus 776,689 (95.3 %); Gammatorquevirus 13,203 (28.0 %) |
| HHV-6B SRR20710641 | none (baseline) | **60.2 %** | 15,800 (knee) | HHV-6B 5,596 |
| | reverse | 53.4 % | 6,412 | HHV-6B 2,405 (0.72 %) |
| | unstranded | 60.2 % | 8,846 | HHV-6B 5,596 (0.88 %) |

- **The 4M-read pattern holds at full depth for pseudoalignment.** Reverse and unstranded recover 46–54 % of reads against 6–9 % for forward on covid (10xv3). The poly-G MW455439.1 artefact falls 33× (x213) and 66× (x216) under reverse, and comes back under unstranded.
- **Correction: the HHV-6B baseline was not run forward.** Its `kallisto bus` call has no strand flag (kallisto 0.52.0, `viralscan_bench`, `-x 10xv2`). It pseudoaligns 60.2 % with exactly the unstranded molecule count (5,596), so it behaved as unstranded. The covid baseline call also has no strand flag (kallisto 0.51.1, `test_viralscan`, `-x 10xv3`) and behaved as forward (6.4 %). Version and chemistry both differ between the two, so which one sets the default is **not separated**. F-020's HHV-6B "forward" row was an explicit `--strand forward`. The ~7–9 % loss is confirmed for the covid baselines only. Tonsil's env was not checked.
- **HHV-6B is strand-sensitive in the viral direction.** Reverse keeps 43 % of the unstranded molecules (2,405 / 5,596), so a reverse-only default would cost HHV-6B sensitivity. This is the R2.6 trade-off, now at full depth.
- **NEW, unvalidated: reverse strand surfaces a large Gammatorquevirus signal on covid.**
  - x213 reverse puts 55,111 molecules on AB303552.1_gene1 and 24,666 on AB303557.1_gene1 (forward: 180 and ~0). Both are TTMDV whole-genome placeholder models. Unstranded keeps it (54,110 / 24,452).
  - Both genomes are already dust-masked (26–28 N), have no homopolymer over 12 nt, and have GC 0.43. So this is **not** the F-019/F-021 homopolymer mechanism.
  - It is sense-strand for a 5′ R2 read, which a real transcript would also be. But 64.9 % of called PBMCs is biologically implausible for a cellular anellovirus transcript (F-022: the tonsil minimap2 screen found 0).
  - **Not validated. Do not report it as infection.** It needs the F-019-style read check (competitive minimap2 against GRCh38 + panel, complexity, position profile, CB-UMI spread). It also needs a rerun on cat42b, whose index the covid runs do not use.
- **Called cells move with strand** (x213: 30,711 reverse vs 45,220 unstranded under emptyDrops), because per-barcode UMI totals change. Viral % of called cells is therefore not comparable across strand modes.

### `--strand auto` check (DEF-02)
Ratios of pseudoalignment rate against unstranded:
- reverse / unstranded: 0.89 (x213 full), 0.85 (x216 full), 0.89 (HHV-6B full; 0.89 at 4M).
- forward / unstranded: 0.12, 0.16, 0.11 (4M).

A fixed τ = 0.8 picks **reverse** for all three 5′ libraries. For a 3′ library it should pick forward, but no 3′ ratio has been measured yet. Caveat: for HHV-6B, reverse is the choice that loses 57 % of viral molecules.
