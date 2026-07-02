# reference_strategy_benchmark

Benchmark comparing viral-detection across **reference strategies** (combined vs
two_step) and **aligners** (ViralScan/kallisto vs STARsolo) on three SRA runs, one per
target virus. Source: `results/reference_strategy_benchmark.tsv` (produced by
`scripts/{prepare,summarize,audit}_reference_strategy.py` + SLURM array). Dataset:
`reference_strategy_refs` / `reference_strategy_fastqs`.

## Status: INCOMPLETE — 4 of 12 rows complete (2026-07-02)

| status | n |
|--------|---|
| complete | 4 |
| blocked | 3 |
| incomplete | 2 |
| running_or_incomplete | 2 |
| failed | 1 |

Only 4 of the planned 12 (design was 18) condition-rows finished. **All four EBV rows
failed/blocked/incomplete**, and every `two_step` viralscan row is blocked. This means
the intended **reference-strategy Selectivity Index** (SI = on-target EBV UMI ÷
off-target HHV-6B+HSV-1 UMI, ranked across aligner×reference) — todo idea 3b — is **not
computable yet**: EBV (the on-target) has no complete rows, and HSV-1 has zero signal
(SI = 0/0). The SI / noisy-channel analyses are blocked pending benchmark completion.

## What the 4 complete rows DO support (two valid comparisons)

Both share the denominator (shared-anchor cells), so counts are comparable within each.

### C1 — HHV-6B (SRR20710641, combined): ViralScan vs STARsolo (aligner)
| aligner | count layer | HHV-6B+ cells | target UMI | off-target cells | off-target UMI |
|---------|-------------|---------------|-----------|------------------|----------------|
| STARsolo | GeneFull.raw | 1295 | 1876 | 0 | 0 |
| ViralScan | per_cell_viral.viral_umi | 1999 | 3217 | 11 | 8 |

**Observation**: ViralScan calls **~54% more HHV-6B-positive cells** (1999 vs 1295) than
STARsolo on the same 3528 anchor cells, but shows a small off-target signal (11 cells /
8 UMI) where STARsolo shows none.
**Caveats**: (i) the two aligners report **different count layers** (`GeneFull.raw` vs
`per_cell_viral.viral_umi`), so this is a detection-sensitivity comparison, not a clean
UMI-for-UMI one; (ii) the extra ViralScan-positive cells and the off-target signal may
partly reflect the known **HHV-6A vs HHV-6B contig ambiguity** (6A cross-mapping), not
pure sensitivity. n=1 condition — preliminary.

### C2 — HSV-1 (SRR8315713, STARsolo): combined vs two_step (reference strategy)
| strategy | HSV-1+ cells | target UMI | off-target |
|----------|--------------|-----------|-----------|
| combined | 0 | 0 | 0 |
| two_step | 0 | 0 | 0 |

**Observation**: HSV-1 is **not detected** in this sample under either strategy — zero
target signal either way, so the two strategies cannot be distinguished here. Consistent
with HSV-1 being absent / the previously documented HSV-1 denominator artifact (there is
no true HSV-1 signal to recover).

## Open questions / next steps
- **Complete the benchmark** (especially the EBV and viralscan/two_step rows) before the
  Selectivity-Index and noisy-channel analyses (todo ideas 3b, 7b) can run.
- Resolve C1's count-layer mismatch (compare on a common layer) and separate true
  HHV-6B sensitivity from HHV-6A cross-mapping before claiming ViralScan is more sensitive.
