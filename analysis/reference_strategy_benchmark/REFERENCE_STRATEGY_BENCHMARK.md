# reference_strategy_benchmark

Benchmark comparing viral-detection across **reference strategies** (combined vs
two_step) and **aligners** (ViralScan/kallisto vs STARsolo) on three SRA runs, one per
target virus. Source: `results/reference_strategy_benchmark.tsv` (produced by
`scripts/{prepare,summarize,audit}_reference_strategy.py` + SLURM array). Dataset:
`reference_strategy_refs` / `reference_strategy_fastqs`.

## Status: COMPLETE — 12 of 12 rows (2026-07-04, run `fresh12b`)

Completed after restoring the scratch-deleted ViralScan index from archive, re-fetching the
scratch-deleted HHV-6B/HSV-1 FASTQs from ENA, and fixing a SLURM `conda activate` PATH-shadowing
bug (see `COMPLETION_PLAN.md`). Canonical results: `results/reference_strategy_benchmark.tsv`.

### Completed 2×2 — target-virus UMI (Selectivity Index = on-target ÷ off-target UMI)

| virus | STARsolo combined | STARsolo two_step | ViralScan combined | ViralScan two_step |
|-------|-------------------|-------------------|--------------------|--------------------|
| EBV    | 79,989 (off 0) | 80,079 (off 0) | 1,479,894 (off 1; SI 1.5M) | 1,434,619 (off 1; SI 1.4M) |
| HHV-6B | 1,876 (off 0)  | 1,876 (off 0)  | 6,944 (off 33; SI 211)     | 6,928 (off 33; SI 211) |
| HSV-1  | 0 (off 0)      | 0 (off 0)      | 49,937 (off 152; SI 328)   | 45,952 (off 141; SI 326) |

**Two findings, both with caveats:**
1. **Aligner axis dominates**: ViralScan recovers far more viral UMI than STARsolo — ~18× (EBV),
   ~3.7× (HHV-6B), and HSV-1 detected *only* by ViralScan (STARsolo = 0 under both strategies).
2. **Reference-strategy axis is minor**: within each aligner, combined ≈ two_step (ViralScan
   combined marginally higher — keeps host-virus ambiguous reads). This is a *smaller* combined-vs-
   two_step gap than the manuscript's 1M-subsample headline (12,255 vs 3,096 = 4×) — reconcile.

**NOT yet a fair comparison (blocks manuscript claims — see COMPLETION_PLAN Step 6):**
- **Count-layer mismatch**: STARsolo = `GeneFull.raw` unique *integer* counts; ViralScan =
  `per_cell_viral.viral_umi` multimap-corrected *fractional* UMI. Part of the aligner gap is the
  count model, not sensitivity.
- **Denominator mismatch**: the `*_fixed` UMI/cell columns are over different barcode universes
  (STARsolo filtered cells vs ViralScan all-barcode), not the shared anchor (EBV 1,908 / HHV-6B
  3,517 / HSV-1 3,307). Re-compute on the shared anchor with a harmonized layer before claiming magnitudes.
- **Off-target = cross-mapping** (HHV-6B→HHV-6A 33 UMI; HSV-1 152) — verify it's not inflating SI.
- **HSV-1 STARsolo = 0** is surprising (the combined index contains HSV-1); confirm STARsolo isn't
  silently dropping the HSV-1 contig before treating "ViralScan detects, STARsolo misses" as real.



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

## Publication Use

**Update 2026-07-04:** the 12-row design is now **complete** (run `fresh12b`;
`summarize_reference_strategy.py` exits 0 with 12/12 complete). Completeness is no longer the
blocker. Re-inclusion in the manuscript is now gated on the **fair-comparison harmonization**
(count-layer parity, shared-anchor denominator, HHV-6A/off-target handling, and confirming the
surprising HSV-1 STARsolo=0) documented in the "Completed 2×2" section above and in
`COMPLETION_PLAN.md` Step 6 — not on whether the rows ran.

## Provenance

Clean-clone-visible provenance for the 2026-06-28 attempt (command manifest, FASTQ
and reference audits, and a per-row `failure_summary.tsv`) is tracked under
[`run_packets/2026-06-28_fresh12/`](run_packets/2026-06-28_fresh12/README.md). The
bulky run tree remains gitignored under `benchmark_runs/`.

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
