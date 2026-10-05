# F-011 — The analytic k-mer model overstates anellovirus detection: P=1.0000 claimed, 0.51 measured

**Date**: 2026-09-27
**Status**: correction recorded; the affected claim is `REF-01`'s parenthetical, not its decision
**Found by**: `CAT-11`, promoting the k-mer capture measurement to a reproducible tool

## What was claimed

`PLAN.md` `REF-01` records, as the evidence for making the expanded anellovirus
panel the default:

> Leave-one-out capture for a *novel* strain: bundled panel 0.04 %
> (P(90 bp fragment captured)=0.024) vs expanded panel 20.5 % (P=1.0000).

## What is measured

`scripts/measure_kmer_capture.py` reproduces every *coverage* number in `REF-01`
exactly, which is how we know the tool is sound:

| Quantity | `REF-01` | Measured |
|---|---|---|
| EBV `NC_007605.1` positive control | 100.0 % over 144,283 31-mers | 144,283 31-mers, coverage 1.0 |
| Bundled 20-genome panel, median coverage | 0.00 % | 0.00 % |
| Bundled panel, zero-coverage genomes | 85.8 % | 85.80 % |
| Betatorquevirus zero-coverage (n=1,542) | 98.4 % | 98.4 % |
| Gamma/Samek/Het/Gyro/Mem zero-coverage | 100 % | 100 % |
| Expanded panel, leave-one-out capture | 20.5 % | 20.55 % |

The one number that does **not** reproduce is the fragment-capture probability:

| | Expanded panel |
|---|---|
| `REF-01` claim | P = 1.0000 |
| Measured, walking real 90 bp windows | **P = 0.5097** |

## Why the two disagree

The claimed value follows the independence model: a 90 bp read holds
`90 - 31 + 1 = 60` k-mer positions, so with per-k-mer capture `c = 0.205`,
`P = 1 - (1 - c)^60 = 1.0000`.

That model assumes shared k-mers are scattered independently along the genome.
They are not — they cluster in conserved blocks, so a read either lands in a
conserved block and matches many k-mers, or lands outside one and matches none.
Reproducing the model against the measurement, per genus:

| Genus | median LOO | analytic P | measured P | overstated by |
|---|---|---|---|---|
| Alphatorquevirus | 0.412 | 1.0000 | 0.788 | 0.212 |
| Betatorquevirus | 0.192 | 1.0000 | 0.498 | 0.502 |
| Gammatorquevirus | 0.198 | 1.0000 | 0.332 | 0.668 |
| Memtorquevirus | 0.098 | 0.9979 | 0.247 | 0.751 |
| Gyrovirus | 0.009 | 0.4137 | 0.064 | 0.350 |

## Why it matters

- **It is corroborated independently.** The 2026-09-27 held-out plant
  (`TONSIL-01`) recovered 73–88 % of reads for a strain with a ~94 % identical
  neighbour and 8–40 % for divergent genera. That is consistent with ~0.5, not
  with 1.0.
- **`REF-01`'s decision still stands.** The expanded panel beats the bundled one
  by every measure (zero-coverage 85.8 % → 0.15 %). Only the "P=1.0000"
  parenthetical is wrong, and it was the optimistic bound rather than the
  reason for the flip.
- **It changes the diversity target.** At P ≈ 0.51 the panel misses about half
  the reads from a strain it does not contain, so "2,042 genomes" is not the
  finish line for anellovirus coverage. `CAT-12` sets the target from this
  measured number.
- **`sensitivity.fragment_capture` is not at fault.** That function models
  *divergence*-driven loss and documents itself as an upper bound on the loss.
  The error is in applying an independence assumption to *observed k-mer
  sharing*, where the positions are correlated.

## Fix

Use the measured `p_fragment` from `scripts/measure_kmer_capture.py`, never the
analytic value, whenever the question is "would we see a strain that is not in
the panel". `CAT-12`'s gate is written against the measured number.
