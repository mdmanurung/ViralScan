# Finding: reference-strategy 2×2 benchmark (complete, 2026-07-04)

Run `fresh12b` — 12/12 complete. Target-virus UMI (Selectivity Index = on-target ÷ off-target):

| virus | STAR·comb | STAR·2step | VS·comb | VS·2step |
|-------|-----------|------------|---------|----------|
| EBV    | 79,989 | 80,079 | 1,479,894 (SI 1.5M) | 1,434,619 (SI 1.4M) |
| HHV-6B | 1,876  | 1,876  | 6,944 (SI 211)      | 6,928 (SI 211) |
| HSV-1  | 0      | 0      | 49,937 (SI 328)     | 45,952 (SI 326) |

- **Aligner axis dominates**: ViralScan >> STARsolo (EBV ~18×, HHV-6B ~3.7×, HSV-1 STAR=0 vs VS ~50k).
- **Reference strategy minor**: combined ≈ two_step within each aligner (smaller gap than the
  manuscript's 4× 1M-subsample headline — reconcile).
- **Confounds (block manuscript claims)**: count-layer (STAR unique int vs VS multimap fractional),
  denominator (`*_fixed` over different barcode universes, not the shared anchor 1908/3517/3307),
  off-target = cross-mapping, and HSV-1 STAR=0 is surprising (index contains HSV-1 — verify).

Canonical: `results/reference_strategy_benchmark.tsv`. Not manuscript-ready until Step-6 harmonization.
