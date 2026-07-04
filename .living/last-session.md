# Last session — 2026-07-03/04 (pub-readiness + reference-strategy 2×2 + multimap profiling)

## Phase 1 — Publication-readiness reconcile + execute (committed; owner-gated remainder)
Host-response honesty fix, docs↔runtime consistency, installed-package CI, provenance. Owner-gated:
manuscript authorship, release (tag/PyPI/Zenodo). See [[decisions]].

## Phase 2 — reference-strategy 2×2 (COMPLETE 12/12) + harmonized
Restored index from archive, re-fetched HHV-6B/HSV-1 FASTQs, fixed SLURM conda-PATH + kallisto
version gotchas, ran all 12 (run fresh12b). Fair-comparison harmonization
(`analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py`): headline collapses — only
**HHV-6B ~1.9×** survives; EBV/HSV-1 are GTF-annotation artifacts; HSV-1 STAR=0 was a regex bug
(fixed + regression test). See [[learnings]] 2026-07-04, [[findings]] reference_strategy_2x2.

## Phase 3 — multimap profiling (RUNNING, detached; ~00:25)
EBV full-depth: 103M BUS records, 81.7% multi-gene ECs; main-pass `itertuples` ≈ 42 min/method.
Recommendation firm: vectorize main pass + EM into sparse EC×gene ops. Final per-method table pending.

## Decisions waiting on user
1. Benchmark → manuscript (HHV-6B-only + caveat, or keep excluded).
2. Push / open PR (all local on claude/multimap-memory-and-showcase; nothing pushed).
3. Implement the multimap speedup (separate from profiling).

**Not pushed.** Plans: docs/superpowers/plans/*, analysis/reference_strategy_benchmark/COMPLETION_PLAN.md.
