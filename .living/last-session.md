# Last session — 2026-07-05 (multimap default → host-conservative; profiling triaged; status review)

## 2026-07-05 — multimap default flip + profiling triage
- **Default multimap method changed `equal` → `host-conservative`** (commits `bb78701`, `7c8ad4d`)
  for specificity in viral detection. NOTE: this reverses the PR 17 decision that had set the
  default to `equal`; docs for both states exist — confirm this is the intended final call.
- Multimap profiling work fully triaged: PLAN.md S9 flipped `[x]`, ANALYSIS_MANIFEST entry set
  `complete`, `register_values.py` writes `numbers.json` (27 values). Redundant wall-time profiler
  is done (no longer running).
- Working tree still has uncommitted `.living/` + PLAN.md + ANALYSIS_MANIFEST + profiling outputs.
  Nothing pushed; still on `claude/multimap-memory-and-showcase`.

---

# Prior session — 2026-07-03/04 (pub-readiness + reference-strategy 2×2 + multimap speedup)

## Phase 1 — Publication-readiness reconcile + execute (committed; owner-gated remainder)
Host-response honesty fix, docs↔runtime consistency, installed-package CI, provenance. Owner-gated:
manuscript authorship, release (tag/PyPI/Zenodo). See [[decisions]].

## Phase 2 — reference-strategy 2×2 (COMPLETE 12/12) + harmonized
Restored index from archive, re-fetched HHV-6B/HSV-1 FASTQs, fixed SLURM conda-PATH + kallisto
gotchas, ran all 12 (fresh12b). Harmonization: headline collapses — only **HHV-6B ~1.9×** survives;
EBV/HSV-1 are GTF-annotation artifacts; HSV-1 STAR=0 was a regex bug (fixed + regression test).
See [[learnings]] 2026-07-04, [[findings]] reference_strategy_2x2.

## Phase 3 — multimap speedup — ✅ DONE (~4× main pass, byte-identical, 567 tests)
All committed on the branch, each gated on golden-equivalence (0.00e+00) + old-vs-new benchmark:
- **Vectorized EM** (sparse mat-vec) — dropped EM from cProfile 86%→~2%.
- **EC-precompute + array iteration** — ~15%.
- **Direct CSR buffer access** (`indptr`/`indices`/`data` + `searchsorted`) — **~2.9×**; kills the
  `_matrix_value` `__getitem__` dispatch that cProfile showed was **86%** of the pass.
- **Corrected the subagent's recommendation**: it flagged "batch fetch `original_counts[cell, gene_list]`"
  (CSR fancy index) — but that goes through the SAME `__getitem__` dispatch and was 1.6× SLOWER
  (measured). The win is bypassing `__getitem__` via raw buffers. See [[learnings]] 2026-07-04.
- Baseline for scale: full-data `equal` = ~5.5 h/method (103M rows) → ~1.4 h after the fix.
- Profiling provenance committed (cProfile outputs + scripts). Redundant wall-time run still burning
  CPU (secondary; can be killed).

## Decisions waiting on user
1. Benchmark → manuscript (HHV-6B-only + caveat, or keep excluded).
2. Push / open PR (all local on claude/multimap-memory-and-showcase; nothing pushed).
3. Kill the redundant wall-time profiler (optional; the actionable cProfile is done).

**Not pushed.** Plans: docs/superpowers/plans/*, analysis/reference_strategy_benchmark/COMPLETION_PLAN.md,
analysis/multimap_profiling/ (profiler scripts+outputs incl. fast_profile.py, interpret_cprofile.py,
cprofile_equal.txt, cprofile_em.txt — subagent-owned; will commit when run completes).
