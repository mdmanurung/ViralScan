# Last session — 2026-07-06 (hostresponse depth-robust close-out)

## 2026-07-06 — Task 1: close-out depth-robust hostresponse module (commit `a1ad4a3`)

Three loose ends from `todo/hostresponse-depth-robust-module.md` folded into a single commit:

- **(a) `evalue_flag` column** — `_per_gene_evalues()` now classifies each stable gene as
  `fragile` / `moderate` / `robust` (E < 1.5 / 1.5–3 / ≥ 3) and writes it to
  `<virus>_depth_diagnostics.csv`.
- **(b) `_panel_depth_adjusted_auc()`** — new helper mirroring `_depth_alone_auc()` that
  augments the stable-gene panel with `log1p(host_depth)` as a covariate, using the identical
  balanced split. Emits `model_auc_depth_adjusted_mean/sd` in `hostresponse_metrics.csv`.
  Decision: report panel+depth (not panel-alone) — see [[decisions]] 2026-07-06.
- **(c) RunConfig wiring** — `hostresponse_label`, `hostresponse_depth_match`,
  `hostresponse_control_mito`, `hostresponse_differential` added to `DEFAULTS`, `RunConfig`
  dataclass, `from_snakemake_config()`, main pipeline parser (`menu.py`), pipeline config dict,
  and Snakemake entry block. Previously CLI-only (only reachable via `viralscan hostresponse`).
- **9 new tests** (5 evalue_flag, 3 depth_adj_auc, 1 RunConfig round-trip); **577 total pass**.
- `todo/hostresponse-depth-robust-module.md` → `Status: done (2026-07-06)`.
- PLAN.md SH1.1 close-out note added; CHANGELOG.md [Unreleased] section updated.
- Gotcha: synthetic depth-proxy test needed `!= "robust"` not `== "fragile"` — see [[learnings]] 2026-07-06.

## 2026-07-06 (continued) — Tasks 2/3: cluster jobs submitted + manuscript covid section written

### Task 2A/2B — STARsolo re-run + TTV read-origin (submitted)
- `diag_viral_read_origin.sh` GENOME repointed to pre-built `combined_GRCh38_2024A_serratus_plus_anellovirus` index, decoupling 2B from the 30-min 2A build step (commit `7191982`).
- Four jobs submitted concurrently:
  - **25149332** — STARsolo build step (`slurm_starsolo_covid.sh build`)
  - **25149333** — TTV read-origin NH-flag test (`diag_viral_read_origin.sh`)
  - **25149334** — Bulk GSE128078 format probe (`--array=0-0 scripts/bulk_viral_scan.sh`)
  - **25149335** — STARsolo sample array (`--array=0-1 --dependency=afterok:25149332`)

### Task 2D — Manuscript covid specificity section (commit `76b769c`)
- Added **Results** subsection: SARS-CoV-2=0 in both covid-era samples + STARsolo confirmation; cell-calling concordance table (ViralScan emptyDrops 30,849 / CellRanger 28,922 / STARsolo EmptyDrops_CR 19,920); GEM-X whitelist mismatch documented.
- Added **Methods** subsection: STARsolo combined-ref run parameters for the cross-check.
- **TTV ~90% deliberately omitted** — held pending job 25149333's NH-flag verdict (F-005 still "under review"). Decision logged in [[decisions]] 2026-07-06.

### Pending (cluster output required)
- **Task 2B verdict**: read log `covid_viralscan/logs/readorigin_25149333.log` — if NH==1 fraction >80%, add TTV paragraph; if NH>1 dominates, close F-005 as host-homology artifact.
- **Task 2C**: feed STARsolo `barcodes.tsv` (job 25149335 output) as `--called-cells-file`; regenerate `covid_viralscan/results/SURVEY_SUMMARY.md`.
- **Task 3A/3B/3C**: inspect bulk probe output (25149334) → submit `--array=0-5` pilot → run `bulk_viral_summarize.py` → flip PLAN P23.op3/op4 `[x]`.

---

# Prior session — 2026-07-05 (multimap default → host-conservative; profiling triaged; status review)

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
