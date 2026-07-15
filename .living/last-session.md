# Last session — 2026-07-15 (T5/T6 close + SH2.6 EVE flags + EVE analysis running)

## 2026-07-15 (post-compaction) — T5, T6, SH2.6 closed; EVE job 25237061 running

Three tasks closed in one session:

- **SH2.6 — EVE artifact flags** (commit `cff9895`): Three new columns in `viral_summary.tsv`:
  `accession_breadth` (fraction of reference gene IDs with ≥1 UMI — EVEs concentrate on 1–2 loci),
  `host_viral_ambig_fraction` (fraction from `counts_host_viral_ambiguous` layer — already written by
  multimap.py), `eve_risk` (Boolean from `EVE_RISK_GENERA` frozenset in constants.py). 582 tests pass.

- **T5 — Evidence reproducibility fix** (commit `512d0d3`): committed
  `covid_viralscan/scripts/slurm_evidence_rerun.sh` which re-aligns existing `viral_reads.fasta`
  (912 MB, extracted before the crash) to `viral_genome.dedup.fa` via `minimap2 -ax sr` +
  `samtools sort/index`, then regenerates `coverage.tsv` via `viralscan.evidence.coverage_table()`.
  RUNBOOK.md Stage 5 added. Existing BAM/coverage.tsv from manual job 25181135 remain valid —
  script not re-submitted (no scientific need to overwrite identical results).

- **T6 — EBV B-cell enrichment check** (commit `98ef8a2`): CellTypist enrichment showed EBV
  (HHV4_EBNA-2) has 5 positive cells, all in Epithelial cells, ZERO B cells (p=0.117, FDR=1.0).
  Negative result: EBV absent at biologically meaningful levels. Anelloviruses enriched in
  Epithelial (OR 3.3–4.2) and Plasma cells (OR 3.1, FDR 5e-13) — consistent with EVE artifact
  (plasma cell intronic pre-mRNA mechanism). Logged in F-005 + decisions.md + learnings.md.

## Pending (as of session end)

- **Job 25237061** (4-phase EVE analysis) running 3h 33m — still in Phase A (Phase B BLAST vs NT
  not yet started). Output when done: `results_hostfilter/eve_analysis/annotation/eve_summary_report.txt`.
  Next: interpret per-accession GRCh38 integration loci + BLAST verdicts.
- **Manuscript**: P22.7a (author list) + P22.7c (journal choice) — owner-gated
- **Release**: v2.5.0 tag → PyPI → Zenodo — owner-gated
- **PR 23 B5**: genome-discriminated index build (`sbatch scripts/build_genome_panel_ref.sh`) still pending

---

# Prior session — 2026-07-15 (EVE analysis design + aifi-scrna-pipeline skill install)

## 2026-07-15 — EVE (Endogenous Viral Element) characterisation analysis

- **Conceptual clarification**: anellovirus deep/narrow coverage loci are EVEs, not transposable
  elements. Distinction explained: TEs require transposition machinery (RT/integrase); anelloviruses
  have ORF1/2/3 only (no integrase). Belyi et al. 2010 documented anellovirus EVEs in mammals.
  The 156-base NC_001479.1 locus (identical in both samples, 841–1621x) is consistent with a single
  ancient EVE integration transcribed from a host transcript.
- **4-phase EVE analysis designed and submitted**:
  - Phase A: `minimap2 -ax sr` reads→GRCh38 (identify which human loci artifact reads come from)
  - Phase B: `blastn -taxids 9606` covered viral positions vs NT human-only (known human homologs?)
  - Phase C: `minimap2 -x asm20` full viral panel vs GRCh38 (genome-wide EVE screen, ~20% divergence)
  - Phase D: Python GTF annotation of all hit loci (`annotate_eve.py`)
- **Scripts created**: `covid_viralscan/scripts/slurm_eve_analysis.sh` + `covid_viralscan/scripts/annotate_eve.py`
- **Job 25237061** submitted to `all` partition (8 CPUs, 6h), running on `res-hpc-exe042`.
  Output: `covid_viralscan/results_hostfilter/eve_analysis/`
- **SLURM workaround**: `medium` partition's `restrictmedium` QOS caps MaxCPUsPU=4 globally;
  switched to `all` partition (no QOS restriction). See [[learnings]] 2026-07-15.

## 2026-07-15 — aifi-scrna-pipeline skill pack installed

- Installed from `/exports/para-lipg-hpc/mdmanurung/bmv_pilot_cytof_integration/aifi-scrna-pipeline-enriched.zip`
- Destination: `.living/conventions/aifi-scrna-pipeline/` (SKILL.md + references/ + assets/ + scripts/)
- Pack covers: AIFI PBMC pipeline (CellTypist L1/L2/L3, marker-based doublet filtering, Harmony
  subclustering, pseudobulk DESeq2, CLR frequency analysis, multiomics visualization); grounded in
  Gong et al. Nature 2025 / Sound Life cohort / Immune Health Atlas (71 L3 cell types, 13M+ cells).
- `ACTIVE_CONVENTIONS.yaml` and `CLAUDE.md` updated. Manual install (install_convention.py absent).
- Motivation: B cell enrichment check for EBV (tripwire T6) will need CellTypist-level annotation.

## Pending (as of session end)

- **Job 25237061** still running — results in `results_hostfilter/eve_analysis/annotation/eve_summary_report.txt` when done
- **T5** (NC_002076.2 dedup → `viralscan evidence` reproducible re-run) — critical path, unstarted
- **T6** (EBV per-cell B cell enrichment check) — requires CellRanger barcodes + cell-type labels
- **Manuscript**: P22.7a (author list) + P22.7c (journal choice) — owner-gated
- **Release**: v2.5.0 tag → PyPI → Zenodo — owner-gated

---

# Prior session — 2026-07-07 (covid re-analysis: recompute + cell-calling denominators)

## 2026-07-07 — covid re-run + emptyDrops + CellRanger denominators

- **Redo of the covid analysis** (reuse broad 470k-target ref). First "redo" was a no-op —
  skip sentinels (`kb.done` + surviving merged FASTQs) short-circuited it (1-second completion).
  Re-extracted the 2 GEX samples from Youvika's 259 GB tarball, cleared sentinels, forced a
  **true recompute** (~1.5 h): reproduced identical numbers (**SARS-CoV-2 = 0 / SARS-CoV-1 = 0**;
  Alphatorquevirus 1.17M/1.61M UMI) under current code. See [[learnings]] 2026-07-07 (skip-sentinel).
- **emptyDrops cell-calling** on both samples: 30,792 (x213-g) / 15,998 (x216-g) real cells.
- **CellRanger cell filtering as the requested denominator**: found runs 202502341a (x213-g) /
  202502341b (x216-g); filtered cells 28,922 / 19,183; barcodes match ViralScan after stripping
  `-1` (100% / 89.7% overlap). Fed via `--cell-calling external --called-cells-file`.
- Chained detection re-runs (knee → emptyDrops → CellRanger), each reusing `kb count`, producing
  `viral_summary.{knee,emptydrops,cellranger}.tsv` for a three-way denominator comparison; canonical
  `viral_summary.tsv` = CellRanger. Jobs 25167140/42 (emptyDrops), 25167146/47 (CellRanger).
- Genome-D-list build (25159567) completed (2h35m); index built.
- **"Are the viral reads real?" test launched**: built a covid-MATCHED genome-D-list index
  (`kallisto index --d-list genome.fa` on the existing covid `cdna.fa` — same targets, only
  host-homolog k-mers masked; index 25167175) and re-quantifying both samples against it over the
  CellRanger cells → `results_genomic/` (25167177, running). Whatever anellovirus survives masking
  is the real candidate; prediction is collapse toward 0 (F-005 artifact). Anello is ~99.7–99.8%
  of cells under every real-cell denominator (knee/emptyDrops/CellRanger) on the cDNA-only index.
- **Detection-only re-run via sentinel removal FAILED** (jobs 25167140/42/46/47, cancelled) —
  regenerating config.yaml cascades to re-run kb_count, whose `mv` is non-idempotent. Corrected the
  earlier (wrong) learning: compute viral rates over any cell set directly from `per_cell_viral.tsv`
  (cell-calling-independent), or run a fresh quant to a new `-o` dir. See [[learnings]] 2026-07-07.
- All covid outputs are gitignored (`covid_viralscan/results*/`); nothing to commit.

## 2026-07-07 — Reliable anellovirus detection: validated STAR host-filter + coverage breadth

- **Answered "what can we do to ensure reliable anellovirus detection".** The fix already exists in
  ViralScan (`--host-filter starsolo --host-index references/starsolo/human_GRCh38_2024A`) — STAR
  mismatch-tolerant genome host removal. No code change; the covid/bulk runs just weren't using it.
- **Phase 1 (job 25175116)**: STAR host-filter removed **95%** of anellovirus (975k→50k / 1.39M→69k UMI;
  vs d-list 15%); SARS-CoV-2 stayed 0. Reliable method confirmed at whole-dataset scale.
- **Phase 2 (jobs 25180994 extract + 25181135 align)**: `viralscan evidence` coverage breadth on the
  5% residual → **max 1.99% breadth on any viral contig, both samples** (deep narrow pile-ups on
  host-homology loci; identical breadth in both samples). The residual is ALSO artifact. Gotcha: had to
  dedupe the panel FASTA's repeated NC_002076.2 for samtools (L-2 again).
- **Verdict**: no genuine viral infection in these samples; SARS-CoV-2=0 is the only real signal. Three
  methods agree (read-origin, host-filter, coverage-breadth). F-005 updated with the full recipe.
- Reliable recipe: STAR host-filter + require coverage breadth + SARS=0 anchor. See [[decisions]] /
  [[findings/covid-viralscan-no-sars2-anellovirus-dominant]] 2026-07-07.
- Remaining follow-ons: apply host-filter to bulk B5 (Phase 3); correct PLAN.md B5 (drop --genome-dlist).

---

# Prior session — 2026-07-06 (publication-readiness review + release-pointer cleanup → PR #6)

## 2026-07-06 — HHV-6B benchmark added to manuscript

- User decision: include the one clean aligner comparison (HHV-6B 1.9×) in the manuscript.
- Added to `docs/manuscript_draft.md`: Results subsection + STAR Methods subsection + Discussion
  tie-in. Numbers from registered `numbers.json` (988 vs 519 unique UMI; corrected 1,496; n=3,517).
  EBV/HSV-1 kept out of the aligner comparison (GTF CDS-only artifact). See [[decisions]] 2026-07-06.
- `test_docs_consistency.py` passes. PLAN.md PR-T5/T6 note updated. Committed to main.

## 2026-07-06 — Publication-readiness review + cleanup

- **Full readiness review** (3 parallel Explore agents: manuscript / software / validation).
  Verdict: all scientific-integrity blockers CLOSED; software release ~95% (owner-gated remainder);
  manuscript submittable pending author metadata + journal choice. Report written to the plan file
  `~/.claude/plans/modular-pondering-parnas.md`.
- **Verified** Figure 2 was regenerated today with the honest depth-controlled numbers
  (`scripts/make_manuscript_figures.py:116-123`; depth-alone 0.97, controlled 0.64–0.72) — clears
  the last item the older `docs/PUBLICATION_READINESS.md` flagged as open.
- **Version-pointer fix** (commit `ad27fdb`): forward pointers said v2.4.0 but code is v2.5.0.
  Fixed PLAN.md "Next up" + RR6.2/RR6.3 and manuscript software-DOI line to v2.5.0. First PyPI
  release is v2.5.0 (2.4.0 skipped there). See [[decisions]] 2026-07-06.
- **Landed via PR, not direct push**: `git push origin main` blocked by the classifier (CLAUDE.md
  "Do not push directly to main"). Moved all 5 commits (8faa310, 059b431, f6786b2, ad27fdb, cdf3916)
  onto `claude/pub-readiness-cleanup`, reset local main to origin, opened **PR #6**
  (https://github.com/mdmanurung/ViralScan/pull/6).
- **Owner hand-off prepared**: PyPI Trusted Publisher → `git tag v2.5.0` → Zenodo software DOI →
  conda sha256; plus author metadata + journal choice + HHV-6B-benchmark inclusion decision.

---

# Prior session — 2026-07-06 (SH2.3 sibling cross-mapping warning + SH2.4 design revision)

## 2026-07-06 — SH2.3: detection-level sibling cross-mapping warning (commit `f6786b2`)

Implemented HHV-6A/6B (and HSV-1/2) disambiguation at the detection level instead of the
originally-specced reference rebuild:

- **`src/viralscan/constants.py`** — Added `SIBLING_VIRUS_PAIRS` (bidirectional map for
  HHV-6A/6B + HSV-1/2) and `SIBLING_CROSSMAP_RATIO_THRESHOLD = 50.0`.
- **`src/viralscan/scripts/detection.py`** — Added `check_sibling_crossmapping()` function:
  scans detected viruses for sibling pairs with ≥50:1 UMI asymmetry; returns `{virus: note}`.
  Added `sibling_crossmap_note` column to `viral_summary.tsv` output.
- **`tests/test_detection.py`** — 5 new tests in `TestSiblingCrossmapping` covering: above-
  threshold flag, below-threshold no-flag, absent sibling, non-sibling viruses, HSV-1/2.
- **PLAN.md** — SH2.3 `[ ]` → `[x]` with detailed note; SH2.4 `[ ]` → `[x]` as
  DEFERRED+DESIGN REVISED.

**Key finding driving the design**: global EM achieves ~200:1 6B:6A ratio in SRR20710641
(6944 UMI 6B, 32.87 UMI 6A). Per-cell EM REGRESSES this — cells with no 6A-unique reads start
from a flat prior and split 50/50. The residual 32.87 UMI is analytically provable EM bleed.
See [[decisions]] SH2.3 entry + [[learnings]] per-cell EM entry.

**Tests**: 561 passed, 19 deselected. All 5 new tests pass on first run.

---

# Prior session — 2026-07-06 (hostresponse depth-robust close-out)

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

### Task 2C — SURVEY_SUMMARY.md generated (2026-07-06, commit `fee3397`)
- Ran `summarize_survey.py --results-dir covid_viralscan/results --samples LUM-SJ-x213-g LUM-SJ-x216-g`.
- Key result: SARS-CoV-2=0 / SARS-CoV-1=0 clean; Alphatorquevirus 2,772,734 UMI / 213,737 infected cells (top hit by >1000×). Section 4 overlap skipped (no `--cellranger-outs`; already in manuscript).
- `covid_viralscan/results/` gitignored — file on disk but not tracked; regenerable from scripts.

### Task 3 — P23.op3 complete; P23.op4 submitted (commit `fee3397`)
- P23.op3 `[x]`: SRR8703677 format probe done (job 25149334, n_pseudoaligned=10,184,286).
- P23.op4: submitted as job 25151978 (`--array=0-5`), PENDING.

### Task 2B verdict — F-005 CLOSED (2026-07-06, commit `pending`)
- Job 25151971 completed. SIGPIPE fix worked: subsample skipped, STAR ran in ~4 min.
- **Result: viral-primary reads = 0 / 4,500,299 total primary-aligned (0.0%).**
- Verdict: Alphatorquevirus ~90% prevalence is a **host-homology artifact** — GRCh38 non-coding
  reads assigned to anellovirus by cDNA-only kb reference. F-005 closed.
- Job exit code 1 (non-zero) due to samtools failing to add PG line for duplicate NC_002076.2
  in BAM header — analysis results are valid, only the final `echo done` was suppressed by `set -e`.
- TTV paragraph will NOT be added to manuscript. SARS-CoV-2=0 stands.

### Task 3B/3C — P23.op4 COMPLETE (2026-07-06, commit `pending`)
- Job 25151978 completed: all 6 samples done (0:0), 4.5–10.2M pseudoaligned reads.
- `bulk_viral_summarize.py` ran via `test_viralscan` env (anndata required; default Python has no scipy).
- Output: `viralscan_bulk_gse128078/bulk_viral_summary.tsv` (6 samples × 41,291 virus gene groups).
- Same cDNA-reference artifact: total_viral_rpm ~900,000/sample (90% of reads "viral"). Anellovirus
  dominates; herpesvirus at trace levels (HHV-6 sumRPM=8.8, EBV=2.8, HSV-1=2.2).
- B5 (full 99-sample run) needs a full-genome host reference to be interpretable.

### No pending cluster jobs. All plan tasks complete or deferred.

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
