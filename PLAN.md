# ViralScan Improvement Plan — Implementation Tracker

This file is the single source of truth for what's done and what still needs
doing. **After every implementation step, tick the relevant checkbox and update
the "Next up" pointer.** Do not mark an item done until the tests pass.

Second-pass audit completed 2026-05-08. All prior PR claims re-verified against
the actual codebase; status corrected where PLAN and code diverged.

Branch: `claude/multimap-memory-and-showcase`
Test command: `PYTHONPATH=src /exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python -m pytest tests/ -q` → 470 passed, 15 deselected (2026-06-25). System python lacks yaml/scipy; use conda env python.

---

## Status legend

- `[x]` — done and verified against the codebase
- `[~]` — partially done / obsolete sub-item
- `[ ]` — not started
- `[!]` — blocked / needs external action

## Next up

→ **Release v2.4.0** — code/tests/gates all green; Phase 6 (RR6.1–6.5) is USER-GATED
  (PR→main, tag→PyPI+ghcr, Zenodo software DOI, bioconda PR). See "Release Readiness".
→ **v2.5 Scientific-Hardening** — Tier 1 (SH1.1–1.5) COMPLETE (2026-07-03): depth-confound
  diagnostics, depth-independent label + depth-matched design, %mito control, whitelist
  preflight, called-cell denominators. Tier 2 tractable done: SH2.1 gene symbols, SH2.2
  genome-wide differential, SH2.5 bulk-claim fix. Remaining SH2.3/2.4 and SH3.1/3.2 are
  DEFERRED with rationale (multi-day, dedicated PRs). See the "v2.5 Scientific-Hardening" section.
→ **PR 23 — Anellovirus into standard combined reference** — code complete (2026-06-24);
  cluster build of anello-augmented `panel.idx` + bulk GSE128078 pilot scan are the remaining
  operational steps (see PR 23 section and "Bulk exploratory scan" below).
→ **Anellovirus reference expansion** — A/B/C/D/E complete + post-review polish applied.
  C.6 (manual Zenodo rebuild with FASTA) deferred until next release. All code/tests done in PR 20.
→ PR 15 Run-context refactor — COMPLETE. S0–S6 showcase findings — all `[x]`.
→ **PR 16 clean-code review Tier 1+2** — bugs and fail-fast hardening — COMPLETE (2026-06-22).
→ **PR 17 rerun-multimap checkpoint** — default changed to `equal`; `viralscan rerun-multimap` added — COMPLETE (2026-06-22).
→ **PR 18 Tier 3 clean-code** — config-key deduplication, main() decomposition, host_filter migration — COMPLETE (2026-06-22).
→ **PR 19 Tier 4 tidy-ups** — EM epsilon guard, inline imports, dead build_multimap_matrix dropped, np.where hoisted, except Exception narrowed — COMPLETE (2026-06-22).
→ **PR 21 hostresponse module** — Luebbert et al. 2026 approach (L2 logistic regression + randomized Lasso stability selection) — COMPLETE (2026-06-23).
→ **PR 21 docs (P21.11)** — user-facing docs for `hostresponse`, `evidence`, `rerun-multimap` — COMPLETE (2026-06-23).
→ **PR 22 publication-readiness** — P22.4 COMPLETE (2026-06-25). P22.6 STARsolo COMPLETE. P22.10 matched-barcode COMPLETE (2026-06-25). P22.5 HSV-1 divergence RESOLVED (2026-06-25). Next: §3.4 host-response numbers; Figure 1–2.
→ **fix(build-ref): Ensembl current_gtf 404** — COMPLETE (2026-07-01). `current_gtf/` symlink removed from Ensembl; switched both URL templates to `release-{N}/fasta/` and `release-{N}/gtf/`; added `_ensembl_release()` helper; added retry to `_list_ensembl_files()`. Unblocks covid_viralscan Stage 2.
→ **covid_viralscan analysis** — Stages 2–4 **COMPLETE + corrected** (2026-07-02/03). Reference build fixed (cDNA-level GTF; dedup guard for dup accession NC_002076.2). Quant hit 3 more bugs, all fixed: SIGPIPE sanity line; gzipped whitelist unreadable by bustools; and — the big one — **wrong 10x whitelist** (GEM-X-5′ chemistry; bundled v3 matched 0.4% of barcodes → all-empty-droplet matrix), fixed by extracting CellRanger's raw barcode universe (job 25140008). **Validated**: 100% CellRanger barcode overlap; emptyDrops 30,849 cells vs CellRanger 28,922 (81% overlap). **Result**: SARS-CoV-2 = 0 (robust); Torque teno ~90% of real cells at ≥5 UMI (vs 37.6% over all barcodes). See `covid_viralscan/RUNBOOK.md` + finding F-005.

→ **Cell-calling + report-both-denominators** — **DONE** (2026-07-03). `cellcalling.py` (external CellRanger/STARsolo list | emptydrops via `emptydrops.R`+DropletUtils | knee | none); `detection.py` viral_summary reports over BOTH called cells (primary) and all barcodes. Fixes the recurring empty-droplet trap. **CLI/config wiring DONE** (2026-07-03): `--cell-calling {knee,emptydrops,external,none}` + `--called-cells-file` in `menu.py`; `RunConfig` fields + `DEFAULTS`; flows through YAML to `detection.py`. Suite 557 passed.

→ **STARsolo combined host+viral (CellRanger-style) on covid** — RUNNING (`slurm_starsolo_covid.sh`; index build 25140485 → run array 25140486). Combined GRCh38+viral+SARS-CoV-2 STAR index, GEM-X chemistry, EmptyDrops_CR cells → cross-check + external called-cell list. Tradeoff: STARsolo counts unique reads only (loses ViralScan's multimap signal; cf. P22.10 STAR 77% vs VS 94% EBV≥1).

---

## Publication readiness (v2.5 release + manuscript honesty) — 2026-07-03

Executed from the reconciled plan
`docs/superpowers/plans/2026-07-03-publication-readiness-99-reconciled.md` on
`claude/multimap-memory-and-showcase`. Original plan Task 1 (package `emptydrops.R`)
and the version cut to 2.5.0 were already done at HEAD and skipped.

- `[x]` **PR-T2 — Installed-package CI + release gate.** `ci.yml` `test`/`integration`
  jobs now `pip install --no-deps -e .` (verified it builds locally; skips the snakemake
  `connection_pool` transitive-dep problem) and drop `PYTHONPATH: src`; smoke test uses the
  installed console script + subcommands. Added an `environment-file` job validating
  `environment.yml`. `release.yml` install-tests the built wheel (asserts `emptydrops.R`)
  before publish. Both workflows parse.
- `[x]` **PR-T3 — Docs/runtime consistency.** New `tests/test_docs_consistency.py` (5 passing).
  Multimap default `host-conservative` → `equal` in quickstart/cli_reference/output_reference
  (recommend `host-conservative` for cross-homology; `em` added to method list); `2.3.0` → `2.5.0`
  container examples; build-ref anellovirus flags + `is_called_cell` + `include_anellovirus`
  documented; `evidence` + `check-whitelist` added to root help; CHANGELOG compare links fixed.
  Dropped the original plan's hallucinated `cell_type_enrichment` API rewrite (api.md was already correct).
- `[x]` **PR-T4 — Distribution completeness.** `MANIFEST.in` refined to keep the unpublished
  manuscript + internal ledgers out of the sdist; conda-recipe sha256 procedure comment already present.
- `[x]` **PR-T5/T6 — Benchmark provenance + exclusion.** Tracked
  `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/` (README, commands,
  audits, truthful per-row `failure_summary.tsv`); added a "Publication Use" exclusion note.
  Six `scripts/*reference_strategy*` helpers committed.
- `[x]` **PR-T7 — Host-response honesty (the integrity gate).** `docs/manuscript_draft.md` §3.4
  now reports the tracked **0.866** headline with same-design depth-alone **0.967** and
  depth-controlled **0.636/0.718** (honest band ~0.64–0.72); heading + Discussion updated.
  Fixed the original plan's cross-design pairing error (it paired the stale 0.845 with 0.967).
  Figure 2 footer caveat added and figure regenerated.
- `[x]` **PR-T8B/T9 — Scope + truth-panel framing.** Manuscript Discussion narrowed (no
  dedicated-tool head-to-head; FPR/FNR not characterized); scaffolds under
  `analysis/dedicated_tool_comparison/` and `analysis/truth_panel/`.
- `[x]` **PR-T10 — Citation hygiene.** VIRTUS2 resolved to a software/repository citation
  (no separate journal article exists). **Author list/affiliations/declarations remain
  owner-gated** (not fabricated).
- `[!]` **Release (Task 12) — owner-gated.** PR→main, tag `v2.5.0`, PyPI/GHCR publish,
  bioconda sha256, Zenodo DOI. Not done here.

---

## Showcase-session findings — 2026-06-21 (full-depth public-data validation)

Building a functionality showcase (`docs/showcase_runbook.md`) and benchmarking against published
studies (`BENCHMARK_COMPARISON.md`) on real public scRNA-seq surfaced one finished optimization and
two real host-filter bugs. Validated on a SLURM full run (HHV-6/EBV/HSV-1 + local skin) against a
combined human+viral index; combined approach reproduces Lareau HHV-6 reactivation (12.6% infected)
and is ~4× more sensitive than the two-step host-first alternative (recovers ambiguous multimapper
viral reads the two-step discards).

- `[x]` **S0 — `multimap.py` memory optimization.** The multimap step OOM-killed full-depth samples
  at 48 GB. Root cause: `pd.read_csv` of `output.bus.txt` loaded the unused `umi` column as tens of
  millions of Python strings, plus a redundant `Int64` recast and 8 defensive sparse-matrix
  `.copy()` in `final_results`. Fix: drop `umi` via `usecols`, categorical barcodes + int32, strip
  on category labels, assign layers without copying. **Validated: full HHV-6 (783k cells) 48 GB-OOM
  → 4.3 GB MaxRSS (~10× reduction); 325 tests pass, behavior unchanged.** Committed `99db0e8`.
- `[x]` **S1 — host-filter non-10x geometry bug.** Closed by S6 (same fix). `host_filter.py` now
  delegates to `evidence.cb_umi_geometry` which handles DROPSEQ (12,8) and raises `ValueError` on
  unknown tech instead of silently falling back to (16,12).
- `[x]` **S2 — `--host-filter` halts after host_filter (feature broken end-to-end).** Root cause:
  conditional `rule host_filter` was defined *before* `rule all` in the Snakefile; Snakemake used it
  as the default target when `host_index` was set. Fix: hoisted `rule all` to always be first;
  `_kb_count_inputs()` now lists filtered FASTQ files (not just the sentinel) as explicit DAG edges
  so the rest of the pipeline runs. Regression guard: `tests/test_snakefile_dag.py` (`TestRuleOrdering`
  asserts first rule == "all"; `TestHostFilterDag` dry-runs with host_index set and asserts all
  rules through `umap` appear). Task 4 re-verified end-to-end.
- `[x]` **S3 — Showcase + benchmark deliverables.** `docs/showcase_runbook.md` (kb-python combined
  workflow, dry-run-validated chemistries 10xv3/10xv2/DROPSEQ) and `BENCHMARK_COMPARISON.md`
  (published-study comparison + combined-vs-two-step; STARsolo + kallisto two-step both = 3096 EBV
  UMI, confirming combined is ~4x more sensitive). Both files committed and dry-run-validated on
  `claude/multimap-memory-and-showcase`. All showcase-session findings S0–S6 now `[x]`.
- `[x]` **S4 — EM multimapper resolution (`--multimap-method em`).** Iterated EM over BUS
  equivalence classes at the gene level (RSEM/kallisto-style; `unique-weighted` is its first
  E-step). bustools count has no single-cell gene-level EM (only `--multimapping` = include-all), so
  this is complementary, not redundant (cf. review PMC7330433). Validated on EBV (= 12014 UMI, ≈
  other methods, since EBV ambiguity is within-virus / total-preserving). 330 tests. Committed
  `9f1963b`.
- `[x]` **S5 — `viralscan evidence` (read-level validation / IGV / BLAST).** New `evidence.py`
  (pure trace+extract) + `scripts/evidence_run.py` + `evidence` subcommand. Pure layer (geometry,
  EC→viral, FASTQ extraction) fully unit-tested. Parse helpers extracted: `_parse_coverage_output`
  and `_parse_blast_output` unit-tested against synthetic tool output without binaries. Live chain
  validated: `tests/integration/test_evidence_chain.py` synthesises a deterministic 2.2 kb viral
  genome + 40 exact-substring reads, then drives the real minimap2→samtools-sort/index→samtools-
  coverage→makeblastdb→blastn wrappers (`align_reads_to_viral`, `coverage_table`, `blast_identity`)
  and asserts BAM+.bai exist, virus_A has > 0 mapped reads with non-zero breadth coverage, and all
  BLAST hits are ≥ 95 % identity to virus_A. Skips gracefully when binaries absent (``have_tools``
  guard). No public data used. End-to-end `viralscan evidence` on a real run-dir remains an
  operational step (needs a kb-python output tree) — not a code gap.
- `[x]` **S6 — apply the evidence-module geometry fix to `host_filter.py`** (closes S1 properly).
  Rewrote `host_filter.py`: deleted `_TECH_PARAMS`/`_cb_umi_lengths()`; `_starsolo_filter` and
  `_kallisto_filter` now call `cb_umi_geometry(technology)` directly; `filter_fastq_pairs` refactored
  as a public, pure function with explicit args (testable without Snakemake); all module-level
  Snakemake bindings guarded under `if "snakemake" in globals():`. Regression guard:
  `tests/test_host_filter.py` (11 tests: geometry resolution incl. DROPSEQ, keep/drop logic with
  synthetic FASTQs, gzipped input, explicit-geometry string, mid-record truncation raises ValueError).
  Post-review fix (code-review pass): `filter_fastq_pairs` now checks `not lines2[0]` (R2 EOF) and
  `not lines1[3] or not lines2[3]` (mid-record truncation) — silent emission of malformed records
  on truncated FASTQ was a latent correctness bug. Also applied the PR-15 `if "snakemake" in
  globals():` guard to `createconfig.py`, which was the last script in `scripts/` missing it.
- `[x]` **S7 — `multimap.py` truncated bus.txt crash (NA in ec/count).** Root cause: `pd.read_csv`
  with `dtype={"ec": "int32", "count": "int32"}` raises `ValueError: Integer column has NA values`
  when the last line of bus.txt is truncated (file written mid-line, missing count field). Affects
  deep samples (e.g., SARS-CoV-2 mock with 13 GB bus.txt). Fix: read without enforced dtype for
  ec/count (pandas reads partial rows as NaN in object/float columns), `dropna`, then cast to int32.
  All 402 tests pass. Committed `cb8de47`.
- `[x]` **S8 — `viralscan` interactive overwrite prompt breaks SLURM resume jobs.** Root cause:
  `check_output()` always calls `input()` when the output dir already exists; SLURM jobs have no
  stdin, so this raises `EOFError` and the job exits immediately as FAIL:viralscan (visible as
  ~25-second completion time in sacct). Surfaced when resume scripts re-ran viralscan on existing
  output dirs after a timeout. Fix: add `--yes` / `-y` flag to the main parser that short-circuits
  the prompt; `check_output()` returns early when `args.yes` is True. All 402 tests pass.

---

## PR 15 — Run-context refactor (architecture deepening)

Root cause: there is no "ViralScan run" module. A run is a loose dict threaded
through Snakemake magic globals, so each worker script re-derives config types,
the kb-python output layout, and file selection locally — and five scripts are
tested through mirror re-implementations that can silently drift. See
`CONTEXT.md` (Run, Run Config, Run Context, Kb Count Outputs).

Locked design decisions:
- Run Context is a **passive** value (`run(ctx)`), not a capability object.
- Config is validated **once on write** (`RunConfig.from_snakemake_config`);
  downstream reads are trusted typed loads (`from_yaml`).
- Current-adata resolution is **by config flag**, not file existence — a missing
  multimap adata when `multimapping` is set is an error, not a silent downgrade.

Success metric: mirror functions in the test suite go 5 → 0.

- `[x]` **Step 1 — `RunConfig`.** New `src/viralscan/runconfig.py`: typed frozen
  dataclass; `from_snakemake_config` (the single coercion+validation checkpoint),
  `from_yaml` (trusted load), `to_yaml`/`to_dict`. `createconfig.py` rewritten to
  delegate. `test_createconfig.py` `_build_cfg` mirror replaced by delegation to
  the real API (~55 assertions now cross the real seam), and the missing
  `"False"`-string falsy cases added — closing the `bool("False") is True` hazard
  at the coercion point. Verified: exact 34-key parity with old output, YAML
  round-trip, 55 createconfig tests + 112 adjacent tests pass (snakemake env).
- `[x]` **Step 2 — `KbCountOutputs`.** New `src/viralscan/kb_outputs.py`: frozen
  dataclass owning the kb-python layout (named paths + `current_adata(multimapping=...)`,
  resolved by flag, no I/O). Repointed `multimap.py` (`define_paths`, multimap
  write), `detection.py`, `umap.py` — no kb-python path literals remain in
  multimap/umap (detection keeps only a user-facing "look in the folder"
  sentence). New `tests/test_kb_outputs.py` (7 tests incl. legacy-fstring parity
  + trailing-sep normalisation). **Correction:** the earlier claim that
  `detection` resolved by file-existence while `umap` used the flag was wrong —
  both already used `if config["multimapping"]`, so this is a pure dedup with
  zero behavior change. Snakefile `mv` block left as-is (it *creates* the layout
  in bash; routing it through Python is out of scope). Full suite: 307 passed.
- `[x]` **Step 3 — `RunContext` + `run(ctx)`.** New `src/viralscan/run_context.py`
  (passive: config dict + `KbCountOutputs`, `from_yaml`/`from_config`). All four
  workers (`analysis`, `multimap`, `detection`, `umap`) are now importable without
  Snakemake — the magic-global wiring runs only under `if "snakemake" in globals():`
  and calls a `run(ctx, …)` entry that sets the run-level state. Extracted pure
  helpers so tests hit production code: `analysis.extract_gene_ids`,
  `multimap.strip_10x_suffix`, `detection.detect_genes`; `umap` lazily imports
  plotly so the module imports in a test env. **Mirror tests eliminated 5 → 0**:
  `test_analysis`/`test_multimap`/`test_umap`/`test_detection` now import the real
  functions (only the deliberate `*_buggy` negative controls remain in
  `test_multimap`). Stale "cannot import directly" docstrings corrected. Full
  suite: 307 passed, ruff clean. **Note:** read-side `config` stays a plain dict
  (helpers use `config.get`); run() sets module globals from ctx — pragmatic
  vs. threading ctx through ~15 `config.get` sites in detection/umap. Migrating
  the read side to typed `RunConfig` is possible future polish. **Polish done:**
  `RunContext.config` is now `RunConfig` (not `dict`). All four workers, plus the
  shared helpers `select_detection_matrix`, `should_write_multimap_evidence`,
  `summarize_multimap_evidence`, and `cell_type_enrichment`, migrated to attribute
  access. Two previously-missing EM tuning fields (`multimap_em_max_iter`,
  `multimap_em_tol`) added to `RunConfig` and `createconfig`. Tests updated.
- `[x]` **Step 4 — `group_genes_by_virus`.** New `src/viralscan/virus_grouping.py`
  with `virus_name_for_gene` + `group_genes_by_virus`, repointed in `detection.py`
  (`_group_viral_genes` deleted, `histogram` inline loop, evidence call) and
  `umap.py` (`gene_to_virus`). **Resolved a real semantic divergence:** the two
  inlined rules disagreed on 151/2692 bundled gene IDs and were both buggy —
  substring over-matched (`EPSTEIN_HHV4_BORF1`→Orf), prefix under-matched
  (`TTV7_gp2`→nothing). Unified on a **boundary-aware** rule (key match when gene
  == key, or starts with key and next char is `_`/digit; longest key wins) — user
  signed off on the output change. New `tests/test_virus_grouping.py` (13 tests)
  guards the divergence cases. Full suite: 320 passed, ruff clean.

---

## Completed work (verified 2026-05-08)

All items below were confirmed present in the codebase.

- PR 1 Hygiene: README fixes, .gitignore, pyproject.toml classifiers/markers/ruff, CHANGELOG, CITATION.cff
- PR 2 Correctness & security: bool flags, subprocess.run, sys.exit/_die, pathlib, whitelist None-check, analysis.py config.get
- PR 3 Partial cleanup: dead code removed, constants.py VIRUS_NAME_MAP, utils.py load_config, umap.py function ordering, unused imports removed
- PR 3 Logging: --verbose/--quiet CLI flags added to menu.py; ANSI print→logging done in Python scripts
- PR 4 Tooling: pre-commit, CI/CD workflows, ruff formatting baseline
- PR 5 Tests: test_ncbi_fetch.py (18), test_cli.py (8 classes), test_createconfig.py (5), test_analysis.py (5), test_errorhandler.py (4), conftest.py; Codecov badge
- PR 6 Reproducibility: environment.yml, Dockerfile, Singularity.def
- PR 7 Docs: full Sphinx docs/ skeleton, .readthedocs.yaml
- PR 9 Partial: type hints on menu.py, utils.py, constants.py
- PR 10 NCBI fetch: ncbi_fetch.py, --ncbi-accession CLI flag, cache, backoff, 18 unit tests
- PR 11 A1–A4: viral_summary.tsv, per_cell_viral.tsv, HTML report (Jinja2), normalized metrics, configurable thresholds
- PR 11 A5: --cell-types flag present in menu.py
- PR 12 Build-ref: build_reference.py, build-ref subcommand, ENSEMBL_SPECIES table, 22 tests
- PR 13: bool normalisation in createconfig.py, umap.py "True" checks removed, Snakefile None→empty check, analysis.py config.get
- PR 14 C2: detection threshold `>` → `>=` in detection.py preprocessing()
- PR 14 C3: gene_id_to_idx dict in multimap.py read_ec() — O(n) list.index replaced
- PR 14 C4: iterrows → itertuples in multimap.py build_multimap_matrix()
- PR 14 C5: stderr merged (2>&1) in Snakefile kb_count shell block
- PR 14 C6: mkdir -p in Snakefile kb_count shell block
- PR 14 C7: redundant f.close() removed from analysis.py
- PR 14 C8: file handle leak fixed in detection.py with-statement
- PR 14 C9: var_names computed once in umap.py (single source)
- PR 14 C10: sc.pp.highly_variable_genes + force-include viral genes before PCA in umap.py
- PR 7 Docs — API reference, README overhaul, vignettes `[x]`
  Rewrote api.md (hand-written), overhauled README (7 sections), wrote 2 vignettes (basic_usage, cell_type_enrichment).

---

## Open tasks — ordered, self-contained, ready to implement

Each task below can be completed independently. All context needed is included.

---

### Task 1 — Fix double-count bug in `umap.py`  `[x]`

**Why:** umap.py lines 330–331 still add `counts_corrected + counts_original`, double-counting
uniquely-mapping reads (same bug as PR 14 C1, which was fixed in detection.py but not umap.py).

**File:** `src/viralscan/scripts/umap.py`

**Find this block (around line 330):**
```python
if "counts_corrected" in adata.layers and "counts_original" in adata.layers:
    adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
```

**Replace with (Option B — only add the extra multimapper share):**
```python
if "counts_corrected" in adata.layers and "counts_original" in adata.layers:
    adata.X = adata.layers["counts_original"] + adata.layers["counts_corrected"]
    # counts_corrected holds only the redistributed multimapper fraction (share per gene
    # when EC maps to >1 gene; 0.0 for unique-mapping ECs), so the sum is correct.
```

Actually confirm the multimap.py share logic first (grep for `share = `). If `share = 0.0`
for unique ECs is already in place (it is — confirmed in audit), the addition is correct as
written in detection.py and just needs the same treatment in umap.py. Verify by reading
`multimap.py` lines around 223–229 before editing.

**Test after:** `PYTHONPATH=src python -m pytest tests/ -q`

**Completed 2026-05-08.** Comment added to umap.py clarifying counts_corrected only carries
multimapper shares (unique ECs are skipped in multimap.py). Regression tests added in
`tests/test_umap.py::TestLayerMergeNoDoubleCount`.

---

### Task 2 — Complete PR 11 A5: cell-type-aware enrichment  `[x]`

**Why:** `--cell-types` flag exists in menu.py (line 265) and is wired through createconfig.py,
but `detection.py` has no code that reads or uses `config["cell_types"]`. The feature is
silently ignored at runtime.

**What to implement in `src/viralscan/scripts/detection.py`:**
1. After `found_genes` is populated in `preprocessing()`, add a new function
   `cell_type_enrichment(adata, found_genes, config)`:
   - Read the barcode→cell_type CSV from `config.get("cell_types")`.
   - For each detected virus, compute per-cell-type viral prevalence (% barcodes ≥1 viral UMI).
   - Run `scipy.stats.fisher_exact` for each cell type vs. all others.
   - Return a DataFrame with columns: virus, cell_type, n_infected, n_total, pct, OR, pvalue, padj (BH).
2. Write output to `{output}/results/cell_type_enrichment.tsv`.
3. Include the table in the existing HTML report (add a section to the Jinja2 template in
   `src/viralscan/templates/`).
4. Skip gracefully if `config.get("cell_types")` is falsy.

**Dependencies already present:** scipy is a transitive dep via scanpy; pandas, jinja2 in requirements.

**Test after:** add a test in `tests/test_cli.py` or a new `tests/test_detection.py` that
mocks the CSV and checks the TSV is written. Then `PYTHONPATH=src python -m pytest tests/ -q`.

**Completed 2026-05-08.** `viralscan/enrichment.py` extracted with `cell_type_enrichment()` and
`write_cell_type_enrichment()`; `detection.py` imports and calls both; HTML report template
includes the cell-type enrichment table; `tests/test_detection.py::TestCellTypeEnrichment`
covers column schema, BH adjustment, and zero-infected-cell edge cases.

---

### Task 3 — Lift magic numbers into config defaults  `[x]`

**Why:** PR 9 magic-number lift is marked `[ ]` deferred. Several hardcoded values in
`umap.py` and `detection.py` should be user-configurable.

**Magic numbers to lift:**

| File | Line (approx) | Variable | Current value | Config key to add |
|---|---|---|---|---|
| `umap.py` | HVG call | `min_mean` | 0.0125 | `hvg_min_mean` |
| `umap.py` | HVG call | `max_mean` | 3 | `hvg_max_mean` |
| `umap.py` | HVG call | `min_disp` | 0.5 | `hvg_min_disp` |
| `umap.py` | neighbors call | `n_neighbors` | 15 (scanpy default) | `umap_n_neighbors` |
| `detection.py` | super-expressor | `se_threshold` | whatever hardcoded value | already exposed as `--se-threshold` — verify it's wired through |

**Steps:**
1. Create `src/viralscan/defaults.py` with a `DEFAULTS` dict of all the above keys and values.
2. In `createconfig.py`, merge `DEFAULTS` under the config YAML so all keys always exist.
3. In `umap.py`, replace hardcoded literals with `config.get("hvg_min_mean", 0.0125)` etc.
4. Expose the keys as optional CLI flags in `menu.py` (use `argparse` defaults that come
   from `DEFAULTS` so help text shows the value).

**Test after:** `PYTHONPATH=src python -m pytest tests/ -q`; add one test in
`test_createconfig.py` asserting DEFAULTS keys are present in the written YAML.

**Completed 2026-05-08.** `src/viralscan/defaults.py` created with DEFAULTS dict; `createconfig.py`
merges DEFAULTS into config YAML; `umap.py` reads all keys via `config.get(...)`; CLI flags
exposed in `menu.py` with defaults from DEFAULTS; `tests/test_createconfig.py` covers all keys.

---

### Task 4 — Add host pre-subtraction option  `[x]`

**Re-verified 2026-06-21 (after S2 fix).** S2 (pipeline halt after host_filter) was a Snakefile
bug not a Task 4 bug. After the S2 fix the pipeline runs end-to-end with `--host-filter`.

**Implemented 2026-05-08.** Two-aligner design, no breaking changes.
- `--host-filter {starsolo,kallisto}` + `--host-index PATH` added to `menu.py`
- `_check_host_filter_tools()` preflight: checks `STAR` or `kallisto`+`bustools`
- `createconfig.py`: writes `host_index`, `host_filter_aligner`, `kb_r1`, `kb_r2` to config YAML;
  `kb_r1`/`kb_r2` point to filtered FASTQs when active, else to original `sample1`/`sample2`
- `scripts/host_filter.py`: new Snakemake script implementing both modes
- `Snakefile`: conditional `host_filter` rule + `_kb_count_inputs()` helper; `kb_count` shell block
  uses `{config[kb_r1]}` / `{config[kb_r2]}` — no conditional logic in the shell
- `docs/faq.md`: new "Reducing false positives" section with STARsolo and kallisto examples

---

### Task 5 — PR 5: integration test skeleton  `[x]`

**Why:** There is no `tests/integration/` directory. The CI matrix has an `integration` mark
registered in `pyproject.toml` but no tests use it.

**What to add:**
1. Create `tests/integration/__init__.py` (empty).
2. Create `tests/integration/test_smoke.py`:
   - One test class `TestSmoke` with a single test `test_cli_help` that runs
     `subprocess.run(["python", "-m", "viralscan.menu", "--help"], check=True, capture_output=True)`
     and asserts returncode == 0. Mark with `@pytest.mark.integration`.
   - One test `test_build_ref_no_kb` (marked `@pytest.mark.integration` and
     `@pytest.mark.network`) that calls `viralscan build-ref --no-kb-ref --host human
     --virus-accessions NC_045512.2 --output /tmp/viralscan_test_ref` and checks that
     `combined.fasta` and `combined.gtf` exist.
3. Update `pyproject.toml` `[tool.pytest.ini_options]` markers to document `integration`.

**Test after:** `PYTHONPATH=src python -m pytest tests/integration/ -m integration -v`

**Completed 2026-05-08.** `tests/integration/__init__.py` and `tests/integration/test_smoke.py` added; `@pytest.mark.integration` and `@pytest.mark.network` gated tests present.

---

### Task 6 — mypy strict mode per-module  `[x]`

**Why:** PR 9 deferred `mypy --strict`. CI currently runs mypy as informational (non-blocking).

**Steps (incremental — do not attempt the whole codebase at once):**
1. Run `PYTHONPATH=src mypy src/viralscan/utils.py src/viralscan/constants.py --strict 2>&1`
   and fix all errors. These two modules are already annotated.
2. Run `PYTHONPATH=src mypy src/viralscan/menu.py --strict 2>&1` and fix errors.
3. Add `[[tool.mypy.overrides]] module = "viralscan.utils" strict = true` etc. to
   `pyproject.toml` so the modules are always checked strictly in CI.
4. Repeat for `ncbi_fetch.py` and `build_reference.py` in separate commits.
5. The Snakefile scripts (`analysis.py`, `detection.py`, `umap.py`, `multimap.py`) use
   snakemake magic globals — exclude them from strict mode using
   `[[tool.mypy.overrides]] module = "viralscan.scripts.*" ignore_errors = true`.

**Test after:** `mypy src/viralscan/utils.py src/viralscan/constants.py --strict` exits 0.

**Completed 2026-05-08.** mypy installed to `.vendor_mypy/`; pyproject.toml `[tool.mypy]` section added with strict overrides for `utils`, `constants`, `menu`, `ncbi_fetch`, `build_reference`; Snakefile scripts excluded with `ignore_errors = true`. All type errors fixed. `PYTHONPATH=src:$PWD/.vendor_mypy python -m mypy src/viralscan` → Success (0 issues, 14 files).

---

### Task 7 — PR 7: rewrite getting_started.ipynb  `[x]`

**Why:** `getting_started.ipynb` has stale cells with errors and largely duplicates the README.
It is the first thing a new user opens.

**What to produce:**
A notebook that can run end-to-end on a small test dataset bundled in `tests/data/` (or
downloaded via a cell that fetches a 100k-read subset of a public SRA accession). Sections:
1. Installation (pip / conda one-liners).
2. Build reference with `viralscan build-ref --no-kb-ref --host human --virus-accessions NC_045512.2`.
3. Run `viralscan run` on the test FASTQ.
4. Inspect `viral_summary.tsv` and `report.html`.
5. UMAP plot.

**Blocked on:** Task 4 (host subtraction) being optional (so the notebook can run without
bowtie2/STAR). The notebook should use `--no-host-subtraction` or just omit `--host-index`.

**Completed 2026-05-08.** All stale cells deleted; 11 new cells with 5-section offline-safe tutorial: Installation, Build Reference, Run ViralScan, Inspect Outputs, UMAP. Uses `RUN_COMMANDS=False` guard flag so the notebook is safe to open without a live ViralScan environment.

---

### Task 8 — PR 8: data unbundling (Zenodo)  `[x]`

**Completed 2026-05-10.** Zenodo DOI: `10.5281/zenodo.20112332`.
- Added `viralscan data fetch` subcommand in `menu.py`.
- Added `viralscan.data_fetch` to resolve Zenodo metadata, download the archive, verify
  the Zenodo checksum plus optional SHA-256, and unpack GTF files to `~/.cache/viralscan/data/`.
- Changed `analysis.py` `obtain_gtf()` to read the cached data directory and raise a clear
  `viralscan data fetch` instruction if the cache is missing.
- Removed `data/*.gtf` from `[tool.setuptools.package-data]` in `pyproject.toml`.
- Updated README and docs with the data-fetch step.

---

### Task 9 — PR 9 remainder: detection/UMAP magic numbers (after Task 3)  `[x]`

Covered by Task 3. This entry is a reminder that Task 3 closes PR 9.

---

### Task 10 — Ambiguity-aware multimapper evidence  `[x]`

**Completed 2026-05-10.** Backward-compatible multimapper evidence feature.
- Added `--multimap-method {equal,host-conservative,unique-weighted}`,
  `--multimap-pseudocount`, and `--multimap-primary-call {legacy,unique-only,confidence}`.
- Preserved default behavior: `equal` method and `legacy` primary call.
- Added `viralscan.multimapping` pure functions for layer construction, detection matrix selection,
  and `multimap_evidence.tsv` summarization.
- Added AnnData diagnostic layers: `counts_multimap_equal`,
  `counts_multimap_host_conservative`, `counts_multimap_unique_weighted`,
  `counts_unique_viral`, `counts_host_viral_ambiguous`, and
  `counts_viral_ambiguous_upper`.
- Added `results/multimap_evidence.tsv` and an HTML report section without changing
  existing `viral_summary.tsv` or `per_cell_viral.tsv` schemas.
- Updated README and docs.

**Remediated 2026-05-10 after review.**
- Host-virus-only equal-split support now reports `low_confidence`, not `ambiguous`.
- `results/multimap_evidence.tsv` is written only when multimapping is enabled.
- Duplicate gene entries in EC mappings preserve legacy equal-split semantics.
- Added `counts_host_viral_selected` diagnostic layer for confidence tiering.

---

---

## Anellovirus reference expansion

ViralScan had only 20 RefSeq anellovirus accessions and a single name-map entry (`"TTV"`).
The Clareau lab's [`clareaulab/anellovirus_reference`](https://github.com/clareaulab/anellovirus_reference)
curates ~2,200 CD-HIT representative human anellovirus genomes plus a rich taxonomy table;
their sequences are NCBI GenBank public records and their CSV is a factual accession/taxonomy
table. We reconcile the two accession sets (~2,042 unique), re-derive sequences from NCBI,
expand the name map to all modern Anelloviridae genera, and bundle everything into the Zenodo
panel (GTF-only fetch → FASTA+GTF+TSV). Reference: `CLAUDE.md` §"What we steal from clareaulab".
clareaulab cited in `docs/reference_panel.md`.

- [x] **0.1** PLAN.md section added; "Next up" pointer updated.
- [x] **0.2** `simple_anello_metadata_V2.csv` inspected: 3545 rows, 2023 CD-HIT representatives,
  columns `Accession, Species, Genus, Family, Virus Name, infer_genus, cdhit_representative`.
  Only 1 accession overlaps with existing 20 ViralScan RefSeq entries. Union = ~2042 unique.
- [x] **A.1** `extras/build_anello_table.py` curation script: reads CSV, keeps representatives.
- [x] **A.2** Extracts existing 20 ViralScan accessions from bundled GTFs; tags both sources.
- [x] **A.3** Union + dedup by bare accession; prefers NC_* RefSeq when both present.
- [x] **A.4** Emits `src/viralscan/data/anellovirus_accessions.tsv`
  (`accession  virus_name  genus  family  source`).
- [x] **A.5** `src/viralscan/anellovirus.py`: `load_accession_table()` + `anello_name_map()` +
  `merged_name_map()`; `importlib.resources` wiring for packaged TSV.
- [x] **B.1** `build_anellovirus_reference()` in `build_reference.py`; reads packaged TSV;
  reuses `ncbi_fetch.fetch_reference()`.
- [x] **B.2** GTF via `_genome_as_transcript_gtf` (whole-genome, `gene_id "{acc}_geneN"`);
  extracted via `_gtf_from_merged_fasta()` helper.
- [x] **B.3** Optional `dustmasker -window 64 -level 30` step (binary is `dustmasker`, not
  `dustmask`); guarded by `shutil.which`; warn + skip if absent.
- [x] **B.4** Optional `cd-hit-est` step; off by default; guarded by `shutil.which`.
- [x] **B.5** `viralscan build-ref --anellovirus` flag in `menu.py`; also `--no-mask` and
  `--cluster`; dispatches to `build_anellovirus_reference()`, skips `--host`/`--virus-accessions`.
- [x] **B.6** Unit tests in `test_build_reference.py` (4 tests, stub `fetch_reference`): FASTA+GTF
  gene_id assertions; mask no-op when `_run_dustmasker` returns False; cluster no-op when
  `_run_cdhit_est` returns False; default accession list loads packaged TSV.
- [x] **C.1** `data_fetch.py:_extract_gtfs` → `_extract_members`: also extracts `.fa/.fasta`
  and `anellovirus_accessions.tsv`.
- [x] **C.2** Manifest + `cache_valid` extended for FASTA/aux file checksums; back-compat kept.
- [x] **C.3** `bundled_anellovirus_fasta()` accessor in `data_fetch.py`.
- [x] **C.4** `menu.py --reference-panel anellovirus` + `build_anellovirus_reference(fasta_path=)`:
  tries bundled FASTA from Zenodo cache first; gracefully falls back to NCBI download if absent.
- [x] **C.5** Tests for C.1–C.4: `TestExtractMembers` (5), `TestCacheValidExtended` (5),
  `TestBundledAnellovirusFasta` (5) — synthetic zip/tar archives with GTFs + FASTA + TSV.
- [ ] **C.6** *(manual)* Rebuild Zenodo archive to include the anellovirus FASTA + accession TSV;
  publish as a new Zenodo version; bump DOI and SHA-256 checksums in `data_fetch.py` manifest
  defaults + `docs/reference_panel.md`. Low priority — `--reference-panel anellovirus` already
  falls back to NCBI download when FASTA is absent from cache.
  Sub-steps:
  - [ ] **C.6a** Run `scripts/build_bundled_panel_ref.py --out $WORKDIR/ref` locally (or via
    `sbatch --wrap`) to produce the anellovirus FASTA + TSV, then bundle into the archive.
  - [ ] **C.6b** Upload new Zenodo version; record new DOI (`10.5281/zenodo.XXXXXXX`).
  - [ ] **C.6c** Patch `data_fetch.py` `_ZENODO_*` constants (URL, checksum, size) and run
    `PYTHONPATH=src python -m pytest tests/test_data_fetch.py -q` to confirm.
- [x] **D.1** `anello_name_map()` in `anellovirus.py`: accession → genus label (rollup).
- [x] **D.2** Anellovirus genus display names added to `VIRUS_NAME_MAP` in `constants.py`.
- [x] **D.3** `merged_name_map()` threaded through `detection.py` and `umap.py`.
- [x] **D.4** Tests: accession bare/versioned/`_geneN` resolve correctly; unmapped falls back.
- [x] **E.1** `docs/reference_panel.md` updated: coverage (20 → ~2 k), `--anellovirus` usage,
  dustmasker/cd-hit-est prerequisites, clareaulab citation + provenance.
- [x] **E.2** Integration test `tests/integration/test_anellovirus_chain.py` (4 tests,
  `@pytest.mark.integration`): synthetic FASTA → `_gtf_from_merged_fasta` → gene_id
  extraction → `virus_name_for_gene` + `group_genes_by_virus` → correct genus labels.
  Full kb-python pipeline is an operational step (needs installed binaries + run-dir).
- [x] **E.3** Full suite green: 366 passed, 15 deselected (2026-06-22, pegasuspy/Python 3.11).
- [x] **E.4** All A/B/D/E rows flipped; "Next up" updated. C rows remain `[ ]` — deferred
  pending a new Zenodo release.

#### Post-review polish (code-review pass, 2026-06-22)
- [x] **R.1** Fixed misleading comment `build_reference.py:343` — now correctly
  documents that the versioned accession is kept (not stripped).
- [x] **R.2** Corrected `--anellovirus` help in `menu.py` — was "Skips
  --host / --virus-accessions" (wrong); now accurately states that `--host` is
  ignored and `--virus-accessions` is treated as an explicit subset.
- [x] **R.3** `_run_dustmasker` and `_run_cdhit_est` now wrap `subprocess.run`
  in `try/except CalledProcessError` → `log.error` + `return False`, matching the
  graceful-degradation pattern already used by the `kb ref` step.
- [x] **R.4** Added `test_empty_fasta_produces_empty_outputs` to
  `TestBuildAnellovirusReference`; suite now 367 passed, 15 deselected.
- [x] **R.5** Fixed `ncbi_fetch.py:_fetch_one` crash on accessions with no CDS features
  (e.g. HM224451.1 in the anellovirus set). `_genbank_to_gtf` now falls back to
  `_whole_genome_gtf_from_fasta` (new helper) when no CDS annotations exist, rather than
  raising `NCBIFetchError`. Surfaced by first live run of `viralscan build-ref --anellovirus`.
  Suite: 367 passed, 15 deselected (2026-06-22).

---

## PR 20 — Anellovirus C.1–C.5: Zenodo FASTA bundling support (2026-06-22)

Extends the Zenodo fetch/cache layer to carry an anellovirus FASTA (and the
accession TSV) alongside the GTF panel, and wires `--reference-panel anellovirus`
to use it when available.

- `[x]` **C.1** `_extract_gtfs` replaced by `_extract_members` in `data_fetch.py`.
  Now categorises members into `{"gtf", "fasta", "tsv"}` and extracts all three
  types from zip and tar.gz archives (only the sentinel filename
  `anellovirus_accessions.tsv` is matched for TSV; any `.fa`/`.fasta` file matches).

- `[x]` **C.2** `cache_valid` extended: checks `fasta` + `fasta_checksum` and
  `tsv` + `tsv_checksum` fields when present in the manifest. Old manifests without
  those fields pass unchanged (backward-compat). `fetch_viral_data` writes the new
  fields when the archive includes a FASTA/TSV.

- `[x]` **C.3** `bundled_anellovirus_fasta(cache_dir)` accessor: reads manifest,
  returns `Path` to the cached FASTA. Raises `ViralScanDataError` with actionable
  messages when the manifest is absent, has no `fasta` key (pre-bundle archives),
  or the file is missing on disk.

- `[x]` **C.4** `build_anellovirus_reference(fasta_path=None)` extended: when a
  `Path` is supplied, logs it and skips the NCBI download step entirely.
  `build_ref_main` (in `build_reference.py`) + `menu.py` wired with
  `--reference-panel anellovirus`: tries `bundled_anellovirus_fasta()`, logs and
  falls back to NCBI download on `ViralScanDataError`.

- `[x]` **C.5** Tests: 15 new tests across 3 classes in `tests/test_data_fetch.py`.
  `TestExtractMembers` (5): zip+tar extraction, `.fasta` suffix, non-matching
  files ignored, wrong-named TSV ignored. `TestCacheValidExtended` (5): valid
  extended manifest, FASTA/TSV checksum mismatch each invalidate, missing FASTA
  file invalidates, old manifest without fasta fields still valid.
  `TestBundledAnellovirusFasta` (5): success path, no manifest, no `fasta` key,
  fasta file missing, `fetch_viral_data` end-to-end with FASTA+TSV in archive.

Verification: `PYTHONPATH=src python -m pytest tests/test_data_fetch.py -q` →
**30 passed**. Full suite: 343 passed, 4 deselected (snakemake env, 2026-06-22).

Note: C.6 (rebuild Zenodo archive with the FASTA included) is a manual step
deferred until the next release. Until then `--reference-panel anellovirus`
gracefully falls back to NCBI download (same behavior as `--anellovirus`).

---

## PR 19 — Tier 4 tidy-ups (2026-06-22)

Batchable maintainability fixes from the clean-code review plan (Tier 4).

- `[x]` **P19.1 — EM epsilon guard** (`multimapping.py`). `if s <= 0.0:` changed to
  `if s <= 1e-12:` so subnormal-weight ECs fall back to equal-split instead of
  only exactly-zero sums.
- `[x]` **P19.2 — Inline imports hoisted to module top**. `import yaml` moved from
  inside `load_config()` to `utils.py` module top. `import base64` / `import datetime`
  moved from inline positions in `_encode_image()` and `generate_html_report()` to
  `detection.py` module top (duplicate inline `import base64` removed).
- `[x]` **P19.3 — Drop dead `build_multimap_matrix`** (`scripts/multimap.py`).
  Production `run()` uses `build_multimap_layers` from `viralscan.multimapping`; the
  old function was never called. Removed. `TestBuildMultimapMatrix` in
  `test_multimap.py` migrated to use `build_multimap_layers(method="equal")` via a
  new `_corrected()` helper — same behavioral assertions, now cross the production seam.
- `[x]` **P19.4 — Hoist `np.where` out of per-cell loop** (`detection.py`).
  `np.where(infected_mask)` was recomputed O(N_infected) times per virus inside the
  per-cell loop. Extracted to `infected_indices` before the loop.
- `[x]` **P19.5 — Narrow bare `except Exception`**. `enrichment.py` catch narrowed to
  `(OSError, pd.errors.ParserError)`; `detection.py` HTML template catch narrowed to
  `TemplateNotFound` (imported alongside `Environment`/`FileSystemLoader`).

Verification: 387 passed, 15 deselected (2026-06-22).

---

## PR 17 — Multimap checkpoint: default → `equal`, `viralscan rerun-multimap` (2026-06-22)

Multimapping is now a resumable checkpoint. Fast runs complete with the default equal-split
method; users can later invoke `viralscan rerun-multimap` to switch algorithms without
re-running the expensive `kb count` pseudoalignment step.

Design insight: `build_multimap_layers()` always pre-stores all three non-EM layers
(`counts_multimap_equal`, `counts_multimap_host_conservative`, `counts_multimap_unique_weighted`)
in every multimap h5ad. Switching between them is a free in-place layer swap; only EM requires
re-processing bus files (which are preserved from the original run).

- `[x]` **P17.1 — Default changed to `equal`** (`src/viralscan/defaults.py`).
  `DEFAULT_MULTIMAP_METHOD` changed `"host-conservative"` → `"equal"`.
  Test `TestBuildMultimapLayers::test_default_method_is_host_conservative` renamed and
  updated to assert `DEFAULTS["multimap_method"] == "equal"`.

- `[x]` **P17.2 — `_swap_multimap_layer(adata_path, new_method)` helper** (`menu.py`).
  Pure testable function. Loads h5ad, overwrites `counts_corrected` from the pre-stored
  layer, updates `uns["multimap_method"]`, writes back. Returns `False` when the target
  layer is absent (older run), signalling caller to fall back to full multimap rerun.

- `[x]` **P17.3 — `viralscan rerun-multimap` subcommand** (`menu.py`).
  Finds all sample subdirs with `log/multimap.done`. For non-EM methods: fast swap via
  `_swap_multimap_layer` (if layers present) — deletes only `detection.done` + `umap.done`,
  snakemake re-runs only those two rules. For EM or absent layers: also deletes `multimap.done`,
  snakemake re-runs from bus file. Updates `config.yaml` `multimap_method` in-place before
  re-invoking snakemake (safe because `create_config.done` still exists).

- `[x]` **P17.4 — Tests** (`tests/test_rerun_multimap.py`, 8 tests).
  `TestSwapMultimapLayer`: swap to each of the three non-EM methods, False on missing layer,
  other layers preserved after swap.
  `TestRerunMultimapParser`: `--help` exits 0, method parsed correctly, unknown method rejected.

Verification: `PYTHONPATH=src python -m pytest tests/ -q` → **375 passed, 15 deselected**.
Files modified: `src/viralscan/defaults.py`, `src/viralscan/menu.py`,
`tests/test_multimapping.py`, `tests/test_rerun_multimap.py` (new).

---

## Verification checklist (run after every task)

```bash
PYTHONPATH=src python -m pytest tests/ -q          # must stay at 223+ passed, 0 failed
PYTHONPATH=src python -m viralscan.menu --help      # smoke — must not crash
ruff check src/ tests/                              # must be clean
PYTHONPATH=src:$PWD/.vendor_mypy python -m mypy src/viralscan  # strict-module check
```

---

## Audit remediation — 2026-05-08 TDD session

All findings from `audits/2026-05-08-full-pipeline.md` were addressed via strict
Red → Green → PLAN.md workflow.  Tests were written *before* the fix was applied.

### Task A10 — Fix §2.2 GTF gene_id parsing  `[x]`

**Module:** `src/viralscan/scripts/analysis.py:obtain_gtf()`
**Severity:** HIGH
**Bug:** `info.split('"')[1]` grabs the first quoted token regardless of attribute
name.  A GTF with attributes in any order (e.g. `source "NCBI"; gene_id "NC_123"`)
returns the wrong accession.
**Fix:** Replaced split with `re.search(r'gene_id "([^"]+)"', info)`.  Added
`if len(cols) < 9: continue` guard for malformed lines.
**Tests:** `tests/test_analysis.py::TestGtfGeneIdAttributeOrder` (4 tests)
**Commit:** `fix(analysis): use regex to extract gene_id attribute by name not position (audit §2.2)`

---

### Task A11 — Fix §3.3 detection_threshold=0 validation  `[x]`

**Module:** `src/viralscan/scripts/createconfig.py`
**Severity:** LOW→WRONG (silent correctness hazard)
**Bug:** `int(cfg_in.get("detection_threshold", 1))` accepts 0 silently, which
causes `total_count >= 0` to always be True — all 195 viruses are "detected".
**Fix:** Added validation block that raises `ValueError` if threshold < 1.
**Tests:** `tests/test_createconfig.py::TestDetectionThresholdValidation` (4 tests)
**Commit:** `fix(createconfig): raise ValueError for detection_threshold < 1 (audit §3.3)`

---

### Task A12 — Fix §3.1 missing random seeds  `[x]`

**Module:** `src/viralscan/scripts/umap.py:viral_neighbor_enrichment()`
**Severity:** MEDIUM (non-reproducible results)
**Bug:** `np.random.permutation(labels)` uses global random state; no seeds on
`sc.pp.pca`, `sc.pp.neighbors`, `sc.tl.umap`.
**Fix:** Added `random_state` parameter to `viral_neighbor_enrichment` with
`np.random.default_rng(random_state)`.  Added `random_state=0` to all scanpy calls.
**Tests:** `tests/test_umap.py::TestViralNeighborEnrichmentReproducibility` (5 tests)
**Commit:** `fix(umap): make permutation test and UMAP reproducible with seeded RNG (audit §3.1)`

---

### Task A13 — Fix §2.3 barcode suffix stripping  `[x]`

**Module:** `src/viralscan/scripts/multimap.py:load_barcodes()` + `normalize_barcodes()`
**Severity:** HIGH (data corruption)
**Bug:** `bc.replace("-1", "")` is a global substitution that corrupts any barcode
with "-1" at a non-trailing position (e.g. "ACGT-1GCTA-1" → "ACGTGCTA" instead
of "ACGT-1GCTA").
**Fix:** Replaced with `bc.removesuffix("-1")` (Python ≥ 3.9) and
`.map(lambda bc: bc.removesuffix("-1"))` for the DataFrame column.
**Tests:** `tests/test_multimap.py::TestLoadBarcodes` (7 tests),
`tests/test_multimap.py::TestNormalizeBarcodes` (2 tests),
`tests/test_multimap.py::TestBuildMultimapMatrix` (3 tests)
**Commit:** `fix(multimap): strip only trailing '-1' suffix using str.removesuffix (audit §2.3)`

---

### Task A14 — Fix §3.2 cache content validation  `[x]`

**Module:** `src/viralscan/scripts/ncbi_fetch.py:_fetch_one()`
**Severity:** MEDIUM (silent use of corrupt/truncated cached files)
**Bug:** `_fetch_one()` only checked `path.exists()` and `st_size == 0`.
The `_checksum()` function was defined but never called, so an interrupted download
that left a non-empty truncated file would be silently reused on the next run.
**Fix:** Added `_cache_valid(path)` helper that checks existence, size, and SHA-256
sidecar file (`.sha256`).  Added `_write_cached(path, content)` helper that writes
both the file and the sidecar.  Replaced all direct `.write_text()` calls in
`_fetch_one()` with `_write_cached()`.
**Tests:** `tests/test_ncbi_fetch.py::TestCacheValidation` (3 tests):
- `test_no_sidecar_triggers_redownload`
- `test_mismatched_sidecar_triggers_redownload`
- `test_valid_sidecar_skips_redownload`
**Commit:** `fix(ncbi_fetch): validate cache with SHA-256 sidecar, re-download on mismatch (audit §3.2)`

---

### Task A15 — §2.1 Detection threshold regression guard  `[x]`

**Module:** `src/viralscan/scripts/detection.py:preprocessing()`
**Severity:** HIGH (would silently break if `>=` was changed to `>`)
**Finding:** The `total_count >= threshold` comparison is correct but unguarded
by any test.  A future refactor changing `>=` to `>` would cause missed detections
at exactly the threshold (the most common edge case).
**Fix:** Tests written; no source change required (code was already correct —
this task adds the regression safety net).
**Tests:** `tests/test_detection.py::TestDetectionThreshold` (7 tests):
- `test_gene_with_zero_counts_never_detected`
- `test_gene_at_threshold_is_detected` (guards the inclusive `>=`)
- `test_gene_below_threshold_excluded`
- `test_host_genes_never_in_found_even_when_high_count`
- `test_unknown_viral_accession_silently_skipped`
- `test_sparse_input_handled_identically`
- `test_total_count_value_is_sum_across_all_cells`
**Commit:** `test(detection): add regression guard for threshold >= filtering (audit §2.1)`

---

### Phase 2 — Integration test skeleton  `[x]`

**Location:** `tests/integration/test_smoke.py`
**Mark:** `@pytest.mark.integration` (excluded from default run)
**Fixture:** Synthetic AnnData built with `anndata` directly (no FASTQs, no network).
**Tests:** `TestEndToEndCountConservation` (4 tests):
- `test_x_equals_corrected_plus_original`
- `test_umi_mass_not_inflated` (grand total conservation)
- `test_gene_totals_match_known_values`
- `test_no_negative_counts_in_x`
**Commit:** `test(integration): add UMI count conservation skeleton (audit §3.4)`

---

### Remediation summary table

| ID  | Severity | Module | Finding | Status | Test class |
|-----|----------|--------|---------|--------|------------|
| §2.1 | HIGH | `detection.py:preprocessing()` | `>=` threshold guard missing | `[x]` DONE | `TestDetectionThreshold` |
| §2.2 | HIGH | `analysis.py:obtain_gtf()` | `gene_id` extracted by position not name | `[x]` DONE | `TestGtfGeneIdAttributeOrder` |
| §2.3 | HIGH | `multimap.py:load_barcodes()` | global `-1` replace corrupts internal substrings | `[x]` DONE | `TestLoadBarcodes`, `TestNormalizeBarcodes` |
| §3.1 | MEDIUM | `umap.py:viral_neighbor_enrichment()` | global RNG → non-reproducible p-values | `[x]` DONE | `TestViralNeighborEnrichmentReproducibility` |
| §3.2 | MEDIUM | `ncbi_fetch.py:_fetch_one()` | truncated cache not detected | `[x]` DONE | `TestCacheValidation` |
| §3.3 | LOW→WRONG | `createconfig.py` | `detection_threshold=0` silently accepted | `[x]` DONE | `TestDetectionThresholdValidation` |
| §3.4 | MEDIUM | `multimap.py` / integration | UMI conservation untested | `[x]` DONE | `TestEndToEndCountConservation` |

All 208 tests pass (`208 passed, 6 deselected`) as of this session.

---

## PR 16 — Clean-code review: Tier 1+2 bug fixes and hardening (2026-06-22)

Three `clean-code-reviewer` agents covered the full codebase (~5,900 lines). This
PR implements the confirmed bugs (Tier 1) and robustness / fail-fast gaps (Tier 2).
Tier 3 (config-key list deduplication, god-function decomposition, `host_filter`
migration) deferred to a follow-up; Tier 4 tidy-ups can ride opportunistically.

Note: two reviewer HIGH findings were **false positives** after verification.
`menu.py:842` appends `os.sep` to every sample output path, so
`f"{config.output}log/found_genes.txt"` in `detection.py` / `f"{config.output}log/analysis.txt"`
in `analysis.py` resolve correctly today. Downgraded to LOW (fragile trailing-sep contract).

- `[x]` **T1.1 — Per-sample runtime is cumulative, not per-sample** (`menu.py`).
  `start = time.time()` was set once before the sample loop; `summary.txt` recorded
  wall-clock since program start for every sample. Fix: removed global `start`;
  added `sample_start = time.time()` at the top of the loop body; changed
  `end - start` → `end - sample_start`.

- `[x]` **T1.2 — EM tuning knobs unreachable from the CLI** (`menu.py`).
  `multimap_em_max_iter` and `multimap_em_tol` existed in `DEFAULTS` and `RunConfig`
  but were absent from both the argparse definition and the `config_args` list passed
  to Snakemake, so `--multimap-method em` always silently used defaults. Fix: added
  `--multimap-em-max-iter` / `--multimap-em-tol` argparse args; added both keys to
  `config_args`.

- `[x]` **T1.3 — Loop-invariant transcript-count I/O** (`menu.py`).
  `_count_lines(transcripts)` and `_count_unique_genes(transcripts)` were called
  every sample iteration despite depending only on the reference. Fix: hoisted both
  above the loop.

- `[x]` **T2.1 — `kb ref` failure swallowed, reports exit 0** (`build_reference.py`).
  `CalledProcessError` was caught, logged, and discarded in both
  `build_combined_reference` and `build_anellovirus_reference`; `build_ref_main`
  then printed "Reference build complete." and exited 0. Fix: both handlers now
  `raise` after logging; `build_ref_main` catches `CalledProcessError` and calls
  `sys.exit(1)`. The intentional "no `kb` on PATH → silent skip" path is unchanged.

- `[x]` **T2.2 — plotly-absent crashes after expensive UMAP compute** (`umap.py`).
  `px` was set to `None` when plotly was absent, but `px.scatter(...)` was called
  unconditionally after the full UMAP + clustering computation. Fix: added an early
  guard at the top of `main()` — `if config.umap and px is None: raise RuntimeError(...)`.

- `[x]` **T2.3 — Raw Python bool interpolated into Snakemake config** (`menu.py`).
  `visual={args.visual}`, `reference={args.reference}`, `umap={args.umap}`,
  `multimapping={args.multimapping}` emitted Python `"True"`/`"False"` strings;
  correct only because `RunConfig._coerce_bool` defends the read side. Fix: added
  `_config_bool(v: bool) -> str` helper emitting `"true"`/`"false"`; replaced all
  four f-string interpolations.

- `[x]` **T2.4 — Sample-ID collision via `Path.name.split("_")[0]`** (`menu.py`).
  Two `--sample1` paths sharing a filename prefix would silently write to the same
  output directory. Fix: added `_sample_id(path)` helper (same derivation, now
  documented); added a pre-loop duplicate-ID guard that calls `_die()` on collision.

Verification: `PYTHONPATH=src python -m pytest tests/ -q` → **367 passed, 15 deselected**.
Files modified: `src/viralscan/menu.py`, `src/viralscan/scripts/build_reference.py`,
`src/viralscan/scripts/umap.py`.

---

## PR 21 — Host-response module: virus-driven gene expression (2026-06-23)

Optional post-pipeline module associating virus presence with host gene expression
via multi-seed L2 logistic regression + randomized Lasso stability selection
(Luebbert et al. 2026 / Meinshausen & Bühlmann 2010). Activated by `--host-h5ad`.

- `[x]` **P21.1 — `pyproject.toml`** — Added `scikit-learn>=1.0` to core deps;
  `enrichment = ["gget>=0.27"]` optional extra; ruff `per-file-ignores` and mypy
  `ignore_errors` entries for `hostresponse.py`.

- `[x]` **P21.2 — `defaults.py`** — Added 4 hostresponse defaults:
  `hostresponse_n_seeds=6`, `hostresponse_n_stab_iter=100`,
  `hostresponse_stab_min_prob=0.6`, `hostresponse_top_n_genes=50`.

- `[x]` **P21.3 — `runconfig.py`** — Added 8 new `RunConfig` fields (`host_h5ad`,
  `hostresponse_n_seeds`, `hostresponse_n_stab_iter`, `hostresponse_use_hvg`,
  `hostresponse_stab_min_prob`, `hostresponse_top_n_genes`, `hostresponse_enrichment`,
  `hostresponse_enrichment_db`) with `or`-based None fallbacks in
  `from_snakemake_config` to handle unspecified optional CLI args.

- `[x]` **P21.4 — `menu.py`** — Added 8 CLI args (`--host-h5ad`, `--hostresponse-n-seeds`,
  `--hostresponse-n-stab-iter`, `--hostresponse-use-hvg`, `--hostresponse-stab-min-prob`,
  `--hostresponse-top-n-genes`, `--enrichment`, `--enrichment-db`) and corresponding
  entries in `_build_config_args`.

- `[x]` **P21.5 — `Snakefile`** — Replaced hardcoded `rule all` inputs with
  `_all_targets(wildcards)` function that conditionally appends
  `log/hostresponse.done` when `config["host_h5ad"]` is set. Added conditional
  `rule hostresponse` after `rule umap`.

- `[x]` **P21.6 — `scripts/hostresponse.py`** (new, ~280 lines) — Key functions:
  `_detect_and_normalize` (raw-count detection heuristic: all-integer + max > 10;
  stores `_raw_depth` before normalization), `_select_features` (HVG or all genes),
  `_balanced_split` (top-50%-depth-filtered 80/20 balanced split),
  `_run_l2_regression` (multi-seed L2 logistic regression returning weights_df +
  metrics), `_run_stability_selection` (randomized Lasso stability selection),
  `_run_enrichment` (gget.enrichr, optional), `run_hostresponse` (per-virus loop).
  Outputs: `<virus>_gene_weights.csv`, `<virus>_stability.csv`,
  `hostresponse_metrics.csv`, optionally `<virus>_enrichment_<db>.csv`.

- `[x]` **P21.7 — `tests/test_hostresponse.py`** (new, 29 tests) — Unit tests for all
  public functions + 6 integration tests for `run_hostresponse`. All 431 tests pass.

- `[x]` **P21.8 — sklearn `penalty='l1'` FutureWarning** (`scripts/hostresponse.py`) —
  Added module-level `_SKLEARN_VER` / `_L1_LR_KWARGS`: sklearn ≥ 1.8 uses
  `solver='saga', l1_ratio=1.0`; older sklearn uses `penalty='l1', solver='liblinear'`.
  Eliminates FutureWarning that would become an error in sklearn 1.10.

- `[x]` **P21.9 — `viralscan hostresponse` standalone subcommand** (`menu.py`) —
  Added `_build_hostresponse_parser` + `_run_hostresponse_subcommand`; registered in
  `create_help()` subparsers and dispatched in `main()`. Takes an existing viralscan
  sample output dir + `--host-h5ad`, loads `RunConfig.from_yaml`, resolves the virus
  h5ad via `KbCountOutputs`, and calls `run_hostresponse()` directly (no Snakemake
  re-run). CLI override flags for all numeric/boolean hostresponse params; falls back
  to config values when not specified.

- `[x]` **P21.10 — `tests/test_hostresponse_subcommand.py`** (new, 20 tests) — Parser tests
  (help, required args, flag defaults, overrides) and dispatch tests (happy path, CLI
  overrides config, error paths for missing config.yaml / analysis.txt / virus h5ad).

Verification: `PYTHONPATH=src python -m pytest tests/ -q` → **451 passed, 15 deselected** (2026-06-23).

- `[x]` **P21.11 — User-facing documentation** — Added `hostresponse` section to
  `docs/cli_reference.md` (flags for both the subcommand and the in-pipeline
  `--host-h5ad` mode), backfilled missing `evidence` and `rerun-multimap`
  sections in the same file, added host-response Feature bullet + User Guide
  subsection to `README.md`, documented `[enrichment]` extra in
  `docs/installation.md`, and added `[Unreleased]` Added/Fixed entries to
  `CHANGELOG.md`.

---

## PR 22 — Publication-readiness gaps (2026-06-23)

Identified by pre-publication audit. Three code items done immediately; four operational items
planned here for tracking.

- `[x]` **P22.1 — EM global-pool caveat** (`src/viralscan/multimapping.py`) — Added
  `**Global-pool design:**` paragraph to `em_gene_abundances` docstring explaining that
  `ec_counts` aggregates multi-gene EC masses across *all* cells before EM; returned `theta`
  is transcriptome-wide (not per-cell); differs from per-cell EM (alevin-fry, STARsolo);
  faster but ignores cell-to-cell abundance variation. Also added `## Limitations` section
  to `README.md` documenting five limitations (FPR/FNR not characterized, global-pool EM,
  untested chemistries, cross-homology inflation, ambient RNA not corrected).

- `[x]` **P22.2 — `viralscan evidence` dispatch test** (`tests/test_evidence_subcommand.py`,
  new file, 18 tests) — `TestEvidenceParser` (15 tests: help, `_subcommand` attribute,
  required args, missing-arg exits, blast/virus/cores/viral-fasta/verbose/quiet flags) +
  `TestEvidenceDispatch` (3 tests: `test_calls_run_evidence`, routing to `run_evidence`,
  non-evidence subcommand does not call `run_evidence`).

- `[x]` **P22.3 — Planted-signal hostresponse test** (`tests/test_hostresponse.py`) —
  Appended `TestHostresponsePlantedSignal.test_planted_genes_rank_high_in_stability`: 200
  cells × 100 genes, 40 virus-positive cells, first 5 genes amplified 10× in positive cells;
  asserts ≥3 of 5 planted genes appear in top-10 by `stab_prob`.

- `[x]` **P22.4 — Full-depth validation** (operational) — SLURM array script written at
  `scripts/slurm_full_depth_validation.sh` (`sbatch --array=0-2`; `-c 8 --mem 32G -t 08:00:00`).
  Covers SRR20710641 HHV-6/10xv3, SRR12682296 EBV/10xv2, SRR8315713 HSV-1/DROPSEQ; ENA-first
  download, viralscan full-depth run, gate check, and `--summarize` helper for
  `BENCHMARK_COMPARISON_full_depth.tsv`. **Submitted 2026-06-24: array job 25082939, summarize
  job 25082940 (afterok). Waiting for cluster results.**
  **Run 2 (2026-06-25, array 25089684): task 0 (HHV-6b) = GATE FAIL (old pattern) but
  viral_summary.tsv written — 1,965/1,292,857 cells (0.152%). Task 1 (EBV) = OUT_OF_MEMORY
  (exit 0:125) at 2h12m — 32 GB insufficient for 61.6M BUS records; estimated peak ~70 GB
  (scales as 4.3 GB × 16.5× from HHV-6 baseline). Task 2 (HSV-1) = OUT_OF_MEMORY at 2h46m
  (36.5M BUS records, ~42 GB estimated). Next: resubmit both EBV + HSV-1 with --mem=128G.**
  **Run 3 (2026-06-25, array 25089827 --mem=128G): COMPLETE.
  Task 1 (EBV, job 25089827_1): COMPLETED at 4h00m, MaxRSS 70.6 GB. Result: 1,252,577 UMI,
  67,254 infected cells (≥1 UMI), 748,518 total, 8.985% infected, 102.3082 UMI/10k;
  ≥10 UMI super-expressors: 2,860 cells (0.382%). GATE PASS.
  Task 2 (HSV-1, job 25089827_2): COMPLETED at 3h33m, MaxRSS 64.4 GB. Result: 48,401 UMI,
  10,455 infected cells, 1,893,827 total, 0.5521% infected, 5.5104 UMI/10k.
  BENCHMARK_COMPARISON.md and docs/manuscript_draft.md Tables 3.2+3.3 updated (2026-06-25).**

- `[x]` **P22.5 — HHV-6 / HSV-1 divergence investigation** (RESOLVED 2026-06-25) —
  Full-depth P22.4 runs confirmed: HHV-6 0.152% (within Lareau range); HSV-1 0.55% (apparent
  divergence from Wyler 2019's 13–19%). Root-cause analysis (P22.5):
  **HSV-1: denominator artifact.** SRR8315713 confirmed = "5 hpi, Rep 2" (NCBI SRA GSM3511326).
  ViralScan reports over 1,893,827 unfiltered barcodes. Over called cells (≥1,000 total UMI):
  18.9% (evonk subsample, 4,571 cells) to 27.1% (full-depth, 4,414 cells) — consistent with
  published 13–19%. Only 1–2 cells exceed 8% viral fraction (Wyler's "high-expressor" threshold),
  indicating this replicate has low lytic burden. SRRs with 95–98% infection (SRR8315729–8315732)
  are later timepoints/higher-MOI conditions in GSE123782.
  **HHV-6: no divergence.** Full-depth result (0.152%) is within Lareau's 0.01–0.3% range.
  Bimodal GMM on log10(HSV-1 UMI+1) over called cells: low component at ~1 UMI (noise),
  high at ~18 UMI (signal); crossover ~2 UMI. At ≥2–5 UMI: 13.5–17.6% infected,
  matching Wyler's 13–19% exactly. Single-UMI counts (~420 cells) are multimapping noise.
  BENCHMARK_COMPARISON.md §Study 3 and Summary Table updated (2026-06-25).

- `[x]` **P22.6 — STARsolo comparison on EBV dataset** (operational) —
  CellRanger binary not available on the cluster. Substituted with STARsolo (STAR 2.7.11b,
  `starsolo` conda env) which is directly equivalent for CB+UMI counting. Scripts written:
  `scripts/slurm_starsolo_ebv_comparison.sh` (genome build + STARsolo run, `sbatch` directly)
  and `scripts/compare_starsolo_viralscan.py` (parse GeneFull filtered matrix, count EBV cells,
  emit comparison TSV). `BENCHMARK_COMPARISON.md` §STARsolo section added with methodology.
  **RESULT 2026-06-25 (job 25089721): 1,909 cells; 1,460 EBV ≥1 UMI (76.48%); 187 ≥10 UMI
  (9.80%). STAR uniquely mapped 86.8%. GeneFull mode captures latent EBV in essentially all LCL
  cells. `starsolo_p22_6/comparison_starsolo_vs_viralscan.tsv` written. BENCHMARK_COMPARISON.md
  and docs/manuscript_draft.md §3.3 updated.**

- `[~]` **P22.7 — Companion manuscript** (operational) — Draft at `docs/manuscript_draft.md`.
  Complete: Abstract, Introduction, full Methods, Results §3.1 (multimapping comparison data
  from BENCHMARK_COMPARISON.md), Results §3.2 benchmark table (all full-depth numbers filled
  2026-06-25: HHV-6b 0.152%, EBV 8.985%, HSV-1 0.5521%), Results §3.3 STARsolo + ViralScan
  full-depth numbers filled (ViralScan 748,518 cells / 67,254 EBV ≥1 UMI; matched-barcode P22.10
  complete), Results §3.4 host-response numbers (2026-06-27: EBV matched anchor, n=1,906;
  corrected multimap >=10 UMI: 1,179 positive / 727 negative; AUROC 0.845 ± 0.032; balanced
  accuracy 0.767 ± 0.030; 15 stable genes), Figure 1–2 generated in `docs/figures/`,
  Discussion, References.
  Commands:
  `NUMBA_CACHE_DIR=/tmp/viralscan_numba_cache MPLCONFIGDIR=/tmp/viralscan_mpl_cache PYTHONPATH=src /exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python scripts/hostresponse_ebv_matched.py --run-dir /exports/para-lipg-hpc/mdmanurung/viralscan_showcase/out_full_depth_wl/lcl_5lines/SRR12682296 --paper-barcodes /exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data/geo_GSE158275/GSM4796271_LCL_777_B958_UMI_barcodes.tsv.gz --output-dir results/hostresponse_ebv_matched --detection-threshold 10 --n-stab-iter 100 --stab-min-prob 0.6 --n-seeds 6`
  `NUMBA_CACHE_DIR=/tmp/viralscan_numba_cache MPLCONFIGDIR=/tmp/viralscan_mpl_cache PYTHONPATH=src /exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python scripts/make_manuscript_figures.py --matched-comparison results/matched_barcode_comparison.tsv --per-gene-comparison results/matched_barcode_comparison_per_gene.tsv --hostresponse-summary results/hostresponse_ebv_matched/hostresponse_summary.tsv --output-dir docs/figures`
  Outputs: `results/hostresponse_ebv_matched/hostresponse_summary.tsv`,
  `results/hostresponse_ebv_matched/hostresponse_metrics.csv`,
  `results/hostresponse_ebv_matched/Epstein-Barr_virus_gene_weights.csv`,
  `results/hostresponse_ebv_matched/Epstein-Barr_virus_stability.csv`,
  `docs/figures/figure1_workflow.{png,pdf}`, `docs/figures/figure2_benchmark.{png,pdf}`.
  **Pending sub-items:**
  - [ ] **P22.7a** — Fill author list + affiliations in `docs/manuscript_draft.md`
    (§ Author contributions; check ORCID for all co-authors).
  - [ ] **P22.7b** — Add GitHub URL (`https://github.com/…/ViralScan`) and data-availability /
    Zenodo DOI statement to the manuscript (Methods §Data Availability).
  - [ ] **P22.7c** — Choose target journal (Bioinformatics App Note / PLOS CompBio /
    GigaScience) and apply its style template; flip P22.7 `[~]` → `[x]` when
    submission-ready.
  Target journals: Bioinformatics Application Note, PLOS Computational Biology, GigaScience.

- `[x]` **P22.8 — mypy clean pass** — Ran mypy 2.1.0 against the 5 strict-mode modules
  (`viralscan.utils`, `viralscan.constants`, `viralscan.menu`, `viralscan.scripts.ncbi_fetch`,
  `viralscan.scripts.build_reference`). Fixed all 33 errors across 6 files:
  (1) `defaults.py`: annotated `DEFAULTS: dict[str, Any]` to fix 21 `runconfig.py` object-vs-typed
  assignment errors; (2) `pyproject.toml`: added `anndata` to `ignore_missing_imports`;
  (3) `menu.py:385`: `KbCountOutputs(rc)` → `KbCountOutputs(Path(rc.output))`;
  (4) `evidence.py`: cast `gzip.open` to `IO[str]`; changed `_parse_coverage_output` /
  `coverage_table` return type to `list[dict[str, str]]` (values are always strings);
  (5) `data_fetch.py:333`: `data_dir / str(fasta_name)` to force `Path`;
  (6) `evidence_run.py`: `_die` → `NoReturn`, added `argparse.Namespace` annotation,
  added `# type: ignore[no-untyped-call]` for untyped multimap helpers.
  Result: `mypy … → Success: no issues found in 5 source files`. 470 tests pass.

- `[x]` **P22.10 — Matched-barcode STARsolo ↔ ViralScan comparison** (COMPLETE 2026-06-25) —
  Re-ran STARsolo (job 25091356, 4h24m) and ViralScan (job 25091357, 4h11m, MaxRSS 70.5 GB)
  with the 10x v2 whitelist (737K). Anchor: 1,906 cells from GSM4796271 (LCL_777_B958).
  **Key results (1,906 shared cells, 0 dropout):**
  - STARsolo: 77.28% ≥1 UMI, 10.23% ≥10 UMI, 3.88% lytic (BZLF1/BRLF1/BHRF1)
  - ViralScan: 93.91% ≥1 UMI, 13.48% ≥10 UMI, 2.73% lytic (VS closer to published 2.2%)
  - Spearman r=0.45, Pearson r=0.42 (n=1,906, p<10⁻⁹⁶)
  - Per-gene: STAR captures LMP-1 (47,220 UMI; 74.4% cells) but misses entire EBNA family (all 0).
    ViralScan recovers EBNA-2 (716 UMI), EBNA-3A (817), EBNA-3B/C (274), EBNA-LP (37) but
    barely detects LMP-1 (99 UMI). Annotation-coverage asymmetry: 16 STAR vs 96 VS EBV entries.
    EBNA failure in STAR consistent with complex Wp/Cp poly-cistronic splicing not captured by GeneFull.
  Results: `results/matched_barcode_comparison{,_per_gene,_ebna_recovery}.tsv`
  §3.3 Table 3.3 updated; matched-barcode interpretation + EBNA recovery paragraph added.

- `[x]` **P22.9 — Publication checklist wrapper** (`scripts/publication_checklist.sh`) —
  Thin SLURM-chaining wrapper that submits the full-depth validation array (`--array=0-2`),
  wires the `--summarize` table-emission step as an `afterok` dependent job (1 CPU / 1 G /
  10 min), and submits the STARsolo EBV comparison independently in parallel — collapsing
  Steps 1-3 of the publication runbook into one command:
  `bash scripts/publication_checklist.sh` (or `--dry-run` to preview).
  Also fixed a pre-existing bug in `slurm_full_depth_validation.sh`: the `--summarize`
  block was positioned after the download/viralscan sections, so calling the script with
  `--summarize` would re-run the array job before emitting the table. Fixed by adding an
  early-dispatch guard right after the array declarations (paths already in scope; no
  kb/snakemake required). Dead duplicate block at the bottom removed.

---

## PR 18 — Tier 3 clean-code: config-key deduplication, main() decomposition, host_filter migration (2026-06-22)

Implements the three Tier 3 maintainability items from the clean-code review
(`/home/mdmanurung/.claude/plans/spicy-herding-whistle.md` findings #8–#10).

- `[x]` **P18.1 — Finding #8: Config-key list triplicated** (`runconfig.py`).
  Added `RunConfig.to_snakemake_config_args() -> list[str]` that derives the
  ``k=v`` list directly from the dataclass fields (bools → `"true"`/`"false"`,
  `None` → `""`, everything else → `str(v)`). This is the single authoritative
  serialisation of `RunConfig` → Snakemake wire format, eliminating the
  hand-maintained parallel list that caused the T1.2 EM-knobs bug.

- `[x]` **P18.2 — Finding #9: `main()` god function** (`menu.py`).
  Extracted two pure, testable helpers:
  - `_build_config_args(args, outs, index, transcripts, f1, s1, s2) -> list[str]` —
    constructs a `RunConfig.from_snakemake_config(...)` from per-sample paths +
    CLI args and returns `.to_snakemake_config_args()`. Replaces the 36-line
    manual list with a 5-line call.
  - `_write_sample_summary(outs, elapsed, n_transcripts, n_genes) -> None` —
    appends the runtime + reference-stat lines to `summary.txt`. The per-sample
    loop body shrank from ~77 lines to ~26 lines.

- `[x]` **P18.3 — Finding #10: `host_filter.py` on the raw-dict config path** (`scripts/host_filter.py`).
  Migrated `main(config: dict, ...)` → `main(config: RunConfig, ...)`. Dict-key
  access (`config["output"]`, `config.get("technology", "10xv3")`, etc.) replaced
  with attribute access (`config.output`, `config.technology`, etc.). Snakemake
  wiring updated from `load_config()` → `RunConfig.from_yaml()`. `load_config`
  import removed. `host_filter.py` is now the last scripts/ module that uses
  the `RunContext`/`RunConfig` pattern — no stale raw-dict consumers remain.

Verification: `PYTHONPATH=src python -m pytest tests/ -q` → **387 passed, 15 deselected**.
Files modified: `src/viralscan/runconfig.py`, `src/viralscan/menu.py`,
`src/viralscan/scripts/host_filter.py`, `tests/test_createconfig.py` (7 new),
`tests/test_cli.py` (5 new).

---

## PR 23 — Anellovirus into the standard combined reference (2026-06-24)

Makes the full 2,022 clareaulab anellovirus accessions part of the default
combined host+viral reference — not behind an opt-in flag.  Both the in-package
CLI (`build-ref`) and the bulk panel builder (`scripts/build_bundled_panel_ref.py`)
are covered.  Key design decisions: anello gene_ids are `{acc}_geneN` (generated
on-the-fly via `_genome_as_transcript_gtf`, same as any un-bundled viral accession);
grouping goes through `merged_name_map()` = `{**VIRUS_NAME_MAP, **anello_name_map()}`;
the 20 RefSeq anello already in the bundled panel are excluded from the anello fetch
(de-dup by `source == "clareaulab"`); >50% NCBI failure aborts, <50% logs and continues.

- `[x]` **P23.1 — `build_combined_reference` default-includes anello.**
  Added `include_anellovirus: bool = True` parameter to `build_combined_reference`
  (`build_reference.py:256`).  When true, loads `anellovirus.load_accession_table()`,
  de-dups against `virus_accessions`, fetches each via `_fetch_one` with fault tolerance
  (>50% failure → RuntimeError; otherwise log + continue), appends to `viral_fasta_path`.
  GTF generation is handled by the existing per-accession `_genome_as_transcript_gtf` step.

- `[x]` **P23.2 — `build_ref_main` passes `include_anellovirus` to combined path.**
  Old `--anellovirus store_true` routed to anello-only build.  New routing: only
  `--reference-panel anellovirus` goes to the anello-only path; `include_anellovirus`
  derives from `getattr(args, "anellovirus", True)` and is passed to `build_combined_reference`.

- `[x]` **P23.3 — `menu.py`: flip `--anellovirus` to `BooleanOptionalAction` default-on.**
  Changed from `action="store_true", default=False` to
  `action=argparse.BooleanOptionalAction, default=True` with updated help text explaining
  `--no-anellovirus` to skip and the `--reference-panel anellovirus` alternative.

- `[x]` **P23.4 — `build_bundled_panel_ref.py`: add Step 4b (anello fetch + GTF).**
  After the curated 195-virus FASTA download, loads clareaulab accessions (2,022),
  fetches each FASTA via `_fetch_one` with same >50% fault-tolerance, generates
  whole-genome GTF via `_genome_as_transcript_gtf`, appends FASTAs + GTFs to
  `combined.fa` / `combined.gtf` in Step 6.  Spot-check for `_geneN` entries added
  to the final verification block.  Fixed "decoy transcriptome" wording → "host
  transcriptome for coexpression".

- `[x]` **P23.5 — `bulk_viral_summarize.py`: `merged_name_map()` integration.**
  Imports `merged_name_map` from `viralscan.anellovirus`; builds `_name_map` once
  before the sample loop; passes `name_map=_name_map` to `virus_name_for_gene` so
  anello `{acc}_geneN` ids group by genus rather than falling through to raw ids.

- `[x]` **P23.6 — `bulk_viral_scan.sh`: fix stale "viral-only / no host decoy" comment.**
  Replaced with accurate framing: combined host+viral index; host reads compete for
  k-mers (no false-positive inflation); counts retained for coexpression analysis.

Verification done (2026-06-24):
- `PYTHONPATH=src python -m pytest tests/test_virus_grouping.py tests/test_build_reference.py tests/test_constants.py` → 51 passed.
- `build_combined_reference` signature confirmed (`include_anellovirus` default=True).
- `virus_name_for_gene("AB026929.1_gene1", name_map=merged_name_map())` → `"Betatorquevirus"` ✓.
- `--anellovirus/--no-anellovirus` present in `menu.py` (BooleanOptionalAction, default=True).

Operational steps remaining (need cluster + network):

- [x] **P23.op1** — Build panel index. **DONE** (2026-07-02, job 25138594,
  `slurm_build_panel_kbref.sh`). `viralscan_panel_ref/ref/panel.idx` 499 MB, `panel.t2g`
  470,533 entries (host ENST 465,769 + anello `_gene` 2,056) — verified via `kallisto inspect`.
  History (2026-07-02): three bugs surfaced building this panel, all now handled.
  (a) Earlier viral-only panel at `viralscan_bulk_gse128078/ref/panel.idx` (2,742 entries,
  0 ENST) was **incomplete** — June 24 build died at the Ensembl `current_gtf` 404
  (fixed 2026-07-01). (b) Re-run (job 25138575) fetched host + 194 curated + 2022 anello and
  wrote `combined.fa` (465,769 ENST + 2216 viral) but died at `kb ref` — `kb` not on PATH.
  (c) **Systemic bug**: `build_bundled_panel_ref.py` (and the `viralscan build-ref` CLI in
  `build_reference.py`) concatenate the Ensembl *chromosomal* host GTF with the *ENST cDNA*
  FASTA → `kb ref` would hang forever (identical to covid Stage 2 bug #1). Fix: regenerated a
  cDNA-level GTF from the existing `combined.fa` (no re-fetch) via `gen_combined_cdna_gtf.py`,
  validated locally (0 chromosomal seqnames, 0 orphan seqnames, 465,769 ENST + 2216 viral),
  and launched `slurm_build_panel_kbref.sh` with a fail-fast seqname guard + keep-first dedup
  fallback + artifact verification. **TODO (source fix)**: make `build_reference.py` /
  `build_bundled_panel_ref.py` emit a cDNA-level host GTF so the native CLI builds host+viral
  correctly — see new PLAN row P23.op1b below.
  Original one-shot `--wrap` recipe (now superseded by the two-stage fetch → kb-ref flow):
  ```bash
  WORKDIR=/exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref
  sbatch --job-name=build_panel_ref --cpus-per-task=2 --mem=16G --time=08:00:00 \
    -o $WORKDIR/build_panel_ref_%j.log -e $WORKDIR/build_panel_ref_%j.err \
    --wrap "NCBI_EMAIL=mikhael.manurung@gmail.com \
            PYTHONPATH=/exports/para-lipg-hpc/mdmanurung/ViralScan/src \
            /exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python \
            /exports/para-lipg-hpc/mdmanurung/ViralScan/scripts/build_bundled_panel_ref.py \
            --out $WORKDIR/ref"
  ```
  Expected runtime: ~4–6 h (downloads ~2,217 FASTA files then `kb ref`).
  Expected outputs: `$WORKDIR/ref/panel.idx`, `$WORKDIR/ref/panel.t2g`, `$WORKDIR/ref/panel.fa`.

- [x] **P23.op1b** — **Source fix for the host-GTF/cDNA-FASTA mismatch (systemic).** DONE
  (2026-07-02). Added `host_cdna_as_gtf()` to `build_reference.py` (seqname = ENST,
  gene_id = the `gene:ENSG…` field, coords 1..len); `build_combined_reference` (the CLI
  path) and `build_bundled_panel_ref.py` Step 6 now emit a cDNA-level host GTF instead of
  concatenating the chromosomal one. Regression test class `TestHostCdnaAsGtf` +
  updated `TestBuildCombinedReference` assert 0 chromosomal/scaffold seqnames leak and
  every GTF seqname is a FASTA header. Full suite green (505 passed, 15 deselected).
  Both `scripts/build_bundled_panel_ref.py` (Step 6) and the native CLI core
  `src/viralscan/scripts/build_reference.py` (`build_combined_reference`) concatenate the raw
  Ensembl *chromosomal* host GTF (seqnames 1/2/X) with the *ENST cDNA* FASTA, then call
  `kb ref combined.fa combined.gtf` — which hangs forever at "Splitting genome" because no
  chromosomal seqname matches a cDNA header. (Discovered 2026-07-02; also the root of covid
  Stage 2 bug #1.) Fix: synthesize a cDNA-level host GTF from the cDNA FASTA headers
  (seqname = ENST, gene_id = the `gene:ENSG…` field, coords 1..len) — the logic already
  exists in `covid_viralscan/scripts/gen_combined_cdna_gtf.py`; lift it into the package.
  Add a regression test asserting the emitted host GTF has 0 chromosomal seqnames and that
  every seqname is a FASTA header. Until this lands, host+viral builds must use the
  `gen_combined_cdna_gtf.py` + `slurm_build_panel_kbref.sh` workaround.

- [x] **P23.op2** — Verify anellovirus entries in `panel.t2g`. **DONE** (2026-07-02):
  `_gene` rows = 2,056 (>0 ✓), ENST rows = 465,769 (host present ✓), total 470,533 (≥200k ✓).
  Also fixed `scripts/bulk_viral_scan.sh` `REFDIR` to point at the new `viralscan_panel_ref/ref`
  panel (was defaulting to the host-less `viralscan_bulk_gse128078/ref` build). Verify command:
  ```bash
  WORKDIR=/exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref
  grep -c "_gene" $WORKDIR/ref/panel.t2g   # must be > 0 (anello _geneN entries)
  grep -c "ENST"  $WORKDIR/ref/panel.t2g   # host transcripts present
  wc -l           $WORKDIR/ref/panel.t2g   # total entries (≥ 200k expected)
  ```

- [ ] **P23.op3** — Format probe: run `kb count -x BULK` on ONE GSE128078 sample to confirm
  output layout before submitting the array. Inspect `counts_unfiltered/` paths and adjust
  `bulk_viral_summarize.py` column parsing if the layout differs:
  ```bash
  # Pick any one SRR accession from scripts/fetch_reference_strategy_fastqs.py or GSE128078
  # then run kb count manually with: -x BULK -i $WORKDIR/ref/panel.idx -g $WORKDIR/ref/panel.t2g
  # and inspect ls -R counts_unfiltered/
  ```

- [ ] **P23.op4** — Pilot array (6 samples, ~1 h each):
  ```bash
  cd /exports/para-lipg-hpc/mdmanurung/ViralScan
  sbatch --array=0-5 scripts/bulk_viral_scan.sh
  # After completion: check sacct -j <JOB> --format=JobID,State,MaxRSS,Elapsed
  # and verify n_pseudoaligned > 0 in at least 3/6 sample logs
  ```

---

## Bulk exploratory scan — GSE128078 (ME/CFS whole-blood, in progress 2026-06-24)

Companion scripting under `scripts/` — does NOT touch the single-cell CLI.
Uses ViralScan's combined host+viral panel with `kb count -x BULK` (kb-python)
for bulk RNA-seq quantification.  After PR 23 the panel now includes the full
~2,022 clareaulab anellovirus accessions alongside the curated 195-virus panel
and Ensembl host cDNA.

Study: GSE128078 / SRP187984 — 99 whole-blood samples, ME/CFS patients + controls
(Illumina HiSeq 2500, paired-end). Scientific goal: exploratory viral-reactivation
screen. FASTQs downloaded from ENA (no sra-tools dependency).

Scripts:
- `scripts/build_bundled_panel_ref.py` — one-time index build (host + curated 195
  + ~2,022 anellovirus; runs `kb ref`). See Step 4b for anello fetch.
- `scripts/bulk_viral_scan.sh` — SLURM array (one task per sample; ENA download
  + `kb count -x BULK`). Pilot: `--array=0-5`; full: `--array=0-98`.
- `scripts/bulk_viral_summarize.py` — aggregates kb count outputs to
  `bulk_viral_summary.tsv` (per-virus RPM, per-sample). Uses `merged_name_map()`
  so anello gene_ids group by genus. Format probe required before first run.
- `scripts/ttv_public_datasets.json` — curated catalog (v1.1) of 16 public
  TTV/anellovirus datasets with SRA accessions, access flags, and download
  recipes. Confirmed open: arze2021_chm (PRJNA679286), wang_oral_srr2037085
  (SRR2037085), tisza2020_elife (PRJNA396064/393166), kraberger2020_brain
  (SRR12450126), devlaminck2013_cell (SRP032345), asct_mngs_2018 (PRJNA504035).

Order of operations:

- [ ] **B1** — Build panel index (= P23.op1). See exact `sbatch --wrap` command in PR 23 section.
  Gate: `ref/panel.idx` and `ref/panel.t2g` must exist and pass the anello grep check (P23.op2).

- [ ] **B2** — Format probe (= P23.op3). Run `kb count -x BULK` on ONE sample; verify
  `counts_unfiltered/` layout matches what `bulk_viral_summarize.py` expects:
  ```bash
  python scripts/bulk_viral_summarize.py --help  # check expected --counts-dir structure
  ls counts_unfiltered/  # should have cells_x_genes.barcodes.txt + cells_x_genes.genes.txt + .mtx
  ```

- [ ] **B3** — Pilot array (6 samples):
  ```bash
  cd /exports/para-lipg-hpc/mdmanurung/ViralScan
  sbatch --array=0-5 scripts/bulk_viral_scan.sh
  # After: check n_pseudoaligned > 0 in 3+ sample run_info.json files
  # and that *.bus files are non-empty (ls -lh <sample>/counts_unfiltered/*.bus)
  ```

- [ ] **B4** — Summarize pilot output:
  ```bash
  PYTHONPATH=src python scripts/bulk_viral_summarize.py \
      --samples-dir /exports/para-lipg-hpc/mdmanurung/viralscan_bulk_gse128078 \
      --t2g /exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref/ref/panel.t2g \
      --output results/bulk_viral_summary.tsv
  # Sanity: grep -c "Alphatorquevirus\|Betatorquevirus\|Gammatorquevirus" results/bulk_viral_summary.tsv
  # Sanity: python -c "import pandas as pd; df=pd.read_csv('results/bulk_viral_summary.tsv', sep='\t'); print(df.shape, df.head())"
  ```

- [ ] **B5** — Full run (99 samples) after pilot is sane:
  ```bash
  sbatch --array=0-98 scripts/bulk_viral_scan.sh
  # Re-run bulk_viral_summarize.py on all 99 samples after completion
  # Expected: ~6–10 h per sample (whole-blood, 30M reads each)
  ```

---

## Release Readiness (2026-07-02) — feature-complete + release before benchmarking

Goal: PyPI + bioconda + container release, full quality-gate hardening. Benchmarking
items (P23.op2–4, B1–B5) and manuscript (P22.7) are **post-release, out of this gate**.

### Phase 0 — verification
- [x] **RR0.1** Baseline: 505 passed, 15 deselected (`test_viralscan` env, 2026-07-02).
- [x] **RR0.2** build-ref combined-reference gate (P23.op1b): fix present + regression-tested;
  independently confirmed via a synthetic `kb ref` smoke on a cDNA-level GTF — completes in
  ~14 s, exit 0, valid t2g (no "Splitting genome" hang). Full real-panel index (P23.op1) is
  post-fix confirmation / benchmarking, SLURM-gated — not a tool-feature blocker.

### Phase 1 — correctness & consistency
- [x] **RR1.1** Default multimap docs reconciled to code (`equal`, per PR 17) + cross-homology
  warning added (CHANGELOG, README ×3). Decision: keep `equal` default + prominent guidance.
- [x] **RR1.2** Canonical repo identity: `pyproject.toml` URLs `emmaevonk`→`mdmanurung`
  (matches README/CITATION). No `emmaevonk` refs remain.
- [x] **RR1.3** `scripts/evidence_run.py` migrated to `RunConfig.from_yaml` (last legacy
  raw-dict script; also fixed the manual `str(run_dir)+"/"` path concat via `os.path.join`).
- [x] **RR1.4** Early `kb` preflight in `build_ref_main` (warns before the long download,
  not after). `evidence` already guards `bustools` inline; `hostresponse` needs no external tools.
- [x] **RR1.5** `RunConfig.from_yaml` normalizes the `output` trailing separator (frozen-safe,
  pre-construction) + `TestFromYamlTrailingSlash` (3 tests). Full suite 508 passed, 15 deselected.

### Phase 2 — hygiene
- [x] **RR2.1** Single-source version: `src/viralscan/__init__.py:__version__ = "2.4.0"`;
  `pyproject.toml` reads it dynamically (`[tool.setuptools.dynamic] version = {attr=…}`).
  Added `viralscan --version` flag.
- [x] **RR2.2** Cut `[Unreleased]`→`[2.4.0] - 2026-07-02`; synced Dockerfile/Singularity/CITATION
  to 2.4.0 (+ CITATION date-released). Authorship left as-is (open, tied to P22.7a).
- [x] **RR2.3** Rebuilt: `python -m build` → viralscan-2.4.0 wheel+sdist; `twine check` PASSED;
  wheel has all modules + anellovirus_accessions.tsv, 0 GTFs; runtime `__version__`=2.4.0.
  (dist/ is gitignored; release.yml builds fresh on tag.)

### Phase 3 — docs
- [x] **RR3.1** Removed stale `getting_started.ipynb` (maintained tutorials live in
  `docs/vignettes/`); updated the CLAUDE.md pitfall note.
- [x] **RR3.2** README: conda install now `pip install .` (note `-e .` for dev); added the
  `connection_pool`/pip caveat. Badges already pointed to `mdmanurung`.
- [x] **RR3.3** `docs/conf.py` `exclude_patterns` now drops manuscript_draft / review-* /
  write-docs-prompt / showcase_runbook from the public build.
- [x] **RR3.4** Added human-facing `CONTRIBUTING.md` (dev env, tests, gates, PR flow).

### Phase 4 — quality gates
- [x] **RR4.4** Ruff ruleset expanded to `E4,E7,E9,F,I,B,UP,SIM` (isort/bugbear/pyupgrade/simplify);
  60+ auto-fixes + `ruff format` (60 files); cosmetic UP007/UP045/SIM117 ignored with rationale;
  substantive findings fixed (SIM115 file handles, B017 broad-raises). `ruff check` clean.
- [x] **RR4.3** mypy: removed `ignore_errors=true` for the 7 Snakemake scripts (+ research
  `reference_strategy`) → now type-checked for real errors (annotation-completeness codes relaxed).
  Fixed 11 latent issues (None-narrowing asserts documenting Snakefile invariants, stale
  `# type: ignore` removals, loop-var type clash). `mypy -p viralscan` clean (29 files).
- [x] **RR4.1** CI `integration` job (micromamba installs kb-python/snakemake/minimap2/
  samtools/blast) runs `pytest -m "integration and not network"` — verified 11 passed locally
  with tools on PATH. Also scoped ruff to the package (`extend-exclude` research/build dirs).
- [x] **RR4.2** Coverage floor `--cov-fail-under=60` in the CI test step (measured 65%).
- [x] **RR4.5** CI `security` job: `bandit --severity-level high` (0 high-sev; subprocess/urllib
  are low/medium and expected) + `pip-audit` (advisory).
- [x] **RR4.6** Added a `research` pytest marker; tagged the two repo-root-script tests
  (hostresponse_ebv_matched, make_manuscript_figures) and excluded them from the default
  hermetic run. Default suite: 503 passed, 20 deselected.

### Phase 5 — packaging
- [x] **RR5.1** bioconda recipe `conda-recipe/meta.yaml` (+ README) — noarch python, PyPI
  source, entry point, run deps (kb-python/kallisto/bustools/star/snakemake-minimal +
  py stack), `viralscan --version`/`--help` tests. Renders/validates. Local `conda build`
  + bioconda PR are **user-gated** (need conda-build + the PyPI sha256).
- [x] **RR5.2** Container: added a `container` job to `release.yml` that builds and pushes
  `ghcr.io/mdmanurung/viralscan:{latest,vX.Y.Z}` on tag (packages:write via GITHUB_TOKEN);
  fixed the release version-check to read the dynamic `__version__`; hardened `.dockerignore`
  (excludes references/benchmark/analysis/data). Image build itself is **user-gated** (no
  container runtime here).

### Phase 6 — release (USER-GATED — everything below needs your action)

Release is **ready**: full suite 503 passed; `ruff check .` / `ruff format --check .`
clean; `mypy -p viralscan` clean; `python -m build` → viralscan-2.4.0, `twine check`
PASSED; `viralscan --version` → 2.4.0; bandit high-sev clean.

- [ ] **RR6.1** Open a PR `claude/multimap-memory-and-showcase` → `main` and confirm CI
  is green (lint + test matrix + new integration + security jobs). Merge.
- [ ] **RR6.2** `git tag v2.4.0 && git push origin v2.4.0` → `release.yml` builds + publishes
  to PyPI (needs the **PyPI Trusted Publisher** configured for project `ViralScan`) and
  builds+pushes the ghcr container.
- [ ] **RR6.3** Post-publish smoke: in a clean env, `pip install ViralScan==2.4.0 && viralscan --version`
  and `viralscan data fetch`.
- [ ] **RR6.4** Archive the GitHub release on Zenodo for a **software DOI** (distinct from the
  data DOI 10.5281/zenodo.20112332); add it to `CITATION.cff` (`identifiers:`) and the README.
- [ ] **RR6.5** (optional) bioconda PR: fill `conda-recipe/meta.yaml` `source.sha256` from the
  PyPI sdist and submit to bioconda-recipes.

---

## v2.5 Scientific-Hardening (2026-07-02) — fold analysis-layer rigor into the package

**Motivation.** The v2.4.0 release covers *quantification* and is release-ready. This
session's EBV / HHV-6B / HSV-1 / covid analyses stress-tested the *scientific-analysis*
layer (`hostresponse`, detection denominators, chemistry handling) and exposed
correctness gaps: every workaround under `analysis/hostresponse_ebv_matched/scripts/`
and `scripts/` is a feature the package lacks. Each item below names the external script
that already implements the methodology, to be folded into the package. Full rationale:
`.living/decisions.md` → "Feature-completeness gap analysis" (2026-07-02); findings
F-001, F-003, F-004, F-005.

Post-release track — does **not** block v2.4.0. Tier 1 is correctness-affecting
(the package's flagship host-response result is depth-confounded today); Tier 2/3 are
capability/robustness enhancements.

**Real-data reconciliation (2026-07-03) — SH1.1–1.3 VERIFIED, not just unit-tested.**
Ran the in-package `run_hostresponse` on the showcase EBV matrices
(`results/hostresponse_ebv_matched/{host_only_matched,ebv_burden_matched}.h5ad`, 1906
cells) three ways; it reproduces the external scripts' AUC pattern to 3 decimals:
`--label raw` → model AUC **0.866**, depth-alone **0.967** (depth beats the model — the
F-001/F-003 confound); `--label cpm` → model **0.672**, depth-alone 0.529; `--depth-match`
→ model **0.680**, depth-alone 0.471 (near chance). The package now yields the honest
AUC ~0.67 with no external script — the SH1.1–1.3 definition-of-done. (`/tmp/reconcile_ebv.py`.)

### Tier 1 — correctness (package can produce *misleading* results today)
- [x] **SH1.1** `hostresponse` depth-robustness reporting — DONE 2026-07-03. Always-on:
  `_depth_alone_auc` (AUC from log-depth ALONE under the identical balanced/top-depth
  split → `depth_alone_auc_mean/sd` in `hostresponse_metrics.csv`), `_per_gene_evalues`
  (depth-adjusted OR + Ding & VanderWeele E-value per stable gene →
  `<virus>_depth_diagnostics.csv`; `n_genes_evalue_ge2`), and a log WARNING when
  depth-alone AUC ≈ model AUC. sklearn-only (no statsmodels dep); C=1.0 keeps E-values
  conservative + stable under quasi-separation. 8 new tests inc. synthetic depth-only
  guard (adjusted OR≈1, E<1.8). Fold from `depth_confounder_check.py`. (F-001/F-003.)
- [x] **SH1.2** Depth-independent label + depth-matched design — DONE 2026-07-03. Opt-in
  (defaults reproduce prior behaviour): `--label {raw,cpm,fraction}` (`_virus_presence_label`:
  cpm/fraction = prevalence-matched top viral-per-host-UMI; cpm≡fraction ranking, host-only
  denominator) and `--depth-match` (`_depth_match_indices`: coarsened-exact depth matching,
  disables the in-split top-depth filter via a new `top_depth_frac` param). `label`/`depth_matched`
  recorded in `hostresponse_metrics.csv`. Wired through `run_hostresponse`, the standalone CLI,
  and `viralscan hostresponse`. 13 new tests. Fold from `depth_matched_reanalysis.py` +
  `cpm_label_crosscheck.py`. Full suite 531 passed.
- [x] **SH1.3** `%mito` control in `hostresponse` — DONE 2026-07-03. On by default
  (`control_mito`, `--no-mito-control` to disable): per-cell %mito (from RAW counts,
  before normalization) added as a covariate to the per-gene E-values, so a mito-QC
  artifact cannot pass as a host-response gene. `_mt_gene_mask` detects MT genes by
  Ensembl ID (13 protein-coding) or `MT-`/`mt-` symbol prefix; leave-one-out drops an
  MT gene from its own %mito covariate (the MT-ND4L self-suppression fix); zero-variance
  %mito is skipped; no-op when no MT genes. `mito_controlled` recorded per virus. 6 new
  tests. Fold from `go_enrichment.py` mito control.
- [x] **SH1.4** Whitelist/chemistry preflight — DONE 2026-07-03. New
  `viralscan/whitelist_preflight.py`: `whitelist_match_rate`/`check_whitelist` sample R1,
  extract the CB (via `evidence.cb_umi_geometry`, so a wrong `--technology` also shows up),
  and report the whitelist match rate. Wired two ways: (a) an automatic best-effort preflight
  in the quant path that WARNs when an explicit `--whitelist` is given and the rate is low
  (never fatal); (b) a standalone `viralscan check-whitelist` diagnostic (exit 1 on mismatch).
  Catches the F-005 silent failure. 11 new tests. Full suite 548 passed.
- [x] **SH1.5** Called-cell denominator — DONE 2026-07-03 (delivered by the parallel
  `cellcalling` commit `317c04a`, verified here). `src/viralscan/scripts/cellcalling.py`
  provides `knee` (dependency-free default), `emptydrops` (DropletUtils via `emptydrops.R`),
  and `external` (CellRanger/STARsolo list) cell-calling; `detection.compute_stats` now
  reports `*_called` stats (n_called_cells, infected_called, pct_infected_called) alongside
  the all-barcode numbers in `viral_summary.tsv`, surfacing the HSV-1 denominator artifact
  by default. 23 cellcalling+detection tests pass. Follow-up (minor): expose the
  `cell_calling` method + emptydrops params as CLI flags (currently config-`getattr`, knee default).

### Tier 2 — capability (needed for the analyses; currently external)
- [x] **SH2.1** Gene-symbol annotation — DONE 2026-07-03. Opt-in `--gene-symbols`
  (network, best-effort): `_map_ensembl_to_symbols` (mygene.info, folded from
  `go_enrichment.py`) + `_add_symbol_column` insert a `symbol` column next to `gene` in
  the gene_weights / stability / depth_diagnostics CSVs. `_looks_like_ensembl` guards
  (skips when genes are already symbols); silent fallback to IDs on any network error.
  Wired through `run_hostresponse`, standalone CLI, and `viralscan hostresponse`. 6 new
  (offline) tests. Full suite 553 passed.
- [x] **SH2.2** Genome-wide depth-adjusted differential test + GO — DONE 2026-07-03. Opt-in
  `--differential`: `_genome_wide_differential` (folded from `go_enrichment.py:partial_assoc`)
  residualizes every gene AND the label on `log(depth)` (+ %mito when on), reports per-gene
  `partial_r`/`p_value`/`fdr`/`direction` over ALL features → `<virus>_differential.csv`
  (`n_differential_fdr05` in metrics). BH-FDR implemented in-package (`_bh_fdr`, no statsmodels
  dep). With `--enrichment`, feeds the FDR<0.05 genes (symbol-mapped when available) to Enrichr.
  Wired through CLI + subcommand. 4 new tests (BH-FDR bounds, recovers true gene while adjusting
  away a depth proxy, integration writes genome-wide table). Full suite 557 passed.
- [ ] **SH2.3** HHV-6A/6B contig-level disambiguation — **DEFERRED (multi-day, reference-level).**
  `virus_grouping` resolves gene-ID *prefixes* (`HUM_HERP6B` vs `HUM_HERP6`), but the real
  ambiguity is reads that pseudo-align equally to the shared 6A/6B contigs — resolving that
  needs per-read alignment evidence against a curated 6A-vs-6B divergent-region model, not a
  naming rule. Correct scope is a reference-build + EM-allocation change, not a quick patch;
  the `equal`/`em` multimap methods already bound the effect. Tracked for a dedicated PR.
- [ ] **SH2.4** Per-cell EM — **DEFERRED (research algorithm).** The current EM resolves
  multimappers over a *global* transcriptome pool; per-cell EM (re-estimating allocation within
  each cell) is a stated manuscript limitation and a genuine algorithm-design task (per-cell
  sparsity, convergence, runtime at 10^5–10^6 cells). Needs its own design + validation PR.
- [x] **SH2.5** BULK mode — DONE 2026-07-03 (removed the unsupported claim). The
  `__init__` docstring said "single-cell/bulk RNA-seq" but bulk is not supported
  (`cb_umi_geometry` has no BULK entry; host-filter/evidence raise `ValueError` without
  a CB/UMI geometry). Corrected the docstring to state single-cell only + why. Full bulk
  support (a no-barcode counting path) is a separate large feature, intentionally deferred.

### Tier 3 — QC / robustness
- [~] **SH3.1** Ambient-RNA / doublet / `%mito` QC — **PARTIAL.** `%mito` is now computed and
  used as a covariate in `hostresponse` (SH1.3), and knee/emptyDrops cell-calling landed
  (SH1.5). Ambient-RNA (SoupX/CellBender) and doublet (Scrublet) correction remain **DEFERRED**
  — each adds a heavy dependency and a pipeline stage; scope as an optional `viralscan qc`
  module in a dedicated PR rather than bolting onto detection.
- [~] **SH3.2** Cell-level BAM for genome-browser inspection — **PARTIAL / DEFERRED.** kallisto+
  bustools emit no BAM, so a full per-cell BAM needs a different aligner path. The intent
  (inspect the reads behind a viral call in IGV) is already served for *viral* reads by
  `viralscan evidence` (minimap2 re-alignment of the traced (CB,UMI) reads → sorted BAM).
  A transcriptome-wide cell-barcoded BAM is out of scope for the pseudoalignment pipeline;
  deferred to a possible STARsolo-backed alternate path.
