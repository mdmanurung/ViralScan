# ViralScan Improvement Plan — Implementation Tracker

This file is the single source of truth for what's done and what still needs
doing. **After every implementation step, tick the relevant checkbox and update
the "Next up" pointer.** Do not mark an item done until the tests pass.

Second-pass audit completed 2026-05-08. All prior PR claims re-verified against
the actual codebase; status corrected where PLAN and code diverged.

Branch: `claude/multimap-memory-and-showcase`
Test command: `PYTHONPATH=src python -m pytest tests/ -q` → 352 passed, 10 deselected (scvi env; 2026-06-21).

---

## Status legend

- `[x]` — done and verified against the codebase
- `[~]` — partially done / obsolete sub-item
- `[ ]` — not started
- `[!]` — blocked / needs external action

## Next up

→ PR 15 Run-context refactor — COMPLETE (steps 1–4 done). Mirror count 5 → 0.
→ **Session 2026-06-21 (showcase validation) findings all closed:** S0 committed `99db0e8`;
  S2 (pipeline halt) fixed and regression-tested; S1/S6 (geometry) unified on
  `evidence.cb_umi_geometry`; S5 parser unit tests added. Task 4 re-verified end-to-end.
→ S5 tail (live evidence-chain validation) complete: `tests/integration/test_evidence_chain.py`
  validates minimap2→BAM→samtools-coverage→BLAST on synthetic reads (no public data, skips when
  binaries absent). S5 flipped to `[x]`.
→ S3 closed: showcase_runbook.md + BENCHMARK_COMPARISON.md already committed; checkbox flipped
  `[~]`→`[x]`. All S0–S6 findings are now `[x]`.

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
