# ViralScan Improvement Plan — Implementation Tracker

This file is the single source of truth for what's done and what still needs
doing. **After every implementation step, tick the relevant checkbox and update
the "Next up" pointer.** Do not mark an item done until the tests pass.

Second-pass audit completed 2026-05-08. All prior PR claims re-verified against
the actual codebase; status corrected where PLAN and code diverged.

Branch: `claude/review-repo-improvements-Sg4Th`
Test command: `PYTHONPATH=src python -m pytest tests/ -q` → 225 passed, 8 deselected.

---

## Status legend

- `[x]` — done and verified against the codebase
- `[~]` — partially done / obsolete sub-item
- `[ ]` — not started
- `[!]` — blocked / needs external action

## Next up

→ All tracked tasks complete.

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
