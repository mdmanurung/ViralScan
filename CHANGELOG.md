# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- **`evalue_flag` column in `<virus>_depth_diagnostics.csv`** — each stable gene is
  now classified as `fragile` (E < 1.5), `moderate` (1.5 ≤ E < 3), or `robust` (E ≥ 3)
  alongside the existing depth-adjusted odds ratio and E-value.
- **`model_auc_depth_adjusted_mean/sd` in `hostresponse_metrics.csv`** — AUC of the
  stable-gene panel logistic with `log1p(host_depth)` added as a covariate (same balanced
  split as the headline model and the depth-alone baseline). Shows how much of the
  stable-gene predictive power survives explicit depth adjustment.
- **`--hostresponse-label`, `--hostresponse-depth-match`, `--hostresponse-control-mito`,
  `--hostresponse-differential`** now exposed on the main `viralscan` pipeline CLI and
  persisted in `RunConfig` / `DEFAULTS` / the Snakemake entry block. Previously these
  knobs were only accessible via `viralscan hostresponse` standalone.

### Changed
- **Default `--multimap-method` is now `host-conservative`** (was `equal`). For viral
  detection on combined host+virus references, specificity is prioritised: host-virus
  ambiguous equivalence-class mass is kept out of primary viral counts by default. Pass
  `--multimap-method equal` for a fast unbiased first pass, or `em` for iterated allocation.

## [2.5.0] - 2026-07-03

### Added — host-response scientific-hardening (v2.5)

- **Depth-confound diagnostics (always on).** `hostresponse` now reports the AUC
  from sequencing depth *alone* (`depth_alone_auc` in `hostresponse_metrics.csv`)
  next to the model AUC, plus per-gene depth-adjusted odds ratios and Ding &
  VanderWeele **E-values** (`<virus>_depth_diagnostics.csv`). Surfaces the fact
  that the raw `counts >= threshold` label tracks library size.
- **Depth-independent labels / design.** `--label {raw,cpm,fraction}`
  (depth-normalized, prevalence-matched positive call) and `--depth-match`
  (coarsened-exact depth-matched cohort) give a depth-independent host-response
  estimate.
- **`%mito` control** (`--mito-control`, default on) — per-cell mitochondrial
  fraction added as a covariate to the per-gene E-values, with leave-one-out for
  mitochondrial genes.
- **`--gene-symbols`** — annotate host-response CSVs with HGNC symbols (mygene.info).
- **`--differential`** — genome-wide, depth/%mito-adjusted differential-expression
  table (`<virus>_differential.csv`; partial correlation, p-value, BH-FDR, direction).
- **`viralscan check-whitelist`** — barcode/whitelist chemistry-mismatch preflight;
  the main run also warns automatically when an explicit `--whitelist` is supplied.
- **Called-cell denominators** — `viral_summary.tsv` reports `pct_infected_called`
  (over knee/emptyDrops/external called cells) alongside the all-barcode rate.

### Changed

- Corrected the package docstring: ViralScan supports single-cell RNA-seq only
  (bulk was never supported).

## [2.4.0] - 2026-07-02

### Fixed
- Release-readiness pass: `evidence_run.py` migrated to the typed `RunConfig`
  (last script on the legacy raw-dict path); `RunConfig.from_yaml` now normalizes
  the `output` trailing separator; early `kb` preflight in `build-ref`; canonical
  GitHub URLs (`mdmanurung/ViralScan`); single-source package version
  (`viralscan.__version__`, read dynamically by `pyproject.toml`).

### Added
- **`viralscan hostresponse` subcommand** — run host-response analysis on a
  completed viralscan output directory without re-running Snakemake. Trains
  per-virus L2 logistic regression models (Luebbert et al. 2026) and runs
  randomized Lasso stability selection to identify stably virus-associated host
  genes. Optional `--enrichment` flag runs pathway enrichment via
  `gget.enrichr` (requires `pip install "viralscan[enrichment]"`). The same
  analysis also runs inline during a full `viralscan` run when `--host-h5ad`
  is supplied.
- **CLI reference docs** for `viralscan evidence`, `viralscan rerun-multimap`,
  and `viralscan hostresponse` — all three subcommands were shipped but
  previously undocumented.
- **`README.md` Limitations section** — documents five known limitations: FPR/FNR
  not yet characterized, global-pool EM (vs per-cell), untested chemistries beyond
  10xv2/v3/Drop-seq, cross-homology inflation for HHV-6/KDM2A, and ambient RNA
  not corrected.
- **EM global-pool caveat in `em_gene_abundances` docstring** — explains that EC
  counts are pooled across all cells before EM; returned theta is transcriptome-wide,
  not per-cell; differs from alevin-fry/STARsolo; faster but ignores cell-to-cell
  abundance variation.
- **`tests/test_evidence_subcommand.py`** (18 tests) — parser tests for all
  `viralscan evidence` flags and dispatch tests confirming `main()` routes
  `_subcommand="evidence"` to `run_evidence()`.
- **Planted-signal test for `run_hostresponse`** (`TestHostresponsePlantedSignal`)
  — synthetic 200-cell matrix with 40 virus-positive cells and 5 genes amplified
  10× in positive cells; asserts ≥3 planted genes appear in top-10 by stability
  probability.

### Fixed
- sklearn ≥ 1.8 `FutureWarning` in host-response stability selection: switched
  to `solver='saga', l1_ratio=1.0` on sklearn ≥ 1.8 (was `penalty='l1'`, now
  deprecated); older sklearn still uses `penalty='l1', solver='liblinear'`.

### Changed
- Multimapping is now a Snakemake checkpoint with a new `viralscan rerun-multimap`
  subcommand to switch allocation methods on an existing run without re-running
  `kb count`. The default `--multimap-method` is `equal` (fastest; enables instant
  in-place layer swaps). `host-conservative` remains available and is **recommended
  when host-virus cross-homology matters** (e.g. HHV-6 / *KDM2A* / *DR1*), where it
  keeps host-ambiguous equivalence-class mass out of primary viral counts.

---

## [2.3.0] - 2026-05-13

### Added
- Sphinx + MyST documentation skeleton under `docs/`; Read the Docs config at `.readthedocs.yaml`.
- `CITATION.cff` for software citation.
- `environment.yml` for reproducible conda environments.
- `Dockerfile` (mamba-based) and `Singularity.def` for containerised HPC runs.
- `tests/test_errorhandler.py` — direct unit tests for every `errorhandler()` branch.
- `tests/conftest.py` — session-level `pyfiglet` stub so unit tests pass without the optional dep.
- Type annotations on all public functions in `menu.py`, `utils.py`, and `constants.py`.
- Codecov integration: coverage XML uploaded in CI; badge added to README.
- `viralscan data fetch` for downloading the external viral annotation panel from Zenodo.
- Ambiguity-aware multimapper evidence outputs and diagnostic AnnData layers.
- Optional host pre-subtraction support before viral quantification.
- Cell-type enrichment tables and HTML report section.

### Fixed
- Comma-separated custom FASTA/GTF reference inputs are now materialized before `kb ref`.
- Comma-separated custom GTF files are all parsed during viral accession discovery.
- CI dependency installation now matches the modules imported by tests and type checks.
- Release/container install metadata now targets the local `2.3.0` package.

---

## [2.2.0] - 2026-05-08

### Added
- **PR 14 Bug fixes & performance**
  - C1: Fixed double-counting of unique reads in multimapping mode (`multimap.py`).
  - C2: Fixed detection threshold off-by-one (`>` → `>=` in `detection.py`).
  - C3: Replaced O(n) `list.index()` scan in EC parsing with O(1) `dict` lookup (`multimap.py`).
  - C4: Replaced `iterrows()` with `itertuples()` for ~3× faster BUS iteration (`multimap.py`).
  - C5: Changed `2> /dev/null` → `2>&1` in Snakefile so "no reads pseudoaligned" is capturable.
  - C6: Changed `mkdir` → `mkdir -p` in Snakefile to survive re-runs.
  - C7: Removed redundant `f.close()` after `with` block in `analysis.py`.
  - C8: Fixed file-handle leak in `detection.py` by wrapping `found_genes_file` in `with`.
  - C9: Consolidated conflicting `var_names` resolution blocks in `umap.py`; added O(1) index.
  - C10: Added `sc.pp.highly_variable_genes()` before PCA in both UMAP branches.
- **PR 11 Interpretation & Reporting**
  - Structured TSV outputs: `viral_summary.tsv` and `per_cell_viral.tsv`.
  - Self-contained HTML report via Jinja2 template with embedded base64 PNG plots.
  - Normalised metrics: viral prevalence (% cells) and viral load per 10k UMI.
  - Configurable thresholds via `--se-threshold`, `--detection-threshold`, `--min-counts`, `--min-genes`.
- **PR 12 Combined Host+Virus Reference Builder**
  - `viralscan build-ref` subcommand.
  - `src/viralscan/scripts/build_reference.py`: Ensembl cDNA download + NCBI viral fetch + `kb ref`.
  - Supported host species table (`ENSEMBL_SPECIES`) in `constants.py`.
- **PR 10 NCBI accession → reference**
  - `src/viralscan/scripts/ncbi_fetch.py`: fetch FASTA + GTF from NCBI by accession.
  - `--ncbi-accession` / `-acc` and `--ncbi-email` CLI flags.
  - Per-accession cache under `~/.cache/viralscan/ncbi/` with SHA256 validation.
  - Accession validation regex; exponential backoff on 429/5xx.
- **PR 5 Tests**
  - `tests/test_ncbi_fetch.py` (18 unit tests).
  - `tests/test_cli.py` (8 test classes, 144 tests total passing).
  - `tests/test_createconfig.py`, `tests/test_analysis.py`, `tests/test_build_reference.py`.
- **PR 4 Tooling**
  - `.pre-commit-config.yaml` (ruff, ruff-format, end-of-file-fixer, etc.).
  - GitHub Actions CI (matrix Python 3.9–3.12 × ubuntu/macos) and release workflow.
- **PR 3 Code cleanup**
  - Removed ~620 lines of dead/commented-out code from `detection.py` and `umap.py`.
  - Moved `VIRUS_NAME_MAP` to `constants.py`; shared `load_config()` in `utils.py`.
  - Centralised logging via `configure_logging()` / `setup_script_logging()`.
- **PR 2 Correctness & security**
  - Fixed `--reference`, `--visual`, `--multimapping` boolean flag bugs (§1.1).
  - Replaced `os.system` and `shell=True` subprocess calls with safe `subprocess.run([...])`.
  - Replaced all bare `exit()` with `sys.exit()` via `_die()` helper.
  - Added `_check_required_tools()` preflight for `kb` and `snakemake`.
- **PR 1 Hygiene**
  - README fixes, expanded `.gitignore`, `pyproject.toml` classifiers and optional deps.

### Changed
- `createconfig.py`: boolean config values are now native Python `bool`, not strings.
- `umap.py`: removed `== "True"` string comparisons for config booleans.
- `analysis.py`: replaced `!= "None"` guard with `config.get('gtf')`.
- Snakefile: whitelist comparison no longer uses literal `"None"`.

---

## [2.1.0] - 2025-01-01

*(No structured changelog was kept before v2.2.0. See git log for history.)*

---

[Unreleased]: https://github.com/mdmanurung/ViralScan/compare/v2.5.0...HEAD
[2.5.0]: https://github.com/mdmanurung/ViralScan/compare/v2.4.0...v2.5.0
[2.4.0]: https://github.com/mdmanurung/ViralScan/compare/v2.3.0...v2.4.0
[2.3.0]: https://github.com/mdmanurung/ViralScan/compare/v2.2.0...v2.3.0
[2.2.0]: https://github.com/mdmanurung/ViralScan/compare/v2.1.0...v2.2.0
[2.1.0]: https://github.com/mdmanurung/ViralScan/releases/tag/v2.1.0
