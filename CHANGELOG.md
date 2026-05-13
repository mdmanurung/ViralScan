# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Changed
- Default multimapper allocation is now `host-conservative`, making combined
  host+virus references the recommended host-aware workflow while preserving
  legacy equal splitting via `--multimap-method equal`.

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

[Unreleased]: https://github.com/mdmanurung/ViralScan/compare/v2.3.0...HEAD
[2.3.0]: https://github.com/mdmanurung/ViralScan/compare/v2.2.0...v2.3.0
[2.2.0]: https://github.com/mdmanurung/ViralScan/compare/v2.1.0...v2.2.0
[2.1.0]: https://github.com/mdmanurung/ViralScan/releases/tag/v2.1.0
