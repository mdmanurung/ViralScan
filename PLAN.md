# ViralScan Improvement Plan — Implementation Tracker

This file tracks progress on the repository improvement plan that was approved
on branch `claude/review-repo-improvements-Sg4Th`. It is the single source of
truth for what's done, in flight, and still pending. **After every
implementation step, update this file** (tick the relevant box, add a one-line
note about what changed, and adjust the "Next" pointer).

The full rationale for each item lives in the original planning document
(`/root/.claude/plans/what-can-be-further-silly-platypus.md`); this file is
the operational checklist.

---

## Status legend

- `[x]` — done and verified
- `[~]` — partially done
- `[ ]` — not started
- `[!]` — blocked / needs decision

## Next up

→ PR 3: Code cleanup (remove commented-out blocks, extract `constants.py`/`utils.py`,
   replace `print` + ANSI with `logging`).

---

## PR 1 — Hygiene quick wins   `[x]`

- [x] README typo `thhe` → `the`
- [x] README "User Guide": de-duplicated index bullets, added 3rd run mode (NCBI accession)
- [x] README: replaced `-reference True` examples with `--reference` flag
- [x] `.gitignore`: expanded with `__pycache__`, `.venv`, `.pytest_cache`, `.ruff_cache`, `.mypy_cache`, `.coverage`, `.idea`, `.vscode`, `.DS_Store`, etc.
- [x] `pyproject.toml`: declared `license = {file = "LICENSE"}`
- [x] `pyproject.toml`: full classifier set (License OSI, Topic Bio-Informatics, Python 3.9–3.12, Dev Status Beta, Audience Science/Research)
- [x] `pyproject.toml`: added `[project.optional-dependencies] dev`
- [x] `pyproject.toml`: added `[tool.pytest.ini_options]` with `network` and `integration` markers
- [x] `pyproject.toml`: added `[tool.ruff]` (line-length 100, target py39)
- [x] `pyproject.toml`: added `pyyaml` and `requests>=2.28` runtime deps
- [ ] Backfill `CHANGELOG.md` (Keep-a-Changelog style)
- [ ] Add `CITATION.cff`

## PR 2 — Correctness & security   `[x]`

- [x] §1.1 `--reference`: `type=bool` → `action='store_true'`
- [x] §1.1 `--visual` / `--multimapping`: `default=True` (string-bug) → `argparse.BooleanOptionalAction`
- [x] §1.2 replaced `os.system('kb ref ...')` with `subprocess.run([...], check=True)`
- [x] §1.2 replaced `subprocess.check_output("wc -l ...", shell=True)` and `cut|sort -u|wc -l` with pure-Python `_count_lines` + `_count_unique_genes`
- [x] §1.3 replaced all 12 bare `exit()` calls with `sys.exit(...)` (via `_die` helper writing to stderr)
- [x] §1.4 fixed validation that printed but didn't stop (now consolidated through `_die`)
- [x] §1.5 inverted reference-flag logic restructured (positive `if args.reference: ... else:`)
- [x] §2.7 added `_check_required_tools()` preflight (`shutil.which` for `kb`, `snakemake`)
- [x] §1.9 `pathlib.Path` / `os.path.join` for output paths
- [ ] §1.6 normalise booleans in `createconfig.py` and remove `== "True"` checks in `umap.py`
- [ ] §1.7 Snakefile: stop comparing whitelist to literal `"None"`
- [ ] §1.8 `analysis.py:45` `!= "None"` → `config.get('gtf')`

## PR 3 — Code cleanup   `[ ]`

- [ ] Remove ~450 lines of commented-out code (`detection.py:275-349`, `umap.py:~526` block)
- [ ] Move duplicated virus-name dict to `src/viralscan/constants.py`; import from `detection.py` and `umap.py`
- [ ] Add `src/viralscan/utils.py` with `load_config(path)` shared helper
- [ ] Replace `print` + ANSI escapes with `logging` (configurable via `--verbose`/`--quiet`)
- [ ] Re-order `umap.py` so `viral_neighbor_enrichment` is defined before its caller
- [ ] `multimap.py:111` `sep=r"\s+"` → `delim_whitespace=True` (silence FutureWarning)
- [ ] Clean up unused `log_done` and duplicate `mkdir` in `createconfig.py`

## PR 4 — Tooling   `[x]`

- [x] `.pre-commit-config.yaml` (ruff, ruff-format, end-of-file-fixer, trailing-whitespace, check-yaml, check-toml, check-merge-conflict, check-added-large-files)
- [x] `.github/workflows/ci.yml` (matrix Python 3.9–3.12 × {ubuntu, macos}; ruff check + ruff format --check + mypy informational + pytest with coverage + CLI smoke test)
- [x] `.github/workflows/release.yml` (tag → build sdist+wheel → PyPI Trusted Publishing, with tag-vs-pyproject version check)
- [x] `pyproject.toml [tool.ruff.lint]`: `select = E4/E7/E9/F`, `E501` ignored, per-file ignores for Snakemake script files
- [x] One-shot `ruff format` over the codebase to establish the baseline
- [x] Fixed 3 py3.9-incompatible nested-quote f-strings in `detection.py`/`multimap.py`/`umap.py` so CI is green from day one

## PR 5 — Tests   `[~]`

- [x] `tests/test_ncbi_fetch.py` (18 unit tests — all passing locally without network)
- [ ] `tests/test_cli.py` (parse `--help`, regression for §1.1 boolean flag bug)
- [ ] `tests/test_errorhandler.py` (each error branch → exit code/message)
- [ ] `tests/test_createconfig.py` (YAML round-trip, booleans)
- [ ] `tests/test_analysis.py` (synthetic GTF fixture)
- [ ] `tests/integration/` smoke test (gated by `pytest -m integration`)
- [ ] Hook coverage reporting + Codecov badge

## PR 6 — Reproducibility   `[ ]`

- [ ] `environment.yml` with pinned conda deps
- [ ] `Dockerfile` (mamba-based, ≤500 MB)
- [ ] Singularity definition for HPC

## PR 7 — Docs   `[ ]`

- [ ] `docs/` skeleton (Sphinx + MyST + RTD theme)
- [ ] Pages: Installation, Quickstart (real worked example), CLI reference, Reference panel curation, Output reference, FAQ
- [ ] Read the Docs hookup
- [ ] Rewrite `getting_started.ipynb` as a real end-to-end notebook (or remove)

## PR 8 — Data unbundling   `[ ]`

- [ ] Move 195 GTFs out of the wheel; ship via Zenodo
- [ ] `viralscan data fetch` subcommand with `~/.cache/viralscan/` + SHA256
- [ ] Drop `[tool.setuptools.package-data]` `data/*.gtf` entry

## PR 9 — Type hints + magic-number config   `[ ]`

- [ ] Type-annotate public functions in `menu.py`, `utils.py`, `constants.py`
- [ ] Lift detection/UMAP magic numbers into config (or `defaults.py`)
- [ ] Enable `mypy --strict` per module incrementally

## PR 10 — NCBI accession → reference (§6 of plan)   `[x]`

- [x] New module `src/viralscan/scripts/ncbi_fetch.py` (efetch + minimal GenBank → GTF)
- [x] `--ncbi-accession` / `-acc` CLI flag, `--ncbi-email` (or `$NCBI_EMAIL`)
- [x] Per-accession cache under `$VIRALSCAN_CACHE` or `~/.cache/viralscan/ncbi/`
- [x] Accession validation regex; rate-limit + exponential backoff retry on 429/5xx
- [x] Mutual exclusion with `--reference` / `-fasta` / `-gtf` enforced in `errorhandler`
- [x] 18 unit tests (validation, location parser, GenBank→GTF, fetch arg validation)
- [ ] Live integration test gated by `@pytest.mark.network`

---

## Verification checklist (run after each PR)

- `python -c "import ast; [ast.parse(open(p).read()) for p in ['src/viralscan/menu.py','src/viralscan/scripts/ncbi_fetch.py']]"` — syntax
- `PYTHONPATH=src python -m pytest tests/ -v` — unit tests
- `PYTHONPATH=src python -m viralscan.menu --help` — CLI parses
- `ruff check .` (once tooling PR is merged)
- `mypy src/viralscan` (once tooling PR is merged)

## How to update this file

Each implementation commit must:
1. Tick the relevant `[ ]` → `[x]` (or `[~]`).
2. Update the "Next up" pointer if the focus has shifted.
3. Note any decisions/blockers under the relevant section.
