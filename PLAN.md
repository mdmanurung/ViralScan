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

→ **PR 5** — Test backfill (`test_cli.py`, `test_createconfig.py`, `test_analysis.py`).
  After PR 5: PR 6 Reproducibility, PR 7 Docs.

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
- [x] §1.6 normalise booleans in `createconfig.py` and remove `== "True"` checks in `umap.py`
- [x] §1.7 Snakefile: stop comparing whitelist to literal `"None"`
- [x] §1.8 `analysis.py:45` `!= "None"` → `config.get('gtf')`

## PR 3 — Code cleanup   `[~]`

- [x] Remove ~620 lines of commented-out code from `detection.py` (super_expressor footer) and `umap.py` (the duplicate `umap()` and the standalone trailing block)
- [x] Move duplicated virus-name dict to `src/viralscan/constants.py`; `detection.py` and `umap.py` now import `VIRUS_NAME_MAP` from there
- [x] Add `src/viralscan/utils.py` with `load_config(path)`; adopted by `analysis.py`, `multimap.py`, `detection.py`, `umap.py`
- [ ] Replace `print` + ANSI escapes with `logging` (configurable via `--verbose`/`--quiet`) — deferred; needs CLI flag plumbing
- [x] Re-order `umap.py` so `viral_neighbor_enrichment` is defined before its caller (now sits above `umap()`)
- [~] `multimap.py` `sep=r"\s+"`: file already uses the raw-string form; the `delim_whitespace=True` fix in the original plan is itself deprecated in pandas 2.2+, so leaving as-is. Marking obsolete.
- [x] Removed unused `file_to_keep` + the duplicate `mkdir` in `createconfig.py`; also dropped the unused `snakefile_dir` reads in `detection.py` and `multimap.py` and the unused `glob`/`yaml` imports in `analysis.py`

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

## PR 11 — Interpretation & Reporting   `[x]`

_Goal: every run produces structured, normalized, publication-ready results._

- [x] **A1 Structured tabular output** — `results/viral_summary.tsv` (per-virus: name, total UMI,
  infected cells, total cells, % infected, UMI-per-10k, clustering p-value) +
  `results/per_cell_viral.tsv` (per-barcode: virus, viral UMI, total UMI, viral fraction)
- [x] **A2 HTML report** — single self-contained `report.html` via Jinja2 template; includes
  run metadata, QC table, per-virus table, all plots embedded as base64 PNG, interpretation
  guidance, and the `viral_neighbor_enrichment` p-value surfaced prominently.
  Add `jinja2` to `pyproject.toml` runtime deps.
- [x] **A3 Normalized metrics** — viral prevalence (% cells ≥1 viral UMI) and viral load per 10k
  UMI computed in `detection.py` and exported to all outputs.
- [x] **A4 Configurable thresholds** — `--se-threshold`, `--detection-threshold`, `--min-counts`,
  `--min-genes` added to CLI and wired through createconfig.py and config YAML.
- [ ] **A5 Cell-type-aware enrichment** *(stretch)* — `--cell-types PATH` (barcode→cell_type CSV);
  per-cell-type viral prevalence table + Fisher's exact test added to report.

## PR 12 — Combined Host+Virus Reference Builder   `[x]`

_Goal: `viralscan build-ref --host human --virus-accessions NC_xxx,... --output ref/`
constructs a kb-ready index combining an Ensembl host transcriptome with NCBI viral sequences._

- [x] **B1 `src/viralscan/scripts/build_reference.py`** — standalone module:
  - `fetch_host_cdna(species, out_dir, cache_dir)` — downloads Ensembl cDNA FASTA + GTF via
    HTTPS FTP mirror; cached + SHA256-verified; supported species table in `constants.py`.
  - `_genome_as_transcript_gtf(fasta_text, accession)` — ports `Viral_GTF_maker.py` logic
    (whole-genome-as-gene) as a pure function suitable for `kb ref`.
  - `build_combined_reference(host_species, virus_accessions, out_dir, ...)` — orchestrates
    Ensembl download → NCBI fetch (`ncbi_fetch.fetch_reference`) → concatenation → `kb ref`.
  - `--no-kb-ref` flag: stop after writing `combined.fasta` + `combined.gtf`.
- [x] **B2 `viralscan build-ref` subcommand** — refactor `menu.py` to use `add_subparsers()`;
  existing behaviour under `viralscan run` (or keep positional default for back-compat);
  new subcommand `build-ref` wires to `build_reference.build_combined_reference`.
- [x] **B3 `ENSEMBL_SPECIES` lookup table** added to `constants.py`.
- [x] **B4 Tests** — `tests/test_build_reference.py`: species lookup, GTF adapter, cache logic,
  mock download (22 new tests; 40 total passing without network);
  `@pytest.mark.network` integration test for SARS-CoV-2.

## PR 14 — Bug fixes & performance (audit 2026-05-08)   `[x]`

_All issues discovered in code audit of 2026-05-08.  Each item is self-contained:
file, exact location, what is wrong, and the required fix._

---

### C1 — Double-counting unique reads in multimapping mode *(critical — statistical)*   `[x]`

**Files:** `scripts/multimap.py` → `build_multimap_matrix()` / `final_results()`;
`scripts/detection.py` → `preprocessing()`

**What is wrong:**
`build_multimap_matrix()` reads the full BUS file — which contains **all** reads
(unique + multimapping).  It redistributes every read into `counts_corrected`.
`final_results()` then also stores the standard `kb count` matrix as
`counts_original`, which already contains the uniquely-mapping reads.
`detection.py` adds the two layers:

```python
# detection.py preprocessing() — WRONG
adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
```

Any read that maps uniquely to one gene appears in **both** layers, so it is
counted twice.  For a typical sample where >90 % of reads are uniquely mapping
this approximately doubles every UMI count, inflating viral load ~2×.

**Fix (choose one — confirm with domain expert before implementing):**

Option A *(replace)*: `counts_corrected` replaces `counts_original`; the
addition in `detection.py` becomes:
```python
adata.X = adata.layers["counts_corrected"]   # correction replaces, not augments
```

Option B *(additive model)*: keep only the **extra** multi-mapped share in
`counts_corrected` by subtracting unique-read contribution from the BUS total.
Unique-mapping ECs have `len(genes_in_ec) == 1`; for those, set `share = 0` so
they are not redistributed into `counts_corrected`, and the addition in
`detection.py` then correctly gives `unique + extra_multimapper_fraction`.

Change `build_multimap_matrix()` in `multimap.py`:
```python
# Only redistribute reads that are genuinely multi-mapping
share = count / len(genes_in_ec) if len(genes_in_ec) > 1 else 0.0
```
And update `umap.py` `main()` analogously — it performs the same addition.

**Also update `umap.py` `main()`** which has the identical pattern:
```python
# umap.py main() — same double-count bug
if "counts_corrected" in adata.layers and "counts_original" in adata.layers:
    adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
```

---

### C2 — Detection threshold off-by-one (`>` should be `>=`)   `[x]`

**File:** `scripts/detection.py` → `preprocessing()`

**What is wrong:**
```python
threshold = config.get("detection_threshold", 1)
if total_count > threshold:          # strict greater-than
```

The CLI help says *"Minimum total viral UMI required to call a virus detected.
Default: 1"*, which implies ≥1.  With `> 1` a gene with exactly 1 UMI is
silently dropped.  Default threshold is 1, so the first detected UMI is always
missed.

**Fix:**
```python
if total_count >= threshold:
```

---

### C3 — O(n_genes) list scan inside EC-parsing loop *(critical — performance)*   `[x]`

**File:** `scripts/multimap.py` → `read_ec()`

**What is wrong:**
```python
gene_indices = [gene_ids.index(gid) for gid in gene_ids_ec if gid in gene_ids]
#                         ↑ list.index() scans up to 60 k entries per call
```

`gene_ids` is a plain Python `list`.  `list.index()` is O(n).  This is called
for every line of the EC file; a human+virus reference has 60 k+ genes and a
typical BUS EC file has hundreds of thousands of lines, making this loop
O(n_genes × n_ec) — potentially tens of billions of comparisons.

**Fix:** build a reverse-lookup dict once before the loop:
```python
# Add before the `with open(ec_file)` block:
gene_id_to_idx: dict[str, int] = {gid: i for i, gid in enumerate(gene_ids)}
```
Then inside the loop:
```python
gene_indices = [gene_id_to_idx[gid] for gid in gene_ids_ec if gid in gene_id_to_idx]
```
Drop the `if gid in gene_ids` membership test (also O(n) on a list) — the dict
lookup handles the missing-key case.

---

### C4 — `iterrows()` on multi-million-row BUS DataFrame *(critical — performance)*   `[x]`

**File:** `scripts/multimap.py` → `build_multimap_matrix()`

**What is wrong:**
```python
for idx, row in bus_df.iterrows():          # iterrows is 10–100× slower than vectorized ops
    bc, ec, count = row["barcode"], row["ec"], row["count"]
```

A 10x PBMC BUS file easily has 10–100 million rows.  `iterrows()` is the
slowest pandas iteration method (returns a full Series per row with dtype
inference overhead).  This can take hours vs. seconds for vectorized
alternatives.

**Fix:** use `itertuples()` for a simple drop-in speedup (5–10×):
```python
for row in bus_df.itertuples(index=False):
    bc, ec, count = row.barcode, row.ec, row.count
    if pd.isna(ec):
        ...
```
Or (preferred, larger speedup): filter the DataFrame to known barcodes and
known ECs first, then use `numpy`/scipy COO matrix construction directly from
the filtered arrays — avoiding Python-level row iteration entirely.

---

### C5 — Snakefile `kb_count`: stderr discarded; error check never fires   `[x]`

**File:** `src/viralscan/Snakefile` → `rule kb_count` (shell block)

**What is wrong:**
```bash
output_log=$(kb count ... 2> /dev/null)   # stderr is thrown away

if echo "$output_log" | grep -q "no reads pseudoaligned"; then   # checks stdout (empty)
    exit 1
fi
```

`kb count` writes **all** progress and error messages — including
"no reads pseudoaligned" — to **stderr**.  By discarding stderr the `grep`
check is checking an empty string and can never trigger, even if alignment
completely fails.

**Fix:**
```bash
output_log=$(kb count ... 2>&1)   # merge stderr into stdout capture
```

---

### C6 — Snakefile `kb_count`: `mkdir` without `-p` fails on re-run   `[x]`

**File:** `src/viralscan/Snakefile` → `rule kb_count` (shell block)

**What is wrong:**
```bash
mkdir {config[output]}kb-python/    # exits non-zero if directory already exists
```

When `--overwrite` is used or the rule is re-run after a partial failure, the
directory exists and `mkdir` fails, aborting the rule.

**Fix:**
```bash
mkdir -p {config[output]}kb-python/
```

---

### C7 — Redundant `f.close()` after `with` block — `analysis.py`   `[x]`

**File:** `scripts/analysis.py` → `obtain_gtf()`

**What is wrong:**
```python
with open(f"{config['output']}log/analysis.txt", "w") as f:
    for v in viral_accessions:
        f.write(v + "\n")
f.close()    # ← f is already closed by the with-statement; this is a no-op
```

Not a correctness bug today, but signals a misunderstanding that could cause
real issues if the pattern is copied.

**Fix:** remove the `f.close()` line.

---

### C8 — File handle leak in `detection.py` when no viral genes found   `[x]`

**File:** `scripts/detection.py` → `main()`

**What is wrong:**
```python
found_genes_file = open(f"{config['output']}log/found_genes.txt", "w")
if len(found_genes_sorted) > 0:
    for g in found_genes_sorted:
        found_genes_file.write(write_to_file)
    found_genes_file.close()    # only closed inside the if-block
# ← file handle leaks if found_genes_sorted is empty
```

**Fix:** use a `with` statement:
```python
with open(f"{config['output']}log/found_genes.txt", "w") as found_genes_file:
    for g in found_genes_sorted:
        ...
        found_genes_file.write(write_to_file)
```

---

### C9 — Inconsistent `var_names` source for viral gene lookup in `umap.py`   `[x]`

**File:** `scripts/umap.py` → `umap()`

**What is wrong:**
Two separate code blocks compute `viral_ids` from different `var_names` sources:
```python
# Block 1 — has_viral check
var_names = adata.raw.var_names if getattr(adata, "raw", None) else adata.var_names
viral_ids = [g for g in found_genes if g in var_names]
has_viral = len(viral_ids) > 0

# ... many lines later ...

# Block 2 — actual viral count extraction
if getattr(adata, "raw", None) and getattr(adata.raw, "X", None) is not None:
    X_counts = adata.raw.X
    var_names = adata.raw.var_names   # ← different condition than block 1
elif "counts" in adata.layers:
    X_counts = adata.layers["counts"]
    var_names = adata.var_names
else:
    X_counts = adata.X
    var_names = adata.var_names

viral_ids = [g for g in found_genes if g in var_names]   # recomputed — could differ
```

If `adata.raw` is set but `adata.raw.X` is `None`, block 1 uses
`adata.raw.var_names` but block 2 falls through to `adata.var_names` — `has_viral`
and the actual extraction use different gene namespaces.

**Fix:** extract `X_counts` and `var_names` once before the `has_viral` check,
reuse the same variables, and compute `viral_ids` only once.

---

### C10 — No HVG selection before PCA in `umap.py` *(statistical / performance)*   `[x]`

**File:** `scripts/umap.py` → `umap()` (both the no-viral and viral branches)

**What is wrong:**
```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# ← highly_variable_genes() call missing
sc.pp.pca(adata, svd_solver="arpack")   # runs on ALL 30k–60k features
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

Running PCA on a host+virus reference matrix with 30 k–60 k features:
1. Makes PCA slow (dominated by O(n_genes²) computation).
2. Statistically poor: the viral genes (tiny fraction of variance) are drowned
   out by thousands of non-informative housekeeping genes.
3. The resulting UMAP neighbourhood graph does not reflect viral expression
   structure.

**Fix (standard scanpy pipeline):**
```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# Force-include viral genes so they survive HVG filter:
adata.var["highly_variable"] |= adata.var_names.isin(list(found_genes))
sc.pp.pca(adata, use_highly_variable=True, svd_solver="arpack")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

---

## PR 13 — Housekeeping finish   `[x]`

- [x] §1.6 Normalize booleans in `createconfig.py` (write Python bool, not string `"True"`);
  remove `== "True"` fallback checks in `umap.py`.
- [x] §1.7 Snakefile: replace `[ "{config[whitelist]}" = "None" ]` with a proper empty-string
  or null check.
- [x] §1.8 `analysis.py` line 45: `!= "None"` → `config.get('gtf')` truthiness check.
- [x] PR 3 logging: replaced all ANSI print calls in Python scripts with logging;
  `--verbose` / `--quiet` CLI flags added to menu.py.
- [ ] PR 5 test backfill: `tests/test_cli.py`, `tests/test_createconfig.py`,
  `tests/test_analysis.py`.

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
