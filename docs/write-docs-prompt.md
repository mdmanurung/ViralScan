# ViralScan Documentation Sprint — One-Off Agent Prompt

> **Usage:** paste this entire file as the first message to a fresh agent session.
> It is self-contained; no additional context is needed.

---

## Mission

You are a technical-writing agent working on **ViralScan v2.3.0** — a
Snakemake-driven Python CLI that quantifies viral load from paired-end FASTQ
samples using `kb-python` (kallisto + bustools).  Your job is to:

1. **Rewrite `docs/api.md`** — full, accurate API reference generated from the
   source docstrings and type hints.
2. **Overhaul `README.md`** — tighten prose, fix inaccuracies, add missing
   sections, make it the definitive landing page.
3. **Write two Jupyter vignettes** in `docs/vignettes/` — executable notebooks
   that demonstrate real end-to-end workflows.
4. **Update `PLAN.md`** checkbox for the docs task when done (see CLAUDE.md
   contract).

Work on the branch `claude/review-repo-improvements-Sg4Th`.  All edits must
pass `PYTHONPATH=src python -m pytest tests/ -q` without new failures.

---

## Codebase snapshot (as of 2026-05-08)

### Package layout

```
src/viralscan/
  __init__.py
  __main__.py
  menu.py          # CLI entry: argparse, validation, orchestration
  utils.py         # load_config(), configure_logging(), setup_script_logging()
  constants.py     # VIRUS_NAME_MAP: dict[str, str]  (195 viruses)
  defaults.py      # DEFAULTS dict consumed by menu.py
  enrichment.py    # cell_type_enrichment(), write_cell_type_enrichment(), _bh_adjust()
  Snakefile        # 6 rules: create_config → kb_count → analysis → multimap → detection → umap
  scripts/
    createconfig.py   # writes per-sample config.yaml
    analysis.py       # parses GTFs, lists viral accessions
    multimap.py       # multimapping correction
    detection.py      # viral detection + visualisations
    umap.py           # UMAP plot
    ncbi_fetch.py     # NCBI E-utilities FASTA+GTF downloader (no Biopython)
  data/*.gtf        # 195 bundled viral reference annotations
  templates/        # Jinja2 HTML report template
```

### Key public API (verbatim signatures — do not invent)

**`viralscan.utils`**
```python
def configure_logging(verbose: bool = False, quiet: bool = False) -> None: ...
def setup_script_logging() -> logging.Logger: ...
def split_comma_paths(value: str | None) -> list[str]: ...
def load_config(path: Union[str, Path]) -> dict[str, Any]: ...
```

**`viralscan.constants`**
```python
VIRUS_NAME_MAP: dict[str, str]   # prefix → human-readable virus name
```

**`viralscan.enrichment`**
```python
def cell_type_enrichment(
    adata: anndata.AnnData,
    virus_obs_col: str,
    cell_type_col: str,
) -> pd.DataFrame: ...

def write_cell_type_enrichment(
    df: pd.DataFrame,
    output_dir: str | Path,
) -> Path: ...

def _bh_adjust(pvals: list[float] | npt.NDArray[np.float64]) -> npt.NDArray[np.float64]: ...
```

**`viralscan.scripts.ncbi_fetch`**
```python
class NCBIFetchError(RuntimeError): ...

def fetch_reference(
    accessions: Iterable[str],
    output_dir: str | Path,
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: Path = DEFAULT_CACHE_DIR,
) -> tuple[Path, Path]: ...   # returns (fasta_path, gtf_path)
```

**`viralscan.menu`** (CLI entry — document subcommands, not internals)
- `viralscan [OPTIONS]` — quantification mode
- `viralscan build-ref [OPTIONS]` — reference builder

### CLI flags (complete, as of v2.3.0)

#### `viralscan` (quantification)

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | required | Output directory |
| `--sample1 PATH` | `-s1` | required | Forward FASTQ (comma-sep for batches) |
| `--sample2 PATH` | `-s2` | required | Reverse FASTQ (comma-sep for batches) |
| `--index PATH` | `-i` | — | Pre-built kallisto index |
| `--transcripts PATH` | `-t` | — | t2g.txt transcript-to-gene map |
| `--reference` | `-ref` | off | Build index from `--fasta` + `--gtf` |
| `--fasta PATH` | `-fasta` | — | FASTA file(s), comma-sep |
| `--gtf PATH` | `-gtf` | GTF file(s), comma-sep | |
| `--ncbi-accession ACC` | `-acc` | — | NCBI RefSeq accession(s), comma-sep |
| `--ncbi-email EMAIL` | | `$NCBI_EMAIL` | NCBI contact email |
| `--technology STRING` | `-x` | `10xv3` | Single-cell technology |
| `--whitelist PATH` | `-w` | bundled | Barcode whitelist |
| `--cores N` | `-c` | `6` | CPU cores |
| `--multimapping` / `--no-multimapping` | `-mm` | on | Multimapping correction |
| `--multimap-method METHOD` | | `host-conservative` | Multimapper allocation |
| `--umap` | `-umap` | off | Generate UMAP |
| `--visual` / `--no-visual` | `-v` | on | Generate visualisations |
| `--detection-threshold N` | | `1` | Min viral UMI to call virus detected |
| `--se-threshold N` | | `10` | UMI count to call a cell a super-expressor |
| `--min-counts N` | | `1000` | Min total UMI per cell (UMAP QC) |
| `--min-genes N` | | `200` | Min detected genes per cell (UMAP QC) |
| `--cell-types PATH` | | — | CSV (barcode,cell_type) for per-type enrichment |
| `--verbose` | | off | Enable DEBUG logging |
| `--quiet` | | off | Suppress INFO |

#### `viralscan build-ref`

| Flag | Default | Description |
|------|---------|-------------|
| `--host SPECIES` | — | Host species (e.g. `human`, `mouse`); see `--list-species` |
| `--virus-accessions ACC [ACC …]` | — | NCBI nucleotide accession(s) |
| `--output PATH` | `viralscan_ref` | Output directory |
| `--ncbi-email EMAIL` | `$NCBI_EMAIL` | NCBI contact email |
| `--list-species` | — | Print all supported host species and exit |

### Output layout

```
output/
├── config.yaml
├── kb-python/
│   ├── output.bus
│   ├── run_info.json
│   └── ...
├── results/
│   ├── viral_summary.tsv       # per-virus: total_umi, infected_cells, pct_infected, umi_per_10k
│   ├── per_cell_viral.tsv      # per-barcode×virus: viral_umi, total_umi, viral_fraction
│   ├── report.html             # self-contained Jinja2 HTML report
│   ├── adata_original.h5ad
│   ├── adata_multimap.h5ad
│   └── plots/
│       └── histogram_<virus>.png
└── logs/
    └── snakemake.log
```

---

## Task 1 — Rewrite `docs/api.md`

Replace the current stub (3 lines of `automodule` directives that produce
near-empty output) with a **hand-written, fully accurate** API reference in
MyST Markdown.

**Requirements:**
- One `## Module` section per public module: `viralscan.utils`,
  `viralscan.constants`, `viralscan.enrichment`,
  `viralscan.scripts.ncbi_fetch`.
- For every public function / class / exception: a `### symbol` heading,
  its full signature in a fenced `python` block, a one-paragraph description,
  and a `Parameters` + `Returns` table (Markdown table, not RST field lists).
- Mark `_bh_adjust` as "internal helper — not part of the public API" but
  document it briefly for contributors.
- Do NOT document Snakemake scripts (`detection.py`, `multimap.py`, etc.)
  directly — they run inside Snakemake and have no importable public API.
  Add a brief note explaining this.
- Keep `{eval-rst}` blocks only where truly needed; prefer plain MyST.

**File to edit:** `docs/api.md`

---

## Task 2 — Overhaul `README.md`

The current README is accurate but thin on context and missing several
important sections added since v2.0.

**Required changes (implement all):**

1. **Badges row** — keep existing CI / codecov / PyPI / License badges; add
   a "Python ≥ 3.9" badge:
   `[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue)](https://pypi.org/project/ViralScan/)`

2. **Introduction** — expand to 3–4 sentences covering:
   - what ViralScan detects (viral load at single-cell resolution)
   - the underlying technology (kallisto | bustools via kb-python)
   - the three operating modes (pre-built index / FASTA+GTF / NCBI accession)
   - LUMC context sentence

3. **Add a "Features" section** (bullet list, after Introduction, before
   Installation):
   - 195 pre-bundled viral reference annotations (GTF)
   - Three flexible reference modes
   - Multimapping correction
   - Per-cell and per-virus summary tables + HTML report
   - Cell-type enrichment analysis (`--cell-types`)
   - NCBI auto-fetch + local cache (`--ncbi-accession`)
   - `viralscan build-ref` combined host+virus reference builder
   - Docker / Singularity containers provided
   - Reproducible via Snakemake

4. **Installation section** — replace the current thin `pip install` block
   with three subsections mirroring `docs/installation.md`:
   - Conda (recommended) — show `conda env create -f environment.yml`
   - pip — show the pip command + note about needing conda for kb/snakemake
   - Container (Docker one-liner)

5. **User Guide** — restructure the current "3 ways to run" prose as three
   clear **Mode** subsections with copy-paste shell blocks.  Add a fourth
   subsection for `viralscan build-ref`.  Remove the "Please Note" prose
   paragraph; fold the information into the relevant mode block.

6. **Output section** — add a brief "## Output" section that lists
   `viral_summary.tsv`, `per_cell_viral.tsv`, `report.html`, and
   `adata_multimap.h5ad` with one-line descriptions.  Link to
   `docs/output_reference.md` for the full schema.

7. **Cite section** — add a "## Citation" section pointing to `CITATION.cff`
   and the LUMC affiliation.

8. **License** — keep as-is.

**File to edit:** `README.md`

---

## Task 3 — Write two Jupyter vignettes

Create the directory `docs/vignettes/` and write **two** notebooks there.
They must be valid `.ipynb` JSON (nbformat 4) and should run cleanly if
`ViralScan` is installed with `kb` and `snakemake` available.  Use
`# [skip-ci]` in the first cell to prevent accidental CI execution.

### Vignette 1 — `basic_usage.ipynb`

**Audience:** first-time user.  **Goal:** run ViralScan end-to-end on a
small public sample and inspect the output.

Cell outline:
1. Markdown intro (what we'll do; prerequisites; public test sample
   `SRR20710651`)
2. `!viralscan --help` — show help
3. Download test FASTQ pair from SRA (use `prefetch` + `fasterq-dump` or
   direct wget of the pre-trimmed subset if available)
4. Build index from NCBI accession (`NC_002021.3` — Influenza A) —
   demonstrate `--ncbi-accession` mode
5. Run ViralScan with `--umap --visual`
6. Load `results/viral_summary.tsv` with `pandas` and display
7. Load `results/adata_multimap.h5ad` with `scanpy` and show `.obs`
8. Display the generated `report.html` path and a saved plot
9. Summary markdown cell

### Vignette 2 — `cell_type_enrichment.ipynb`

**Audience:** analyst who already ran ViralScan.  **Goal:** demonstrate the
`--cell-types` workflow and `viralscan.enrichment` Python API.

Cell outline:
1. Markdown intro
2. Load `results/adata_multimap.h5ad` (use the output from Vignette 1 or a
   synthetic AnnData if the file doesn't exist)
3. Assign mock cell types from a CSV — show how to prepare the
   `barcode,cell_type` CSV
4. Import `from viralscan.enrichment import cell_type_enrichment,
   write_cell_type_enrichment`
5. Run `cell_type_enrichment()` and display the result DataFrame
6. Run `write_cell_type_enrichment()` and confirm the TSV was written
7. Bar chart of enrichment p-values with `matplotlib`
8. Summary / interpretation markdown cell

**Format rules for both notebooks:**
- `metadata.kernelspec` → `python3`
- No hard-coded absolute paths; use `pathlib.Path` and relative paths
- Each code cell ≤ 30 lines
- Markdown cells use ATX headings (`##`)
- Do not import `viralscan.menu` directly (it triggers pyfiglet + argparse
  at import time)

---

## Task 4 — Update PLAN.md

After all three writing tasks are done, open `PLAN.md` and:
- Find the open task (or add a new row) for "PR 7 Docs — API reference,
  README overhaul, vignettes"
- Flip its checkbox to `[x]`
- Add a one-line note under the relevant PR section:
  `Rewrote api.md (hand-written), overhauled README (7 sections), wrote 2 vignettes (basic_usage, cell_type_enrichment).`

---

## Constraints and style rules

- **No hallucination:** every function signature, parameter name, default
  value, and output column name must be taken verbatim from the source code
  snapshot above.  Do not invent flags or options.
- **No shell=True** in any code cells that call subprocesses.
- **Booleans:** never `type=bool`; use `argparse.BooleanOptionalAction` or
  `action='store_true'` — follow existing code style.
- **Paths:** `pathlib.Path`, never manual string concatenation.
- **README length:** aim for ~150–200 lines — concise but complete.
- **API reference length:** aim for ~200–250 lines — all public symbols
  documented, no padding.
- Commit message style: `docs: overhaul README, api.md, add vignettes`

---

## Verification checklist (run before marking done)

```bash
# Tests still pass
PYTHONPATH=src python -m pytest tests/ -q

# Notebooks are valid JSON
python -c "import json, pathlib; [json.loads(p.read_text()) for p in pathlib.Path('docs/vignettes').glob('*.ipynb')]"

# README renders (check links exist)
python -c "
import re, pathlib
readme = pathlib.Path('README.md').read_text()
links = re.findall(r'\[.*?\]\((.*?)\)', readme)
for l in links:
    if not l.startswith('http') and not l.startswith('#'):
        assert pathlib.Path(l).exists(), f'Broken link: {l}'
print('All local README links resolve.')
"

# PLAN.md checkbox flipped
grep -E '\[x\].*vignettes|api\.md|README overhaul' PLAN.md
```

---

*End of prompt. Begin with Task 1 (`docs/api.md`), then Task 2
(`README.md`), then Task 3 (vignettes), then Task 4 (PLAN.md).*
