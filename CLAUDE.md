# CLAUDE.md

Guidance for Claude Code (and other Claude agents) when working in this
repository.

## What this repo is

**ViralScan** is a Snakemake-driven Python bioinformatics CLI that quantifies
viral load from paired-end FASTQ samples using `kb-python` (kallisto +
bustools). The package is exposed as the `viralscan` command (entry point in
`src/viralscan/menu.py`).

Layout:

```
src/viralscan/
  menu.py                  # CLI entry, argparse, validation, top-level orchestration
  Snakefile                # 6 rules: create_config → kb_count → analysis → multimap → detection → umap
  scripts/
    createconfig.py        # writes the per-sample config.yaml consumed by the rules
    analysis.py            # parses GTFs and lists viral accessions
    multimap.py            # multimapping correction
    detection.py           # viral detection + visualizations
    umap.py                # UMAP plot
    ncbi_fetch.py          # download FASTA + GTF from NCBI by accession (no Biopython)
  data/*.gtf               # 195 bundled viral reference annotations
tests/                     # pytest suite (run with PYTHONPATH=src)
PLAN.md                    # implementation tracker — keep this current!
```

## The PLAN.md contract  (IMPORTANT — non-negotiable)

`PLAN.md` is the authoritative checklist for the in-flight repo improvement
plan. **ALWAYS update `PLAN.md` whenever you finish an implementation step,
in the same commit as the implementation itself.** This is mandatory, not
optional — treat it as part of the definition of "done" for any task in this
repo. Do not consider an implementation complete (do not commit, do not push,
do not open the PR) until `PLAN.md` has been updated.

Each implementation commit must:

1. Flip the relevant checkbox: `[ ]` → `[x]` (or `[~]` for partial, `[!]` for
   blocked).
2. Update the "Next up" pointer at the top if the focus has changed.
3. Add a short note under the relevant PR section if you made a non-obvious
   decision or hit a blocker.

If a piece of work isn't in `PLAN.md` yet, add a row before starting it. Do
not silently skip this — the user relies on `PLAN.md` to see progress
between sessions.

The locked v3 rationale lives in
`docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md`.
`PLAN.md` is the operational tracker.

## Working on the codebase

### Branch

All v3 in-flight work happens on `codex/viralscan-v3`. Do not push directly to
`main`.

### Running the tests

The local environment cannot reliably `pip install -e .` because a transitive
snakemake dep (`connection_pool`) fails to build with the system setuptools.
Use `PYTHONPATH` instead:

```
PYTHONPATH=src python -m pytest tests/ -v
PYTHONPATH=src python -c "from viralscan import menu; menu.create_help()"  # smoke test
```

Network-hitting tests are gated by `@pytest.mark.network` (see
`pyproject.toml`); run them explicitly with `pytest -m network`.

### External tools

`viralscan` invokes the `kb` (kb-python) and `snakemake` binaries at runtime.
`menu.py:_check_required_tools()` does a `shutil.which` preflight. When
adding new shell-outs, prefer `subprocess.run([...], check=True)` (list form,
no `shell=True`) — this is enforced by §1.2 of the plan.

### Style

- No bare `exit()` — use `sys.exit(<code>)` (or `_die()` in `menu.py`).
- No `shell=True` with user-controlled paths.
- Booleans on the CLI use `action='store_true'` or
  `action=argparse.BooleanOptionalAction`, never `type=bool`.
- Path operations: `pathlib.Path` or `os.path.join`, never manual string
  concatenation / slash trimming.
- Tests live in `tests/` and run from the repo root.

### Commit hygiene

Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `chore:`) are
preferred but not yet enforced. When you complete a PLAN.md row, mention it
in the commit body, e.g. `feat(ncbi): add accession-based reference fetch
(closes PLAN PR 10)`.

## Common pitfalls

- The legacy `createconfig.py` writes some YAML values as strings (e.g.
  `"True"`/`"False"`). Downstream `umap.py` has both `if config["umap"]:` and
  `if config["umap"] == "True":` checks. PLAN §1.6 fixes this — until then,
  be careful when adding new boolean config keys.
- The stale `getting_started.ipynb` was removed (2026-07-02) and
  `basic_usage.ipynb` was superseded (2026-07-20) by an 8-vignette suite under
  `docs/vignettes/` (index: `docs/vignettes/README.md`; design:
  `VIGNETTES_PLAN.md`). Six execute in CI on synthetic/committed data; the two
  `[skip-ci]` ones (`quickstart_fastq_to_viral_load`, `building_a_reference`)
  build an index / run `kb count`. Keep CI notebooks runnable — pass `RunConfig`
  (not a dict) to `cell_type_enrichment()`.
- The 195 GTFs in `src/viralscan/data/` are 84 % of the package size. PLAN §3.6
  / PR 8 will move these to Zenodo; do not add more without discussing.


## Installed Convention Packs

- **skill-bridge** — See `.living/conventions/skill-bridge/analysis-conventions.md`

- **bioinformatics** — See `.living/conventions/bioinformatics/analysis-conventions.md`

- **idea-generator** — See `.living/conventions/idea-generator/analysis-conventions.md`

- **report-generator** — See `.living/conventions/report-generator/analysis-conventions.md`

- **robust-analysis** — See `.living/conventions/robust-analysis/analysis-conventions.md`

- **aifi-scrna-pipeline** — See `.living/conventions/aifi-scrna-pipeline/SKILL.md`

## Mycelium living-repo layer

This repo is now a mycelium "living repository." A `.living/` memory layer
records decisions, learnings, and findings across sessions; the convention
packs above guide analysis, reporting, and review; SessionStart/PostToolUse/Stop
hooks (in `.claude/settings.local.json`) enforce a post-action logging protocol.

- **`.living/decisions.md` / `.living/learnings.md` / `.living/findings/`** —
  append here after significant work (the hooks will remind you). These capture
  *why* and *what was learned*; they complement, they do not replace, `PLAN.md`.
- **`PLAN.md` remains the authoritative in-flight tracker** (see "The PLAN.md
  contract" above). It is non-negotiable and takes precedence for
  implementation-step tracking. Mycelium's `.living/` layer is additive context,
  not a second source of truth for the plan.
- **`.living/INDEX.md`** — knowledge map; skim before decisions in a known area.
- **Skills**: `/mycelium:ingest` (register data), `/mycelium:analyze` (run an
  analysis under conventions), `/mycelium:report`, `/mycelium:review` (6-agent
  scientific review), `/mycelium:ideas`, `/mycelium:core` (crystallize / todo-idea).
- **Structure note**: the minimal scaffold was chosen — the standard mycelium
  `algorithms/`, `analysis/`, `reference_material/` top-level dirs were pruned
  because ViralScan already has `src/`, `scripts/`+`results/`, and `references/`.
  Only `data/` (ingest target), `todo/`, `skillpacks/`, and `.living/` were kept,
  so `validate_structure.py`'s top-level-dir check is intentionally relaxed.
