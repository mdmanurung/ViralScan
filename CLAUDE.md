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

The full rationale (why each item exists, line-number references, etc.) lives
in the original planning doc at
`/root/.claude/plans/what-can-be-further-silly-platypus.md`. `PLAN.md` is the
operational tracker.

## Working on the codebase

### Branch

All in-flight work happens on `claude/review-repo-improvements-Sg4Th`. Do not
push directly to `main`.

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
- `getting_started.ipynb` has stale errors and largely duplicates the README;
  prefer updating the README until the docs site exists.
- The 195 GTFs in `src/viralscan/data/` are 84 % of the package size. PLAN §3.6
  / PR 8 will move these to Zenodo; do not add more without discussing.
