# ViralScan Publication-Readiness Plan — RECONCILED against HEAD (2026-07-03)

This supersedes `2026-07-03-publication-readiness-99.md`. Every task was verified
against the actual repository state on branch `claude/multimap-memory-and-showcase`.
Items already complete are marked DONE and skipped; items where the original plan
was factually wrong are corrected with the reason.

**Execution decision:** stay on `claude/multimap-memory-and-showcase` (per CLAUDE.md).
Do **not** use the original Task 0 worktree step — it would drop the uncommitted
manuscript edits and untracked scripts. Update `PLAN.md` alongside the work
(CLAUDE.md contract). Stop before release push/tag/publish (owner-gated).

---

## Already DONE at HEAD (skip)

- **Task 1 (package `emptydrops.R`)** — `pyproject.toml:68` already has
  `"viralscan.scripts" = ["emptydrops.R"]`; `tests/test_cellcalling.py:20` already
  has the resource test. Commit `e6b1186` landed it. **Skip entirely.**
- **Task 0 Step 3/4 (execution log)** — `docs/PUBLICATION_READINESS_EXECUTION.md`
  already exists. Append to it, don't recreate.
- **Version cut to 2.5.0** — `__init__.py`, `pyproject.toml`, `CITATION.cff`,
  `docs/conf.py`, `docs/cli_reference.md` header already at 2.5.0.
- **Docs depth-confound honesty** — `docs/cli_reference.md` §"Depth-confound
  controls" and `docs/output_reference.md` §hostresponse already explain the
  depth confound (AUC 0.97 vs 0.87, honest ~0.67). Only the **manuscript** is stale.

## Corrections to the original plan (do NOT apply verbatim)

- **Task 3 Step 6 (`cell_type_enrichment` API) — DROP.** The original rewrite to a
  `SimpleNamespace(config=..., stats=..., per_cell=...)` signature is hallucinated.
  The real signature is `cell_type_enrichment(adata, group_by_virus: dict, cfg: dict)`
  and `docs/api.md:30-88` already documents it correctly. Leave it.
- **Task 7 AUROC pairing — CORRECTED.** The original pairs raw **0.845** (all-cell
  design) with depth-alone **0.967** (a *different*, balanced+depth-filtered design)
  — apples-to-oranges. The tracked summary `hostresponse_summary.tsv` already reports
  `auc_mean = 0.866`, and `depth_confounder.txt:14` states depth-alone 0.967 pairs
  with the host-gene headline **0.866 in the same design**. Fix: update the manuscript
  headline to the tracked **0.866 ± 0.036** (also fixes a text/figure mismatch), then
  the 0.967 pairing is same-design and honest. Depth-controlled: CPM 0.636,
  depth-matched 0.718; honest band ~0.64–0.72.
- **Task 11 Step 6 grep — RESCOPE.** Do not forbid `AUROC 0.845 …` (Task 7 no longer
  keeps it — it's replaced by 0.866) and exempt `docs/showcase_runbook.md` +
  `docs/PUBLICATION_READINESS*.md` version guards from the `2.3.0` sweep.
- **CHANGELOG (Task 3 Step 9)** — a `[2.5.0]` section already exists; there is **no**
  `[2.4.0]` release. Only fix the bottom compare links to match existing sections
  (`[Unreleased]` → v2.5.0…HEAD; add `[2.5.0]` → v2.3.0…v2.5.0). Do not invent a 2.4.0 link.
- **Task 4 Step 4 (CITATION.cff DOI)** — no fake DOI present; no change until Zenodo.

---

## Tasks to execute

### T2 — Installed-package CI + release install-test
`.github/workflows/ci.yml`, `.github/workflows/release.yml`.
- `test` job: add `python -m pip install --no-deps -e .` (skips the snakemake
  `connection_pool` transitive-dep build that motivated `PYTHONPATH`), drop
  `PYTHONPATH: src`, switch smoke test to the installed console script + subcommands.
- `integration` job: add `pip install --no-deps -e .`, drop `PYTHONPATH: src`.
- `release.yml` build job: add a wheel install-test in a venv that asserts
  `viralscan --version/--help` and `emptydrops.R` is packaged, before publish.
- Verify `--no-deps -e .` actually builds locally before committing.

### T3 — Docs consistency (verified subset only)
- `tests/test_docs_consistency.py`: version/default/schema guards, adjusted to
  reality (cli_reference already 2.5.0; `is_called_cell` now documented).
- `docs/installation.md:59-67`: `2.3.0` → `2.5.0` (docker/singularity).
- Multimap default `host-conservative` → `equal` in `docs/quickstart.md:130,173`,
  `docs/cli_reference.md:21-23,61`, `docs/output_reference.md:156`; recommend
  `host-conservative` for host-virus cross-homology. Add `em` to the method list.
- `docs/cli_reference.md` build-ref table: add real flags `--anellovirus/--no-anellovirus`
  (default on), `--reference-panel anellovirus`, `--no-mask`, `--cluster` (verified in menu.py:152-184).
- `docs/api.md` `build_combined_reference`: document `include_anellovirus=True` (real, build_reference.py:360).
- `docs/output_reference.md` per_cell_viral.tsv: add `is_called_cell` (real, detection.py:425/470).
- `src/viralscan/menu.py:743-748`: add `evidence` + `check-whitelist` to the root help block.
- `CHANGELOG.md` bottom compare links (see correction above).

### T4 — Distribution completeness
- `MANIFEST.in`: LICENSE/README/CHANGELOG/CITATION/env/containers + package data;
  scope docs narrowly (do not ship the manuscript/superpowers plans in the sdist).
- `conda-recipe/meta.yaml`: add the "fill sha256 after PyPI" procedure comment.

### T5 — Promote benchmark provenance out of ignored dirs
- Create `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/`
  (README, failure_summary.tsv, copies of commands.jsonl/fastq_audit.tsv/reference_audit.tsv).
- Correct provenance language in the data manifests + `reference_manifest.json`.
- The six untracked `scripts/*reference_strategy*` helpers get committed with T5.

### T6 — Exclude incomplete reference-strategy benchmark from manuscript claims
- Manuscript does not cite it (verified). Add a "Publication Use" exclusion note to
  `analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md`; note the
  exclusion in `docs/PUBLICATION_READINESS.md`.

### T7 — Correct host-response manuscript text + Figure 2 (the honesty gate)
- `docs/manuscript_draft.md`: heading `:85`; rewrite `:87-91` with the 0.866 headline,
  same-design depth-alone 0.967, depth quintile trend, depth-controlled 0.636/0.718,
  honest band ~0.64–0.72; add a depth limitation to Discussion `:101`.
- `scripts/make_manuscript_figures.py`: add a depth-confound caveat line to the Fig 2 footer.
- Regenerate `docs/figures/figure2_benchmark.{png,pdf}`.

### T8B — Narrow manuscript scope (dedicated-tool comparison deferred)
- Cannot run Venus/ViralTrack here. Add the scope-narrowing + explicit
  no-dedicated-comparison limitation to `docs/manuscript_draft.md`; note the deferral
  in `docs/PUBLICATION_READINESS.md`. (Optionally scaffold `analysis/dedicated_tool_comparison/` as future work.)

### T9 — Truth-panel specificity/sensitivity framing
- `analysis/truth_panel/README.md`; ensure README + manuscript state FPR/FNR are not
  yet formally characterized.

### T10 — Citation hygiene (author metadata is owner-gated)
- Resolve the VIRTUS2 reference (verify a DOI or fall back to VIRTUS). **Author list,
  affiliations, lead contact, acknowledgments, contributions, declarations remain the
  owner's to supply** — leave placeholders, flag clearly. Do not fabricate.

### T11 — Verification
- `compileall`, unit tests (conda env python), `python -m build` + `twine check`,
  wheel install smoke, docs build, rescoped stale-claim grep. Append results to the
  execution log.

### PLAN.md
- Add a "Publication readiness (v2.5 release + manuscript honesty)" section tracking T2–T11.

## Owner-gated / NOT done here
- Task 12 release: PR→main, tag `v2.5.0`, PyPI publish, GHCR, bioconda sha256, Zenodo DOI.
- Manuscript authorship block + declarations.
- Dedicated-tool head-to-head benchmark (T8A) and FPR/FNR truth panel execution.
