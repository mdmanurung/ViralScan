# Handoff (RECONSTRUCTED) - optional scAutoQC in ViralScan, Codex thread 2026-10-19

**Provenance:** rebuilt on 2026-10-19 from Codex rollout `~/.codex/sessions/2026/10/19/rollout-2026-10-19T00-59-20-01a1513d-b41e-7ac1-bcdd-3c2f3796e5c6.jsonl`
(user messages, assistant messages, and the thread's own HANDOFF.md write at line 3758/3764), then cross-checked against the repo.
**This file authorizes no commit, push, install, or job submission.** Branch `codex/drop-py39` @ `5e551c9`; all thread work is uncommitted.

## Goal
Optional (`--scautoqc`) SCTK scAutoQC flags added in place to ViralScan H5ADs and UMAP outputs. Cells are only flagged, never removed.
Earlier in the thread: read the pachterlab LSCHWCP_2023 notebooks (see Dead ends / notebook work) and grounded the QC design empirically.

## Next action
Read-only: inspect the SRR8315713 (HSV1) FASTQs, merged reference and config, then draft an isolated run manifest (outputs outside `results/`, no scheduler).
Put the scope question (Open items 1) to the user before any execution.

## State
- Verified now: `PYTHONPATH=src /tmp/viralscan-scautoqc-env/bin/python -m pytest tests/test_scautoqc.py tests/test_run_safety.py -k scautoqc` gives **3 passed, 27 deselected** (rollout shows the same at lines 3028/3337/3395/3530). The 3 tests are 2 in `test_scautoqc.py` (in-place matrix preservation; CLI opt-in and run-config) plus 1 in `test_run_safety.py` (fingerprint only when enabled).
- Rollout-only receipts (not re-run): Ruff, format, `git diff --check` and `environment.yml` parse OK; disposable Snakemake `scautoqc` + `umap` rules passed on HHV6B GSM6416024 (5,310 cells: 3,787 pass / 1,523 fail) and `scautoqc` on EBV (2,647: 2,006 / 641); count matrix and obs/var indices identical before and after. Detection receipts there were fixtures, so this is wiring evidence, not detection validation.
- Live Scanpy BioMart + `mitochondrial_genes` query worked on a temp HHV6B copy: 35,085 of 47,259 features got symbols, 13 mito; cache written to `results/scautoqc_gene_annotations.tsv`. `pybiomart` was a missing import until added to extras.
- Uncommitted, thread-owned (confirmed by `git status`): new `src/viralscan/scautoqc.py` (241 lines), `src/viralscan/scripts/scautoqc.py` (70), `tests/test_scautoqc.py` (51); edits to `Snakefile` (+15: opt-in rule + `scautoqc_done` into detection inputs), `cli/parsers.py` (+9 `--scautoqc`), `menu.py` (+1), `runconfig.py` (+2), `scripts/umap.py` (+42/-6), `environment.yml` (+7/-1: python-igraph, leidenalg, pip sctk==0.2.0 + pybiomart), `pyproject.toml` (full extra: `scanpy[leiden]`, sctk==0.2.0, pybiomart). Also touched: `run_safety.py` (+1, fingerprint skip when False), `tests/test_run_safety.py` (+21), docs `cli_reference.md`, `installation.md`, `output_reference.md`.
- Scratch (exist): `/tmp/viralscan_scautoqc_{pilot,current_test,livequery,snakemake}`, plus env `/tmp/viralscan-scautoqc-env` (py3.11, scanpy 1.11.5, sctk 0.2.0, plotly, snakemake; no kb/kallisto/bustools).
- Empirical pilot (rollout only; exploratory): all 13 HHV6B Track A samples, 85,399 cells: default `^HB` passes 22,694 cells (26.6% of viral molecules), curated globin passes 41,467 (70.2%); EBV `combined_off` 3 samples, 9,755 cells: 4,833 to 7,322 pass. Default `^HB` also hits HBEGF/HBP1/HBS1L (the repo test now asserts HBEGF is not flagged).

## Locked decisions (user words, 2026-10-18/19)
- "note to just flag the cells for now, do not remove the poor QC cells"
- "include the qc workflow integral in UMAP outputs and initial h5ad construction ... add flags to either do it or not, so its optional ... design the function such that the outputs can modify objects in place so not starting from scratch"
- "for resulting single-cell dataset, flag QC cells using scautoqc as in teichlab ... automatic_qc"; use scanpy `queries.biomart_annotations` / `mitochondrial_genes` for symbols and mito/ribo/Hb sets.
- Thread conclusion (its own, not the user's): flags are review metadata; marker definitions (globin, ribosomal) are a sensitivity choice, not validated biology.

## Dead ends / cautions
- SCTK 0.2.0 lacks the tutorial's `multi_resolution_cluster_qc`; only cellwise + generic clustering ran. SCTK 0.2.0 raises TypeError when every cell fails (guarded at `scautoqc.py:180`).
- `from __future__` in a Snakemake `script:` module is a SyntaxError (the thread hit it once; `scripts/scautoqc.py` is clean now). `src/viralscan/scautoqc.py` does use it, which is fine since it is imported, not a rule script.
- Env `viralscan_test_full` was reported absent/broken in-thread (Numba import error elsewhere); I could not verify that. Do not rely on it; the `/tmp` env passed.
- Thread deleted the prior HANDOFF.md (a perf-delivery handoff) and wrote its own; that file has since been replaced again by another session. Do not look for the thread's HANDOFF.md.
- Notebook work (LSCHWCP) was analysis only, no repo edits: host masking is a sensitivity axis; validate with real positive controls; keep read-competition ties/no-hits; study thresholds are not transferable; use group-disjoint host-response folds. Top suggestion: a controlled host-filter-strategy comparison with sequence-matched positive controls. `viralscan_lshwcp_spec` backlog is stale versus PLAN.md.

## Read first
- `src/viralscan/scautoqc.py:54` (`add_scautoqc`), `:11` (`query_gene_annotations`), `:180` (all-fail guard)
- `src/viralscan/scripts/scautoqc.py:13` (run), `src/viralscan/Snakefile:252` (rule), `:133` (detection input)
- `src/viralscan/scripts/umap.py:194` (QC call) and `:214` (legacy count/gene cell filter skipped when `--scautoqc`)
- `tests/test_scautoqc.py:1`, `tests/test_run_safety.py` (`test_scautoqc_is_fingerprinted_only_when_enabled`)

## Open items
1. **Needs a user decision:** real FASTQ-to-UMAP validation on SRR8315713 (`benchmark_inputs/reference_strategy/SRR8315713/`, 1.6 + 4.3 GB gz; local FASTQ sets run 5-9 GB per pair). Needs `kb`/kallisto/bustools (none on PATH now) in an isolated pinned env, and a reference (candidate `viralscan_showcase/.../fullrun/refs/merged/`). Decide scope, CPU/RAM, and who installs.
2. Scope risk: `scautoqc` runs after multimap, before cell calling, so it may process all barcodes. Confirm intended cell set and runtime.
3. Behavior change to review: with `--scautoqc`, `umap.py:214` skips the legacy n_counts/n_genes filter, so standard UMAP/HVG also runs on low-count cells.
4. `PLAN.md` has no scautoqc row (grep: none). CLAUDE.md requires a PLAN.md row and update before any commit; add it first.
5. Clusterwise QC (multi-resolution) and ribosomal-pattern sensitivity across samples remain untested; `twostep` / `combined_artefact` EBV runs not compared.

## Verification status of this reconstruction
- Verified against repo: file list, line counts, diff content, docs edits, scratch dirs, SRR8315713 FASTQ sizes, absence of kb/kallisto/bustools, 3-test pass (re-run).
- Not verifiable here: HHV6B/EBV cell counts and retention percentages, Snakemake rule runs, live BioMart result (rollout and `/tmp` summaries only), Ruff/format results, and the state of `viralscan_test_full`.
