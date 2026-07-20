---
orphan: true
---

# ViralScan vignette suite — plan & status

Vignettes are the **reproducible public face** of ViralScan. Every notebook is
designed so an outside reader on a laptop can either run it end-to-end, or is
told exactly why not and how to obtain the data. No institutional paths, no
private `evonk` index (contrast `docs/showcase_runbook.md`, which is internal).

Each vignette is grounded in one of the manuscript's result narratives
(`docs/manuscript_draft.md`) so the tutorials double as reproducible evidence
for the paper — including the caveats, not just the headline numbers.

## Data strategy (tiers)

| Tier | Source | Ships? | CI-executed? |
|------|--------|--------|--------------|
| Synthetic | AnnData/FASTQ built in-notebook | code only | yes (fast, deterministic) |
| Panel fetch | `viralscan data fetch` (Zenodo, SHA-256) | downloads | network-gated |
| Tracked outputs | committed `results/hostresponse_ebv_matched/*.csv/.txt` | yes | yes |
| Reproduce-paper | 1M-read EBV subsample (SRR12682296) via ENA/SRA | no (external) | no — `[skip-ci]` |

Mechanics/API demos use synthetic or tracked data and **execute in CI**. Anything
that builds a human index or runs full `kb count` is `# [skip-ci]` with a
copy-paste public download block. Reproduce-paper vignettes use the **1M-read
subsample** the manuscript uses, never full 112M-read depth.

## Conventions

- Jupyter under `docs/vignettes/`, matching the existing two + Sphinx (`docs/conf.py`).
- Header block: the paper claim reproduced + minimum `viralscan` version
  (multimap flags require ≥ 2.3.0; current 2.6.0).
- Keep consistent with `test_docs_consistency.py` (CLI flags / output columns).
- Pass `RunConfig(...)`, **not** a dict, to `cell_type_enrichment()` (the old
  enrichment vignette passed a dict and was silently broken under `[skip-ci]`).

## Catalog & status

Legend: [ ] todo · [~] in progress · [x] done · CI = executes in CI · SKIP = `[skip-ci]`

### Tier A — core on-ramp
- [x] **V1 quickstart_fastq_to_viral_load** — default quant + `check-whitelist`; reads `viral_summary.tsv`/`report.html`. Data: 1M EBV subsample. SKIP. *(revamp `basic_usage.ipynb`)*
- [x] **V2 building_a_reference** — `data fetch` (Zenodo) + `build-ref`; toy build executes, human build SKIP.

### Tier B — differentiators
- [x] **V3 multimapping_correction** — `--multimap-method`, `rerun-multimap`, `--multimap-primary-call`; synthetic EM layers. CI. ⚠ depends on the uncommitted `--multimap-primary-call` matrix diff — land that first or pin to released behavior. Min ver 2.3.0.
- [x] **V4 cell_calling_denominators** — `--cell-calling knee|emptydrops|external`; `pct_infected` vs `pct_infected_called`. Synthetic. CI.
- [x] **V5 cell_type_enrichment** — Fisher OR + BH; `--cell-types`. Synthetic. CI. *(upgrade existing; fix dict→RunConfig)*

### Tier C — trust & advanced
- [x] **V6 specificity_true_negative** — SARS-CoV-2=0 interpretation + `check-whitelist` chemistry preflight + cell-calling concordance. Whitelist demo CI; real run SKIP.
- [x] **V7 qc_and_read_evidence** — `evidence` subcommand, sibling cross-mapping flags, `accession_breadth`, `host_viral_ambig_fraction`; EVE/host-homology caveat callout. Synthetic. CI.
- [x] **V8 host_response_depth_control** — `hostresponse`; reads committed result CSVs; teaches raw AUROC 0.87 → 0.64–0.72 under depth control. CI (tracked data).

## Cross-cutting
- [x] `docs/vignettes/README.md` index mapping vignette → subcommand → paper claim.
- [x] Wire into `docs/index.md` / Sphinx toctree.
- [x] Shared "data setup" block with public download commands (no `/exports`).
