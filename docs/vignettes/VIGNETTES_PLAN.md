---
orphan: true
---

# ViralScan vignette suite — plan & status

Vignettes are the **reproducible public face** of ViralScan. Every notebook is
designed so an outside reader on a laptop can either run it end-to-end, or is
told exactly why not and how to obtain the data. No institutional paths, no
private `evonk` index (contrast `docs/showcase_runbook.md`, which is internal).

Each vignette teaches one v3 feature on synthetic data, with a negative or
ambiguous example. Real-library numbers from earlier drafts were removed (DOC-06);
benchmark claims live in the claims registry, not in notebooks.

## Data strategy (tiers)

| Tier | Source | Ships? | CI-executed? |
|------|--------|--------|--------------|
| Synthetic | AnnData/FASTQ built in-notebook | code only | yes (fast, deterministic) |
| Panel fetch | `viralscan data fetch` (Zenodo, SHA-256) | downloads | network-gated |
| Public-data demo | 1M-read EBV subsample (SRR12682296) via ENA/SRA | no (external) | no — `[skip-ci]` |

Mechanics/API demos use synthetic data and **execute in CI**. Anything
that builds a human index or runs full `kb count` is `# [skip-ci]` with a
copy-paste public download block. Public-data vignettes use the **1M-read
subsample**, never full 112M-read depth.

## Conventions

- Jupyter under `docs/vignettes/`, matching the existing two + Sphinx (`docs/conf.py`).
- Header block: the paper claim reproduced + minimum `viralscan` version
  (v3 only: `em-global`/`em-cell`, `viral_molecules_total_est` columns).
- Keep consistent with `test_docs_consistency.py` (CLI flags / output columns).
- Pass `RunConfig(...)`, **not** a dict, to `cell_type_enrichment()` (the old
  enrichment vignette passed a dict and was silently broken under `[skip-ci]`).

## Catalog & status

Legend: [ ] todo · [~] in progress · [x] done · CI = executes in CI · SKIP = `[skip-ci]`

### Tier A — core on-ramp
- [x] **V1 quickstart_fastq_to_viral_load** — default quant + `check-whitelist`; reads `viral_summary.tsv`/`report.html`. Data: 1M EBV subsample. SKIP. *(revamp `basic_usage.ipynb`)*
- [x] **V2 building_a_reference** — `data fetch` (Zenodo) + `build-ref`; toy build executes, human build SKIP.

### Tier B — differentiators
- [x] **V3 multimapping_correction** — `--multimap-method` (em-global/em-cell), `rerun-multimap`, `selected-method`; synthetic EM, plus ambiguous (rare host) and no-op negative cases. CI.
- [x] **V4 cell_calling_denominators** — `--cell-calling external` executed, `emptydrops` described (needs R), `knee` sensitivity-only; `pct_infected` vs `pct_infected_called`; wrong-list failure. Synthetic. CI.
- [x] **V5 cell_type_enrichment** — Fisher OR + BH; `--cell-types`. Synthetic. CI. *(upgrade existing; fix dict→RunConfig)*

### Tier C — trust & advanced
- [x] **V6 specificity_true_negative** — zero-count interpretation + `check-whitelist` chemistry preflight + zero-over-real-cells check. Synthetic. CI; serves as the negative example.
- [x] **V7 qc_and_read_evidence** — `evidence` subcommand, sibling cross-mapping flags (ambiguous EM-bleed example plus no-flag negatives), `accession_breadth`, `host_viral_ambig_fraction`; EVE/host-homology caveat callout. Synthetic. CI.
- [x] **V8 host_response_depth_control** — `hostresponse` label/matching helpers on a synthetic no-host-biology cohort; raw depth-alone AUROC vs cpm/depth-matched, plus a weak-signal ambiguous case. CI.

## Cross-cutting
- [x] `docs/vignettes/README.md` index mapping vignette → subcommand.
- [x] Wire into `docs/index.md` / Sphinx toctree.
- [x] Shared "data setup" block with public download commands (no `/exports`).
