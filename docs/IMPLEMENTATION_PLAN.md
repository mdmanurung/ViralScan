# ViralScan roadmap — sequenced implementation plan

Back-to-back build order for every item in [`docs/ROADMAP.md`](ROADMAP.md),
dependency-ordered into phases. Each item carries an **autonomy tag** because
"implement all of it" spans three very different kinds of work:

- 🟢 **AUTO** — I can implement + test end-to-end now (code + synthetic tests).
- 🟡 **COMPUTE** — I can write the code/harness, but *running it to completion*
  needs a cluster, large refs, or many GB of real data.
- 🔴 **OWNER** — needs a human action or decision (release upload, DOI, ship-scope,
  installing third-party tools, signing key).

Global guardrails (every feature, no exceptions): work on a feature branch off
`main`; keep the full suite green; **golden-gate** anything touching
`build_multimap_layers` or detection numerics (all-layers rtol=1e-9); update
`PLAN.md` in the same commit (repo contract); don't break the CLI/output
contract; add tests before calling an item done.

## Dependency graph (the only hard orderings)

```
B3 ─→ B4                  A1 ─→ A5
A1, A3 ─→ F6              F1(merge) ─→ F2(release)
```
Everything else is independent and ordered by risk/value, not necessity.

---

## Milestone 0 — Ship checkpoint (do first; runs in parallel with Phase 1)

Get the current 25-commit `pub-readiness-hygiene` branch out the door so features
land on a released base, not a giant unmerged branch.

| # | Item | Tag | Notes |
|---|------|-----|-------|
| 0.1 | **F1** open PR #7 → main | 🟢 | I can open it now. |
| 0.2 | **F3** functional-script path hygiene (~14 `evonk`-path SLURM scripts) | 🟢/🔴 | Code is AUTO; the *ship-scope* decision (env-var-ize vs exclude) is OWNER. |
| 0.3 | **F4** stop tracking `ref/10x…whitelist.txt` + `analysis/**/*.pdf` | 🟢 | gitignore + fetch-doc. |
| 0.4 | **F1 merge → F2 release** | 🔴 | Merge, tag v2.6.0, PyPI/conda, Zenodo DOI — owner. Signing-key fix first. |

**Gate:** PR green in CI; branch merged. (F2 release deferred to Milestone 7 if you'd
rather ship once, with all features.)

---

## Phase 1 — Robustness & CI hardening (protect everything downstream) — ~2–3 d

| # | Item | Tag | Deliverable / gate |
|---|------|-----|--------------------|
| 1.1 | **D1** non-unique `var_names` safety | 🟢 | Validate/dedup var_names at load (or first-occurrence-safe `matrix_for_genes`); test with duplicate accessions → no crash. |
| 1.2 | **D2** EVE `gi\|…\|ref\|NC_…\|` sseqid unwrap | 🟢 | small regex fix + test. |
| 1.3 | **D3** `host_viral_ambig_fraction` >1.0 clamp/doc | 🟢 | clamp to [0,1] or document; test the edge. |
| 1.4 | **F5** CI gates: nbmake (6 CI vignettes) + deptry/bare-venv import | 🟢 | `.github` workflow; catches doc-rot + undeclared deps (both bit us this session). |

**Gate:** new CI jobs green; suite green. *Why first: these prevent silent breakage of
the features that follow (the enrichment vignette silently broke once; CI installs
`--no-deps` so a missing dep hides).* 

---

## Phase 2 — Host-response completion — ~2–3 d

| # | Item | Tag | Deliverable / gate |
|---|------|-----|--------------------|
| 2.1 | **C1** gene-symbol annotation + pathway enrichment | 🟢 | Ensembl→symbol (pybiomart/gget) + gseapy/enrichr on the **depth-controlled** labels; surface symbols in `hostresponse` outputs + HTML. Test on the tracked matched-EBV CSVs. |
| 2.2 | **C2** extra depth-control designs (propensity/E-values in report) | 🟢 | optional; extends existing depth diagnostics. |

**Gate:** `hostresponse` emits gene symbols + an enrichment table; tests green.

---

## Phase 3 — Evidence positional features (the requested cluster) — ~4–6 d

| # | Item | Tag | Deliverable / gate |
|---|------|-----|--------------------|
| 3.1 | **A1** read-start distribution + `--dedup {umi,markdup,none}` | 🟢 | `evidence.read_start_distribution()` + CLI + `read_start_profile.tsv` + plot. Full plan in [`todo/read-start-distribution.md`](../todo/read-start-distribution.md). Unit tests on synthetic SAM. |
| 3.2 | **A3** cell-level viral BAM with CB/UB tags | 🟢 | carry `(CB,UMI)` into BAM tags; IGV group-by-CB recipe. |
| 3.3 | **A5** subgenomic-RNA junction detection (deps A1) | 🟢 | leader–body junction pile-up caller on the read-start profile; coronavirus test fixture. |
| 3.4 | **F6** IGV + read-start vignette (deps A1, A3) | 🟢 | new `docs/vignettes/` notebook, CI-runnable on synthetic BAM. |

**Gate:** `viralscan evidence --read-start-profile` produces the TSV+plot; A5 flags a
planted junction; vignette executes.

---

## Phase 4 — Reference & specificity — ~3–5 d code (+ compute)

| # | Item | Tag | Deliverable / gate |
|---|------|-----|--------------------|
| 4.1 | **G1** reference-annotation provenance in outputs | 🟢 | record + surface the exact viral accessions/GTF version used (manuscript argues for explicit annotation reporting). |
| 4.2 | **A4** combined-genome D-list specificity mode | 🟡 | `build-ref --genome-dlist genome.fa` code is AUTO + a small-genome test; the *real* human-genome D-list build is COMPUTE (~64 GB / ~8 h on a cluster). |
| 4.3 | **G2** mouse host support (`build-ref --host mouse`) | 🟢 | Ensembl mouse cDNA path + test. |

**Gate:** D-list masks a planted host-homologous k-mer on a toy genome; mouse ref builds.

---

## Phase 5 — Validation & benchmarking — 🟡🔴 (code AUTO, running gated)

The highest paper-credibility work, but it needs a cluster, real samples, and
third-party tools. I write the harnesses; running them to completion is COMPUTE/OWNER.

| # | Item | Tag | Deliverable / gate |
|---|------|-----|--------------------|
| 5.1 | **B3** complete the reference-strategy benchmark (4/12 rows) | 🟡 | re-run the failed/blocked SLURM conditions → analyzable `reference_strategy_benchmark.tsv`. Needs the cluster + FASTQs. |
| 5.2 | **B4** Selectivity Index + noisy-channel bits (deps B3) | 🟢 | analysis on B3 output — AUTO once B3 data exists. |
| 5.3 | **B1** head-to-head vs Venus + ViralTrack | 🟡🔴 | harness AUTO; installing the tools + running on EBV/HHV-6B/HSV-1 is OWNER/COMPUTE. |
| 5.4 | **B2** gold-standard truth panel (planted-read simulation) | 🟡 | AUTO: the spike-in simulator + FP/FN/ROC/LoD scoring code + a small synthetic demo. COMPUTE: running at scale on real host libraries. *Biggest credibility lever.* |

**Gate:** B2 simulator recovers a known planted rate on synthetic data; B1/B3 tables
populate when run on the cluster.

---

## Phase 6 — Performance & scale — ~3–6 d

| # | Item | Tag | Deliverable / gate |
|---|------|-----|--------------------|
| 6.1 | **E2** chunked/streaming BUS aggregation | 🟢 | cap peak RSS on deep samples; golden-gate identical output. |
| 6.2 | **E1** numba the collapse loop | 🟢 | ragged CSR arrays + `@njit`; golden-gate; benchmark another 5–10×. |
| 6.3 | **A2** per-cluster / per-cell EM | 🟢🟡 | `--multimap-scope {global,cluster,cell}`; code AUTO, but validating it *helps* (vs the global default) needs the deep samples → COMPUTE. Research-grade; benchmark before default-on. |

**Gate:** E1/E2 golden-identical + faster/leaner; A2 matches global scope when `scope=global`.

---

## Milestone 7 — Release (owner) — 🔴

**F2** cut the tag, publish PyPI/conda, mint the Zenodo DOI for the GTF panel, update the
manuscript data-availability. Fix the SSH signing key first so the release commits are
signed.

---

## Critical path & honest scope

- **What I can build back-to-back autonomously (🟢):** Phase 1 → 2 → 3 → most of 4 →
  5.2/5.4-code → 6. That's the bulk of the *code*, each with tests and golden gates,
  landing as a sequence of small PRs. Realistic order-of-magnitude: several focused
  sessions per phase.
- **What blocks (🟡🔴):** the D-list build (A4), the benchmarks (B1/B3), the at-scale
  truth-panel run (B2), and every release/ship-scope/DOI/signing action. I'll write the
  code and stop at the "run this on the cluster" / "you decide ship scope" boundary.
- **Recommended real order:** ship (M0) → harden (P1) → the features you asked for (P3,
  with C1 from P2 folded in) → specificity (P4) → then the resource-gated validation and
  perf as compute frees up.

## How I'll drive it
One feature per branch/commit, PLAN.md + `.living/` updated each step, full suite +
golden green before moving on, and I'll pause at each 🟡/🔴 boundary and hand you the
exact command/decision needed rather than fake completion.
