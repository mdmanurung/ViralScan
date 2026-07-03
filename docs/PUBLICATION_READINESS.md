# ViralScan — Publication Readiness Assessment

_Assessed 2026-07-03 (v2.4.0 + v2.5 scientific-hardening on `claude/multimap-memory-and-showcase`)._

## Verdict — two separate tracks

Readiness splits cleanly into two efforts with **very different timelines**. Conflating
them hides a days-vs-months difference.

| Track | Readiness | Gating work |
|-------|-----------|-------------|
| **A. Software release** (PyPI / bioconda / container) | **~85% — days of mechanical cleanup** | repo hygiene, version sync, DOIs, conda sha256 |
| **B. Methods manuscript** (Cell Reports Methods) | **Not ready — one integrity fix + one comparison are the long pole** | correct the confounded AUC in the artifact; add a head-to-head vs a dedicated tool |

**The one thing that genuinely blocks an *honest* submission** is the depth-confounded
host-response AUROC still printed in the manuscript (§3.4 / Figure 2). The fix already
exists in the code (v2.5); it just has to propagate to the manuscript artifact and figure.

---

## Track A — Software release (near-ready)

The package itself is in good shape: MIT `LICENSE` present, 557 tests passing, CI matrix
(8 OS×Py combos) + lint/mypy/security/integration jobs, Sphinx/RTD docs, CLI reference,
vignettes, `environment.yml` with pinned external tools, and the 195 GTFs externalized to
Zenodo. The remaining blockers are mechanical:

- [ ] **Repo hygiene / `.gitignore` gaps (do first).** Untracked HPC artifacts sit in the
  repo root (`Log.out` ~380 KB, `vs_val_*.err/out`) and `.claude/`, `.specify/`,
  `benchmark_runs/`, `covid_viralscan/{data,results,viralscan_ref}/`, `logs/`,
  `docs/figures/` are unignored — a naive `git add .` before a public push would ship AI
  session logs and job files. Gitignore these before any public release.
- [ ] **`docs/conf.py` version is hardcoded `2.3.0`** → read it dynamically from
  `viralscan.__version__` (RTD currently advertises the wrong version on every page).
- [ ] **Commit `docs/figures/`** (workflow + benchmark PNG/PDF) — referenced by the docs
  build but untracked.
- [ ] **Stale `dist/`** holds a 2.3.0 wheel; rebuild on release (the release workflow does
  this on tag, but the local dir misleads).
- [ ] **Cut the version.** The branch already contains the v2.5 features but is versioned
  2.4.0 with them under CHANGELOG `[Unreleased]`. Decide 2.4.0-vs-2.5.0, bump
  `__version__`, move the changelog section, sync `CITATION.cff`/Dockerfile/Singularity.
- [ ] **bioconda:** fill `conda-recipe/meta.yaml` `sha256` (currently all-zeros placeholder)
  after the PyPI release.
- [ ] **Software DOI:** archive the tagged release on Zenodo and add the DOI to
  `CITATION.cff` (currently only a GitHub URL). Distinct from the data DOI.

**JOSS note (only if you also pursue a JOSS software paper):** JOSS additionally requires a
`paper.md` in JOSS format and a Zenodo software DOI. These are **not** gates for the Cell
Reports Methods track — ignore unless JOSS is a target.

Estimated effort: **~1 focused day**, none of it scientific.

---

## Track B — Methods manuscript (blocked)

The draft (`docs/manuscript_draft.md`, ~72% complete) is well-written: complete Intro,
Results §3.1–3.4, Discussion, and STAR Methods. The gaps are one integrity fix, one
scientific addition, and submission mechanics.

### B1. Integrity blocker — the confounded AUROC (must fix)

`docs/manuscript_draft.md:89` reports the host-response result as **"AUROC 0.845 ± 0.032"**,
and the Figure 2 footer reports the same regime. This is the **raw `>=10 UMI` label**, which
this session proved is depth-confounded: on the same data, **sequencing depth alone predicts
EBV status at AUC 0.967** — higher than the 0.866 model. The honest, depth-independent value
is **~0.67** (`--label cpm` → 0.672; `--depth-match` → 0.680; depth-alone drops to ~0.5).

Aggravating: the Discussion limitations (`:101`) list EM / BAM / Ensembl-IDs but **do not
mention the depth confound at all** — so 0.845 is presented *uncaveated*, and the section
heading (`:85`, "leakage-controlled") implies a rigor the reported run does not have (gene
leakage is controlled; depth is not).

**Mitigating:** the number appears *only* at §3.4 + the Figure 2 footer — **not** in the
abstract, highlights, or eTOC — so the revision is contained.

Fix (code already produces the honest numbers — see `analysis/hostresponse_ebv_matched/scripts/reconcile_package_hostresponse.py`):
1. Re-run `viralscan hostresponse … --label cpm` (and/or `--depth-match --differential
   --gene-symbols`) to regenerate `results/hostresponse_ebv_matched/hostresponse_metrics.csv`.
2. Update `:89` to the ~0.67 figure, **report the depth-alone baseline (0.97) alongside it**,
   and revise the `:85` heading.
3. Regenerate `docs/figures/figure2_benchmark.*` from the new CSV.
4. Replace the Ensembl-ID stable-feature list (`:91`) with **gene symbols** (v2.5
   `--gene-symbols`) and add one sentence of biological interpretation.

### B2. The real methods-paper gate — a head-to-head comparison

A tool paper needs at least one head-to-head comparison against a **dedicated** viral
scRNA-seq quantifier. Right now the only comparison is against **STARsolo used as a general
aligner** (one virus, one dataset, mismatched annotations: 16 vs 96 EBV loci), and the
"3.64× sensitivity" claim is measured against ViralScan's *own* alternate modes. Add a
comparison against **Venus (Luebbert et al., which the paper already cites)** or ViralTrack
on at least one shared dataset. This is close to mandatory and is the long pole.

### B3. Strengthens the paper (reviewers will push, not disqualifying)

- **FPR/FNR against known truth.** The README itself notes these are uncharacterized. The
  three published-study reproductions (HHV-6B 0.15%, HSV-1 13–18% over called cells, EBV) are
  *genuine external validation* — not circular — but a spike-in / simulation with ground-truth
  single-cell infection status would answer "what's the single-UMI false-positive rate?"
- **Reference-strategy benchmark is 4/12 complete** (EBV STARsolo rows fail on a FASTQ quality
  error; the three ViralScan two-step rows are blocked by `kb` not on PATH; HSV-1 STARsolo
  detects 0 cells). Either complete it or cut it from the paper — as-is it is not presentable.
- **Reproducibility gap:** the headline benchmarks used a private "evonk" Serratus index, not
  a `viralscan build-ref` output. A reviewer cannot reproduce Figure 2 from the documented
  build command. Re-run the key benchmark from a `build-ref` reference, or document the exact
  index and publish it.

### B4. Submission mechanics (user actions)

- **Authorship** — populate the author list, affiliations, ORCIDs, lead contact,
  acknowledgments, CRediT contributions, and declaration of interests (all placeholders).
  _This is the authors' call, not an automated edit._
- **Verify citations** — the VIRTUS2 reference (`:202`) is flagged as having no Crossref
  record; resolve to a real DOI, fall back to the original VIRTUS citation, or drop the
  "VIRTUS2-like" framing.
- **Verify the Zenodo data DOI** (`10.5281/zenodo.20112332`) is live.
- **Apply the Cell Reports Methods template** (structured abstract, STAR format, figure
  legends section, reference style) — PLAN P22.7c.

---

## What is genuinely solid (don't redo)

- EM recovers **3.64×** more EBV UMI than unique-only counting (internal, reproducible).
- Full-depth reproductions within published ranges: **HHV-6B 0.15%** (Lareau range),
  **HSV-1 13–18%** over called cells (Wyler 13–19%).
- **Called-cell vs all-barcode denominator** teaching point (HSV-1 0.55% → 13–18%) — a real
  methodological contribution, now in the tool (v2.5 SH1.5).
- STARsolo↔ViralScan concordance (Spearman r=0.45) with the interpretable LMP-1/EBNA
  annotation-dependence finding.
- The host-response depth-confound analysis itself — now the *honest* result and a
  reproducible in-package capability (v2.5 SH1.1–1.3, verified on real data to 3 decimals).

---

## Prioritized path to submission

1. **(days, mechanical)** Track A hygiene + version cut → software is releasable.
2. **(hours, code exists)** B1: propagate the honest AUC + gene symbols to §3.4 and Figure 2.
   *This is the honesty gate — do it regardless of everything else.*
3. **(weeks)** B2: one head-to-head vs Venus/ViralTrack.
4. **(weeks, optional-but-expected)** B3: FPR/FNR characterization; finish or cut the
   reference-strategy benchmark; re-run headline benchmarks from a public reference.
5. **(user)** B4: authorship, citation verification, journal template.

**Bottom line:** the *software* can be published within days. The *manuscript* is one
must-fix correction away from honest, and one head-to-head comparison away from competitive
for a methods venue. The correction is cheap because the fix already lives in the code.
