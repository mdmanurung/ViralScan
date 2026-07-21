# ViralScan — Publication Readiness (delta review, 2026-07-17)

_Independent multi-agent delta review (7 dimensions, adversarially verified) against
`docs/PUBLICATION_READINESS.md` (2026-07-03). Supersedes that doc where they conflict;
see "Corrections to the 2026-07-03 baseline" at the end. A companion code-craft review
lives at `.claude/reviews/main/summary.md` (Linus-style)._

## Verdict — the two-track asymmetry holds, and it's the headline

| Track | State | Gating work |
|-------|-------|-------------|
| **A. Software release** | **Near-ready — days.** CI green (582 tests), Trusted-Publisher release flow, containers, docs all in place. | ~3 concrete fixes before tagging (below), then user-gated tag→publish + post-tag DOIs. |
| **B. Journal manuscript** | **Gated on one venue decision that cascades — plus 2 venue-independent blockers.** Scientifically honest core. | Fix the ethics/data-availability + authorship blockers regardless of venue; then pick a venue (that choice sets the severity of everything else). |

Do not flatten these into one list — engineering is days-away and mechanical; the journal
track has real blockers and a decision that changes what "blocker" even means.

---

## What is done well (don't redo)

**Engineering**
- **Config-boolean-string hazard is fully closed.** `runconfig._coerce_bool` + explicit
  `_TRUE/_FALSE_STRINGS` eliminated the old `"True"`/`"False"` YAML pitfall. Verified: no
  module still does `== "True"`. (The CLAUDE.md "known pitfall" no longer exists.)
- **`RunConfig` / `RunContext` / `KbCountOutputs`** are exemplary seams: single coercion
  checkpoint, a testable frozen context that killed fake-Snakemake mocks, and the kb output
  path-convention centralized as properties. `to_snakemake_config_args()` iterates
  `dataclasses.fields`, so a new field can't silently drop from the wire format.
- **Release engineering:** dynamic version single-sourcing with a tag↔`__init__.py` check,
  PyPI Trusted Publishing (OIDC, no stored token), correct `MANIFEST.in`/package-data globs,
  tight `.dockerignore`, structurally-correct bioconda recipe.
- **CI:** 8-combo OS×Py matrix (fail-fast off) + separate lint/security/integration/conda-env
  jobs; disciplined network gating; audit-anchored regression tests; a `test_docs_consistency`
  guard against version/default drift.
- **Determinism:** fixed seeds (`DEFAULT_SEEDS=[0,1,10,42,100,1234]`) throughout host-response.
- **Subprocess discipline:** list-form `subprocess`, no `shell=True`, no bare `exit()` (`_die`),
  `BooleanOptionalAction` for CLI bools.

**Manuscript / science**
- **AUROC depth-confound self-audit is exemplary honesty** — headline 0.866 reported *with*
  the depth-alone baseline (0.967) and the depth-controlled band (0.64–0.72). (This closes the
  2026-07-03 "B1 integrity blocker" — it is done.)
- **HSV-1 denominator+threshold artifact** rigorously decomposed (0.55% all-barcode → 13–18%
  over called cells); **called-cell denominator** is a genuine methodological contribution.
- **SARS-CoV-2 = 0** cross-validated; **anellovirus treated as a host-homology limitation,
  not signal**; EBV/HSV-1 aligner-comparison exclusion (GTF CDS-only) clearly explained so the
  single clean HHV-6B comparison doesn't read as cherry-picking; **3.64× EM** internally
  consistent and correctly scoped.
- **Reproducibility positives:** figures regenerate from committed TSVs; seeds pinned; the three
  public benchmark datasets carry full GEO/SRA accessions; Zenodo data DOI confirmed live.
- **Unusually honest Discussion** that explicitly disclaims the head-to-head and FPR/FNR gaps.

---

## Track A — Engineering: fix before tagging

**Must-fix (before `git tag`):**

1. **Strip institutional paths from *shipped* package files.** `src/viralscan/reference_strategy.py:24`
   ships `STAR_BIN = "/exports/…/mdmanurung/conda/envs/starsolo/bin/STAR"`, an `os.chdir("/exports/…")`,
   and private `fastq_root`/`genome_fasta` defaults — and it is **never imported at runtime**
   (verified). It is benchmark scaffolding: move it out of `src/viralscan/` (to `analysis/` or
   `scripts/`). Both this review and the Linus code review independently flagged it as the #1 code
   issue. _(Also remove the guarded private R-lib path in `scripts/emptydrops.R:34`. Verified: the
   `dir.exists()` guard makes it functionally harmless off-HPC — this is optics/hygiene, not a
   functional break — but it should not ship.)_
2. **Reconcile CHANGELOG + version, then tag.** No git tag exists, and HEAD is ~20 feature commits
   past the `v2.5.0` bump (EVE flags, sibling cross-mapping, `--genome-dlist`, and a **behavior-
   changing default**: multimap-method `equal`→`host-conservative`), all sitting in CHANGELOG
   `[Unreleased]`. Tagging HEAD as `v2.5.0` would ship behavior the v2.5.0 notes don't describe.
   Bump to **2.6.0**, close `[Unreleased]` as `[2.6.0]`, sync `CITATION.cff` / `Dockerfile` /
   `Singularity.def` / `conda-recipe/meta.yaml`, then tag.
3. **Add `--no-deps` to the Dockerfile pip install.** CI uses the `--no-deps` workaround for the
   `connection_pool`/setuptools build failure; the Dockerfile does not, so a released container
   risks failing to build. (Verified the CI workaround; apply the same to the container.)

**Should-fix (low effort, batch with the above):**
- **Declare `anndata`** in `pyproject.toml` (it's an *eager* top-level import in `scripts/multimap.py`,
  currently satisfied only transitively via scanpy; CI's `--no-deps` install is blind to this).
  Sync `environment.yml` too (missing `anndata` **and** `scikit-learn`).
- **Governance / .gitignore:** `covid_viralscan/results_genomic/` and `results_hostfilter/` are
  **untracked** right now — one `git add .` from committing human-subjects clinical outputs. Add
  explicit ignores. (See the manuscript IRB blocker below — same underlying data.)

**Post-tag / polish (not release-gating):**
- Confirm the **PyPI Trusted Publisher is actually configured** for this repo/environment, or the
  release job fails on tag push (not verifiable from the repo — a pre-tag checklist item).
- Pin a **Snakemake ceiling** (`>=7,<8` or validate v8 — v7/v8 have breaking API changes; `--use-conda`
  is deprecated in SM8).
- 195 GTFs are still git-tracked (**downgraded to minor** — verified the wheel *and* sdist exclude
  them; only clone weight is affected); `conda-recipe/README.md` stale at 2.4.0; `LICENSE` year (2025)
  + holder ("ViralScan", not a legal person); CITATION.cff lacks ORCID + a Zenodo software DOI (post-tag);
  container images tag-pinned not digest-pinned; no arm64; Python 3.13 absent from the matrix.

**Test-coverage gaps (matter for a "does it actually run?" reviewer):**
- Integration tests install `kb-python` but **never invoke `kb count`** on real FASTQs — the core
  quantification has **zero automated end-to-end coverage** (verified). Add one tiny real-FASTQ
  smoke run.
- `evidence_run.py`'s orchestrator (`run_evidence()`) is mocked out — never executed in tests.

**Docs gaps (block a clean tag / first impression):**
- **EVE artifact output columns undocumented** in `output_reference.md` and CHANGELOG.
- **Four depth-confound `hostresponse` flags** now on the main CLI are **undocumented** in
  `cli_reference.md`.
- CHANGELOG `[2.5.0]` reads as a tagged release that never happened (fixed by must-fix #2).
- No `hostresponse` vignette; `quickstart.md` never mentions `hostresponse` (venue-conditional).

---

## Track B — Journal: two venue-independent blockers, then a decision

**Venue-independent blockers (fix regardless of target):**

1. **Data Availability + ethics** (`manuscript_draft.md:140`). **GSE210063 is mislabeled
   "COVID-era clinical"** — it is HHV-6B CAR-T (Lareau 2023; the table at line 170 labels it
   correctly, so the paragraph contradicts the table). And the actual COVID libraries
   (LUM-SJ-x213-g / x216-g) are called "unpublished" with **no accession and no IRB/consent/ethics
   statement anywhere** (grep-confirmed zero ethics text). No bioinformatics venue accepts human
   clinical data without an ethics disclosure. One paragraph rewrite fixes both.
2. **Authorship** — every author field, affiliation, ORCID, CRediT role, lead contact, and
   declaration-of-interests is a placeholder. User-gated; cannot submit anywhere without it. (Also
   reconcile the cross-artifact author inconsistency: `pyproject.toml`/`CITATION.cff` list only
   "Emma Vonk" with a personal email; the git committer is `mdmanurung`; the LICENSE holder is
   "ViralScan"; the manuscript list is a third placeholder.)

**The master gate — pick a venue (it cascades):**
The manuscript header targets **Cell Reports Methods** and is formatted with STAR Methods + eTOC,
but PLAN P22.7c lists the **application-note tier** (Bioinformatics App Note / PLOS Comp Biol /
GigaScience). These are irreconcilable without a decision, and the decision sets the severity of
almost everything else:

- **Methods venue (Cell Reports Methods):** a **head-to-head vs Venus or ViralTrack** is a near-
  blocker (months of work), and the private benchmark index becomes a hard reproducibility blocker.
- **App-note / resource (PLOS Comp Biol, GigaScience, Bioinformatics App Note):** the head-to-head
  is an already-written *disclosed limitation* (Discussion line 124), not a blocker. GigaScience
  additionally requires the **benchmark index to be publicly reproducible**; PLOS Comp Biol is the
  **lowest-format-friction** path given the current STAR-Methods-shaped draft.

**Recommendation:** target **PLOS Computational Biology** (or Bioinformatics App Note). Reserve Cell
Reports Methods only if you're willing to add the dedicated-tool head-to-head.

**Rigor/honesty majors to fix before submission (venue-independent unless noted):**
- **HHV-6B threshold-mixing:** the manuscript compares ViralScan's 0.152% (**≥1 UMI**) against
  Lareau's 0.2%, but Lareau's 0.2% is a **≥10 UMI super-expressor** rate — apples-to-oranges.
  Match thresholds or state the mismatch explicitly.
- **"Under 2 h on an 8-core node"** runtime claim (abstract) is **unsupported** by any documented
  benchmark in the repo — back it with a timed run or drop the number.
- **eTOC/Summary frames EM as the headline feature, but the shipped default is `host-conservative`**
  — feature/default mismatch a reviewer will catch.
- **Citations:** Venus (ref 10) points to a bioRxiv DOI — reconcile against the published version;
  **ViralTrack is named in the Discussion but has no reference entry** — add it.
- **Benchmark index reproducibility** (private "evonk" Serratus index; Figure 2 not reproducible
  from a documented `build-ref` command) — **venue-conditional** (blocker at methods venue; major at
  app-note; publish the index or re-run from a public one). **Software DOI** absent from CITATION.cff
  (post-tag).

---

## Cross-cutting themes (one fix resolves several findings)
- **Hardcoded institutional paths + private index** — appears in 4 dimensions (reference_strategy.py,
  emptydrops.R, benchmark SLURM scripts, Figure 2 index). One cleanup pass + one reproducibility
  decision.
- **Version / CHANGELOG / DOI drift** — packaging, docs, reproducibility, manuscript all point at the
  same root: HEAD past 2.5.0, no tag, "DOI pending." Reconcile CHANGELOG → tag → fill DOIs, once.
- **The Data Availability paragraph** (line 140) is the *same* paragraph flagged by two dimensions;
  one rewrite fixes the GSE mislabel and the IRB gap.
- **Stale internal pointers** — "557 tests" (→582), PLAN's 0.845 AUROC pointer, and a "TTV read-origin
  test pending" HTML comment at manuscript line 4 (F-005 is closed; prose at line 122 is already
  correct). One doc-sync pass.

## Coverage gaps the dimensions didn't own (worth a look)
- **Community-health files** (CODE_OF_CONDUCT, issue/PR templates) — expected for a healthy public
  release; CONTRIBUTING exists, the rest don't.
- **Data licensing of bundled third-party content** — MIT covers the code, not the 195 GTFs, the
  Serratus-derived index, or Ensembl/GRCh38-derived content; resource venues want explicit data-
  licensing statements.
- **Human-subjects governance** for the untracked COVID clinical outputs (see Track A governance item).

---

## Corrections to the 2026-07-03 baseline (now stale)
- **"557 tests passing" → 582** (20 network-deselected).
- **"B1 AUROC integrity blocker: must-fix" → DONE.** Manuscript reports 0.866 caveated with depth
  baselines; the Discussion names the confound. No longer open.
- **"Pre-tag cleanup DONE; only tag/publish + DOIs remain" → not fully accurate.** HEAD drifted ~20
  commits past 2.5.0 (behavior-changing default in `[Unreleased]`); institutional paths still ship;
  no tag exists.
- **VIRTUS2 citation → resolved** (software-repo form). **Zenodo data DOI → confirmed live.**
- **Framing around "Cell Reports Methods" → stale/undecided.** PLAN P22.7c targets the app-note tier;
  the venue is the open master decision, not a settled assumption.

_Method: 7 parallel dimension reviewers grounded as a delta vs the 2026-07-03 doc, each finding's
blocker/major claims independently re-checked against current source (adversarial verify), plus a
completeness critic. Severities here reflect the post-verification adjustments (e.g. emptydrops.R
libpath and the 195 git-tracked GTFs were both downgraded to minor on verification)._
