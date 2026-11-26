# ViralScan roadmap to 3.0.0 stable

Status as of 2026-11-17. Status words come from `PLAN.md`, the reproducibility repo's tracker and its gate dashboard. This page makes no analysis claim, so it cites no figure. Rows named `LIKE-THIS-01` are `PLAN.md` rows.

## 1. Where we are

**Merged to `main`:**
- Run contract: one defaults table (`CONTRACT-DEFAULTS-01`), one fingerprint rule with `--resume` refused across versions (`CONTRACT-FP-01`), version constants and a stricter manifest schema (`CONTRACT-VERSIONS-01`), a visible EM non-convergence warning (`CONTRACT-DIAG-01`).
- Interfaces: deprecation helper, exit-code table, API stability section and release policy (Phase 2); packaging fixes, the CHANGELOG breaking block, a migration section and `docs/limitations.md` (Phase 3).
- Safety net: golden digests of the tiny run's deterministic outputs (`GOLDEN-01`).

**Open:** `TABLE-SCHEMA-01` (column registry for the 12 default tables). Older pull requests from 2026-11-04: the G3 amendment schema (#195), the version bump for outside laboratories (#197) and the 2026-11-04 decisions (#206).

**Not started:** consolidation Phases B to F, the first public tag, outside-laboratory testing, G3 to G5, the stable release.

## 2. Release line

`3.0.0b1` for the outside laboratories, then `3.0.0rc1`, then `3.0.0`. `3.0.0rc1` was never published, so no public version exists yet.

| Tag | Meaning | Gate to cut it |
|---|---|---|
| `3.0.0b1` | Laboratories test the install, the first run and the API. The CLI and outputs can still change on their feedback. No validation claim. | Consolidation Phases A to F merged; owner steps in section 5 done; tester packet ready. |
| `3.0.0rc1` | Contract frozen; fixes only. Defaults provisional until G5. | Laboratory feedback closed (`RC-05`); default decision rows approved from G5 output; `defaults_status` set to `final`. |
| `3.0.0` | Stable and citable. | G3, G4, G5 and G6 closed; DOIs; archives; manuscript audited. |

A laboratory-driven rename is allowed during `b1`. After `rc1`, anything breaking is 4.0.

## 3. Immediate fixes

1. **`TABLE-SCHEMA-01` ship scope.** Done: `table_columns.json` is listed in `config/public_ship_scope.json`. A new packaged file always needs that entry; `tests/test_ship_scope_tracked_tree.py` catches it only with `PYTHONPATH=.`.
2. **Older pull requests** (decided by the maintainer, one at a time):
   - #206 (2026-11-04 decisions): its content is already on `main` (no `strong` tier in `evidence_tiers.py`; the `GRILL-DECISIONS-01` row exists). Recommend closing as superseded.
   - #197 (`3.0.0b1` version bump): conflicts with `main`. Keep as the model for the beta bump and redo it after the consolidation.
   - #195 (G3 amendment schema): mergeable. Keep until its companion pull request in the reproducibility repo lands.
3. **One re-pin of the reproducibility repo.** `docs/output_reference.md`, `README.md`, `docs/faq.md` and `docs/v3_counting_contract.md` are hash-pinned there, so five of its tests fail by design against a newer `main`. After the Phase B and C documentation merges: set `PACKAGE_BASE_SHA`, refresh `pkg/`, update the inventory rows, keep claim markers byte-identical, run both validators. This is a commit in the other repository and needs the maintainer's approval.

## 4. Consolidation (`API-CONSOLIDATE-01`, about 3 to 4 weeks of agent time)

Goal: every function has a purpose, the parts fit, nothing is duplicated, and a new user needs few commands. Each unit is one branch and pull request, has a `PLAN.md` row, and must leave the golden digests unchanged unless its row says otherwise.

### Phase A: safety net
| Unit | State | Acceptance |
|---|---|---|
| `GOLDEN-01` | merged | deterministic outputs of the tiny run pinned |
| `TABLE-SCHEMA-01` | in review | 12 default tables pinned and documented |
| `TABLE-SCHEMA-02` | in review | the same for tables written by optional arms and post-hoc verbs (`molecule_support`, `gene_program_*`, `strand_qc`, evidence tables), each from its own small run; 11 tables pinned under `optional_tables` |

### Phase B: remove what has no purpose (about 2 days)
- `DEAD-01`: delete, or move to tests, the 33 production functions reached only from tests (for example `evidence.extract_viral_reads`, `detection.build_sensitivity_table`, `multimapping.resolve_cb_umi_molecules`, the test shims in `scripts/host_filter.py` and `scripts/anello_align.py`). First check `tests/_frozen_oracles` and the reproducibility repo for importers.
- `CURATION-MOVE-01`: move `hpv_genes`, `repeat_merges` and the maintainer-only gate and reconcile functions (about 600 lines) out of the shipped package; update the ship-scope file and the maintainer scripts.
- `evidence_tiers` and the unwired half of `tier_measures` stay as documented staged code. Wiring them changes calls and needs its own decision (`D-TIER`).

### Phase C: one implementation per helper (about 3 to 4 days)
- `ONLIST-01` (bug fix): one on-list resolver that lets `ChemistryError` propagate. Two arms catch it as `ValueError` today and silently run without barcode correction when `ngs_tools` is missing.
- `BARCODE-01` (maintainer decision): one `-N` GEM-group rule. The six copies agree on `-1` and differ on `-2` and above. Record each copy's output first; switch only if they agree on real inputs; otherwise the maintainer picks the rule (default proposal: `cellcalling._strip_suffix`).
- `SAM-01`: one SAM layer (CIGAR, `AS:i:`, read number, `samtools view` stream) replacing five copies. Merge the two "representative read" rankings only if the golden run shows identical picks.
- `SEQIO-01`: one FASTA reader and one GTF reader replace six and eleven copies; reuse `utils.open_maybe_gzip` and `run_safety.sha256_file`.
- `TOOLCHECK-01`: one host-filter tool check and one STARsolo barcode-argument builder.
- `WHITELIST-FOLD-01`: fold `whitelist_preflight` into `chemistry`.

### Phase D: one data flow through a run (about 4 to 6 days)
- `VALIDATE-ONCE-01`: validate once on `RunConfig` (three passes today: `errorhandler`, a probe `RunConfig`, the real one).
- `CONFIG-KEYS-01`: flags for the four config-only keys (`emptydrops_fdr`, `emptydrops_lower`, `knee_min_umi`, `cell_caller_rscript`) at their current defaults; delete the unreachable fields (`overwrite`, `f1`, `reference`) after confirming no reader.
- `CONFIG-LOAD-01`: one config loader; `with_overrides` for the rerun verbs.
- `WORKER-SHAPE-01`: one worker shape (`run(ctx)`, no module globals). Low-risk workers first; `detection`, `multimap`, `umap`, `gene_programs` last, or documented as deferred. Remove about 50 `getattr` fallbacks on a typed `RunConfig`.
- `ERRORS-01`: a `ViralScanError` base for the 14 exception classes and one catch in `menu._dispatch`; exit codes unchanged.
- `LAYERS-01`: break the function-level import cycles that matter; add `viral_summary.tsv` to `V3Outputs` (11 literal paths today).

### Phase E: a small, consistent CLI (breaking, about 3 to 4 days)
- `RUNDIR-01`: every post-hoc verb takes `--run-dir <run root>` and an optional `--sample`; one resolver; one freshness guard.
- `CMDGROUP-01`: groups `check {chemistry,run,env}`, `discover {run,collapse,states}`, `rerun {multimap,programs,hostresponse}`, `build {ref,star-index,translated-index}`; `host-loss` under a hidden `dev` group. Old names stay for one MINOR with a deprecation pointer.
- `FLAGS-01`: `run --help` shows about 15 essential and common flags and `--help-all` the rest. Remove flags with no effect (`--f1`, `--multimap-primary-call`, `--se-threshold`, `--positive-control-scope virus_key`). Consistent names (`--programme-*`, molecules versus UMI, `--host-filter` versus `--host-subtract`, one called-cells name, one list syntax). Renames are deprecated aliases.

### Phase F: first run in five minutes (about 2 to 3 days)
- `FIRSTRUN-01`: tchem in `environment.yml`, `pixi.toml`, the Dockerfile and `doctor`; quickstart and README use `--kb-params auto`; fix the README `hostresponse` example and the vignette's mixed reference modes.
- `HAPPYPATH-01`: one tested page from install to the report, and a `build-ref` preset for the shipped panel (the bundled-cache wording of the 2026-11-04 decision stays).
- `OPENFIRST-01`: one answer to "what do I open first".
- `COLNAMES-01` (maintainer decision): one column vocabulary and a short core view of the 63-column `viral_summary.tsv`; the rename list is approved as its own row before any edit.
- `API-FACADE-01` (done 2026-11-26, `codex/api-facade-01`): `viralscan.open_run(path)` with `.summary()`, `.cells()` and `.adata()`. The vignettes still import `scripts.detection` and `scripts.cellcalling` for compute helpers (`compute_stats`, `external_cells`, `check_sibling_crossmapping`) that run on synthetic matrices, not on a run directory; moving those behind a stable name is a separate unit.
- `DEMO-01`: `viralscan demo`, a bundled tiny dataset and a one-command check. It is the tester packet for the laboratories.
- Small and independent: a Python 3.13 CI job (or an explicit "unsupported" note); the two `[skip-ci]` vignettes on a schedule in `notebooks-heavy.yml`.

## 5. Beta (`3.0.0b1`) checklist

Agent work: Phases A to F merged and green on `main` including the integration tests; `CHANGELOG`, `docs/migration.md`, `docs/release_policy.md` and `docs/cli_reference.md` consistent; version bump at the six `scripts/check_versions.py` sites (basis: #197); `image-audit.yml` green or exceptions listed; a dry-run tag on TestPyPI; the tester packet.

Owner-only steps:

| Step | Row | Why it blocks |
|---|---|---|
| Authors, order, ORCIDs, contacts, licence holder | `REL-15` | citation order cannot change after the first public tag |
| Protected `testpypi` and `pypi` environments, Trusted Publishers, GHCR permissions | `RC-02` | the publish job cannot run |
| Sign off 8 OSV and 8 grype exceptions | `REL-11` | image audit gate |
| Delete the stale `v3.0.0rc1` tag; push the beta tag on a green `main` SHA | `RC-03` | tag rules |
| Recruit three outside laboratories covering the supported chemistries | `RC-04` | `rc1` gate |

## 6. Laboratory testing and `3.0.0rc1`

- `RC-04`: three laboratories run a clean install, the demo, `validate-run` and one real run each.
- `RC-05`: every correctness, data-loss, installation or documentation failure is tracked to closure and all gates rerun. Laboratory-driven renames land in this window only.
- Freeze: contract constants, schemas, the table registry and the CLI are unchanged from here. `defaults_status` stays `provisional` until the defaults freeze in section 7.

## 7. Science gates (reproducibility repository; the long pole)

These decide whether the defaults are validated.

1. **G3, preregistration.** Merge the SCI-05 round fixes, write the freeze sidecar, then the maintainer's signature and the freeze commit. Open maintainer-owned values: `seeds.generation`, the factor levels, the holdout salt digest. The protocol is `status: draft` today.
2. **G4, references.** Byte rebuild of the panel, the GRCh38 D-list, the calibrated holdout safeguard, the `REF-11` archive with a working DOI (owner), and the reproducibility repo's `pkg/` checkout matched to `PACKAGE_BASE_SHA`.
3. **G5a and G5b.** Truth panel generated and scored; public positives and comparators run on the beta build; uncertainty reported; frozen bundle. An orthogonal anellovirus-positive sample (`REF-10`) comes from a collaborator, or the claim stays screening-only.
4. **Defaults freeze (`DEF-FREEZE-01`).** The maintainer approves each row below from G5 output; none is a silent default. The proposal for each row is the current package default.

| Row | Choice | Current default | Closing evidence |
|---|---|---|---|
| `D-INF` | `infected_min_molecules` | 2.0 | per-cell sensitivity and false-positive rate at 1, 2 and 3 molecules |
| `D-MINCELL` | `detection_min_cells`, `detection_threshold` | 1, 1 | host-only and host-homology negatives |
| `D-MM` | multimap allocation | `host-conservative` | planted ambiguous molecules |
| `D-CC` | `emptydrops_lower`, `emptydrops_fdr` | 100, 0.01 | called-cell sweep |
| `D-PROG` | `programme_productive_min_molecules`, `programme_mixed_veto` | 1, `catalogue` | held-out lytic control |
| `D-TIER` | evidence-tier cutoffs (wire or keep staged) | `uncalibrated` | G5 training negatives |
| `D-HOM` | `REF-08` homology threshold | 95% shared 31-mers | diverged-homologue arm |
| seeds | `seeds.generation`, `emptydrops_seed` | none asked yet | maintainer supplies; no proposed value |

If any default changes, ship a new release candidate on that build and rerun the gates on it. Then `defaults_status` becomes `final` (it is outside the fingerprint, so resume is unaffected).

5. **G6a and G6, documentation and manuscript.** Make `validate_claim_registry.py --coverage` a failing check, put the manuscript in scope, regenerate `docs/manuscript_draft.md` from artifacts, run the independent audits. Today 3 of 10 registry claims are validated and the manuscript is a placeholder.

## 8. Stable `3.0.0`

`STB-01` prove the release-candidate fixes did not change frozen outputs; `STB-02` reserve the software and benchmark DOIs (owner); `STB-03` final metadata, Production/Stable classifier, tag, wheel, sdist, OCI image with `latest`, SIF, locks, checksums, SBOM, attestations; `STB-04` put the PyPI sdist SHA-256 into the Bioconda recipe, lint, build and submit; `STB-05` archive, verify every DOI and link, update `CITATION.cff`, submit the manuscript.

## 9. After 3.0.0

- **Patch line:** `OPS-01` to `OPS-05` (support policy, reference-update cadence, false-positive and false-negative triage, 90-day installation monitoring, truth-panel refresh).
- **3.1 candidates (additive):** an EM-convergence flag column in `viral_summary.tsv` (declined for 3.0); host-barcode overlap reporting once a threshold is defined; 5' chemistry support (needs its own validation; no accuracy claim today); wiring the evidence tiers (needs G5 training negatives); more Python versions and platforms; removal of aliases deprecated in 3.0 (no earlier than 3.1.0).
- **Research backlog:** the appendix below.

## 10. Decisions only the maintainer can make

| Decision | Needed by | Proposal |
|---|---|---|
| `BARCODE-01` rule for `-2` GEM groups | Phase C | keep `_strip_suffix` unless the six copies disagree on real inputs |
| `COLNAMES-01` rename list | Phase F | propose after Phase E; approve as its own row |
| Default rows `D-INF` to `D-HOM`, seeds | G5 | current defaults, decided from G5 output |
| `REL-15` author list and order | before `b1` | owner supplies |
| Pull requests #195, #197, #206 | now | keep #195 and #197; close #206 as superseded |
| Re-pin of the reproducibility repository | after Phase B and C documentation | one re-pin, approved as a commit |

## 11. Critical path, sizing and risks

Agent work in order: immediate fixes (1 day), Phases A and B (about 3 days), C and D (about 8 to 10 days, units in parallel), E (3 to 4 days), F (2 to 3 days), beta preparation (1 to 2 days): about 3 to 4 weeks before `b1` when reviews and CI are quick. External and review time is not estimated: owner steps, laboratory turnaround, the G3 review and signature, the G4 archive, G5 compute and audit. The science gates can start as soon as G3 is signed and do not depend on the consolidation, except that G5 must run on the beta build, which is why the contract and CLI freeze before the laboratories.

Risks:
1. A refactor silently changes a result. Mitigation: the golden digests and the integration suite on every pull request.
2. CI-only failures (ship scope, Python 3.10 and macOS, the `-W` docs build) after a green local run. Mitigation: run the full chunked suite with `PYTHONPATH=src:.` and read CI before merging.
3. Laboratory feedback after `b1` forces more renames. Accepted until `rc1`.
4. G5 contradicts a default. Handled by the new-candidate rule in section 7.
5. Concurrent work in the shared checkout. Use one git worktree per unit.

## 12. How progress is tracked

- Each unit has a `PLAN.md` row, updated in the same commit as the code. The parent rows are `RC2-PLAN-01` (read it as the beta line) and `API-CONSOLIDATE-01`.
- This page carries the roadmap; the gate dashboard in `PLAN.md` stays the authority for G0 to G8.
- Every pull request states its expected output change (none for refactors) and merges only after full CI, including the integration tests, is green.

---

## Appendix: research and feature backlog (2026-07, unchanged)

The items below predate the 3.0 release plan. Their sequencing section is historical and does not bind the release line above.

> Path placeholders: `$COLLAB_ROOT` = collaborator root (institutional filesystem paths, replaced by placeholders).

Future work, grounded in (M) the manuscript's stated limitations
(`docs/manuscript_draft.md` Discussion), (R) code-review findings from this
session, and (T) existing `todo/TODOLIST.md` items. Each item notes rationale,
approach, effort, dependencies, and grounding tag. This is a planning document,
not a commitment — research items are marked exploratory.

Tiers: **P0** publication/release-blocking · **P1** high-value near-term ·
**P2** medium · **P3** long-term / research.

---

### A. Scientific & methods

#### A1 · Read-start distribution + PCR-dup handling in `evidence` — **P2**
Per-position read-start histogram along the viral genome (3′ bias, subgenomic-RNA
junctions, EVE hotspots, coverage evenness). Adds `read_start_distribution()` +
`--read-start-profile`/`--dedup {umi,markdup,none}` to `viralscan evidence`.
Full plan: [`todo/read-start-distribution.md`](../todo/read-start-distribution.md).
*Grounding:* user request; extends the paper method (Chen et al. Nature 2024). *~1–2 d.*

#### A2 · Per-cluster / per-cell EM allocation — **P3 (exploratory)**
The EM estimates **one transcriptome-wide** abundance vector (`em_gene_abundances`,
global pool) and allocates every cell's multimappers with it. In highly-infected
cultures, per-cluster (EM within Leiden clusters) or per-cell (alevin-fry-style)
allocation may be preferable. *Approach:* optional `--multimap-scope
{global,cluster,cell}`; cluster scope reuses the existing scanpy clustering, runs the
existing vectorised EM per group. Watch memory/runtime — global is the memory-cheap
default. *Grounding:* (M) "estimates one transcriptome-wide abundance vector … per-cell
allocation may be preferable." *Large; benchmark against the global default before
shipping.*

#### A3 · Cell-level viral BAM with CB/UB tags — **P2**
`evidence` extracts reads to FASTA then realigns; produce a first-class per-cell BAM
with `CB`/`UB` tags so users can load viral reads in IGV **grouped by cell**. *Approach:*
carry the `(CB,UMI)` from the read name into BAM tags during/after
`align_reads_to_viral`; document an IGV "group by CB" recipe. *Grounding:* (M) "does not
produce cell-level BAM files … genome-browser inspection remains a separate alignment
step." *~1–2 d.*

#### A4 · Combined-genome specificity mode / genome D-list — **P1**
The cDNA-only host reference cannot absorb reads from GRCh38 **non-coding** regions
that resemble viral references → spurious viral signal (the anellovirus / TTV
host-homology artifact: ~90% "viral" under cDNA-only vs 0/4.5M on a combined genome).
*Approach:* official `--genome-dlist genome.fa` in `build-ref` (mask host-homologous
k-mers via kallisto D-list; partially designed per PLAN 2026-07-06) **and/or** a
documented "genomic specificity control" workflow (orthogonal combined-genome STAR/kb
run) surfaced as a first-class QC. *Grounding:* (M) cDNA-only host-homology limitation;
"libraries with abundant non-coding reads … should include an orthogonal alignment."
*Medium; the D-list build is compute-heavy (~64 GB / 8 h).*

#### A5 · Subgenomic-RNA junction detection (coronaviruses) — **P3 (exploratory)**
Once A1 lands, detect canonical leader–body junction pile-ups as a
genuine-transcription signature (vs genomic contamination). *Depends on A1.* *Research.*

---

### B. Validation & benchmarking

#### B1 · Head-to-head vs dedicated viral scRNA tools — **P1**
The manuscript explicitly did **not** benchmark against Venus (Luebbert et al. 2024) or
ViralTrack (Bost et al. 2020) — the STARsolo comparison tests a general splice-aware
aligner only. *Approach:* extend the reference-strategy harness to run Venus/ViralTrack
on the same EBV/HHV-6B/HSV-1 samples; compare sensitivity, specificity, gene attribution,
runtime, memory. *Grounding:* (M) "did not complete a head-to-head comparison against
dedicated viral single-cell tools." *Large; strengthens the paper materially.*

#### B2 · Gold-standard truth panel (planted-read simulation) — **P1**
No ground truth today — validation uses published rate ranges + matched comparisons.
*Approach:* spike known viral reads into a host library at known per-cell rates →
formal false-positive / false-negative rates, ROC, limit-of-detection. *Grounding:*
(M) "Formal FP/FN rates require negative controls and planted viral-read simulations …
left as a release-gated validation extension." *Large; the single biggest credibility
lever for the tool.*

#### B3 · Complete the reference-strategy benchmark — **P1 (T)**
Only 4/12 rows done; all EBV rows failed, viralscan/two-step conditions blocked. Re-run
the failed SLURM conditions so `results/reference_strategy_benchmark.tsv` is analyzable.
*Grounding:* (T) TODOLIST open. *Unblocks B4/B5.*

#### B4 · Selectivity Index + noisy-channel bit-cost — **P2 (T, blocked on B3)**
SI = on-target ÷ off-target viral UMI per (aligner × reference) → ranked reference
recommendation; per-pipeline I(true;measured) in bits (HHV-6A/6B BSC, HSV-1 denominator
artifact). *Grounding:* (T) TODOLIST, blocked on B3.

---

### C. Host-response module

#### C1 · Gene-symbol annotation + pathway enrichment — **P2**
Host-response reports Ensembl feature IDs when the matrix lacks symbols, limiting
biological interpretation of the depth-robust gene set. *Approach:* map Ensembl→symbol
(gget/pybiomart), run enrichment (gseapy/gget enrichr — partially wired via
`viralscan[enrichment]`) on the **depth-controlled** labels, not the raw one. *Grounding:*
(M) "reports Ensembl feature IDs … limiting biological interpretation until annotation
and enrichment are added." *~1–2 d.*

#### C2 · Additional depth-control designs — **P3**
The module already reports a depth-alone baseline, CPM-normalised and depth-matched
labels. Extend with propensity/covariate-adjusted designs and per-gene E-value reporting
surfaced in the HTML report. *Grounding:* (M) depth-confound discussion. *Exploratory.*

---

### D. Robustness / QC / correctness (review findings)

#### D1 · Non-unique `var_names` handling — **P2 (R)**
`matrix_for_genes` uses `var_names.get_indexer`, which raises on duplicate accessions;
the codebase elsewhere anticipates non-unique names (`super_expressor` calls
`var_names_make_unique()`). *Approach:* validate/dedup var_names once at load, or make
`matrix_for_genes` first-occurrence-safe. Not a regression (old code also raised), but a
latent crash on duplicate reference accessions. *Grounding:* (R) primary-call review.

#### D2 · EVE annotate robustness — **P3 (R)**
`_is_chromosome_subject` uses `re.match` anchored at start; old-style `gi|…|ref|NC_…|`
BLAST sseqids downgrade to `subject_not_chromosome` (safe failure, modern BLAST -outfmt6
unaffected). Add a `gi|`/`ref|` unwrap if legacy BLAST DBs are in scope. *Grounding:* (R)
EVE review. *Low.*

#### D3 · `host_viral_ambig_fraction` >1.0 edge — **P3 (R)**
Under some data shapes the ambiguity numerator (separate layer) can exceed the adata.X
denominator → fraction >1.0. Pre-existing; document or clamp. *Grounding:* (R) detection
review. *Low.*

---

### E. Performance & scale

#### E1 · numba the collapse loop — **P3**
The multimap loop is now ~3× faster (collapse) + memory-lean. For very deep samples, a
numba-jitted inner loop over the collapsed rows could give another 5–10×. *Approach:*
ragged CSR-style arrays for `ec_info`, `@njit` the per-record scatter. *Grounding:*
profiling; the loop is the runtime ceiling for 100M-record samples. *Medium; only if
deep-sample runtime becomes a user complaint.*

#### E2 · Chunked/streaming BUS processing — **P3**
`build_multimap_layers` materialises the full `bus_df` (~5 GB RSS on deep samples).
Process the BUS file in chunks (aggregate `(cell,ec)` incrementally) to cap peak memory.
*Grounding:* profiling (peak-RSS on deep samples). *Medium.*

---

### F. Usability / packaging / release

#### F1 · Open PR #7 → main — **P0 (T)**
The `claude/pub-readiness-hygiene` branch (25+ commits) needs its PR merged. Note the
git SSH commit-signing key mismatch (`user.signingkey`) — recent commits may be unsigned.
*Grounding:* (T).

#### F2 · Tag v2.6.0 + publish (PyPI/conda) + Zenodo DOI — **P0**
Version is reconciled to 2.7.0; cut the tag, publish the wheel, mint a Zenodo DOI for the
GTF panel (already on Zenodo) and cite it in the manuscript data-availability statement.
*Needs owner (release action).* *Grounding:* release hygiene.

#### F3 · Functional-script path hygiene — **P1 (T, R)**
~14 SLURM/infra scripts (`covid_viralscan/scripts/*.sh`, `scripts/*.sh`) and
`reference_manifest.json` embed `VS_CONDA_ENV=$COLLAB_ROOT/...` (third-party path)
as a default. Env-var-ize per script or exclude from the shipped artifact — decide ship
scope first. *Grounding:* (T, R) publication-hygiene branch goal. *Medium; per-script,
breaking-risk.*

#### F4 · Stop tracking downloadable/generated data — **P2 (R)**
`ref/10x_version2_whitelist.txt` (3.8 MB, standard 10x file, referenced by 2 SLURM
scripts) and tracked `analysis/**/*.pdf` build artifacts. Fetch-instead-of-track the
whitelist; gitignore generated outputs. *Grounding:* (R) bloat analysis.

#### F5 · CI gates: notebook execution + declared-deps — **P2**
Add an nbmake gate for the 6 CI-runnable vignettes (so `[skip-ci]` docs don't rot — the
enrichment vignette silently broke once) and a `deptry`/bare-venv import check to catch
undeclared deps despite CI's `--no-deps` install. *Grounding:* (R) vignette + CI-no-deps
learnings. *~1 d.*

#### F6 · Docs: IGV recipe + read-start vignette — **P2**
Once A1/A3 land, add a "read-start profile & genome-browser inspection" vignette and an
IGV group-by-CB recipe to the 8-vignette suite. *Depends on A1, A3.*

---

### G. Reference data & annotation

#### G1 · Reference annotation quality — **P2**
The EBV LMP-1/EBNA gene-attribution divergence vs STARsolo traces to **annotation**, not
the aligner (Serratus 96 EBV entries vs STARsolo's 16). Curate/version the viral GTF
panel; surface the reference annotation provenance in outputs (the manuscript argues for
"explicit reporting of viral reference annotations"). *Grounding:* (M) annotation-vs-model
discussion. *Medium.*

#### G2 · Mouse / non-human host support — **P3**
Confirm/extend `build ref --host mouse` (Ensembl) for non-human single-cell viral studies.
*Grounding:* generality. *Medium.*

---

### Suggested sequencing

1. **Ship the branch:** F1 (PR), F3 (path hygiene), F2 (tag/DOI) — publication-blocking.
2. **Strengthen the paper:** B2 (truth panel) and B1 (dedicated-tool benchmark) are the
   two highest-credibility additions; B3 unblocks the reference analyses.
3. **Feature depth:** A1 (read-start, requested) → A3 (cell BAM) → A4 (genome specificity).
4. **Polish:** C1 (gene symbols), F5 (CI gates), D1 (var_names robustness).
5. **Research:** A2 (per-cell EM), A5 (sgRNA junctions), E1/E2 (perf at scale).
