# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

## [2026-07-02] Reference-strategy Selectivity Index blocked — benchmark is only 4/12 complete

**Context**: Next-task attempt at the reference-strategy Selectivity Index (idea 3b) on `results/reference_strategy_benchmark.tsv`.

**Decision**: Did NOT build the SI/channel framework — inspection showed only 4/12 rows complete (all EBV rows failed; every viralscan/two_step blocked; HSV-1 has zero signal). Forcing an SI ranking on 4 heterogeneous rows would over-interpret sparse data. Instead registered the benchmark as an analysis, documented the 2 comparisons the complete rows DO support (HHV-6B: ViralScan ~54% more positive cells than STARsolo, with count-layer + HHV-6A/6B caveats; HSV-1 undetected under either strategy), and flagged the SI (3b) / noisy-channel (7b) todos as BLOCKED pending completion. Added a high-priority "complete the benchmark" todo.

**Consequences**: The reference-strategy thread needs the SLURM benchmark re-run (EBV especially) before quantitative analysis. The HHV-6B aligner observation is n=1 and confounded by count-layer differences — not a firm finding.

**Tags**: reference-strategy, benchmark, incomplete-data, blocked, honesty

## [2026-07-02] Intermediate-attractor hypothesis not supported — host state is a continuum

**Context**: The 4-persona-convergent hypothesis (specificity gap = a discrete primed/intermediate host-cell state) was the top open question. Tested it depth-aware (idea 1a/4a/4b), since the gap is largely a depth artifact.

**Decision**: Ran `intermediate_state_test.py` (UMAP + HDBSCAN + KMeans silhouette + forced k=3, on all cells AND the depth-matched subset). Verdict: **not supported** — HDBSCAN finds no clusters, silhouette ~0.15 (weak) even at matched depth → the host transcriptome is a continuum, not discrete attractors. A naive k=3 split of all cells just stratifies by depth. The only genuinely distinct, depth-independent population is a small high-viral-burden lytic tail (~5.7%, median viral fraction 0.20). Logged as finding F-004; two UMAP figures saved.

**Consequences**: Closes the convergent hypothesis negatively — a good example of a cross-lens-appealing idea that the data (with depth control) do not support. Downstream state-modeling should use a continuum + rare-lytic framing, not tristable attractors. Did not add to the report (the report is about the classifier); kept as an analysis finding. Follow-up: is the lytic tail host-transcriptionally distinct (BZLF1), or only in viral burden?

**Tags**: cell-state, continuum, clustering, ebv, lytic, hypothesis-not-supported, depth-control

## [2026-07-01] CPM cross-check + gene characterization corroborate a real, coherent EBV signal

**Context**: Follow-up to the depth-matched re-analysis — cross-check the ~0.72 with a second de-confounding method and identify the depth-robust genes.

**Decision**: (a) CPM-label cross-check (`cpm_label_crosscheck.py`): a depth-normalized EBV label (per-10k-host-UMI, corr with depth −0.09) drops depth-alone AUC to 0.53 and gives host-gene AUC 0.64. Two independent methods (depth-matched 0.72; CPM 0.64) agree the signal is real; the CPM estimate is slightly lower because the 15 genes were selected on the raw label (partial OOD test) — reported honestly as effect ~0.64–0.72, not a single point. (b) Mapped the 5 depth-robust genes via mygene.info: LTA/ING1/MT-ND4L up, EVI2B/MACROD2 down — biologically coherent with EBV programs (LTA is an NF-κB/LMP1 target; the strongest mechanistic hit). Flagged MT-ND4L's mitochondrial-QC caveat. Formal GO enrichment not run (n=5 too small; gseapy absent) — narrative characterization instead.

**Consequences**: Report (depth section + appendix), F-001 (added CPM refines-evidence + gene identities), analysis doc updated (6pp, 0 drift). Remaining: re-select a depth-robust panel on a depth-independent label then run powered GO; MT-fraction control for MT-ND4L; cross-line generalization.

**Tags**: confounding, cpm-normalization, gene-annotation, ebv, nf-kb, corroboration

## [2026-07-01] Depth-matched re-analysis resolves the EBV signal to a real AUC ~0.72

**Context**: The depth-confounder check showed the AUC 0.866 was depth-inflated. Ran the principled fix (idea 5a follow-up) to get the honest effect size.

**Decision**: Used a depth-matched case-control design (coarsened exact matching on host-depth quantile bins) rather than a CPM-threshold label — it removes the confound by construction and needs no arbitrary new threshold. Result: 1006 cells (503/class), EBV+ vs EBV− depth medians 15,522 vs 15,718 (Mann-Whitney p=0.99), depth-alone AUC 0.48 (chance). At matched depth the host genes still discriminate at **AUC 0.72** (15 genes) / 0.68 (5 depth-robust genes).

**Consequences**: The finding is resolved, not killed — a real, depth-independent host-response signal of AUC ~0.72 exists; the raw-count label inflated it to 0.866. Updated report (abstract + depth section + conclusion), F-001 (contradicted→refined with the matched-set evidence), analysis doc, manifest. This "truth in between" is the honest landing after the review→ideas→confounder-check→re-analysis chain. Cross-validation of ~0.72 with a CPM label + biological characterization of the 5 robust genes remain as open todos.

**Tags**: confounding, sequencing-depth, ebv, case-control-matching, resolution, causal-inference

## [2026-07-01] Depth-confounder check overturns the EBV headline; revised report + F-001 to "contradicted"

**Context**: Ran the high-priority depth-confounder todo (idea 5a) against the AUC 0.866 host-response headline. First verified `_raw_depth` is host-only (the todo's key risk) — it is.

**Decision**: Reported the result honestly as a MAJOR caveat. Evidence: EBV+ rate 27.5%→96.3% across host-depth quintiles; depth-alone AUC 0.803 (all cells) and 0.967 within the headline's balanced+depth-filtered design (beating the 0.866 model); only 5/15 stable genes survive depth adjustment (E≥2). Verified the counterintuitive 0.967 (advisor flag): the top-50% filter widens the between-class depth gap (EBV+ median 34,981 vs EBV− 17,868, 3% overlap), so it amplifies rather than removes the confound.

**Framing (per advisor)**: "substantially depth-confounded and dominated by a depth-driven label," NOT "the finding is fake" — host features are depth-normalized and ~5 genes carry a depth-robust signal. Revised F-001 to `contradicted`, added methodological finding F-003 (raw-count thresholds confound with depth), softened the report title/abstract/conclusion, added a dedicated depth section (report now 6pp).

**Consequences**: Root cause is the ≥10-raw-UMI label; the balancing/depth-filtering does not fix it. Principled fix (named, not run): depth-normalized EBV label (CPM/fraction), depth-matched case-control, or depth in the label definition, then re-estimate on the depth-robust genes.

**Tags**: confounding, sequencing-depth, ebv, review-followup, causal-inference, surprising-result

## [2026-07-01] Ideation session: 4 lenses converge on the specificity gap as a real cell state

**Context**: Phase 6 (`mycelium:ideas`) — 7-persona brainstorm on the EBV host-response finding + reference benchmark; 14 ideas.

**Decision**: Prioritized and promoted 6 ideas to `todo/`. The two high-priority items (depth-confounder E-value check; intermediate-attractor test) were chosen because a convergent hypothesis emerged: **four independent personas (evolutionary biology, stem-cell biology ×2, causal inference) independently flagged the classifier's specificity gap (0.744 < 0.822) as possibly a real "primed/intermediate" cell state rather than noise.** Cross-lens convergence is treated as a stronger signal than any single idea.

**Consequences**: The two high-priority todos are designed to be run together — the intermediate-attractor cluster must be checked against the depth confounder so a cluster isn't just a library-size artifact. Full idea set retained in `analysis/ideas/2026-07-01-cross-disciplinary-brainstorm/`.

**Tags**: ideas, brainstorm, ebv, specificity-gap, convergence, todo

## [2026-07-01] Fixed HVG feature-selection leakage; corrected metrics went UP, not down

**Context**: The mycelium review (F1) found HVG selection was fit on all 1906 cells before the train/test split — feature-selection leakage that, in principle, inflates held-out AUC/MCC. User asked to fix + re-run + regenerate the report.

**Decision**: Moved HVG selection inside the CV fold (`_hvg_mask` on training cells only) via an opt-in `use_hvg=True` path in `_run_l2_regression` (default path kept bit-identical so 30 tests stay green). `run_hostresponse` now passes the full gene matrix for leakage-free per-fold selection; stability selection stays on the HVG subset (descriptive; F6 unchanged). Also fixed F2 (module `detection_threshold` default 1→10), F3 (report "seed 42"→"seed 0"), F5 (MCC rationale), F7 (numbers.json provenance via a committed register script).

**Result (verified, not a bug — but within noise)**: corrected metrics differ from the pre-fix values by **less than one seed-SD on every metric** (AUC 0.845→0.866 vs sd≈0.036; MCC 0.539→0.570 vs sd≈0.077; specificity 0.710→0.744 vs sd≈0.053). Because HVG selection is *unsupervised* (variance-based, never uses labels), this leak class is expected to be negligible — and it was. Do NOT assert a causal mechanism for a sub-SD shift (an earlier draft did; corrected per advisor). The pipeline change is confirmed correct (per-fold path ran, 15-gene stable set identical, train/test partition unchanged, mask/fit see only training cells); the metric change itself is not distinguishable from seed noise (the pipeline retains only mean/sd, not per-seed values, so a paired test isn't available).

**Consequences**: Report headline is AUC 0.866 / MCC 0.570 (the methodologically clean pipeline), framed as "unchanged within noise" vs the leaked version. gene_weights.csv now lists only genes selected in ≥1 fold (~3.8k). My pre-fix "probably lower" prediction was wrong; the honest conclusion is "it was nothing, within noise."

**Tags**: leakage, hvg, cross-validation, mcc, review-fix, surprising-result

## [2026-07-01] Comprehensive report built; scitexlintr unavailable, review consolidated

**Context**: Phase 4 (`mycelium:report`) for hostresponse_ebv_matched. Planning brief: comprehensive shape, Tier B, AUC headline + MCC, baseline = chance + GSE158275.

**Decision**: Built a comprehensive LaTeX report with SciVal-wrapped manifest values, one figure (metrics bar, sha256-fingerprinted), compiled to a 5-page PDF. Framed the quantitative baseline as **chance** (AUC 0.5 / MCC 0) and treated **GSE158275 (SoRelle et al., eLife 2021)** as data source + motivation, NOT a numeric baseline — because that study's barcodes define our positive/negative labels, so it cannot be an independent comparator.

**Deviations from the full report flow** (both forced/pragmatic):
- `scitexlintr` install (external git package) was blocked by the sandbox. Substituted a manual snapshot-vs-manifest drift check (34 SciVal uses, 0 drift) covering scitexlintr's load-bearing snapshot-mismatch rule.
- Phases 4–6 blind reviewers consolidated into a single fresh-context sub-agent (read only .tex + manifest) instead of three sequential loops. It returned 4 findings (lbfgs gloss, changelog framing, MCC chance baseline in abstract, specificity precision); all patched.

**Consequences**: `siunitx` also missing from the TeX install (removed; unused). Specificity precision (0.710 vs 0.71) fixed via a manifest `display` override to stay drift-clean. Full audit trail in `reports/.compile-log.md` + `.manifest.json`.

**Tags**: report, latex, scitexlintr, blind-review, baseline, gse158275

## [2026-07-01] Add MCC by extending the analysis code and re-running (not hand-computing)

**Context**: Phase 4 report planning brief asked to feature MCC, which the analysis did not compute. The report convention forbids typing numbers not produced by the code.

**Decision**: Added `matthews_corrcoef` to `_run_l2_regression` in `src/viralscan/scripts/hostresponse.py`, re-ran the identical pipeline (fixed seeds) on the saved matched h5ads in the `bioenv` conda env, validated reproduction (thresholded metrics bit-identical; AUC to 3 dp; stable gene set identical), then `register_value`'d MCC 0.539±0.059. Adopted the re-run outputs wholesale for internal consistency.

**Alternatives considered**:
- Approximate MCC from aggregate sensitivity/specificity + prevalence — rejected: Jensen's-inequality bias vs the per-seed mean, and it's a non-code-grounded number.
- Re-run the full wrapper from the external run_dir — rejected: that dir is gone; `write_hostresponse_summary` regenerates the summary from the output dir alone.

**Consequences**: `hostresponse.py` gains one metric; all small outputs regenerated from one run. AUC headline moved 0.84456→0.84487 (unchanged at 0.845). 30 module tests pass. The external showcase run_dir is no longer available, so full re-preparation from FASTQ is not reproducible locally.

**Tags**: analyze, mcc, reproducibility, register-value, report

## [2026-07-01] Register (not re-run) the EBV host-response analysis; waive scilintr FPs

**Context**: Phase 3 (`mycelium:analyze`) on an analysis (`hostresponse_ebv_matched`) that was already run with complete outputs (AUC 0.845, CV mean±sd). User asked for the quicker path.

**Decision**: Documented and registered the existing analysis under `analysis/hostresponse_ebv_matched/` (doc + ANALYSIS_MANIFEST + `register_value` for 10 headline numbers) rather than re-running heavy compute. Ran scilintr on the script; resolved its 2 `unchecked-cache` findings with structured waivers (they are input-location resolution in `_h5ad_path`, not output caching — a false positive), not behavioral edits.

**Alternatives considered**:
- Re-run the full pipeline under robust-analysis — rejected: slow, needs the external showcase run dir; outputs already exist.
- Behaviorally "fix" the scilintr cache findings — rejected: would change a correct input-path helper; a waiver is the right call.

**Rationale**: Fastest correct registration; preserves reproducibility (script + numbers.json) without altering validated results.

**Consequences**: Re-introduced `analysis/` (now non-empty, so it earns its place). `register_value` `computed_at` shows `<stdin>` because it was invoked via a heredoc, not the script itself — cosmetic. A threshold sensitivity sweep remains a recommended follow-up (added to todo).

**Tags**: analyze, scilintr, register-value, ebv, robust-analysis

## [2026-07-01] Register reference set in place rather than moving into data/raw/

**Context**: The reference-strategy reference set is ~64 GB, immutable, already built in place (`references/` + out-of-repo archive paths), and already has machine-readable provenance (`reference_manifest.json`) and a SHA256 audit (`reference_audit.tsv`).

**Decision**: Ingested it as dataset `reference_strategy_refs` **in place** — `data/raw/reference_strategy_refs/` holds only a pointer/large-file doc; `data/metadata/reference_strategy_refs/` holds schema/provenance/summary that cite the existing manifest + audit as the authoritative checksum source. Gitignored the bulk data. Committed `reference_manifest.json` and `reference_audit.tsv` as the provenance records.

**Alternatives considered**:
- Copy/move data into `data/raw/` — rejected: 64 GB, would duplicate and break in-place aligner paths.
- Recompute SHA256s into `provenance.md` — rejected: `reference_audit.tsv` already has them; duplication invites drift.

**Rationale**: Honors mycelium's large-file convention (gitignore + document + checksums) while reusing existing prior art (manifest, audit) instead of regenerating it.

**Consequences**: DATA_MANIFEST `raw_path` points to a doc, not data. Anyone rebuilding must follow `build_commands` in the manifest and verify against the audit. If either filesystem root moves, manifest paths break.

**Tags**: ingest, large-files, references, provenance, dry

## [2026-07-01] Minimal mycelium scaffold; PLAN.md stays authoritative

**Context**: ViralScan is a mature repo with its own layout (`src/`, `scripts/`+`results/`, `references/`) and a strict, pre-existing `PLAN.md` contract for tracking in-flight work. Adopting mycelium risked (a) cluttering the tree with unused standard scaffold dirs and (b) creating a second, competing source of truth for the plan.

**Decision**: Ran `mycelium:core init` with a **minimal** scaffold — pruned the standard `algorithms/`, `analysis/`, `reference_material/` top-level dirs (kept only `data/`, `todo/`, `skillpacks/`, `.living/`). Installed core packs (robust-analysis, report-generator, idea-generator) + domain packs bioinformatics and skill-bridge. Documented in `CLAUDE.md` that **`PLAN.md` remains authoritative**; `.living/` is additive context only.

**Alternatives considered**:
- Full scaffold as-is — rejected: adds empty dirs duplicating existing structure.
- Map mycelium manifests onto existing dirs — rejected for now: more upfront work; can revisit.

**Rationale**: Keeps the tree clean and avoids a two-tracker conflict while still enabling the `.living/` memory layer, convention packs, and enforcement hooks.

**Consequences**: `validate_structure.py` reports 6 errors for the pruned dirs — this is **intentional and expected**, not a real failure. The SessionStart health hook does NOT check those dirs, so no per-session nagging. `.claude/settings.local.json` hooks were kept user-local (absolute paths), not committed.

**Tags**: mycelium, repo-structure, tooling, plan-md, scaffold
