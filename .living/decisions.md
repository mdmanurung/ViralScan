# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

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
