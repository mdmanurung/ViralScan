# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

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
