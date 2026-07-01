# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

## [2026-07-01] Minimal mycelium scaffold; PLAN.md stays authoritative

**Context**: ViralScan is a mature repo with its own layout (`src/`, `scripts/`+`results/`, `references/`) and a strict, pre-existing `PLAN.md` contract for tracking in-flight work. Adopting mycelium risked (a) cluttering the tree with unused standard scaffold dirs and (b) creating a second, competing source of truth for the plan.

**Decision**: Ran `mycelium:core init` with a **minimal** scaffold — pruned the standard `algorithms/`, `analysis/`, `reference_material/` top-level dirs (kept only `data/`, `todo/`, `skillpacks/`, `.living/`). Installed core packs (robust-analysis, report-generator, idea-generator) + domain packs bioinformatics and skill-bridge. Documented in `CLAUDE.md` that **`PLAN.md` remains authoritative**; `.living/` is additive context only.

**Alternatives considered**:
- Full scaffold as-is — rejected: adds empty dirs duplicating existing structure.
- Map mycelium manifests onto existing dirs — rejected for now: more upfront work; can revisit.

**Rationale**: Keeps the tree clean and avoids a two-tracker conflict while still enabling the `.living/` memory layer, convention packs, and enforcement hooks.

**Consequences**: `validate_structure.py` reports 6 errors for the pruned dirs — this is **intentional and expected**, not a real failure. The SessionStart health hook does NOT check those dirs, so no per-session nagging. `.claude/settings.local.json` hooks were kept user-local (absolute paths), not committed.

**Tags**: mycelium, repo-structure, tooling, plan-md, scaffold
