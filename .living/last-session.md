# Last session — 2026-07-03 (publication-readiness reconcile + execute)

**Goal**: Review, then reconcile-and-execute the agent-authored publication-readiness plan
(`docs/superpowers/plans/2026-07-03-publication-readiness-99.md`).

**What happened**:
- Reviewed the plan (durable review at the plan-file path), then wrote a reconciled version
  (`...-reconciled.md`) after verifying every task against HEAD.
- Executed everything not owner-gated: 7 commits on `claude/multimap-memory-and-showcase`.

**Key changes**:
- **Manuscript §3.4 host-response honesty (integrity gate)** — headline now the tracked
  **0.866**, paired same-design with depth-alone **0.967**, depth-controlled **0.636/0.718**
  (honest ~0.64–0.72). Figure 2 footer caveat + figure regenerated. Fixed the original plan's
  cross-design 0.845↔0.967 pairing error.
- **Docs↔runtime** — multimap default `equal`, 2.5.0 examples, build-ref anellovirus flags,
  `is_called_cell`, `include_anellovirus`, root-help subcommands, CHANGELOG links; new
  `tests/test_docs_consistency.py`.
- **CI/release** — installed-package testing (`pip install --no-deps -e .`), env-file job,
  wheel install-test before publish.
- **Provenance** — tracked reference-strategy run packet + truthful `failure_summary.tsv`;
  benchmark excluded from manuscript claims. MANIFEST keeps the manuscript out of the sdist.
- **Scope** — manuscript narrowed (no dedicated-tool head-to-head; FPR/FNR uncharacterized);
  scaffolds under `analysis/dedicated_tool_comparison/` + `analysis/truth_panel/`.

**Verification**: 562 tests pass (20 deselected); compileall clean; wheel+sdist build; twine
check PASSED; wheel has `emptydrops.R`; sdist excludes manuscript.

**Owner-gated / next**:
- Manuscript authorship + declarations (7 `[to be supplied]` placeholders remain).
- Release: PR→main, tag `v2.5.0`, PyPI/GHCR, bioconda sha256, Zenodo software DOI.
- Dedicated-tool head-to-head benchmark + FPR/FNR truth panel (scaffolds in place).

**Not pushed.** See [[decisions]] (2026-07-03 reconcile entry) and [[learnings]]
(verify-agent-plans-against-head).
