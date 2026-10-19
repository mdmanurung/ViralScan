# Shifted-log normalization qualification

User requirement: use the linked Single-cell best practices procedure, normalize_total(target_sum=None), then natural log1p; no fixed TP10K. User next selected implementation after the requirement was recorded. Baseline6bd1652; isolated codex/shifted-log-normalization.

Four call sites changed in two production modules: host-only PCA, both UMAP branches and raw host-response preprocessing. Inputs/count layers, raw-depth covariates, preprocessed-input handling, thresholds, labels and downstream fitting code remain unchanged. Manifest records the conditional raw-input method, scope and Scanpy version; logs identify the transformation.

Local234regressions pass, including13new cases: independent nonempty-cell formula over dense/CSR and integer/float32/float64, empty-row preservation, both UMAP branches and host-only graph reach, raw label counts, fitted host-response features and immutable input files. Two no-HVG fixtures now use genuinely constant host expression; all skip assertions and thresholds retained. Historical weight golden unchanged, tested with explicitly preprocessed legacy features plus original raw-count layer.

Before source change,11new cases failed on main (red_baseline.log). Initial zero-cell formula assumption failed for CSR: installed Scanpy1.11.5 uses target30 for CSR vs40 dense on cell totals[0,40,20,100]. Exact native API retained; no representation-independent zero-cell median claim. Formula equivalence is asserted on nonempty inputs; zero-cell tests assert empty/finite rows and raw preservation. Initial failure receipts retained.

Paired synthetic fits use identical inputs/paths/seeds, each original and candidate repeated twice. All raw outputs are deterministic within each arm; input file hashes unchanged. A/B changes five outputs, including the manifest and fitted tables. In this bounded random fixture AUC0.5→0.375, weight absolute sum9.158124→9.323273, stable genes8→5. These are expected scientific changes, not improvement evidence; no production calibration or biological validity claim.

Native documented CLI21passed/no skips, covering doctor, kb reference build, counts, visual artifacts, exact read evidence and validate-run. Governance21passed; Ruff/format277files and mypy60source files pass. Full PR CI remains the merge gate.

Public documentation is pinned to its immutable implementation commit. Preserve provenance/shifted-log-2026-10-19 after rebase, together with all three existing perf provenance tags. Raw evidence stays in this research packet and is not shipped. No jobs, data downloads or previous analysis-output rewrites.
