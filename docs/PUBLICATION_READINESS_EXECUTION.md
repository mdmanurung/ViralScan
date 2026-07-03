# Publication Readiness Execution Log

Started: 2026-07-03
Branch: claude/multimap-memory-and-showcase
Plan: docs/superpowers/plans/2026-07-03-publication-readiness-99.md
Reconciled: docs/superpowers/plans/2026-07-03-publication-readiness-99-reconciled.md

## Gates

- [x] Wheel includes every runtime resource. (`pyproject.toml` package-data already had
      `emptydrops.R`; local `pip install --no-deps -e .` confirmed it is importable next to
      `cellcalling.py`; release workflow now asserts it.)
- [x] CI tests installed package, not only PYTHONPATH. (`ci.yml` `test`/`integration` now
      `pip install --no-deps -e .`, drop `PYTHONPATH`, smoke-test the console script.)
- [x] Release workflow install-tests wheel before publish. (`release.yml` build job.)
- [x] Public docs agree with runtime defaults. (`equal` default; 2.5.0 examples;
      `tests/test_docs_consistency.py` passes.)
- [x] Manuscript host-response AUROC is deconfounded / explicitly labelled confounded.
      (§3.4 reports 0.866 headline + depth-alone 0.967 + depth-controlled 0.636/0.718;
      Figure 2 footer caveat; figure regenerated.)
- [x] Incomplete benchmark claims are removed or benchmark is complete. (Reference-strategy
      benchmark excluded from manuscript claims; provenance tracked in a run packet.)
- [x] Clean clone can regenerate every cited lightweight result. (Benchmark helper scripts +
      run-packet provenance tracked.)

## Final verification (2026-07-03)

- `python -m compileall -q src tests scripts`: PASS
- `PYTHONPATH=src pytest tests/ -q`: **562 passed, 20 deselected** (6m22s)
- `pip install --no-deps -e .` (fresh venv): PASS — console script + `emptydrops.R` importable
- `python -m build`: PASS — wheel `viralscan-2.5.0-py3-none-any.whl` + sdist
- `twine check dist/*`: **PASSED** (wheel + sdist)
- Wheel contains `viralscan/scripts/emptydrops.R`: PASS
- sdist ships CITATION/CHANGELOG and **excludes** the manuscript: PASS
- Manuscript stale `AUROC 0.845` removed; honest 0.866/0.967/0.636/0.718 present: PASS
- `.github/workflows/{ci,release}.yml` parse as YAML: PASS

## Owner-gated remainder (not done here)

- Manuscript authorship/affiliations/lead-contact/acknowledgments/contributions/declarations
  (7 `[... to be supplied]` placeholders remain by design).
- Release: PR→main, tag `v2.5.0`, PyPI/GHCR publish, bioconda sha256, Zenodo software DOI.
- Dedicated-tool head-to-head benchmark and FPR/FNR truth panel execution (scaffolds tracked).
