---
description: >
  Convert a ViralScan scientific audit into an ordered TDD implementation plan.
  For each audit finding, writes a failing pytest first (red), then implements
  the fix (green), then updates PLAN.md. Output: new PLAN.md tasks + test stubs
  committed together. Use after running viralscan-scientific-audit.prompt.md.
name: "Audit → TDD Implementation Plan"
argument-hint: "Path to audit file, e.g. audits/2026-05-08-full-pipeline.md"
agent: "agent"
model: ['Claude Sonnet 4.5 (copilot)', 'Claude Opus 4.5 (copilot)']
tools: [codebase, search, usages, findTestFiles, problems, runCommands, terminalLastCommand, terminalSelection]
---

You are converting a ViralScan scientific audit report into a **test-driven
implementation plan**. You work in strict Red → Green → PLAN.md order:
write a failing test first, implement the fix, mark the task done.

## Inputs

- **Audit file** (argument, or discover the most recent `audits/*.md`).
- **PLAN.md** — the operational tracker you must keep current.
- **Source tree** under `src/viralscan/` and tests under `tests/`.
- **Test command:** `PYTHONPATH=src python -m pytest tests/ -v`

## Operating Rules

1. **Never** skip the failing-test step. A fix without a prior failing test is
   not TDD and must be treated as incomplete.
2. **Tests must validate scientific correctness**, not just execution. Distinguish:
   - ✅ `assert corrected_adata.X.sum() == raw_adata.X.sum()` (count conservation)
   - ❌ `assert corrected_adata is not None` (execution check only)
3. Follow the repo style: `PYTHONPATH=src`, `pytest.mark.network` for NCBI calls,
   `tmp_path` fixtures, no `shell=True`, `sys.exit` not bare `exit()`.
4. Update `PLAN.md` after every task: flip `[ ]` → `[x]`, update "Next up".
5. Commit message format: `fix(<module>): <one-line description> (audit <section>)`.

---

## Phase 0 — Triage and Ordering

Read the audit file fully. Extract every finding from:
- **Section 2** (High-Risk / Definite Errors)
- **Section 3** (Medium-Risk / Likely Problems)
- **Section 4** (Unclear Assumptions — include only those with a concrete
  verifiable behavior, skip pure domain-review items)

For each finding, record:
```
ID       : <audit section number, e.g. 2.1>
Severity : High | Medium | Low
Module   : <file:function>
Claim    : <one-sentence scientific or correctness claim to verify>
Test kind: unit | integration | property
```

**Ordering rules (apply in order):**
1. Any finding where existing code provably produces wrong numbers → first.
2. Reproducibility breaks (missing random seeds) → second, because they block
   validation of all other fixes.
3. Silent data-corruption paths (barcode mangling, cache reuse) → third.
4. Validation gaps (threshold ranges, GTF format checks) → fourth.
5. Coverage gaps (no integration test) → last.

Print the ordered triage table before starting any implementation.

---

## Phase 1 — For Each Task, in Order

Repeat the following cycle for every finding in the triage table.

### Step 1 · Understand the current behavior

- Read the relevant source file(s) and the existing test file(s).
- Identify the *exact* line(s) the audit flagged.
- Confirm the bug/gap is still present in the current codebase before proceeding.
  If the audit was wrong or the issue is already fixed, mark it `[x]` in PLAN.md
  with note "Pre-existing fix confirmed YYYY-MM-DD" and move on.

### Step 2 · Write the failing test (RED)

Create or extend a file in `tests/` following the pattern:
- `test_<module>.py` for unit tests.
- `tests/integration/test_pipeline_<scope>.py` for integration tests.
- `tests/property/test_<module>_properties.py` for property/conservation tests.

**Test template (adapt as needed):**

```python
# tests/test_<module>.py  — added for audit finding <ID>
import pytest
# ... imports ...

class Test<FindingName>:
    """Audit <ID>: <one-line description of the scientific claim>."""

    def test_<specific_behavior>(self, <fixtures>):
        """
        GIVEN: <setup — synthetic input that exposes the bug>
        WHEN:  <the function under test is called>
        THEN:  <the scientifically correct assertion>

        Regression for: audits/<audit-file>.md §<ID>
        """
        # Arrange — create minimal synthetic data that exposes the issue
        ...

        # Act
        result = function_under_test(...)

        # Assert — scientific claim, not just execution
        assert ..., "<human-readable failure message>"
```

**Required assertion patterns by finding type:**

| Finding type                    | Required assertion                                          |
|---------------------------------|-------------------------------------------------------------|
| UMI count conservation          | `total_after == total_before` (float tolerance ≤ 1e-6)     |
| Double-count bug                | `result.X.sum() < raw.X.sum() * 1.01` (not inflated)       |
| Detection threshold FP          | known-zero gene not in `found_genes` at threshold=1         |
| GTF gene_id extraction          | correct accession extracted from adversarial GTF line       |
| Barcode suffix stripping        | barcode with mid-string "-1" unchanged except trailing      |
| Cache content validation        | truncated cached file triggers re-download                  |
| Random seed reproducibility     | two runs with same seed produce identical UMAP coordinates  |
| Threshold boundary              | threshold=0 raises `ValueError`; threshold=1 accepted       |

Run the test and **confirm it fails** before proceeding:
```
PYTHONPATH=src python -m pytest tests/test_<module>.py::Test<Finding>::test_<behavior> -v
```
Paste the failure output into your working notes. If the test passes already,
the bug was already fixed — update PLAN.md accordingly and skip to the next task.

### Step 3 · Implement the fix (GREEN)

Make the **minimal change** that makes the failing test pass.

- Do not refactor unrelated code.
- Do not add docstrings or comments to code you did not change.
- Prefer `pathlib.Path` over string paths; list-form `subprocess.run`; no `shell=True`.
- If the fix requires a config-schema change, also update `createconfig.py` validation
  and add a test in `test_createconfig.py`.

Run the full suite after each fix:
```
PYTHONPATH=src python -m pytest tests/ -v
```
All previously passing tests must still pass. If any regress, fix them before moving on.

### Step 4 · Update PLAN.md (DONE)

Add or update the task block in `PLAN.md`:

```markdown
### Task <N> — <Short title>  `[x]`

**Audit finding:** <audit file>§<ID>
**Severity:** High | Medium | Low
**Fix:** <one-sentence description of the change>
**Test:** `tests/test_<module>.py::Test<Finding>::test_<behavior>`
**Regression:** <what property the test permanently guards>
```

Update the "Next up" pointer at the top of `PLAN.md` to the next open task.

### Step 5 · Commit

```
git add <changed source files> <changed test files> PLAN.md
git commit -m "fix(<module>): <description> (audit §<ID>)"
```

Do not batch multiple findings into one commit unless they are
inseparable (e.g., a shared helper used by two fixes).

---

## Phase 2 — Integration Test Skeleton

After all individual fixes are done, add one integration smoke-test that runs
the full multimapping → detection → umap chain on **synthetic minimal data**
(no real FASTQs needed; use pre-built `.h5ad` fixtures):

```python
# tests/integration/test_pipeline_counts.py
class TestEndToEndCountConservation:
    """Full pipeline: UMI counts must be conserved from multimap → detection → umap."""

    def test_total_umi_conservation(self, minimal_h5ad_fixture):
        """
        GIVEN a synthetic h5ad with known total UMI count T,
        WHEN  the multimap correction, detection, and umap steps are applied,
        THEN  adata.X.sum() at each stage never exceeds T * (1 + floating_tolerance).
        """
        ...
```

This test must be runnable in CI without network access or real FASTQs.
Use `pytest.mark.integration` and add that marker to `pyproject.toml`.

---

## Phase 3 — Final PLAN.md Summary

After all tasks and the integration test are committed, append a summary block
to `PLAN.md`:

```markdown
---
## Audit remediation — 2026-05-08 full-pipeline (completed YYYY-MM-DD)

| Finding | Severity | Status | Test |
|---------|----------|--------|------|
| §2.1 detection threshold / no FDR | High   | [x] | `test_detection.py::TestDetectionThreshold` |
| §2.2 GTF gene_id parsing           | High   | [x] | `test_analysis.py::TestObtainGtf`          |
| §2.3 barcode suffix stripping      | High   | [x] | `test_multimap.py::TestLoadBarcodes`       |
| §3.1 missing random seeds          | Medium | [x] | `test_umap.py::TestReproducibility`        |
| §3.2 cache content validation      | Medium | [x] | `test_ncbi_fetch.py::TestCacheValidation`  |
| §3.3 threshold not positive        | Low    | [x] | `test_createconfig.py::TestThresholdValidation` |
```

Fill in actual task numbers and test identifiers from the work done.

---

## Constraints

- **No** new external dependencies unless strictly necessary and approved.
- **No** changes to public CLI flags or config key names without a note in PLAN.md.
- **No** `toarray()` calls on matrices larger than a defined safe size constant
  without a guard (see `constants.py`).
- Boolean config values must be native Python `bool`, never the string `"True"`.
- All new code must pass `ruff check` (configured in `pyproject.toml`).

---

## Definition of Done

The audit-to-plan conversion is complete when:

1. Every High-Risk and Medium-Risk finding has a named, passing pytest test that
   verifies the scientific claim (not just execution).
2. Every Low-Risk finding has at minimum a `pytest.warns` or `pytest.raises` guard.
3. `PYTHONPATH=src python -m pytest tests/ -v` passes with zero failures and no
   new skips compared to the baseline.
4. `PLAN.md` is updated with all new task blocks marked `[x]` and the remediation
   summary table.
5. All changes are on branch `claude/review-repo-improvements-Sg4Th`.
