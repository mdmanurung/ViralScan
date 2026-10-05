# Package completion plan: finish the software before more experiments (2026-10-04)

Status: **proposed**. User directive (2026-10-04): prioritise finishing
package development over new experiments. This plan supersedes the ordering in
`2026-10-04-anello-remaining-work.md`:
- its A1 (REL-16) and A3 (evidence measures) are absorbed here;
- its A2 (census on counted reads), Track B (`ANDET-09f`) and Track C
  experiments are **deferred**.

## Definition of "package finished"

The tracker's own gates:
- **G0**: governance.
- **G1**: software. "All count invariants, schema checks, safety scenarios,
  rerun consistency and the full tiny workflow are green. No known correctness
  or data-loss defect may remain."
- **G2**: distribution. "Clean lock-created env plus Docker and Apptainer run
  the same tiny workflow with identical validated counts; installable
  schemas/assets; security/licence reports clean; exact-SHA release workflow
  dry-run verified."
- **G6a**: docs, minus the numeric claims, which need G5 by design.

Milestone tag: `3.0.0.devN` after G1 (`DEF-09`), with `defaults_status =
provisional` in `run_manifest.json`. That stamp lets the package ship while
experiment-dependent defaults stay open: DEF-02 strand, ANDET-07 panel default,
ANDET-09 on/off.

`3.0.0rc1` stays gated on G5 (science) as the tracker defines it. A
software-only RC would need the gate definitions changed (decision D3).

## Baseline measured 2026-10-04 (HEAD c8bc8c7)

| Gate command | State |
|---|---|
| `pytest tests/ -q` | **green**: 1,681 passed, 113 deselected |
| `ruff check .` | **red**: 18 errors. 11 in untracked `extras/vendor_sources/`, 2 in `extras/`, 5 in `src/` + `tests/` |
| `ruff format --check .` | **red**: 26 files. ruff 0.16 now also formats Markdown under `.living/`, `audits/`, `todo/`; about 12 are tracked `src/`/`tests/`/`extras/` files |
| `pytest -m "integration and not network"` | **hangs, plus 1 fail and 2 skips** (details below) |
| `check_data_governance.py` | green |
| G0 archive build | blocked: no PyPA `build` frontend (old "Do now") |
| Sphinx (G6a) | not installed in any env |
| Docker | no daemon on the cluster |
| Apptainer | available as module `container/apptainer/1.3.3` |
| Lock tools | `uv`, `pixi` in `~/.local/bin`; an untracked `uv.lock` in the repo root |

Integration details:
- `tests/integration/test_tiny_end_to_end.py:108` builds its index with PATH
  (conda) `kallisto index`, then `kb count` runs kb's bundled kallisto on it.
  That spins at about 760 % CPU indefinitely (killed after 20 min). This is the
  REL-16 version-lock trap, so **REL-16 is a G1 blocker, not only an evidence
  bug**.
- `test_anellovirus_chain.py::test_build_anellovirus_reference_produces_labelable_gtf`
  **fails**. The build's low-complexity gate rejects NC_002076.2 (74/74
  k-mers) because `dustmasker` is absent, and the test does not skip when it
  is missing.
- `test_evidence_chain.py:86` and `test_exact_lineage.py:94` **skip**: no
  minimap2, blastn or makeblastdb in `viralscan_bench`.

---

## M0 — Reconcile the tracker (≈ 2 h, one docs commit)
Stale or duplicate rows found by the 2026-10-04 survey:
- CAT-05 is `[ ]` at line 2856 and `[x]` at line 3285: merge them.
- CAT-27 (EBER audit) is done in practice (CAT-32 outcome, line 3646): close it
  with evidence.
- MECH-D and DEF-02 describe the same chemistry module: keep one row.
- SW-22 folds into MECH-C.
- Line ~1640: the DEF-02 sub-row says TTMDV is "unvalidated", but EXPL-TTMDV
  is `[x]`.
- ANDET-09f's "09a still needs the cat42b index build" is stale (built Oct 3).
- ANDET-06 contradicts REF-01's `[x]` "now default": fix REF-01's wording.
- CAT-15 cites old PROG-07 numbers.
- CAT-12's SARS-CoV-2/influenza bullets are superseded (CAT-21/22); CAT-20's
  ICTV refresh was retracted.
- Two "Blocked" notes are under the wrong heading (lines 1970 and 2180).
- Closed rows with open residue (GOV-06, SW-10, SW-20, MECH-A): move the
  residue to its owner row.
- Replace the stale "Do now" block (line 340) with a pointer to this plan.

**Exit:** `PLAN.md` has no row contradicting another; the "Next action" pointer
names M1.

## M1 — G1 green, then tag `3.0.0.devN` (≈ 3–5 days)
1. **Lint (≈ 30 min).**
   - Add `extend-exclude` for `.living/`, `audits/`, `todo/`,
     `extras/vendor_sources/` and `*.md`. These are notes and vendored code,
     not package source.
   - `ruff format` the tracked code, and fix the 7 real findings
     (`anello_align.py` UP035/UP031, three tests, `extras/cat09_*`).
   - Exit: `ruff check . && ruff format --check .`.
2. **`REL-16`: resolve kallisto/bustools through kb_python (≈ 3 h).**
   - One resolver built on `kb_python.config.get_kallisto_binary_path()` /
     `get_bustools_binary_path()`.
   - Route `evidence.py:392-437` and `scripts/multimap.py:142-174` through it.
     The multimap bustools runs on **every default run**; the conda bustools
     differs by size, at the same version.
   - The integration tests build indexes with the **resolved** kallisto
     (`test_tiny_end_to_end.py:108`, `test_exact_lineage.py:33`).
   - Record path and md5 in `evidence_manifest.json`, and in
     `run_manifest.json` under a key **excluded from `run_fingerprint`**
     (`run_safety.py:110`), so `--resume` survives.
   - Exit: `test_tiny_end_to_end` completes; a `--resume` fingerprint test
     passes.
3. **A test environment with the full tool set (≈ 1 h, needs D1).**
   - minimap2, BLAST+ (dustmasker, blastn, makeblastdb) and R (for the
     EmptyDrops paths) from bioconda/conda-forge.
   - `test_anellovirus_chain` gains the same `have_tools` skip as its
     siblings, so a missing tool skips instead of failing.
   - Exit: the integration suite runs with 0 skips for tool reasons.
4. **Correctness-defect triage (≈ ½ day to triage, then 1–3 days of fixes).**
   - Each candidate is checked against the code before it becomes a G1 blocker.
   - Confirmed defects get fixed in M1; features and refactors move to M4.

   | Row | Candidate defect | Check |
   |---|---|---|
   | REF-13 | duplicate `NC_002076.2` header crashed the July evidence job at samtools | do the current build paths emit it? |
   | MECH-C (+SW-13, SW-22) | `rerun-multimap` rewrites `config.yaml`, so mtime triggers a full rerun; the config round-trip loses types | reproduce with the tiny fixture |
   | CAT-37 | GTF cache has no format-version stamp, so stale GTFs are reused | read the cache-key code |
   | CAT-41, ANDET-05 | wrong display names / two genus names per genome in outputs | grep the catalogue |
   | PROG-14 | unresolved programme markers pass silently | read the resolver |
   | ANELLO-14 | integration test expects old gene IDs | run it |
   | DEF-06 | counting-contract property tests for within-sibling allocation are missing | test gap, not a defect |
   | MECH-F | off-list barcode drops not counted in the audit | today's TSO audit shows the drop is real and silent |
   | SW-07 | host–virus-ambiguous boundary count undefined | scope check |
5. **Gate run and tag (≈ 1 h).**
   - Run all four G1 commands green on one SHA.
   - `DEF-09`: tag `3.0.0.devN`, stamp `defaults_status=provisional`.
   - Record the evidence in PLAN.
   - Pushing the tag is the user's call (D2).

**Exit:** the G1 block (PLAN ~line 816) passes in full; `G1` flips to `[x]` in
the gate dashboard.

## M2 — G2 distribution (≈ 1–2 weeks, part of it CI-bound)
1. **Dependency tiers (`REL-01`, `REL-02`).**
   - Snakemake and the workflow dependencies move out of the mandatory pip
     dependencies.
   - Python 3.11 linux-64 is declared the canonical full-workflow toolchain.
   - `doctor --profile pip` narrows to match.
2. **Locks (`REL-03`, `CAT-25`, `DEF-08`).**
   - Python lock with `uv`.
   - External-tool lock (kb, kallisto, bustools, STAR, samtools, minimap2,
     BLAST+, snakemake) with `pixi` or conda-lock.
   - An exact tool-version manifest; a snakemake 9 dry-run in CI.
   - Needs D1.
3. **Archives and clean installs (`G0` close-out, `REL-04`, `REL-05`,
   `CAT-13`).**
   - `python -m build` with the prescribed member check.
   - A wheel/sdist size gate.
   - `doctor --profile pip` and `validate-run` on a packaged v3 fixture, from
     both clean installs.
4. **Containers and parity (`REL-06`–`REL-10`).**
   - Docker built in CI from the lock and the tested wheel.
   - OCI digest, then an Apptainer SIF built on the cluster from that digest.
   - The same tiny workflow in conda, Docker and Apptainer must give identical
     counts and hashes, written to `analysis/release_parity/parity_report.json`.
   - Needs D2.
5. **Supply chain and release dry-run (`REL-11`–`REL-14`).**
   - Audit the locked env; OCI scan, SBOM, licence report.
   - SHA-pinned Actions; attestations; protected environments.
   - An exact-SHA release workflow dry-run without publishing.
   - Needs D2 and repo-admin access; `REL-15` is the user's metadata.

**Exit:** the G2 gate text (PLAN ~line 1896) is satisfied and the parity report
is committed.

## M3 — G6a docs, non-numeric (≈ 3–4 days, parallel with M2)
1. A Sphinx env (needs D1). `sphinx -W` and linkcheck green.
2. `DOC-01`/`DOC-04`: pip-tier vs full-tier install text; data-fetch wording.
   `ANDET-06`: the REF-01 default claim, the host_filter STAR-defaults claim, and
   the cli_reference ranking.
3. `DOC-07`: the 6 CI notebooks run in CI; fix the 3 that fail on column drift;
   schedule the full notebooks.
4. `DOC-09`: claim-registry schema, validator, stale-hash detection, coverage
   check. Parser/default parity must be exact.
5. `OPS-01`/`OPS-02`: support, deprecation and security policy; the
   reference-update cadence.

Out of scope until G5: `DOC-08`, `DOC-10` (numeric claims). `DOC-05` needs
`REF-11` (Zenodo DOI, external).

**Exit:** the G6a command block passes except the claim-coverage items that
depend on G5.

## M4 — Package features and refactors (scope decision D4)
Not G1 defects. Recommendation: before rc, land only what changes the **output
schema** (breaking it after rc costs users) plus the anellovirus evidence
columns.
- **Recommended before rc:**
  - `ANELLO-PRIOR.3` evidence wiring (`complex_body_fraction`,
    `reagent_fraction` in `alignment_qc.tsv`);
  - `ANDET-03` (`claim_scope` column);
  - `ANDET-01` (accession breadth);
  - `CAT-18` (tier column, if two-index is chosen);
  - `MECH-B` (virus-level detection; changes `viral_summary` semantics).
- **After rc or never:**
  - `SW-09` (CLI split), `MECH-E` (reference pipeline), `SW-03` (assignment
    evidence file);
  - `CAT-19`, `CAT-22`, `PROG-12`/`PROG-13`, `DEF-01`, `DEF-04`, `DEF-05`,
    `TONSIL-02`;
  - `CAT-04`, `CAT-06`, `CAT-16`, `CAT-28`, `CAT-29`, `CAT-36`, `CAT-38`,
    `CAT-39`.

## Deferred — experiments (do not start until M1–M3 close)
- `ANDET-09f`, plus the plant regeneration and filter sweep.
- The census on the counted reads (A2 of the anello plan).
- `ANELLO-PRIOR.1`, `ANELLO-PRIOR.2`.
- Everything classed EXPERIMENT in WP3–WP7 and WP9: `SCI-03..05`, `REF-06`,
  `REF-08`, `REF-09`, `HOST-05`, `SENS-06`, `HPV-11`, `ANELLO-12`/`-13`,
  `ANDET-08`, `CAT-07`/`-08`/`-12`/`-14`/`-15`/`-23`/`-26`/`-40`,
  `PROG-15`/`-16`/`-18`, `VAL-04`/`-05`/`-08..10`, `RUN-01..05`, `CMP-00..06`,
  `RES-01..05`, `MS-01..05`.
- `VAL-01..03`, `VAL-06`, `VAL-07` (truth-panel generator and scorer) are code,
  but they exist only for G5. Deferred by default; the user can pull them
  forward.

## Decisions needed (they gate execution, not this plan)
- **D1 — network installs.** The standing "no PyPI installs" rule blocks:
  - G0's `build` frontend;
  - Sphinx;
  - the test tool env (minimap2, BLAST+, R);
  - the lock resolution.

  Ask: one approval for dedicated envs (`viralscan_test_full`, `viralscan_docs`,
  `viralscan_build`) from conda-forge/bioconda/PyPI. The alternative is offline
  artifacts supplied by the user.
- **D2 — push, CI and registries.**
  - Docker parity and the release dry-run only run in GitHub Actions.
  - GHCR and protected environments need repo-admin.
  - Pushing is the user's call.
- **D3 — release target.** Confirm that "package finished" = G0 + G1 (devN tag)
  + G2 + G6a non-numeric, with rc1 still gated on G5. The alternative is a
  software-only rc1, which needs the gate definitions amended.
- **D4 — M4 scope.** Schema-changing features before rc (recommended), or all
  of M4, or none.
- **External (unchanged):**
  - `REL-15` (authors, licence holder, ORCID, …);
  - `REF-11` (Zenodo DOI 10.5281/zenodo.20112332 currently returns 404; it
    blocks REL-05, DOC-05, G2 and GOV-06);
  - `RC-02`;
  - `scancel 25651950`.

## Order
| Week | Work |
|---|---|
| 1 | M0 → M1.1–M1.3 → M1.4 triage and fixes → M1.5 devN tag |
| 2–3 | M2.1–M2.3 locally; M2.4–M2.5 once D2 is granted; M3 alongside |
| then | M4 per D4 → rc1 preparation stays gated on G5 (deferred experiments) |
