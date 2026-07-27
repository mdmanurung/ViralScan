# Review — codex/viralscan-v3 vs main — 2026-07-27

**Scope**: branch vs main, focused on this session's work
**Files reviewed**: 8 code files + `analysis/v3_validation/protocol.yaml` (branch total: 445 files, +46506/−15166)
**Sub-agents run**: 6 of 6

Two Major findings, seven Minor. One reported Major was dropped in synthesis as a
false positive; a second was downgraded to Minor on verification — see "Corrections
made during verification".

## Key decisions in this analysis

- **Allocation and resampling unit is the whole biological sample** — for both the
  train/holdout split and the bootstrap, explicitly rejecting cells, molecules,
  reads, and technical repeats. The correct anti-pseudoreplication choice.
- **Cell-calling anchor is computed once per dataset/library and reused across every
  workflow row** — from a host-only matrix, so it carries no infection outcome.
  See F2, which breaks its reproducibility.
- **Threshold calibration is a closed pre-declared grid, training-only, under a hard
  zero-false-call constraint** with a one-standard-error conservative tie-breaker.
  This is the standard mitigation for optimal-cutpoint bias, applied correctly.
- **The v3-vs-2.2.0 improvement claim is restricted to one axis** — molecule
  precision/recall/F1 on `counts_unique` vs `counts_original` — because 2.2.0's
  reported product mixes UMI-deduplicated integers with read-weighted fractions.
  General-superiority claims are prohibited.
- **`E1`'s count invariant is scored only against tools that declare it** (today,
  ViralScan alone). Correct — scoring comparators against a claim they never made
  would make a v3 safety hypothesis permanently unpassable. Watch the framing
  downstream.
- **Cell calling now fails closed** (`SW-11`), with `--cell-calling none` as the
  explicit opt-out. See F6 for the one cost.
- **Frozen-input identity is established once by a full-stream audit and re-checked
  at run time by byte count only** (`audit_fastq_pair.py` streams SHA-256 and MD5;
  `verify_frozen_fastq` compares `st_size`). A deliberate cost trade-off, honestly
  documented. See F1 for the residual gap.
- **Ledger tamper-evidence is split across two modules** by git-dependency:
  `check_ledger_append_only.py` anchors on git history, `validate_v3_protocol.py`
  must run without git. Principled, not accidental duplication.

## Questions for the analyst

- **Was `seeds.cell_calling` meant to reach emptyDrops, or is the R-level seed
  considered out of scope for reproducibility?** F2 assumes the former, because the
  protocol names `seeds.cell_calling` as the seed source for the shared anchor.
- **Between the freeze-time audit and a run, is the storage holding the frozen
  FASTQs immutable?** F1's severity depends entirely on this. On a shared
  filesystem with restores and symlinks, a size-preserving content swap is
  plausible; on write-once archival storage it is close to theoretical.
- **When `VAL-01` freezes the abundance and homology level counts, will you recompute
  the per-stratum holdout sample count?** At the illustrative 63 strata it is ~1.2
  per stratum against your own 3-sample floor for an interval (F4).
- **Is `EBV/KSHV` intended as a discrimination-gradient anchor, or will it be read as
  a claim that they are siblings?** The protocol declares it; the shipped
  `SIBLING_VIRUS_PAIRS` constant disagrees (F5).
- **Should `PLAN.md`'s dated evidence log be treated as append-only history or as
  current state?** F7 only matters under the second reading.

## Findings

### Bioinformatics

#### Major

##### F2. The frozen cell-calling seed never reaches emptyDrops
`src/viralscan/scripts/cellcalling.py:208`
```python
        return emptydrops_cells(
            obs,
            mdir,
            rscript=getattr(config, "cell_caller_rscript", "Rscript"),
            fdr=float(getattr(config, "emptydrops_fdr", 0.01)),
            lower=float(getattr(config, "emptydrops_lower", 100)),
            niters=int(getattr(config, "emptydrops_niters", 10000)),
        )
```
**Why it matters here**: the plumbing is complete on both ends and severed in the
middle. `emptydrops.R:65` calls `set.seed(seed)`, reading `args[[6]]`;
`emptydrops_cells` appends `str(seed)` as the sixth argv element
(`cellcalling.py:148`). But `call_cells` passes `rscript`, `fdr`, `lower`, and
`niters` and stops — so emptyDrops always runs at the signature default
`seed=100` (`cellcalling.py:135`), and the protocol's frozen
`seeds.cell_calling: 20260727002` reaches nothing. `config` has no seed field for
it to read; the `call_cells` docstring lists
`emptydrops_fdr/emptydrops_lower/emptydrops_niters` and no seed. The shared
outcome-independent cell anchor — reused across every workflow row, and the
denominator under every cell-level metric including the one axis the v3-vs-2.2.0
comparison is allowed to use — is therefore not reproducible from the frozen
protocol as the protocol promises. This is also default-parameter smuggling: a
scientific constant living in a signature rather than in the declared
configuration.
**Fix**: add an `emptydrops_seed` config field wired from `seeds.cell_calling`,
pass it through `call_cells`, and drop the signature default so an unset seed is
an error rather than silently 100.

#### Minor

##### F5. `SIBLING_VIRUS_PAIRS` omits the EBV/KSHV pair the protocol declares
`src/viralscan/constants.py:156`
```python
SIBLING_VIRUS_PAIRS: dict[str, str] = {
    "Human herpesvirus 6": "Human herpesvirus 6b",
    "Human herpesvirus 6b": "Human herpesvirus 6",
    "Human herpesvirus 1": "Human herpesvirus 2",
    "Human herpesvirus 2": "Human herpesvirus 1",
}
```
**Why it matters here**: `protocol.yaml:515` declares three sibling pairs for
`H4`/`E5`/`D13`/`D24`, the third being `[EBV, KSHV]`. The shipped constant
consumed by `check_sibling_crossmapping()` has two. The protocol registers a public
EBV+KSHV co-infection dataset (`kshv_ebv_gse154900`), so the gap is reachable.
**Fix**: add the pair to the constant, or state in the protocol that the third pair
is scored analytically rather than by the shipped cross-mapping check.

##### F6. `Rscript` is absent from the CLI preflight, so an R-less run fails late
`src/viralscan/menu.py:38`
```python
REQUIRED_TOOLS = ("kb", "snakemake")
```
**Why it matters here**: `cell_calling=auto` resolves to `emptydrops` whenever no
external list is supplied (`cellcalling.py:188`) — the documented default path.
With `SW-11` making cell calling fail closed, an environment without
`Rscript`/DropletUtils now aborts the run *after* `kb_count`, `analysis`, and
`multimap` have completed, rather than in `_check_required_tools`
(`menu.py:1429`). Failing closed is right; failing closed late is expensive.
**Fix**: add `Rscript` to the preflight when the resolved method is `emptydrops`.

### LLM coding antipatterns

#### Major

##### F3. `check_git_sha_fields` has fail-open branches that skip verification silently
`scripts/check_ledger_append_only.py` (`observed is None` branch, and the
`shown.returncode != 0` continue)
```python
            if observed is None:
                # The section did not exist at the base commit ...
                pass
```
**Why it matters here**: when the claimed `digest_scope` names a section absent from
the base commit, or `git show` cannot read the path, the record's before-digest is
accepted unverified — a check that reports success without having checked. This is
the same class as the `except Exception` removed from `cellcalling.py` this session.
Two independent reviewers found it; it matches the project's own open `R11-F2`.
**Fix**: distinguish "genuinely genesis" (before == after) from "could not verify",
and make the latter an error.

### Data pipeline & leakage

#### Minor

##### F1. Frozen FASTQ identity is audited once and never re-verified at run time
`scripts/run_fresh_control.py:49`
```python
    if path.stat().st_size != expected_bytes:
        raise FreshControlError(f"stored byte count drifted: {path}")
    if re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None:
        raise FreshControlError(f"invalid audited storage SHA-256: {path}")
```
**Why it matters here**: the SHA-256 is real — `scripts/audit_fastq_pair.py:63`
streams both storage and content digests, and `freeze_control_inputs.py` records
them. `verify_frozen_fastq`'s docstring ("Verify runtime size against an identity
established by a full-stream audit") and its `runtime_size_verified` field are both
accurate about what happens. The residual gap is the window *between* the audit and
each run: a size-preserving content change — a stale restore, a repointed symlink,
the wrong sample of equal length — is not detected at run time, and
`prepare_fresh_controls.py:151` re-checks the same way when building the task
manifest, even though it already imports `hashlib` and hashes the cache manifest at
line 111. Attempt 3 will read these files weeks after the audit.
**Fix**: recompute the digest at manifest-build time at minimum, keeping the size
check as a cheap pre-filter; hashing a multi-GB FASTQ is small against a `kb count`
run. If the run-time check stays size-only by design, rename
`audited_storage_sha256` to something that cannot read as "verified just now".
**Also**: `test_fresh_control_refuses_fastq_storage_size_drift` passes a fabricated
`"a"*64` digest and only exercises the size branch. That is currently all it *can*
do; it should gain a real hash-mismatch case alongside the fix.

### Statistics & causal inference

#### Minor

##### F4. Per-stratum holdout counts likely fall below the protocol's own interval floor
`analysis/v3_validation/protocol.yaml` (`partitions.minimum_samples_rationale`,
`calibration.uncertainty.minimum_samples_for_interval: 3`)
**Why it matters here**: at the illustrative 63 strata, 252 samples and 76 in
holdout is ~1.2 holdout samples per stratum, against a stated 3-sample minimum for
a BCa interval. Most per-stratum holdout estimates would be points with
not-estimable intervals, weakening the per-stratum reporting `E3`/`E4` promise. No
acceptance rule is invalidated — neither endpoint carries a pass/fail threshold.
**Fix**: recompute once `VAL-01`/`REF-08` freeze the level counts; if it stays near
1, report per-stratum results at a coarser aggregation, as
`uncertainty.stratified_by` already does by dropping homology.

### Documentation & schema fidelity

#### Minor

##### F7. `PLAN.md` dated evidence entries state sections are frozen; all four are pending
`PLAN.md:707`
```
  `failure_and_deviation_reporting` are frozen, schema-validated sections. The
```
**Why it matters here**: `partitions`, `calibration`, `workflow_matrix`, and
`failure_and_deviation_reporting` are all `status: pending` in the live protocol,
having been downgraded after SCI-05 round 1. The work-package rows are correct and
say "**`pending`**, not frozen", so the authoritative status is right; only the
dated evidence log reads as current. Under the project's own convention that log is
append-only history.
**Fix**: append a superseding line to the log rather than editing the historical
entries.

##### F8. "the last two returning no blocker" invites reading the reviews as passed
`PLAN.md:23`
**Why it matters here**: accurate — rounds 10 and 11 each returned no blocker — but
both verdicts were `does-not-pass` with open Majors. One reviewer read it as a
contradiction, which is evidence the phrasing misleads.
**Fix**: say "no blocker, verdict still does-not-pass, Majors open".

##### F9. `detect_cells` gained a parameter without a docstring update
`src/viralscan/scripts/detection.py:329`
```python
def detect_cells(adata, found_genes, summary, viral_count_matrix=None):
```
**Why it matters here**: `viral_count_matrix` changes which count layer the function
reads — it feeds `resolve_count_matrix` on the next line. The Parameters block stops
at `summary`. Sibling functions in the same diff had their docstrings updated.
**Fix**: document the parameter and its default behaviour.

## Corrections made during verification

- **F1 downgraded from Major to Minor.** The first pass reported that
  `verify_frozen_fastq` "never calls `hashlib`" and that `audited_storage_sha256`
  named a guarantee the code did not deliver. Both were wrong on checking: `hashlib`
  is imported and used at `run_fresh_control.py:232`, and the FASTQ digests come
  from a real full-stream audit in `audit_fastq_pair.py`. The docstring and
  `runtime_size_verified` describe the size-only run-time check accurately. What
  remains is the narrower and genuine point above.
- **Dropped entirely: "nothing enforces the two schema copies stay identical"**
  (code-quality, Major). `tests/test_validation_protocol.py:752`
  `test_packaged_schema_matches_the_canonical_schema` asserts byte-identity, and
  both copies were verified identical.

## What was checked but is fine

- **Statistics & causal inference**: no p-value family, so no multiple-comparison
  gap; cutpoint selection correctly mitigated; the one-axis restriction on the
  v2.2.0 comparison is a genuine like-for-like choice.
- **Data pipeline & leakage**: leakage prohibitions cover reference, D-list, and
  feature-universe construction; the cell anchor is host-only and outcome-free;
  nothing compares across the mismatched 2.2.0 unit boundary.
  `prepare_fresh_controls.py:111` does hash the viral-data cache manifest and refuse
  on drift.
- **Bioinformatics**: emptyDrops parameters otherwise match the frozen contract;
  knee is correctly confined to sensitivity-only; sample structure (2000 cells, 25k
  reads, 1 technical replicate, 10% mixed) is realistic and justified; no surviving
  HIV assumption anywhere.
- **LLM antipatterns**: the two `except Exception` blocks in the comparison scripts
  convert to explicit failure status with nonzero exit — the opposite of the removed
  antipattern. No hallucinated APIs found.
- **Doc & schema fidelity**: schema copies byte-identical; all 24 hardcoded
  denominator tuples and 11 audit-artifact column sets match the protocol exactly;
  ledger docstrings are carefully self-limiting rather than overclaiming.
- **Code quality**: no bare `exit()`, no `shell=True` with user paths, no
  string-concatenated paths, no institutional absolute paths introduced by this diff.

## Notes

- **F2 and F3 are both "did this ever run" questions**, and neither would surface in
  a passing test suite. F2's seed never reaches R; F3's branch reports success
  without checking. F1 shares the shape but not the severity — there the check that
  does not run is a *deliberate* omission, documented as such.
- **Eleven rounds of SCI-05 could not have found F2.** Those rounds read
  `protocol.yaml`. A frozen constant that no call site passes is invisible from the
  protocol side and invisible from the output side — the run succeeds and produces
  plausible cells. Every entry under `seeds:` deserves a "does this reach its
  consumer" test; see the tripwire audit.
- **Pre-existing institutional absolute paths** exist in ~45 tracked files outside
  this diff (e.g. `scripts/hostresponse_ebv_matched.py:22`). Out of scope here, but
  it contradicts the `CLAUDE.md` rule and deserves its own pass.
