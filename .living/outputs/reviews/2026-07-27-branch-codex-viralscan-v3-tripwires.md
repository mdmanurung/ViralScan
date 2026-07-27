# Tripwire audit & plan — codex/viralscan-v3 vs main — 2026-07-27

**Static review**: `.living/outputs/reviews/2026-07-27-branch-codex-viralscan-v3.md`
**Mode**: audit (this document describes what would be tested; no code was run
and nothing was modified)
**Hooks the project has**: none of the four standard ones, but three
project-specific analogs — see "Hooks" below

## What this audit does

Names the failure modes on this branch worth testing *behaviorally* rather than
by reading code, says what each test would actually do, and points at the ones
runnable today with no setup.

Three of the seven are runnable today, and all three would **fail on the current
branch** — they are the executable form of static findings F1, F2, and F3. A
tripwire that fails is doing its job; it converts "I read the code and the run-time
check is size-only" into "I swapped the file and nothing noticed".

The standard starter four (`missing-counts-file`, `missing-metadata-sample`,
`label-permutation`, `toy-contrast-direction`) map poorly here. This branch's
scientific boundaries are not a counts/metadata join — they are *frozen input
identity*, *frozen parameter reach*, *append-only evidence*, and *fail-closed cell
calling*. The tripwires below are the project-specific equivalents; tripwire5 is
the closest thing to `label-permutation`, and tripwire7 is the closest thing to
`toy-contrast-direction`.

## Tests that would apply here

### tripwire1. The frozen-input-really-checked check
**Watches for**: a frozen FASTQ whose *contents* changed after the freeze-time
audit while its byte count stayed the same — a stale restore, a repointed symlink,
the wrong sample of equal length.
**How it'd work**: copy a frozen FASTQ, flip one byte in place (size preserved),
point a control-input manifest at the copy with the *original* SHA-256, then call
`verify_frozen_fastq` and `prepare_tasks`. Both must raise.
**Slug**: `frozen-input-content-swap` (fault-injection category)
**Related static finding**: F1
**Today**: **runnable, and it fails.** The digests themselves are real —
`scripts/audit_fastq_pair.py:63` streams SHA-256 and MD5 over both mates at freeze
time. But nothing recomputes them afterwards: `verify_frozen_fastq`
(`run_fresh_control.py:49`) compares `st_size` and regex-matches the recorded
digest's *shape*, and `prepare_fresh_controls.py:151` does the same when building
the task manifest. So a post-audit same-size swap passes at both points.
**Starter check** (runnable today, no setup): a ~25-line pytest using `tmp_path` —
no real FASTQ needed, two files of equal length suffice. This is also the fix's
regression test, so writing it is not throwaway work.
**Note**: the existing `test_fresh_control_refuses_fastq_storage_size_drift` passes
a fabricated `"a"*64` digest and only exercises the size branch. That is currently
all it can do — no code path can detect a content mismatch — so it should gain a
real hash-mismatch case alongside the fix.
**Scope limit**: this tripwire checks that a swap *would be caught*, not that none
has happened. For the latter, recompute all ten digests once against the manifest;
nothing has re-verified them since the audit.

### tripwire2. The frozen-seed-actually-reaches-R check
**Watches for**: a seed declared and frozen in `protocol.yaml` that never reaches
the process it is supposed to make reproducible.
**How it'd work**: two forms, cheap and expensive.
*Cheap*: capture the `Rscript` argv/stdin that `emptydrops_cells` builds and
assert the value from `seeds.cell_calling` appears in it.
*Expensive*: run cell calling twice on one fixture with two different configured
seeds and assert the called-barcode sets differ, then twice with the same seed and
assert they are identical. The second half is the real reproducibility claim.
**Slug**: `frozen-seed-reaches-consumer` (known-answer category)
**Related static finding**: F2
**Today**: **the cheap form is runnable and would fail.** The seed does reach argv
when supplied — `emptydrops_cells` appends `str(seed)` as element six
(`cellcalling.py:148`) and `emptydrops.R:65` calls `set.seed(seed)` — but
`call_cells:208` passes `rscript`, `fdr`, `lower`, and `niters` and stops, so the R
step always runs at the signature default `100` regardless of the frozen value
`20260727002`. The expensive form needs `Rscript` + DropletUtils and a small matrix
fixture — which is `SW-10` in WP1.
**Starter check**: assert the seed is *plumbed*, not that it *works* — mock the
subprocess and check the value in the command. Ten lines, catches the entire
defect.
**Generalization worth doing once**: every entry under `seeds:` in `protocol.yaml`
deserves this test. A declared seed no call site passes is not frozen, it is
decorative, and nothing in the current gate distinguishes the two.

### tripwire3. The ledger-tamper-really-caught check
**Watches for**: an append-only evidence chain that reports success without having
verified anything.
**How it'd work**: four perturbations of a committed `deviations.yaml`, each
expecting a nonzero exit from `scripts/check_ledger_append_only.py` —
(a) mutate an existing record's `after` digest;
(b) delete a record from the middle;
(c) reorder two records so `record_sha256` no longer chains;
(d) **the fail-open probe** — set a record's `digest_scope` to a section that does
not exist at the base commit, and separately make `git show` fail for the scoped
path.
**Slug**: `ledger-append-only-tamper` (fault-injection category)
**Related static finding**: F3
**Today**: **runnable.** (a)–(c) should pass; **(d) is expected to fail** —
`check_git_sha_fields` treats a section absent at the base commit as acceptable
(`if observed is None: pass`) and `continue`s past a failed `git show`, so both
accept the before-digest unverified. This is the executable form of R11-F2.
**Starter check**: (a)–(c) need a scratch git repo, ~40 lines with
`subprocess.run(["git", "init"], ...)` in `tmp_path`. (d) is three more lines once
that harness exists.
**Why it matters that this one is behavioral**: this session already shipped one
ledger fix that was a **no-op** — `field_corrections` was written but
`_ledger_chain_for_section` never read it. It survived review because the four
findings beside it got tests and it did not. A tamper tripwire would have caught
it; reading the diff did not.

### tripwire4. The cell-calling-fails-closed check
**Watches for**: a cell-calling failure that silently degrades to "every barcode is
a cell", which does not lose the result — it changes the denominator under every
reported viral rate.
**How it'd work**: five injected failures, each expecting `CellCallingError` and a
nonzero exit rather than a full-barcode mask —
`Rscript` absent from `PATH`; DropletUtils not installed; an external barcode list
that is empty; one whose canonical barcodes collide; one that matches zero
observed barcodes. Then the control: `--cell-calling none` must succeed and report
over all barcodes, because that is the deliberate opt-out.
**Slug**: `cell-calling-fault-injection` (fault-injection category)
**Related static finding**: none — this is the behavioral confirmation of the SW-11
fix, plus the late-failure cost noted in F6.
**Today**: the four `external_cells` cases are runnable now against `call_cells`
directly. The two R-environment cases need either an `Rscript` stub on `PATH` or
`SW-10`'s end-to-end fixture.
**Starter check**: the missing-`Rscript` case is a one-line `monkeypatch.setenv`
on `PATH` — and it doubles as the test for F6, since it should ideally fail in
preflight rather than after `kb_count`, `analysis`, and `multimap` have run.

### tripwire5. The calibration-never-sees-holdout check
**Watches for**: threshold calibration, reference construction, D-list
construction, or feature-universe construction drawing on holdout information —
the leakage firewall `partitions.leakage_prohibitions` promises.
**How it'd work**: the label-permutation shape, adapted. Permute holdout
membership among samples, rerun calibration, and assert the selected threshold and
the constructed feature universe are **bit-identical** across permutations. If
either moves, holdout information reached a training-only step.
**Slug**: `partition-permutation` (metamorphic category)
**Related static finding**: none static — this is the behavioral counterpart to the
F17 amendment, which added the prohibition as *text*.
**Today**: **cannot run.** No SCI-04 executor exists — that is precisely what the
`execution_gate_unimplemented` blocker records. There is nothing to permute
against yet.
**Starter check**: none honest. Do not fake one; a permutation test against a
pipeline that does not exist would pass vacuously, which is worse than not running
it. Write this the same week the executor lands, before the first real calibration
run.

### tripwire6. The tracker-claims-match-the-protocol check
**Watches for**: `PLAN.md` asserting a protocol section is frozen when the live
`protocol.yaml` says `pending` — the exact drift F7 found.
**How it'd work**: parse `protocol.yaml` for every section's `status`, grep
`PLAN.md` for sentences naming those sections alongside "frozen" / "pending", and
flag disagreements. Run it as a test, so the tracker cannot drift silently between
sessions.
**Slug**: `tracker-status-freshness` (freshness category)
**Related static findings**: F7, F8
**Today**: **runnable, ~30 lines**, and it would flag the `PLAN.md:707` evidence
entry. Needs one wrinkle: the dated evidence log is append-only history by project
convention, so the check must either scope itself to the work-package rows or
accept a dated entry carrying an explicit superseding line.
**Why here specifically**: `PLAN.md` is the contract this repo says is
authoritative, and it is maintained by hand across sessions by an agent that also
edits the thing it describes. That is the standard setup for silent drift.

### tripwire7. The molecule-vs-read-unit check
**Watches for**: a comparison that rescales across the v2.2.0/v3 unit boundary —
summing UMI-deduplicated integers with read-weighted fractional multimap shares
and reporting the total as UMIs.
**How it'd work**: the toy-known-answer shape. Build a tiny fixture with a
hand-computable answer: one barcode, two viral genes, a known number of reads
collapsing to a known number of UMIs, and one multimapping equivalence class.
Assert v3's `counts_unique` equals the hand-computed molecule count, and assert the
comparison harness refuses to place v2.2.0's mixed-unit total on the same axis.
**Slug**: `unit-boundary-known-answer` (known-answer category)
**Related static finding**: none — this is the behavioral guard on the decision to
restrict the v3-vs-2.2.0 claim to one axis.
**Today**: partially runnable — the v3 half needs `SW-10`'s fixture; the refusal
half can be tested against the comparison scripts now.
**Starter check**: hand-compute the expected molecule count for a five-read toy
case and assert it, before any real comparison runs. This is the cheapest possible
defense against the finding recorded in F-007, and it is the kind of number that
becomes unfalsifiable once real data is flowing.

## Hooks this project would need to run the rest

None of the four standard hooks are present. But this repo has three
project-specific analogs that are closer than the count suggests, and naming them
matters because the gap is smaller than "0 of 4" implies:

- **Per-task status JSON** (`scripts/run_fresh_control.py` writes `run_id`,
  `stage`, `attempt_id`, `scientific_parameter_hash`, `exit_code`) — this *is* a
  checkpoint log, one record per terminal stage rather than one line per boundary.
  It is what tripwire1 and tripwire4 would assert against.
- **The hash-chained deviation ledger** — a stronger evidence trail than a drop
  ledger, and the substrate tripwire3 perturbs.
- **`validate_v3_protocol.py --phase`** — a gate that already knows how to refuse.
  The nearest thing to `--stop-after`.

What is genuinely missing:

- **A runnable end-to-end fixture** (`SW-10`, WP1). It is the single blocker on
  tripwire2's expensive form, tripwire4's two R cases, and tripwire7's v3 half.
  One tiny fixture unlocks most of this document.
- **An SCI-04 executor.** Nothing to permute for tripwire5.
- **A label/partition declaration.** Once the executor exists, tripwire5 needs to
  know which field carries holdout membership.
- **Per-boundary checkpoint emission inside a single run.** The status records are
  per-task; a within-run log would let a tripwire assert *where* a pipeline
  stopped, not just that it did.

## What this audit does NOT cover

- Whether the frozen FASTQs on storage are currently correct. These tests check
  that a swap *would be detected*, not that no swap has happened. If you want the
  second, recompute all ten digests once against the manifest — worth doing after
  F1 is fixed, since nothing has ever verified them.
- Numerical reproducibility against a prior run — that is a snapshot test.
- Reproducibility of manuscript figures — a CI concern (`G6`).
- Runtime and resource regressions, including whether the new 384 GiB highmem tier
  is actually sufficient. That is an empirical question for attempt 3.
