# Handoff — ViralScan resource policy after PR #76

**Date:** 2026-10-19 · **Branch:** shared root `codex/drop-py39` @ `5e551c9` · **Status:** PR #76 merged; CPU-budget proposal awaits user approval.

## Goal
Finish B11/B30 resource correctness with explicitly approved total CPU-budget semantics and native regression evidence.

## Next action
Re-present the pending CPU-budget proposal below and obtain the user's compatibility decision before changing workers or scheduler declarations. The handoff request is not approval of that proposal.

## Pending proposal (verbatim)
Next is **CPU-budget enforcement (B11/B30)**. Currently, `--cores 6` permits kb to request eight workers, while CD-HIT requests all visible CPUs.

I recommend:

1. Make `--cores` the total CPU budget across workflow and reference-building entrypoints.
2. Bound native workers and scheduler reservations to that budget.
3. Requalify locally at 1/6/8/32 cores, comparing native outputs and failure behavior against the frozen baseline.

The alternative is preserving fixed workers and rejecting insufficient budgets—including today’s six-core default.

Approve the bounded-worker approach?

No worker changes made yet. Changed thread counts may change output bytes; any differences require a separate compatibility decision.

## State
- PRs #73/#74/#75 and #76 delivered; PR #76 rebase merge is `6bd1652779dc64e5ba68cbcdafe0b322ecbcdb6c`, verified against remote main this session.
- Persistent main is `../vs_work/main`; it was clean at the merge verification. Shared root is intentionally older and dirty: 22 tracked edits plus untracked research work; never reset, clean, stash, or stage wholesale.
- PR #76 passed all 12 CI checks: Linux 5739 passed/57 skipped and macOS 5737 passed/59 skipped per Python 3.10–3.12; native CI 87 passed/23 skipped.
- Local native 54 passed/no skips; scheduler/failure 34 passed; paired DAG 32 cases/64 runs passed; 83 runtime files matched the tested freeze.
- Both owned resource worktrees and their local/remote feature branch are already removed after archive verification. Do not repeat merge or cleanup.
- No jobs launched during this continuation. Current unrelated scheduler/process state is unverified.
- Root `HANDOFF.md` was absent when this handoff was written; earlier contents were carried forward from this conversation.

## Locked decisions
- Explicit legacy profile is delivered: kb count/ref `-t 8`, reference homology minimap2 `-t 3`; CD-HIT `-T 0 -M 8000`, numerical environment, scientific settings and scheduler declarations were preserved.
- Changed effective workers and total-admission policy need a separate compatibility decision; the last user action requested a handoff, not approval.
- Exactness uses raw A/A and A/B artifacts; never relax comparators, normalize away differences, or regenerate historical goldens to obtain a pass.
- Rebase-only PR delivery; preserve unrelated shared-root edits and concurrent `src/viralscan/scripts/umap.py` work.
- Preserve non-release tags `provenance/perf-wave1-2026-10-18`, `provenance/perf-wave2-2026-10-18`, and `provenance/perf-legacy-workers-2026-10-19`.
- No SLURM jobs, large builds/downloads, paid calls, or retained-data deletion are authorized by this completed resource phase.

## Dead ends — do not redo
- PR #76 is already merged; an earlier automatic approval rejection and pending status were superseded by live remote verification and the user's merge request.
- Adding Snakemake thread declarations alone can change inherited OMP/BLAS/NumExpr workers; native argv equality alone does not establish unchanged execution.
- Requested worker counts differ from transient sampled thread peaks; tiny fixtures do not prove sustained production utilization or whole-run gains.

## Read first
- `../vs_work/main/docs/plans/2026-10-19-performance-resource-policy.md:1` — current worker inventory, compatibility boundary and required matrix.
- `../vs_work/main/docs/plans/2026-10-17-remaining-performance-candidates.md:1035` — B11/B30 remaining scope and stop conditions.
- `../vs_work/main/docs/plans/2026-10-19-legacy-worker-qualification.md:1` — qualified legacy behavior and evidence limits.
- `analysis/computational_bottlenecks_2026-10-17/followups/20261019-084842-resource-legacy/FINAL_DELIVERY.json:1` — final delivery receipts; adjacent `qualification_status.json` reconciles merge/cleanup state.
- `.living/archive/2026-10-18-performance-landing/resource-legacy/RECOVERY.md:1` — verified worktree archives and complete Git bundle restoration.
- `.mycelium/last-session.md` — Mycelium continuation context; reconcile with current receipts.
- `.living/decisions.md` and `.living/learnings.md` — standing project decisions; find by topic heading, not positional IDs.

## Open items
- B11/B30: decide budget semantics before implementation; production concurrency and utilization qualification remain open.
- G02: production attribution and whole-run timing remain unmeasured; no whole-run speedup claim.
- NORM-SHIFT-01: user requires `normalize_total(target_sum=None)` then `log1p`, replacing fixed TP10K; separate scientific follow-up with raw counts/provenance and downstream validation.
- `analysis/normalization_shifted_log_2026-10-19/` is now untracked in the shared root; its ownership and progress were not inspected. Reconcile before touching it.
- Track A and other unrelated analyses: use their dedicated handoffs, including `HANDOFF_trackA_2026-10-19.md`; this handoff does not authorize their execution or deletion.
