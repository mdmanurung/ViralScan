# Handoff — ViralScan audit fixes: land PR #71, clean worktrees, then AUD3-08

**Date:** 2026-10-15 · **Branch:** PR branch `codex/aud3-integration` @ 2a5a1a2 (worktree `vs_work/aud3-integration`); root checkout is `codex/drop-py39` @ 5e551c9, shared with Codex and another Claude session, so never edit or checkout there · **Status:** PR #71 open, CI re-running (8 pass, 4 pending at last look); nothing merged yet

(`HANDOFF.md` in this directory belongs to a parallel Track A session; not touched.)

## Goal
Merge PR #71 (11 verified audit-fix units, AUD3-01a…09b), remove my worktrees and branches, then implement AUD3-08 (resume guard). Done = PR merged with green CI, worktrees gone, AUD3-08 PR open.

## Next action
`gh pr checks 71` (repo `mdmanurung/ViralScan`). If all pass: `gh pr merge 71 --merge` (earlier PRs used merge commits; never `--admin`). If a job fails: read `gh run view <id> --log-failed`; the first run failed only `tests/test_hostresponse_warnings.py::test_outputs_unchanged_for_fixed_seed` (platform tolerance, fixed in 2a5a1a2).

## State
- Uncommitted in root checkout: Codex/`.living/`/Track A files, not mine. My notes are untracked in `biomni_bundle/plans3/verification/` (PR_DESCRIPTION.md, report, matrix, plan, ~40 repro scripts).
- Tests: full suite on the integration branch at 0d6370e = 4240 passed, 15 skipped (baseline main 4121); ruff, mypy, inventory, governance, protocol validators pass.
- Running: only GitHub CI on #71. All agents finished.
- Mine to remove after merge (12 + 1 worktrees): `vs_work/aud3-{01a,01b,02,03,04,05,06a,06b,07,09a,09b,integration}` and `vs_pinned/87c9e35` under `/exports/archive/hg-funcgenom-research/mdmanurung/`; local branches `codex/aud3-*` (11 unit branches + `codex/aud3-integration`, also on origin). Leave other `vs_pinned/*` and `/tmp/...scratchpad/wt_*` worktrees (other sessions). Delete a branch only after checking its tip is an ancestor of `origin/main`.
- Stale 0-byte `.git/worktrees/<name>/index.lock` appears after test runs; check `ps` for a live git process, then remove (own worktrees only).

## Locked decisions (user, 2026-10-15)
- Frozen protocol definitions never change (`analysis/v3_validation/protocol.yaml` harmonization block); only warnings, NA-with-reason, docs.
- Per-unit `codex/aud3-*` branches off `origin/main`, one worktree + agent each; push/merge of #71 explicitly authorised by the user.
- AUD3-08: option (a), resume-time guard as its own small unit after the PR: for every sample with `log/kb.done`, require its counts file, else stop ("delete log/kb.done or use --overwrite").
- Accepted: digest-only re-pin of claim `v3-reference-homology`; `--cv auto` as the pipeline default; AUD3-05 flags only near-exact copies (threshold is REF-08's).
- Governance pins (inventory, claims registry) re-pinned once, last.

## Dead ends — do not redo
- System `python3` is 3.6; use the conda env `viralscan_test_full` python.
- pytest `--basetemp` dir must exist, else hundreds of fixture errors.
- `ruff check <file>` alone reports I001 on `detection.py`; CI runs `ruff check .`, which passes.
- Running many full suites at once makes the R cell-calling test flaky (`AUD3-FLAKY`).
- The audit's own repro scripts hard-code `/workspace/vs_audit/...`; use the patched copies in `repro/`.

## Read first
- `PLAN.md` rows `AUD3-*` (on the PR branch) — decisions and open rows.
- `analysis/audit_2026-10-verification/verification_report.md` (on the PR branch) — outcome table per finding.
- `src/viralscan/run_safety.py:407-449` (resume branch), `src/viralscan/Snakefile:136-140` (`kb_count` outputs), `menu.py:~1297,~1371` — AUD3-08.
- `biomni_bundle/plans3/verification/repro/P7-01b_resume_after_purge.py` — reproduces silent reuse after a purge.
- `.mycelium/run/claude/1abbef7e-ed3e-43e2-9d5f-e6fd737c4347/last-session.md` — earlier handoff for this task (partly superseded).

## Open items
- `AUD3-08` guard — start after #71 merges, off the new `origin/main`.
- `AUD3-02-OPEN` (viral-free UMAP null still above nominal), `AUD3-FLAKY`, `AUD3-G3-EVIDENCE`, `AUD3-REF08-DIVERGED` — logged in PLAN, not scheduled.
