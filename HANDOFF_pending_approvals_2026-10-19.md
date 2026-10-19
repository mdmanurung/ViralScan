# Pending approvals and ownership — 2026-10-19 (written by Claude, Track A session)

Authorises nothing. Each row needs the user. Shared root `codex/drop-py39` @ 5e551c9 is dirty: never reset/clean/stash/stage wholesale. Clean main: `../vs_work/main` @ 6bd1652.

## Needs a user decision
| # | Item | Owner | Decision needed |
|---|------|-------|-----------------|
| 1 | NORM-SHIFT-01 (`normalize_total(target_sum=None)` then `log1p`) | Codex 01a15016, `/tmp/viralscan-shifted-log-candidate`, branch `codex/shifted-log-normalization` c528ab3 | push + open PR; Track A keeps pinned CP10k (user approved), shifted-log sensitivity run after merge |
| 2 | scAutoQC (`--scautoqc`) | Codex 01a1513d; uncommitted in shared root (`src/viralscan/scautoqc.py`, `scripts/scautoqc.py`, `tests/test_scautoqc.py`, edits to Snakefile/parsers/menu/runconfig/scripts/umap.py/environment.yml/pyproject.toml) | commit it; scope of real FASTQ-to-UMAP validation on SRR8315713 (needs kb/kallisto/bustools install + 5-9 GB run). Handoff reconstructed in `HANDOFF_scautoqc_2026-10-19.md` |
| 3 | CPU budget B11/B30 | Codex 01a15412, `/tmp/viralscan-cpu-budget` (branch `codex/total-cpu-budget`, clean) | user said "approve" 14:11 (per Codex rollout); byte-difference findings need a separate decision. `HANDOFF.md` still says pending (stale) |
| 4 | TRACKA-HOST-02 scaling pilot (80 -> ~300 genes, 200 cells, 1 CPU) | Codex 01a15152, `analysis/trackA_host_programs/`, handoff only in `/tmp/handoff-1Dvp4X.md` | new SLURM submission or reuse of allocation 25833654; full-cohort fitting NOT authorised |
| 5 | Perf remainder (B21/B13/B28/B16, C00, PR) | Codex, `/tmp/viralscan-perf-remaining-candidates` (`codex/perf-remaining-candidates`) | continue timing/native/C00; whole-run proposal `docs/plans/2026-10-19-performance-whole-run-proposal.md`; no SLURM without approval |
| 6 | Track A FASTQ/BUS deletion (hsv1, hhv6b, ebv, hnscc) | Claude | dry-run first, then explicit OK per dataset; hhv6b may be dry-run alone once its anello classify+evidence finish |
| 7 | Git delivery of Track A (report, scripts, `.living`, PLAN.md rows) | Claude | one commit via a separate worktree after the report is done; user reviews and pushes |

## Unassigned prototypes (stay with Codex unless handed over)
- Viral read cleanup/assembly pilot VR-01..05: `analysis/viral_read_cleanup/`, stopped before real data (7 Ruff findings, repeat-offset regression); resume in a fresh attempt dir, not v3.
- Translated-protein search: `analysis/translated_search/`, reference built, known-positive pilot partial.

## Housekeeping (not mine, left alone)
- Wedged Codex sandbox `find` over /exports (pid 1621907, state D, ~5 h).
- Stray half-checked-out `/tmp/viralscan-cpu-budget-candidate`.

## Added after the scAutoQC reconstruction (read-only verified)
- scAutoQC has no PLAN.md row; CLAUDE.md requires one before any commit.
- With `--scautoqc`, `src/viralscan/scripts/umap.py:214` skips the legacy count/gene cell filter, so the standard UMAP/HVG also run on low-count cells: decide before commit.
- The thread also touched `run_safety.py`, `tests/test_run_safety.py` and three docs files (not listed above). Focused tests re-run: 3 passed. Reconstructed handoff: `HANDOFF_scautoqc_2026-10-19.md`.
