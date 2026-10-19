# Handoff — Track A: anello v3 + real libraries, hostresponse controls, deletions held

**Date:** 2026-10-19 (evening) · **Branch:** shared root `codex/drop-py39` @ `5e551c9` (shared with Codex; never edit/checkout/stage wholesale) · **Status:** ~290 anello jobs and 2 hostresponse jobs queued; nothing committed.
Written here, not `HANDOFF.md`: root `HANDOFF.md` is the Codex performance/CPU-budget thread and must not be overwritten (it was already clobbered once).
`R=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/trackA_round1`, `P=$R/anello_settle_2026-10-19`, `M=/exports/archive/hg-funcgenom-research/mdmanurung`. Env `$M/conda/envs/viralscan_test_full/bin` on PATH; pinned code `$M/vs_pinned/97a64f6/src`. Tonsil cancelled.

## Goal
Report ViralScan as recovering the previously reported viruses (and what else), settle the anellovirus question under the frozen PREREG, finish hostresponse controls, then (user-gated, per dataset) delete FASTQs/BUS.

## Next action
When `squeue -u $USER | grep -c tA_anello` shows the v3 chain done (classify jobs 25842488/90/93/95/98/500), score v3 and read the gates:
`cd $R; P=$PWD/anello_settle_2026-10-19; for ds in hnscc_gse164690 hsv1 ebv; do for s in shipped relaxed; do python scripts/trackA_anello_score.py $P/plant_v3_tail5/$ds $P/plant_v3_tail5/$ds/$s/classify $P/plant_v3_tail5/$ds/$s/score.tsv; done; done`
Gates: STAR exact@100 >= 0.90, candidates >= 0.50 (`$P/PREREG.md`, `$P/PREREG_ADDENDUM_1.md`). Then give the user a Track A status table (standing rule).

## State
- v2 plant scores (rescored with CB-correction): exact@100 STAR shipped/relaxed hhv6b 0.952/0.965 PASS; hnscc 0.707/0.778, ebv 0.728/0.802, hsv1 0.808/0.857 FAIL (poly-A tail geometry: pinned planter tails 1/3 of 3' reads). Old raw-key scores: `plant/<ds>/<setting>/score_rawkey_v2.tsv`.
- Queued: v3 plant chain `$P/plant_v3_chain_jobs.txt` (5% tail, out `$P/plant_v3_tail5/`); hhv6b real libs `$P/real_hhv6b_jobs.txt` (26 STAR + classify); hsv1/ebv/hnscc real libs `$P/real_other_jobs.txt` (111). Outputs: `anello_align[_relaxed]_2026-10-19/<ds>/<lib>/`, `anello_classify[_relaxed]_2026-10-19/<ds>/<lib>/`. Not yet checked that the first real classify reads the real BAM path correctly.
- Hostresponse done: ebv kNN at caliper 0.5 AUC 0.966 (0.3: 0.960; PC-only 0.86 -> upper bound; old files `*.caliper0.3.*`). hnscc 0.787, hpv16 0.907 (kNN files in `hostresponse_knn_2026-10-19/`).
- Hostresponse queued: `trackA_hostresponse_xlib.py` jobs 25842522 (hsv1 GSM3511335 vs 37C uninfected GSM3511328,GSM3511329, has batch null) and 25842523 (hhv6b GSM6416024 vs D38 GSM6416026, no null possible). Output `hostresponse_xlib_2026-10-19/`.
- Running elsewhere: agent `ab6f235e…` prepares rc1 blockers in worktree `$M/vs_work/rc1-prep` (branch `codex/rc1-blocker-prep`, uncommitted: Action SHA pins, REF-11 bundle script + draft Zenodo metadata, PANEL-01 WP5 sbatch, REL-15/06 checklist). Review its report before anything else touches that worktree.
- Uncommitted by me: `scripts/trackA_anello_{score,plant}.py` edits (`plant.py.v2bak` kept), `.living/{learnings,decisions}.md` entries, `HANDOFF_pending_approvals_2026-10-19.md`, `HANDOFF_scautoqc_2026-10-19.md`. Disk 94%.

## Locked decisions (user approved the grilling rounds 2026-10-19)
- No PREREG threshold changes; v3 re-plant only via addendum 1 (sha 5f75b6a7…), v2 kept as pessimistic bound; rate chosen after measuring real tails (ebv 0.4%, hnscc 0.1%, hsv1 2.6% trailing A>=10).
- Run real arrays for hsv1/ebv/hnscc now; no throttle (scheduler fair-share); remote blastn on candidates ONLY after showing the user the read/library count first.
- hostresponse: Track A keeps pinned CP10k (shifted-log is a post-merge sensitivity); hsv1/hhv6b via cross-library controls; headline = adjusted-over-batch-null AUC for hsv1, raw AUC + batch caveat for hhv6b (same-donor GSM6416023 is infected, 282k HHV-6B UMI), kNN AUC with diagnostics for hnscc/ebv/hpv16; "association, not a virus effect".
- Track A is a separate internal report, not in the rc1 manuscript (manuscript = frozen v3 artifacts only). Priority: Track A -> rc1-blocker prep -> bookkeeping for done Codex threads -> park perf remainder, VR-01..05, translated search until after rc1.
- Deletion needs explicit OK per dataset; dry-run hhv6b alone once its anello classify + evidence finish; keep every twostep `host_filtered` and `lineage_bus` until anello is done; hpv16 FASTQs already gone.
- Git: one commit via a separate worktree after the report is done (report folder, scripts, `.living`, PLAN.md rows incl. ANDET-09f note); the user reviews and pushes. Claude owns Track A only; Codex-owned threads are listed in `HANDOFF_pending_approvals_2026-10-19.md`.
- NORM-SHIFT-01 PR: user (not Claude) pushes. scAutoQC: keep legacy cell filter, validate on a ~1M-pair subset (user installs toolchain). CPU budget: any raw-byte difference stops and returns to the user. TRACKA-HOST-02 scaling pilot authorised, run by its Codex owner.

## Dead ends — do not redo
- Scorer keyed on raw planted CB (fixed; STARsolo rewrites CB to a whitelist neighbour). hsv1 Drop-seq has no whitelist -> exact keys.
- `--label cpm` for hostresponse (wrong cells); raw counts into hostresponse (prep script writes log1p CP10k).
- k>=3 matching for ebv (383 controls < 316 x 3); a looser caliper only adds positives.
- EM rerun for ebv GSM4796271 (raw `output.bus` deleted); `hostresponse` CLI inside run dirs (rewrites `run_complete.json`).
- Unquoted paths / `python` not on PATH in sbatch wrappers; never use `->` in agent briefs (clobbered 10 FASTQs once).

## Read first
- `$P/PREREG.md`, `$P/PREREG_ADDENDUM_1.md` — rules, gates, tail-rate addendum.
- `$R/scripts/trackA_anello_{plant,score,classify,align}.py`, `trackA_hostresponse_{knn,xlib,prep}.py`.
- `analysis/trackA_report_2026-10-18/FINDINGS_2026-10-18.md`, `completeness_critic.md`.
- `.living/learnings.md` and `decisions.md` entries dated 2026-10-19; `HANDOFF_pending_approvals_2026-10-19.md`; `.mycelium/last-session.md`.

## Open items
- After v3 scoring: if a dataset still fails, verdict is NOT RESOLVABLE for it; else classify real libs -> blastn count to user -> approved remote blastn -> evidence for the 5 largest kallisto anello calls (`$R/evidence_anello_2026-10-19/`) -> verdicts + FINDINGS + `.living` finding + PLAN.md `ANDET-09f` note.
- Check older jobs: `tA_advance`, validate-run arrays, combined_artefact_L, summaries/provenance chain, evidence-extra array 25836013, scCLR UMAPs 25836018 (only `tA_advance` was still queued).
- Report fixes: counting convention (`tiers.tsv` only), re-verify paper quotes, section-8 hostresponse numbers, per-molecule assignment counts vs EM pilot, PLAN.md row for the upstream `--multimap-molecule-assignments` patch.
- Deletion dry-runs (`delete_dataset.py`/`delete_bus.py`): re-run `check_ena_readcounts.py` first; snapshot `provenance/` + `concat_plan.tsv` before deleting.
