# Dataset run plan: every dataset not yet fully analysed (2026-10-05)

Status: **proposed**. Nothing here has been submitted. Tracker rows: PLAN.md
WP6A (`RUN-01`…`RUN-05`), WP6C (new: `DSR-01`…`DSR-08`), WP5 (`VAL-01`…`VAL-10`),
`TONSIL-02`, `PROG-16`, `PROG-18`, `HPV-11`, `DEF-07`/`CAT-30`.

Standing constraints:
- **Package first** (user, 2026-10-04, `2026-10-04-package-completion-plan.md`).
  Experiments wait until M1–M3 close, unless the user explicitly releases a row.
- **Exploratory vs reportable.** Any run before the G4 reference freeze and the
  G3 defaults freeze is *exploratory*: it uses cat42d and pre-freeze defaults,
  and its numbers must be rerun before they are reported (R2.4).
- **Never seen before freeze:** `pbmc_10x_healthy_v3` (VAL-04 external
  negative) and any healthy-PBMC library that may become a
  `synthetic_factorial` background ("none already seen", R2.2).
- **User owns infra:** downloads, installs, pushes and CI are the user's.
- **Latest ViralScan version for every dataset** (user, 2026-10-05). See
  "Version rule" below. No result from an older version counts as analysed.

## Version rule

- **Latest today:** `3.0.0.dev1`, commit `12ceac7`
  (`git describe` = `v3.0.0.dev1-31-g12ceac7`, 2026-10-05).
- **Every existing run is `3.0.0.dev0`** (from `results/reference_provenance.json`):
  EBV, HHV-6B and HSV-1 (all arms under `ebv_latest_ref_2026-09-27/runs/`),
  HPV16 (`runs`, `runs_cat42`, `runs_cat42b`, `runs_cat42d`), and COVID
  (`results_v3_strand_reverse`, `results_v3_strand_unstranded`). They predate
  DEF-01, DEV-022…024 and the `3.0.0.dev1` bump, so **every dataset needs a
  rerun** to count as analysed on the latest version.
- **The version string is not enough.** `3.0.0.dev1` already covers 31 commits,
  and runs record only `viralscan_version`, not the commit. Each job script
  therefore writes `git -C $REPO describe --always --dirty > code_sha.txt` into
  the run dir and refuses to start on a `-dirty` tree (the COVID runs already
  did this). Adding the commit to `run_manifest.json` is tracked as `DSR-08`.
- **One round, one commit.** All datasets in a round run on the same pinned
  commit, recorded in the table below. If the code changes mid-round, either
  finish the round on the pinned commit or restart it; never mix commits.
- **Reportable numbers** come only from the G3/G4-frozen release, so this
  round is exploratory and its version is recorded for comparison.

## Overview

`Last version` = ViralScan version of the most recent run. `Target` = the
commit this round must run on (latest at launch; today `3.0.0.dev1 @ 12ceac7`).

| # | Dataset | Last version | Status | Blocker | Next action | When |
|---|---|---|---|---|---|---|
| 1 | SRR12682296 EBV (10x v2) | 3.0.0.dev0 | old ref, knee | none | `RUN-02` rerun on cat42d + emptydrops | now (exploratory) |
| 2 | SRR20710641 HHV-6B (10x 5′) | 3.0.0.dev0 | old ref, knee; strand unresolved; v3 two-step failed | geometry + strand | measure geometry, then `RUN-01` | now (exploratory) |
| 3 | SRR8315713 HSV-1 (Drop-seq) | 3.0.0.dev0 | old ref, knee; kb allowlist drops 28 % (F-018) | MECH-D/SW-21 decision | `RUN-03` after SW-21 | after SW-21 |
| 4 | GSE189670 HPV16 rafts (10x v3) | 3.0.0.dev0 | cat42d + emptydrops | none | rerun on target; read-validate HPV118/HPV29 (`DSR-03`) | now |
| 5 | COVID x213/x216 (10x 5′) | 3.0.0.dev0 | anellovirus artefact; reverse run on old ref | DEF-07 wording (user) | rerun cat42d + reverse + read filter (`DSR-04`) | now (exploratory) |
| 6 | SFL tonsil x223 (10x 5′ R2-only) | never run in ViralScan | minimap2 screen only (F-010) | `TONSIL-02` geometry | 1 M-read strand/geometry test | now |
| 7 | GSE190558 KSHV latent/lytic (10x v3) | never run | not downloaded | lane manifest frozen first | user downloads; freeze manifest | after manifest |
| 8 | GSE154900 KSHV+EBV BC-1 PEL (10x v3) | never run | not downloaded | read-cap decision frozen first | decide cap; user downloads | after cap decision |
| 9 | GSE164690 HNSCC HPV± (10x 3′ v2, ~1 TB) | never run | scoped only | download size | one-lane pilot HN18 vs HN01 (`DSR-06`) | after package |
| 10 | E-MTAB-13687 tonsil atlas (2.1 TB) | never run | scoped only | size; role | negative-control subset (`DSR-07`) | after package |
| 11 | `pbmc_10x_healthy_v3` | never run (correct) | — | G3 + G4 | `VAL-04` on the frozen release | after freeze only |
| 12 | 4 synthetic datasets | not generated | — | VAL-01 → G3 → G4 | build generator (`VAL-01/02/03`) | after G3 |
| 13 | reagent / empty-droplet controls | — | unselected | no source | select before freeze or report unsupported | before freeze |
| 14 | anellovirus orthogonal positive | — | none exists | no sample | keep searching; screening-only claim until found | open |

Unusable, no run plan:
- **GSE208653** (cervical, 518 GB): no per-sample HPV truth.
- **CELLxGENE HPV collection**: processed matrices only, no viral features, no raw reads.
- **CELLxGENE CMV / Sound Life (GEO GSE271896, dbGaP phs003841)**: raw FASTQ controlled; GEO has H5 only; latent CMV ≈ 0 per sample.

## Shared run spec

Every run below fills these fields. Defaults unless the dataset block says otherwise.

| Field | Default |
|---|---|
| ViralScan version | the round's pinned commit (today `3.0.0.dev1 @ 12ceac7`); `code_sha.txt` in every run dir; clean tree only |
| Reference | cat42d: `/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42d/build/{panel.idx,panel.t2g}` + viral-only GTF extracted once from its `combined.gtf` (`DSR-01`) |
| Host strategy | combined (no `--host-filter`); STAR two-step as a second arm only where named |
| Cell calling | `--cell-calling emptydrops` (knee removed, A7) |
| Read filter | off; `--read-filter artefact` as a second arm where named |
| Strand | measured, not assumed (`--strand` per 1 M-read test) for 5′ libraries |
| Geometry (`-x`) | measured: R1 length + `viralscan check-whitelist` match rate |
| Cores | 4 per job (F-024: CPU efficiency 25–40 % at 8–16) |
| Cost | combined ~1.3 core-h, two-step ~5.6 core-h per 50 M pairs (F-024) |
| Inputs | checksums verified against `source_urls.tsv`; public positives are `recorded-reverify` |
| Off-target check | every call other than the expected virus with ≥ 3 molecules gets read validation (`DSR-02`) |
| Outputs | run dir + `run_manifest.json`, `sacct` line, and a one-line result in `.living/findings/` |

## Standard arm matrix (`DSR-10`, user 2026-10-06)

Every dataset gets every arm "whenever possible". This overrides "second arm
only where named" in the table above. Common to every arm: cat42d, the DSR-01
GTF, measured `-x`/strand and emptydrops. A cell is `todo`, a job id, `done`
or `n/a: <reason>`. Never leave a cell blank.

| # | Dataset | combined, filter off | combined, `artefact` | two-step (virus-only cat42d, `DSR-11`) | `DSR-02` off-target |
|---|---|---|---|---|---|
| 1 | EBV SRR12682296 | todo | todo | todo | todo |
| 2 | HHV-6B SRR20710641 | todo | todo | todo | todo |
| 3 | HSV-1 SRR8315713 | todo | todo | todo | todo |
| 4 | HPV16 GSE189670 (×2 rafts) | todo | todo | todo | todo |
| 5 | COVID x213 / x216 | todo | todo | todo | todo |
| 6 | SFL tonsil x223 | pilot 25702321 | after pilot | after `DSR-11` | todo |
| 7 | KSHV GSE190558 | todo (manifest first) | todo | todo | todo |
| 8 | KSHV+EBV GSE154900 | todo (cap first) | todo | todo | todo |
| 9 | HNSCC GSE164690 pilot | todo | todo | todo | todo |
| 10 | Tonsil atlas E-MTAB-13687 | todo | todo | todo | todo |

## Per-dataset blocks

### 1. EBV SRR12682296 — `RUN-02`
- Done: combined, two-step, corrected runs on `ebv_latest_ref_2026-09-27/full_panel`, knee. HHV-6B call = telomere repeats (F-025).
- Missing: cat42d reference, emptydrops, read-filter arm.
- Run: `-x 10xv2`, combined, two arms (`--read-filter` off / artefact). 127 M pairs → ~3.3 core-h per arm.
- Accept: EBV dominant; latent programme (LMP-1, EBNA) present; no off-target virus survives `DSR-02`.
- Tie-in: confirms whether the cat42d GRCh38 D-list removes the telomere sink (F-025).

### 2. HHV-6B SRR20710641 — `RUN-01`, `DSR-05`
- Done: combined (forward), 5′ forward/reverse/unstranded variants, knee. HHV-6A residual under read validation (`DSR-05`, job 25701668).
- Missing: measured geometry (protocol says 10x v3; the 5′ runs used `-x 10xv2`); strand choice; two-step (job 25652882_1 failed, cause unknown).
- Steps:
  1. Read the 25652882_1 log for the failure cause.
  2. 1 M-read test: R1 length, `check-whitelist` against v2 and v3 lists, host mapping rate under forward / reverse / unstranded.
  3. Full run with the measured `-x` and `--strand`, combined + two-step. 115 M pairs → ~3 core-h combined, ~13 core-h two-step.
- Accept: HHV-6B dominant; HHV-6A flagged as bleed or explained by `DSR-02`.

### 3. HSV-1 SRR8315713 — `RUN-03`, `PROG-16`
- Done: combined + two-step on the old ref; F-018 shows kb's data-derived allowlist drops 28 % of HSV-1 molecules.
- Blocker: SW-21 / MECH-D decision on Drop-seq barcode correction.
- Run after SW-21: `-x dropseq`, combined, emptydrops. 100 M pairs → ~2.6 core-h.
- Accept: HSV-1 productive programme (PROG-16: no latent calls); HHV-2 stays flagged as bleed.

### 4. HPV16 GSE189670 — `HPV-11`, `DSR-03`
- Done: cat42, cat42b, cat42d runs with emptydrops. On cat42d the Gammatorquevirus poly-A sink is gone; HPV16 7,726 vs 8.
- Open: HPV118 (18 / 12 molecules) and HPV29 (11 / 13) appear in **both** the infected and the parental raft, so they are not HPV16 cross-mapping alone.
- Run: rerun both rafts on the target commit (cat42d, emptydrops, `-x 10xv3`), then `DSR-02` read validation on HPV118/HPV29 in both rafts.
- Accept: each off-target HPV is explained (host, low complexity, HPV16 shared k-mers) or kept as a real finding with read evidence.

### 5. COVID x213 / x216 — `DEF-07`, `CAT-30`, `DSR-04`
- Done: forward, reverse, unstranded runs on the old ref; SARS-CoV-2 = 0; anellovirus = poly-G artefact (F-019), still dominant under reverse.
- Run: cat42d, `--strand reverse` (F-020), two arms read filter off / artefact, emptydrops.
- Accept: SARS-CoV-2 stays 0; anellovirus collapses on cat42d + artefact, or the survivors go to `DSR-02`.
- Blocker for publishing anything: DEF-07 retraction wording needs the user's sign-off.

### 6. SFL tonsil x223 — `TONSIL-02`
- 2026-10-06, released by the user (`DSR-12`). R1 is 28 bp and R2 90 bp; 1,605,890,600 pairs; cellranger reports 95 % valid barcodes. The whitelist is cellranger's raw barcode universe (3,491,354 barcodes). Code is pinned at `vs_pinned/4346dc8`.
- Pilot: job 25702321, first 4 M pairs, `viralscan_work/sfl_tonsil/pilot/`.
- Done: TONSIL-01 minimap2 screen, 0 TTV / 0 HPV (F-010).
- Steps:
  1. 1 M-read subsample: `-x 0,0,16:0,16,28:1,0,0`, cellranger barcodes as `-w`, host mapping under forward / reverse / unstranded (`--strand` now exists).
  2. `viralscan check-whitelist`.
  3. Full native run, cat42d, measured strand, emptydrops.
- Accept: native counts agree with TONSIL-01 (0 TTV, 0 HPV after validation); any call goes to `DSR-02` and per-donor attribution via cellhashr singlets.

### 7. KSHV GSE190558 — `RUN-04`, `PROG-16`
- 4 GSMs × 4 lane SRRs (SRR17180386–SRR17180401), iSLK.219 latent vs induced.
- Before any run: freeze the lane-concatenation manifest (URLs, md5, lane → GSM) in the protocol asset `GSE190558_fastq_manifest`.
- User step: download (ENA `filereport` for sizes first).
- Run: `-x 10xv3`, combined, emptydrops, one run per GSM after lane concatenation.
- Accept: induced > uninduced lytic programme; latent genes (LANA, vFLIP, vCyclin) in both.

### 8. KSHV + EBV GSE154900 — `RUN-05`, `PROG-18`
- 2 GSMs, 4 SRRs (SRR12287751–SRR12287754), BC-1 PEL (KSHV + EBV co-infected).
- Before any run: decide and freeze full depth vs read cap in `GSE154900_fastq_manifest`.
- User step: download.
- Run: `-x 10xv3`, combined, emptydrops, per GSM.
- Accept: KSHV and EBV both detected, mostly latent with a lytic minority (PROG-18); co-infection per cell reported, not assumed.

### 9. HNSCC GSE164690 — `DSR-06`
- 6 HPV+ / 12 HPV− tumours, 10x 3′ v2, about 1 TB.
- Pilot: one lane of HN18 CD45− (HPV+) vs HN01 (HPV−). User downloads the two lanes only.
- Run: `-x 10xv2`, combined, emptydrops.
- Accept: HPV16 (or the clinically typed HPV) in HN18 tumour cells, ~0 in HN01. A full run only if the pilot passes.

### 10. Tonsil atlas E-MTAB-13687 — `DSR-07`
- Role: specificity control (EBV+ cells expected ≈ 0–2 per sample).
- Pilot: every scRNA library is a hashed donor pool of ~100 GB, so the pilot takes one run each: the lowest-accession run of scRNA Samples 1, 2 and 3 (ERR13027010, ERR13027122, ERR13027027), chosen from metadata alone. Downloaded under `DSR-09`.
- Run: geometry measured per library, cat42d, emptydrops.
- Accept: no probable/strong call beyond rare EBV; every call goes to `DSR-02`.

### 11. `pbmc_10x_healthy_v3` — `VAL-04`
- Do not run before G3 and G4. Resolve FASTQ URLs and hashes now (metadata only, no quantification).

### 12. Synthetic datasets — `VAL-01`…`VAL-10`
- `synthetic_factorial`, `synthetic_host_only`, `synthetic_mixed_host_virus`, `synthetic_host_homology`.
- Chain: VAL-01 generator (design in `2026-10-05-val01-generator-design.md`) → G3 sign-off → G4 reference freeze → VAL-08 golden tiny panel → VAL-10 full panels.
- Open input: `rerun-multimap` cost is unmeasured and decides whether the defaults grid fits the 11 k core-h budget (grid-cost report, 2026-10-05). Time one `rerun-multimap em-cell` before VAL-10 sizing.

### 13–14. Optional controls
- Reagent and empty-droplet controls: select before freeze or report the rows as unsupported.
- Anellovirus orthogonal positive: no sample exists; anellovirus results stay screening-only.

## New tooling item

`DSR-02` — one parameterised read-validation script replacing the two hand-edited
`run.sh` copies (`hhv6b_in_ebv_2026-10-05/`, `hhv6a_in_hhv6b_2026-10-05/`).
Inputs: sample FASTQs, run dir, target transcript set, CB/UMI lengths. Steps:
EC lookup → bus records → read extraction (UMI exact, CB ≤ 1 mismatch) →
competitive minimap2 against GRCh38 + panel → per-read verdict (viral-best,
host-best, tie, repeat class).

`DSR-08` — record the git commit (`git describe --always --dirty`) in
`run_manifest.json` and `reference_provenance.json`, next to
`viralscan_version`, so a run's exact code is recoverable without a side file.

## Run log (fill per run)

| Dataset | Run dir | ViralScan version | Commit | Reference | Date | Result |
|---|---|---|---|---|---|---|
| | | | | | | |

## Suggested order (once experiments are released)

1. `DSR-08` (commit in the manifest) and `DSR-02` script, then `DSR-05` (HHV-6A result) and `DSR-03` (HPV16 rerun + HPV118/29 check).
2. `RUN-01` geometry/strand test and `TONSIL-02` 1 M-read test: small, unblock two datasets.
3. `RUN-02`, `DSR-04` reruns on cat42d: ~10 core-h total.
4. GSE190558 / GSE154900 manifests and cap decision, then the user downloads.
5. `DSR-06`, `DSR-07` pilots.
