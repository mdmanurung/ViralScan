# Protocol pending-section drafts for the G3 sign-off (Q83)

**Status (2026-10-13, later):** the user accepted all 25 USER CHOICE recommendations (section 10). The
pending-phase pastes are applied, re-based onto the post-#46 protocol, as `DEV-028`
(`analysis/v3_validation/deviations.yaml`), which lists each UC and whether it moved protocol text.
Not pasted from this draft: the `g3-freeze` stanza and the sidecar file (G3 is the user's, after SCI-05
round 13); the section 6 outcome definition (DEV-027 already defines it) and its last clause (it contradicts
Q86). The body below is the original proposal, left as written for provenance.
**Original status:** proposal. Nothing here was applied: `analysis/v3_validation/protocol.yaml`,
`deviations.yaml` and both schema copies were untouched.
**Base:** `origin/main` @ `88f02a1`, protocol v0.3.0, ledger DEV-000..DEV-026.
**Authority:** Q83 (user, 2026-10-13): an agent drafts the six `pending` sections from recorded
decisions and marks every choice that needs the user; the user signs the filled protocol at G3.
**Owner of the next step:** the user (choices), then whoever holds the open protocol edit
(DEV-027, digest scope, workflow coverage) to paste the blocks after that edit lands.

## 0. Read this first

1. **Five of the six sections already have complete bodies.** `partitions`, `calibration`,
   `defaults_selection`, `workflow_matrix` and `failure_and_deviation_reporting` are written
   (`protocol.yaml:528`, `735`, `933`, `1936`, `2292`) and only carry `status: pending`.
   What stops them freezing is not missing text. It is (a) the freeze flip and its co-flips
   (section 1), (b) values no human chose (the `A`/`N` rows below), and (c) round-12 findings.
   Re-pasting 800 lines of YAML would add nothing, so each section below gives a **delta**.
2. **`independent_review` has no body and the schema forbids one.** The top level is
   `additionalProperties: false` and has no such key (`schemas/v3/validation_protocol.schema.json`).
   The only protocol entry is `planned_freeze_sections.independent_review` (`protocol.yaml:2451`).
   Its content (reviewers, findings, digests, Git SHA) belongs in the "SCI-05 freeze sidecar"
   (`analysis/v3_validation/README.md:50`), for which no path or format is recorded.
   `schemas/v3/protocol_lock.schema.json` is the Track B discovery-study lock, **not** this sidecar.
3. **Every section schema is closed** (`additionalProperties: false`). New keys need a schema
   edit, which is out of scope. All blocks below only use existing keys.
4. **Q53-Q83 are not on `origin/main`.** They are on `origin/codex/trackb-gates-q53-q70`
   (cited as `trackb:` below). **Q71-Q76 are not in any ref of this repo**
   (`git log --all -G "Q7[1-6]([^0-9]|$)"` finds only two unrelated GTF/data commits and this
   branch; `HANDOFF.md`, `HANDOFF_def02_2026-10-03.md` and `SESSION_RESUME.md` in the shared
   checkout do not contain them either). They may be in the study repo or session notes. If they
   bear on these sections, they are missing here.
   **REG-05 does not exist** in any ref (`git log --all -G "REG-05"` and a grep of `PLAN.md`,
   `docs`, `.living`, `analysis` are empty), so the split decisions here come from D1/D1.1
   (D:1633-1643), VAL-01 D2-D6 (D:1645-1657) and the 2026-07-27 SCI-03 write-up (D:1058-1102).
   **VAL-04** is PLAN:5255 (the presumed-negative PBMC control); it raises UC-P4.
5. **Hashes are not written into the blocks.** Digests depend on `frozen_at`, on the other
   agent's pending edit, and on each other; decisions.md records three hand-copied digests that
   were wrong (`.living/decisions.md:1165`). Blocks carry `contract_sha256: RECOMPUTE`;
   section 9 gives the recompute command. All validation below used the validator's own functions.
6. **Top-level `status` stays `draft` after G3.** The training-phase checks
   (`validate_v3_protocol.py:1371-1413`) also need REL-03 pins, REF-09, VAL-09 and `data_hashes`,
   none of which the six sections can close. See UC-I4.

Provenance classes used in every source table:

| Class | Meaning |
|---|---|
| **U** | user decision, recorded (decisions.md / PLAN WP1E table / Q-number) |
| **A** | agent or reviewer-round decision, recorded, **no user sign-off found** |
| **N** | no recorded source: becomes a `USER CHOICE` |

Line numbers are for `origin/main` @ `88f02a1` unless prefixed `trackb:`.
`D` = `.living/decisions.md`, `P` = `PLAN.md`, `Y` = `analysis/v3_validation/protocol.yaml`.

## 1. What G3 flips, and what the validator enforces

| Flip | Enforced by |
|---|---|
| `partitions.status: frozen` + `frozen_at` | `_validate_sci03` (`validate_v3_protocol.py:103-146`) |
| ... needs `seeds.root/split/cell_calling/evidence_sampling/bootstrap` frozen | `:150-155` (already frozen, `Y:1641-1676`) |
| ... needs `viral_abundance`, `host_virus_homology`, `chemistry` frozen with levels | `:189-198` |
| `calibration.status: frozen` needs `viral_abundance` frozen with at least 7 levels | `:204-217`, `MINIMUM_ABUNDANCE_LEVELS=7` |
| `workflow_matrix` / `failure_and_deviation_reporting` frozen + `frozen_at` | `_validate_sci04` (`:250-265`) |
| each flipped section's planned entry flips in step | `:122-129`, `:256-263` |
| each frozen section needs a ledger chain ending at its current digest | `validate_amendment_ledger` (`:1192-1217`) |
| `defaults_selection`, `execution_readiness`, frozen factors sit under `frozen_inputs` | `frozen_inputs_sha256` (`:1003-1023`) |
| ... so changing them needs a DEV record **and** a `recorded_deviations` entry | `:1179-1191`, `:1192-1217` |

Co-flips for the G3 commit (all in block `g3-freeze` below, plus two list edits the block
cannot express): flip the 8 pending factors; delete the closed blockers `partitions_metrics`,
`workflow_rows`, `independent_review_round_2`, `independent_review` from
`execution_readiness.training_blockers` (`Y:2343-2388`); append one ledger record per frozen
section and one for `frozen_inputs`. `workflow_matrix.environment_pinning` stays `pending`
(REL-03), as the validator permits at the draft phase.

## 2. Partitions

**Schema requirements** (`$defs/partitions`): required `status, resolution_task, unit,
stratification_factors, algorithm, sort_key, holdout_fraction, leakage_prohibitions,
outcome_blinding, contract_sha256`; `holdout_fraction` in (0,1); `minimum_samples_per_stratum >= 4`;
`generated_sample_structure` closed (`cells_per_sample, reads_per_cell, technical_replicates_per_sample,
rationale, file_layout`); `frozen_at` format `date`. Cross-checks: section 1.

**Delta:** none forced. No recorded decision changes any value in this section. Only the freeze
flip (block `g3-freeze`) and the confirmations UC-P1..P3.

**Sources**

| # | Value | Where | Class | Source |
|---|---|---|---|---|
| P1 | unit = background donor library (real background), generated sample (synthetic) | `Y:532-542` | U | D:1633-1641 (D1, D1.1) |
| P2 | at least 10 donor libraries per chemistry (3 in holdout at 0.3) | `Y:541-542` | U | D:1641 |
| P3 | per-dataset reading of D1.1 (real background: library; synthetic: sample) | `Y:538-542` | A | D:1661 ("B2 reading") |
| P4 | stratify by viral_abundance, host_virus_homology, chemistry | `Y:543-546` | A | `Y:601-620`; SCI-03 round 1; no user row |
| P5 | `reads_per_cell: 25000` | `Y:549` | U | D:1647 (D5) |
| P6 | `cells_per_sample: 2000`, `technical_replicates_per_sample: 1` | `Y:548,550` | A | round 8 R8-F4 (`sci05_review_round8.yaml:261`); `Y:551-565` |
| P7 | sibling pairs HHV-6A/6B, HSV-1/2, EBV/KSHV | `Y:577-580` | A | `sci05_review_round8.yaml:70,479`; `Y:581-597` |
| P8 | `one_member_only` sibling regime, seeded draw | `Y:591-592` | U | D:1652 (D2) |
| P9 | `minimum_samples_per_stratum: 4` | `Y:621` | A | R7-F2 (`sci05_review_round7.yaml:147`) |
| P10 | largest-remainder apportionment, cap `ceil(N*f)`, ties by stratum key | `Y:650-669` | A | D:1071-1077; round 1 |
| P11 | byte-wise ASCII sort key | `Y:676-681` | A | round 1; `Y:676-681` |
| P12 | single-sample stratum goes to training and is reported | `Y:660-662,691-695` | A | D:1074-1076 |
| P13 | `holdout_fraction: 0.3` | `Y:682` | A | D:1078-1081; P:6308-6309 (no user row) -> UC-P1 |
| P14 | `holdout_evaluations_allowed: 1` | `Y:732` | A | D:1078-1081 -> UC-P1 |
| P15 | only `training-and-holdout` datasets are split; context-only rows exempt | `Y:683-695` | U | D:1847 (grill 2) |
| P16 | "template/locus" = host and challenge sources; viral genomes shared | `Y:698-704` | U | D:1653 (D4) |
| P17 | reference, D-list and feature universe built from training-eligible information only | `Y:714-720` | A | round 1 F-series; `Y:2460-2466` |
| P18 | blinding: owner-only truth outside repo, sha256 only, opaque read names, shuffled order, unblinding needs written sign-off | `Y:721-731` | U | D:1654 (D6); P:2000 (R3.4) |
| P19 | `seeds.split: 20260727001` | `Y:1659` | A | frozen 2026-07-27, D:1058 |
| P20 | `pbmc_10x_healthy_v3` is excluded as a planting background and is the external-evaluation negative | `Y:1180-1184,1298-1305` | U, **conflicting** | D:1639 (D1, 2026-10-05) says excluded; P:27 (grill 2, 2026-10-12) says "the D1 background and T4 negative" -> UC-P4 |

20 of 20 rows have a recorded source (U: 8, A: 12; P13 and P14 are A-only and go to UC-P1; P20's two
sources disagree and go to UC-P4).

```yaml pending-edit
# Partitions: no edit is forced by any recorded decision. Deliberately empty.
{}
```

**USER CHOICE**

- **UC-P1 holdout fraction 0.3 and exactly one holdout evaluation.** Only the agent's 2026-07-27
  design write-up records them (D:1078). Options: (A) keep 0.3 / 1 (10 libraries per chemistry
  gives exactly 3 held out, the arithmetic at `Y:541-542`); (B) 0.2 (2 libraries); (C) 0.4
  (4 libraries, fewer training libraries for tuning). **Recommend A.**
- **UC-P2 generated sample structure.** 25,000 reads per cell is yours (D5). 2,000 cells and one
  technical replicate came from review round 8. Options: (A) keep; (B) 1,000 cells (halves the
  50 M-pair sample cost used in F-024, but a 1 % infected fraction then gives about 10 cells);
  (C) two replicates (the bootstrap forbids treating them as independent, so only cost rises).
  **Recommend A.**
- **UC-P3 sibling pairs.** Three pairs were fixed by review rounds 7-8, not by you. Options:
  (A) keep; (B) drop EBV/KSHV (the most distant pair); (C) add a fourth pair. **Recommend A**; the
  `Y:581-597` rationale (a sibling metric with only near-identical pairs cannot show where
  discrimination starts) is sound.
- **UC-P4 is `pbmc_10x_healthy_v3` a planting background?** The protocol and D:1639 exclude it
  because it is the VAL-04 external-evaluation negative (`Y:1180-1184`). The grill-2 summary
  (P:27) calls it "the D1 background and T4 negative". If "D1" there means the VAL-01 real-PBMC
  background, the same donor library would be both a planted-and-twinned background (so inside
  the partition) and an external negative, a leakage path `partitions.leakage_prohibitions` does
  not cover. If "D1" means the M-04 real-PBMC mode (D:1825), there is no conflict. Options:
  (A) keep it external-only, as the protocol says, and source the 10x v3 backgrounds from other
  donors (the >= 10 libraries per chemistry rule then needs 10 other 10x v3 libraries);
  (B) make it one background library and drop its external-evaluation role (changes a
  `required` dataset and its endpoints' denominators); (C) use it for both with the twin arm
  alone scored as the negative (not an independent sample). **Recommend A**; please confirm what
  grill 2 meant.

## 3. Calibration

**Schema requirements** (`$defs/calibration`): required `status, resolution_task, threshold_search
{data, grid, objective, prohibited}, tie_breaker, metrics[{id, endpoint_id, definition,
denominator_rule}], limit_of_detection{endpoint_id, estimator, response, reported_quantity,
extrapolation}, uncertainty{resampling_unit, prohibited_units, replicates, interval, stratified_by},
contract_sha256`. Validator: every metric and the LOD name a real endpoint (`:157-174`); `viral_abundance`
frozen with at least 7 levels before freeze (`:204-217`).

**Delta:** none forced by a recorded decision. Two **conflicts** between recorded user decisions
and the draft are routed to UC-C3, UC-C5 and UC-C6 instead of being silently rewritten.

**Sources**

| # | Value | Where | Class | Source |
|---|---|---|---|---|
| C1 | threshold grid: molecules {1,2,3,5,10} x features {1,2,3} x support {0.5,0.75,0.9} | `Y:743-747` | A | D:1082-1087; round-12 M3: no code reads it -> UC-C1 |
| C2 | objective: maximise F1 | `Y:748-752` | U | P:1998 (R3.2 "SCI-03 F1") |
| C3 | hard constraint: zero probable-or-strong calls in host-only and homology negatives | `Y:748-752` | U | P:1989 (R2.3) |
| C4 | no eligible grid point = calibration fails, no tier | `Y:758-765` | A | D:1082-1087; conflicts with R3.2 -> UC-C6 |
| C5 | 1-SE rule, most conservative point | `Y:766-771` | U | P:1987 (R2.1); ordering rule A: D:1088-1091 |
| C6 | SE from the resampling-unit bootstrap | `Y:772-779` | A | round 1 |
| C7 | metrics M1-M9 and their denominators | `Y:780-858` | A | endpoints `Y:1089-1159`; rounds 1-5 |
| C8 | LOD: probit, one fit per virus and chemistry, LOD95, no extrapolation | `Y:859-892` | A | D:1092-1096; R9-F1/DEV-016 (`sci05_review_round10.yaml:126`) -> UC-C3 |
| C9 | at least 7 abundance levels, bracketing the knee | `Y:874-881` | A | R5-F4 (`sci05_review_round5.yaml:218`); `validate_v3_protocol.py:22-26` |
| C10 | abundance levels 1,3,10,30,100,300,1000 | `Y:1546` | U | D:1646 (D3) |
| C11 | 2,000 bootstrap replicates, BCa 95 % | `Y:907-910` | A | P:6305-6313 (2026-07-27); no user row -> UC-C2 |
| C12 | at least 3 resampling units for an interval; else not estimable | `Y:911-919` | A | `sci05_review.yaml:111-116` (round 1) |
| C13 | resampling unit = donor library, cluster bootstrap | `Y:894-906` | U | D:1641 (D1.1) |
| C14 | stratify by chemistry and viral_abundance where carried | `Y:920-927` | A | round 1 |

14 of 14 rows have a recorded source (U: 5, A: 9).

**Finding for UC-C1.** `grep -rli probable src scripts` returns only `scripts/v3_score.py`. No
shipped code assigns `probable` or `strong`; `docs/output_reference.md:760` says they "require
calibrated read-level evidence and are never assigned from molecule counts alone". The grid
therefore tunes a tier rule that has no implementation, and R2.3's "zero probable/strong calls"
constraint binds defaults that cannot yet emit those tiers.

```yaml pending-edit
# Calibration: no edit is forced. Deliberately empty.
{}
```

**USER CHOICE**

- **UC-C1 the 45-point grid.** Agent judgement (D:1082-1087), unconsumed by code, and the tier
  producer it tunes does not exist (finding above). Options: (A) keep as preregistration and
  require the tier implementation to expose exactly these three knobs before G3 (the grid is a
  contract on that code); (B) delete the grid until the tier code exists and freeze calibration
  later; (C) replace with knobs read from the implemented tier rule. **Recommend A**, because
  the grid must be fixed before outcomes either way, and B delays G3 on unwritten code.
- **UC-C2 bootstrap size.** 2,000 replicates, BCa, 95 %: agent only. Options: (A) keep;
  (B) 1,000; (C) 10,000 (cost is negligible at this sample count, but the number is arbitrary).
  **Recommend A.**
- **UC-C3 LOD grain.** Draft: one probit per virus within chemistry (round 9-10). You recorded
  "LOD95 per genus reported, not pass/fail" (grill 2, P:27). Options: (A) keep per-virus for the
  synthetic estimator and report a genus roll-up descriptively outside the frozen estimator;
  (B) change the frozen estimator to per-genus (pools viruses that share no detection curve,
  the defect round 9 found); (C) both frozen. **Recommend A**; grill 2's sentence reads as a
  reporting rule for real data, not for the synthetic probit. Please confirm.
- **UC-C5 claim-gate minimums and scored T2 rows (cross-cutting, also `workflow_matrix`).**
  Recorded in grill 2 (D:1847, 1850) and absent from the protocol: at least 3 donor-distinct
  libraries for any T2/T3 sensitivity claim, at least 10 independent T4 libraries before
  "specificity assessed", HPV16 raft pair and iSLK.219 become scored T2 rows. Round-12 B1 found
  the protocol has no T2 class (`analysis_role` and `truth_status` have no value for it, and
  `harmonization.cell_universe.public_positive_role` is a schema `const`), and Q80 settled
  only V03 (trackb D:1874). Options: (A) keep these in the claim registry and Track B lock,
  protocol unchanged, T2 rows `context-only`; (B) add E11/D26 and a `design-known` truth status
  (schema + validator + DEV-027 scope, a large change); (C) drop the scored-T2 plan. **Recommend A**
  for G3; revisit B only if a T2 claim is meant to be preregistered.
- **UC-C6 no-eligible-grid-point policy.** `Y:758-765`: calibration fails and no tier is reported.
  Your R3.2 (P:1998) says "no passing option -> most conservative, reported 'screening only'" for
  every objective including "SCI-03 F1 (thresholds, tiers)". The protocol text (2026-07-27)
  predates R3.2 (2026-09-29). Options: (A) keep "no calibrated tier" for tiers, "screening only"
  for defaults (a zero-false-call violation makes any tier untrustworthy); (B) apply R3.2 to tiers:
  adopt the most conservative point (10 molecules, 3 features, 0.9) labelled screening only;
  (C) two-level: B but only if the constraint is violated by fewer than a stated count (needs a
  number, so not recommended). **Recommend A.**
- **UC-C4 acceptance scope.** `H2`/`E2` fail on any probable-or-strong call in synthetic exact
  negatives (`Y:1046-1052`, `1106-1111`). Grill 2 says "no numeric specificity pass/fail" (D:1850).
  Options: (A) confirm grill 2 means presumed-negative real data only, so the exact-negative
  zero-call rule stays (recommended; the calibration hard constraint cannot be relaxed without
  defeating its purpose); (B) extend grill 2 to the synthetic negatives and demote E2 to
  "report the upper bound with n". Endpoints are hashed by `harmonization.dependent_fields_sha256`,
  which the schema pins as a `const`, and again under `frozen_inputs`, so B means both schema
  copies plus two ledger scopes (the round-12 B1 mechanism), a large change. **Recommend A.**

## 4. Defaults selection

**Schema requirements** (`$defs/defaultsSelection`): required `status, resolution_task,
hard_constraint, tie_breaker, defaults[{id, objective, conservative_direction, grid}], grid_rule`.
When `status: frozen`, every default needs a non-empty `grid`. `frozen_at` is a plain string here (no
`date` format). There is no `contract_sha256`: this section is digested under `frozen_inputs`
(`validate_v3_protocol.py:1016-1018`), so any edit needs a DEV record.

**Delta:** one candidate-structure paragraph in `grid_rule` (UC-D1, recommended option A). The
grids themselves are complete (DEV-021..023) and are not touched.

**Sources**

| # | Value | Where | Class | Source |
|---|---|---|---|---|
| D1 | frozen at G3, before any tuning; no tuning before G4 | `Y:937-939` | U | P:1987, 1990 (R2.1, R2.4) |
| D2 | training partition only; already-seen data excluded | `Y:940-942` | U | P:1988 (R2.2) |
| D3 | hard constraint: zero probable and strong calls in training negatives | `Y:943-946` | U | P:1989 (R2.3) |
| D4 | 1-SE conservative tie-break, SE from the partitions.unit bootstrap | `Y:947-952` | U | P:1987 (R2.1) |
| D5 | no passing option: most conservative, "screening only" | `Y:953-956` | U | P:1998 (R3.2) |
| D6 | per-chemistry defaults allowed, selected on that chemistry's units | `Y:957-960` | U | P:1993 (R2.7) |
| D7 | 20,000 core-h; 55 / 25 / 10 / 10 % split | `Y:961-965` | U | P:1983 (Q10d), 1997 (R3.1) |
| D8 | host_strategy grid [combined, exact-fragment-star-two-step], recall | `Y:967-974` | U | P:1976 (Q7); D:1669-1670; deviations.yaml:723 (DEV-021); objective P:1998 |
| D9 | star_alignment grid [pinned, star-default]; star-default conservative | `Y:975-989` | U | D:1679-1682 (confirmed 2026-10-05); deviations.yaml:747 (DEV-022) |
| D10 | read_artefact_filter grid [off, artefact]; CLI default off | `Y:990-1006` | U | D:1688-1690; deviations.yaml:770 (DEV-023); definition A: D:1696-1700 (DEV-024) |
| D11 | multimap grid host-conservative, unique-weighted, em-global, em-cell; `equal` context only | `Y:1007-1019` | U | D:1674; objective P:1998 (R3.2) |
| D12 | EM max-iter 100, tolerance 1e-6 | `Y:1012-1013` | A | shipped defaults; D:1674 says "at their defaults" |
| D13 | cell_calling grid [emptydrops], fdr 0.01, lower 100, niters 10000, not tuned | `Y:1020-1030` | U | D:1673; values frozen in `Y:132-139` |
| D14 | within-sibling D3 variant only after the counting contract is amended | `Y:1016-1019` | U | P:1991 (R2.5); D:1674 |
| D15 | `conservative_direction` of each default | `Y:969,977,992,1009,1022` | A | star-default confirmed U (D:1682); the others agent |

15 of 15 rows have a recorded source (U: 13, A: 2).

```yaml pending-edit
defaults_selection:
  grid_rule: >-
    Each grid must hold explicit candidate values before G3; a default with an
    empty grid cannot be frozen and cannot be tuned. DEV-021 fills host_strategy,
    multimap_allocation and cell_calling, DEV-022 fills star_alignment, and
    DEV-023 fills read_artefact_filter. Candidate structure: host_strategy,
    star_alignment and read_artefact_filter share one objective, recall of
    planted truth molecules, and are tuned jointly. One candidate is one
    combination of the three. star_alignment applies only under
    exact-fragment-star-two-step, so the candidates are the 2 combined-host-plus-virus
    combinations (read_artefact_filter off or artefact) and the 4 two-step
    combinations (2 star_alignment points by 2 read_artefact_filter points), 6 in
    all. multimap_allocation is scored on every candidate's kb count output at no
    quantification cost, and cell_calling has one point. Ties inside the 1-SE
    band are broken by host_strategy first, then star_alignment, then
    read_artefact_filter, each toward its conservative_direction.
```

**USER CHOICE**

- **UC-D1 candidate structure.** The section says each default has its own grid and objective but
  never says whether candidates are combinations. `star_alignment` only exists under two-step, and
  the filter changes the STAR input, so they interact. Options: (A) joint over the three
  recall-objective defaults, 6 candidates, multimap scored on all of them (block above);
  (B) coordinate descent in the listed order (cheaper, order-dependent, a second preregistered
  rule needed); (C) full cross including the 4 multimap modes as candidates (24 candidates; the
  multimap objective differs, so it would need its own tie-break). **Recommend A.**
- **UC-D2 grid cost against Q10d.** P:128 lists "check the grid cost against Q10d" as still open.
  Arithmetic from recorded figures, not a decision: F-024 gives 1.3 core-h (combined) and 5.6
  (two-step) per 50 M pairs; the filter pass is about 1.15 plus a 1.3 kb count (`Y:1004-1006`).
  Candidates A cost 1.3 + 2.45 + 5.6 + 5.6 + 2 x 2.45 = **about 19.9 core-h per sample**, up to about
  **52** if each of the 18 extra multimap reruns cost a full combined run (1.8; unmeasured, F-024).
  Shared defaults tune once on 10x v3 (R3.1). Taking the real-background factorial only,
  7 training libraries x 28 samples (14 cells, planted plus twin, `Y:640-645`) = 196 samples, i.e.
  **about 3.9k core-h (up to 10.2k)** against the shared share of 11k (55 % of 20k). The upper bound is
  tight and excludes the synthetic-background datasets. Options: (A) accept, measure the multimap
  rerun cost on the pilot panel first and re-check; (B) restrict the shared tuning set to a stated
  subset of libraries; (C) raise Q10d. **Recommend A.**
- **UC-D3 R2.5 within-sibling D3 variant.** Options: (A) defer to a later release (the counting
  contract and its conservation claim must be amended first, `Y:1016-1019`); (B) amend now and add it
  as a fifth multimap point. **Recommend A.**
- **UC-D4 which defaults are chemistry-specific.** R2.7 allows them and R3.1 budgets 25 % plus
  10 % confirmation, but no default is declared chemistry-specific. Options: (A) none at G3: tune
  once on 10x v3, confirm on the other two chemistries, promote only by amendment; (B) preregister
  `read_artefact_filter` per chemistry (the TSO reagent shape is 10x-specific); (C) every default
  per chemistry. **Recommend A.**

## 5. Workflow matrix

**Schema requirements** (`$defs/workflowMatrix`): required `status, resolution_task,
primary_comparison_rules, required_row_record_fields, environment_pinning{status, resolution_task,
rule}, workflows[{id, tool, mode, role, reference_id, dataset_ids}], no_silent_algorithm_patch,
expected_row_count, contract_sha256`. Validator: every `native-published` arm is claimed by exactly
one `matched-accession-index` arm of the same tool; `expected_row_count` equals the enumerated rows
(90 today); at the training phase every row needs `tool_version` and `container_digest`
(`:287-388`).

**Scope limit.** I do not touch `workflows[].dataset_ids`. Workflow coverage (round-12 M5, the
zero-row `synthetic_mixed_host_virus`) and V03 belong to the open protocol edit.

**Delta:** one sentence on version policy in `environment_pinning.rule` (UC-W2). Pin values are data
for REL-03, not choices.

**Sources**

| # | Value | Where | Class | Source |
|---|---|---|---|---|
| W1 | rows enumerated in the matrix, expected 90 (13 workflows x 6 datasets + 12) | `Y:2007-2193,2289` | A | D:1104-1131 (SCI-04), D:1128 (count check) |
| W2 | both host strategies are candidates: ViralScan combined and two-step | `Y:2008-2037` | U | P:1976 (Q7) |
| W3 | comparators STARsolo, traditional subtraction (with/without CB-UMI), Venus, Viral-Track, VIRTUS | `Y:2038-2166` | A | D:1104-1131, D:1133-1170 |
| W4 | every comparator has a native and a matched-accession arm | `Y:2086-2193` | U | P:1999 (R3.3 "both index arms"); fairness fix D:1133-1170 (A) |
| W5 | ViralScan 2.2.0 is a comparator, both arms, `counts_original` layer | `Y:2167-2233` | A | D:1241-1248 (you asked; the agent decided); round 5 |
| W6 | primary comparison rules (same FASTQ sha256, same release, anchor, audited intersection, unique-only) | `Y:1940-1956` | A | harmonization (frozen, SCI-02) |
| W7 | 17 required row-record fields | `Y:1966-1983` | A | P:6333 |
| W8 | kallisto two-step excluded with revisit condition | `Y:2245-2252` | A | D:1117-1122 |
| W9 | native-versus-matched claim scope, breadth and chemistry-confound disclosures | `Y:1957-1965,2266-2288` | A | rounds 2-5 |
| W10 | comparator compute is a separate budget of about 10k core-h | (not in the section) | U | P:1999 (R3.3) |
| W11 | `environment_pinning` stays pending until REL-03; the matrix freezes which rows exist | `Y:1984-2006` | A | D:1109-1116 -> UC-W1 |
| W12 | `execution_gate` declared, `status: unimplemented`, never called by the CLI | `Y:1992-2006` | A | D:1172-1217 |

12 of 12 rows have a recorded source (U: 3, A: 9).

```yaml pending-edit
workflow_matrix:
  environment_pinning:
    rule: >-
      Every workflow must carry an exact tool version and container digest before
      any outcome-generating run. The matrix below freezes which rows exist and
      what must be pinned; REL-03 supplies the values. A row whose environment is
      unpinned may not execute. Version policy: each comparator is pinned at its
      latest tagged release on the REL-03 pin date, never a branch head, unless
      that release changed the published defaults the row measures, in which case
      the release named in the tool's publication is pinned and the substitution
      is recorded in a deviation record. ViralScan 2.2.0 is fixed.
```

**USER CHOICE**

- **UC-W1 freeze the matrix before the pins exist.** `tool_version` and `container_digest` are
  inside each workflow item, so they are inside the frozen digest and every later fill is an
  amendment (DEV record, `outcome_triggered: false`). The agent chose this on 2026-07-27
  (D:1109-1116). Options: (A) freeze at G3, fill pins by one REL-03 amendment (recommended;
  matches the recorded design and the validator's training-phase pin gate); (B) hold the matrix
  until REL-03 so it freezes once with pins, delaying G3 on packaging; (C) split the pin fields out
  of the digest (schema change). **Recommend A.**
- **UC-W2 comparator version policy.** No recorded decision. Options: (A) latest tagged release at
  the pin date (block above); (B) the version in each tool's paper; (C) chosen per tool in a DEV
  record. **Recommend A**, so comparators are measured as users run them today.
- **UC-W3 the standard-arm-set rule against matrix breadth.** Your 2026-10-06 rule: the same
  analyses on every dataset whenever possible, an arm that cannot run is recorded n/a (D:1714-1725).
  `W2` (two-step, ViralScan's own arm) runs on 6 of 12 datasets and has no row for
  `pbmc_10x_healthy_v3`, the control datasets or the KSHV and anellovirus sets (`Y:2026-2037`).
  The comparator rows are deliberately six. This overlaps the open workflow-coverage edit; I only
  flag it. Options: (A) extend W2 to every dataset where it can run and record n/a reasons for the
  rest; (B) keep six and state in the matrix that the rule applies to exploratory rounds only;
  (C) leave. **Recommend A** (about 5.6 core-h per 50 M pairs per row, F-024).

## 6. Failure and deviation reporting

**Schema requirements** (`$defs/failureAndDeviationReporting`): required `status, resolution_task,
row_failure_policy_ref` (const `harmonization.row_failure_policy`), `required_failure_fields,
retry_lineage, deviation_record{required_fields, disclosure_rule}, amendment_rule, prohibited,
contract_sha256`. The row statuses and comparability semantics live in the frozen
`harmonization.row_failure_policy` (`Y:459-491`), not here.

**Delta:** a definition of "outcome" appended to `deviation_record.disclosure_rule` (UC-F1,
recommended option A). Round-12 M6 found "outcome" undefined and every DEV-020..026 record
blanket `outcome_triggered: false` (report at the SCI-05 round-12 scratchpad, not committed). Q79 (trackb D:1873) assigns M1-M5 and the Minors to agents and does **not** mention M6.

**Sources**

| # | Value | Where | Class | Source |
|---|---|---|---|---|
| F1 | failure semantics inherited from `harmonization.row_failure_policy` | `Y:2296,459-491` | A | SCI-02, frozen |
| F2 | 8 required failure fields | `Y:2297-2305` | A | `workflow_failures.tsv` columns `Y:497-498`; P:6292-6298 |
| F3 | every attempt is retained, a retry appends | `Y:2306-2309` | A | `Y:470-472` |
| F4 | 11 required deviation fields incl. `digest_scope`, `git_sha` | `Y:2310-2322` | A | D:1133-1217 (rails) |
| F5 | outcome-triggered deviations are disclosed as such | `Y:2323-2326` | A | precedent D:1033-1040 (the 1.1.0 highmem amendment) |
| F6 | frozen-section change = new digest + DEV record + full rerun of affected rows | `Y:2327-2330` | A | D:1172-1217 |
| F7 | 5 prohibitions (no row deletion, no replace-by-retry, no imputed accuracy, no unlabelled success subset, no unrecorded amendment) | `Y:2331-2336` | A | exclusion X5 `Y:1931-1934`; D:1133-1170 |
| F8 | resource-escalation precedent for exploratory rows: resubmit OOM/TIMEOUT once at 2x, cap 480 G, flags/code/reference locked | (not in the section) | U | D:1734-1739 (DSR-13), D:1766-1769 -> UC-F2 |

8 of 8 rows have a recorded source (U: 1, A: 7).

```yaml pending-edit
failure_and_deviation_reporting:
  deviation_record:
    disclosure_rule: >-
      A deviation prompted by observing an outcome is recorded with
      outcome_triggered true and disclosed as such in every report that depends
      on the affected rows. An outcome is any viral call, viral molecule or cell
      count, evidence tier, or per-virus metric produced by any tool on a protocol
      dataset. Host-only, reference-only, read-structure, chemistry and cost
      measurements are design measurements, not outcomes. A dataset on which an
      outcome was produced before the freeze commit is context-only (R2.2). Every
      such run is listed in a hash-chained explore ledger whose sha256 is
      recorded in the freeze sidecar, and a deviation record may describe no
      workflow row as outcome-free if a ledger lists a run of its dataset.
```

**USER CHOICE**

- **UC-F1 definition of "outcome" and the pre-freeze run inventory (M6).** Options: (A) the
  definition above, with the explore ledger (`scripts/explore_ledger.py`, P:5618) hashed into the
  sidecar and any pre-freeze outcome run making its dataset context-only (follows R2.2);
  (B) any run of any tool on a protocol dataset counts (stricter; PRE-01 QC and cost runs would
  count, and V03/CAL-11 could never be scored); (C) leave undefined and reword the DEV texts.
  **Recommend A.** It also needs the rule "no run of any V03 library before the commit that
  registers V03" in the study repo's CAL-11 gating.
- **UC-F2 retries for protocol rows.** The protocol retains every attempt but sets no limit or
  escalation rule. The exploratory-round rule is yours but scoped to exploratory runs (D:1734-1739).
  Options: (A) adopt it for protocol rows: only OOM/TIMEOUT may be resubmitted, once at 2x memory
  or time (cap 480 G), no change to flags, code, reference or seed, all attempts retained; (B) no
  automatic retry, a failed row stays failed; (C) unlimited retries with lineage. **Recommend A**,
  but it needs one new sentence in `retry_lineage` once you confirm (not drafted: no recorded
  decision applies it to scored rows).

## 7. Independent review

**Schema requirements.** No top-level key. The one entry is
`planned_freeze_sections.independent_review` = `{status, resolution_task, freeze_required,
description}` (`$defs/plannedSection`). `freeze.freeze_manifest` (string or null) is the sidecar
pointer, required only when the top-level `status` is `frozen`. **No schema exists for the
sidecar itself**, so the block below can only be parse-checked and digest-checked, not
schema-validated.

**What is recorded**

| # | Value | Where | Class | Source |
|---|---|---|---|---|
| I1 | review history: rounds 1-11 recorded, round 12 (independent agent) report written | `analysis/v3_validation/sci05_review*.yaml` | A | round-12 report in the scratchpad (not committed) |
| I2 | verdict flips only when a round passes clean with zero open Majors | P:2548-2551 | A | P:2548-2551 -> UC-I1 |
| I3 | G3 = schema-valid, no unresolved finding, hashed, no outcome job before the freeze commit | P:2559-2561 | A | P:2559-2561 |
| I4 | user signs at protocol freeze, defaults freeze, holdout unblinding | P:2001 | U | P:2001 (R3.5) |
| I5 | the user signs G3; round 12 = independent agent + Biomni | trackb `P:2569` | U | trackb D:1861-1868 (Q53, Q66) |
| I6 | agents fix M1-M5 and the Minors, then round 13; Biomni runs on round 13 | trackb D:1873 | U | Q79, Q82 |
| I7 | G3 freezes right after the user's sign-off | D:1850 | U | grill 2 |
| I8 | sidecar holds reviewer identity, findings, resolutions, protocol digest, Git SHA | `Y:2451-2455`, `README.md:50` | A | `Y:2451-2455` |
| I9 | a lock cites the sha256 of the protocol file bytes, not a re-serialisation | trackb D:1864 | U | Q59 (for the Track B lock; adopted here by analogy) |

9 of 9 rows have a recorded source (U: 5, A: 4).

```yaml pending-edit
planned_freeze_sections:
  independent_review:
    description: >-
      Reviewer identity, findings, resolutions, protocol and section contract
      digests, Git SHA, explore-ledger sha256 and the user's G3 sign-off, in the
      freeze sidecar named by freeze.freeze_manifest.
```

Sidecar skeleton (proposal for the file `freeze.freeze_manifest` would point to; field names follow
the `sci05_review_roundN.yaml` headers; every value below is a placeholder to fill at G3):

```yaml sidecar
schema_version: "1.0.0"
protocol_id: viralscan-v3-validation
gate: G3
frozen_at: "YYYY-MM-DD"
protocol_path: analysis/v3_validation/protocol.yaml
protocol_file_sha256: "sha256sum of the committed bytes"
ledger_sha256: "sha256sum of deviations.yaml"
git_sha: "base commit the freeze was written against"
section_contract_sha256:
  partitions: RECOMPUTE
  calibration: RECOMPUTE
  workflow_matrix: RECOMPUTE
  failure_and_deviation_reporting: RECOMPUTE
  frozen_inputs: RECOMPUTE
defaults_selection_sha256: "canonical sha256 of the section (it has no digest field)"
explore_ledger:
  path: "owner-held path"
  sha256: "sha256sum"
review:
  round: 13
  record: analysis/v3_validation/sci05_review_round13.yaml
  verdict: pass
  open_blockers: 0
  open_majors: 0
  reviewers:
    - {id: R-agent, identity: "independent reviewer agent, cold read", verdict: pass}
    - {id: R-biomni, identity: "Biomni external review (Q82)", verdict: pass}
user_signoff:
  name: mdmanurung
  date: "YYYY-MM-DD"
  statement: "I have read the filled protocol and sign G3."
```

**USER CHOICE**

- **UC-I1 pass bar.** Recorded bar: a round that passes clean with zero open Majors (P:2548-2551).
  Options: (A) zero open Blockers and Majors, Minors dispositioned (recommended, the recorded bar);
  (B) A plus zero open Minors; (C) allow a Major to be accepted by you with a written reason in
  the sidecar. **Recommend A.**
- **UC-I2 reviewers and tie-break.** Recorded: independent agent each round, Biomni on round 13
  (Q82). Not recorded: what happens when they disagree. Options: (A) any open Blocker or Major
  from either reviewer blocks G3 (recommended); (B) the agent review gates and Biomni is advisory;
  (C) you adjudicate disagreements case by case. **Recommend A.**
- **UC-I3 sidecar path and format.** Options: (A) `analysis/v3_validation/freeze_sidecar.yaml`,
  committed, `protocol_file_sha256` over the file bytes (Q59 convention), set as
  `freeze.freeze_manifest` (recommended); (B) owner-only outside the repo (nothing committed pins
  it); (C) fold it into `sci05_review_round13.yaml`. **Recommend A.** A schema for it is a separate
  task.
- **UC-I4 one freeze or two.** The top-level `status: frozen` needs every training blocker closed
  (REL-03 pins, REF-09, VAL-09, `data_hashes`, `seeds.generation`), which cannot happen by G3.
  Options: (A) two steps: G3 freezes the six sections, factors and the sidecar with top-level
  `status: draft`, and the top-level flip is a later ledgered commit once the blockers close
  (recommended; matches D:1109-1116 and the validator, which accepts a frozen section under a
  draft top level); (B) one freeze after all blockers close, so G3 waits on REL-03 and the
  reference work. **Recommend A.**
- **UC-I5 G3 date and sign-off wording.** Yours; the date is the `frozen_at` of every section and
  the sidecar.

## 8. Factor levels (prerequisite of the partitions and calibration freeze)

Not one of the six sections, but `partitions` needs three factors frozen and `calibration` needs
`viral_abundance` frozen, and `Y:1537-1538` says each factor freezes at G3 "after independent review
and user sign-off". The levels you have not individually signed are the agent's VAL-01 §4
proposals.

| Factor | Levels | Class | Source |
|---|---|---|---|
| viral_abundance | 1, 3, 10, 30, 100, 300, 1000 | U | D:1646 (D3) |
| sibling_virus_similarity | disjoint_cells, co_infected, one_member_only | U | D:1652 (D2) |
| host_virus_homology | none, repeat_homology | U (at your instruction) | D:1660 (section F) |
| infected_cell_fraction | 0.01, 0.04, 0.15 | A | `2026-10-05-def00-protocol-amendment.md:42` (B10, "D3 consequence") |
| low_complexity | 0, 0.01, 0.05 | A | `2026-10-05-def00-protocol-amendment.md:43` (B11) |
| pcr_duplication | 2, 5, 10 | A | `2026-10-05-def00-protocol-amendment.md:44` (B12) |
| cb_umi_collisions | 0, 0.005, 0.02 | A | `2026-10-05-def00-protocol-amendment.md:45` (B13) |
| ambient_index_hopping | 0, 0.01, 0.05 | A | `2026-10-05-def00-protocol-amendment.md:46` (B14) |
| chemistry | 10xv2, 10xv3, drop-seq | already frozen | `Y:1629-1634` |

- **UC-X1 sign off the five agent-proposed factor levels.** Options: (A) accept as written;
  (B) change levels (the panel size and the 42/252 stratum arithmetic move with them);
  (C) freeze only the three stratification factors at G3 and the rest at G4. **Recommend A.**

## 9. The G3 stanza, and how this was verified

```yaml g3-freeze
partitions:
  status: frozen
  frozen_at: "YYYY-MM-DD"
  contract_sha256: RECOMPUTE
calibration:
  status: frozen
  frozen_at: "YYYY-MM-DD"
  contract_sha256: RECOMPUTE
workflow_matrix:
  status: frozen
  frozen_at: "YYYY-MM-DD"
  contract_sha256: RECOMPUTE
failure_and_deviation_reporting:
  status: frozen
  frozen_at: "YYYY-MM-DD"
  contract_sha256: RECOMPUTE
defaults_selection:
  status: frozen
  frozen_at: "YYYY-MM-DD"
planned_freeze_sections:
  partitions: {status: frozen}
  calibration_metrics_lod: {status: frozen}
  workflow_matrix: {status: frozen}
  failure_and_deviation_reporting: {status: frozen}
  independent_review: {status: frozen}
freeze:
  freeze_manifest: analysis/v3_validation/freeze_sidecar.yaml
factors:   # merge by id
  - {id: viral_abundance, status: frozen, pending_reason: null}
  - {id: infected_cell_fraction, status: frozen, pending_reason: null}
  - {id: host_virus_homology, status: frozen, pending_reason: null}
  - {id: sibling_virus_similarity, status: frozen, pending_reason: null}
  - {id: low_complexity, status: frozen, pending_reason: null}
  - {id: pcr_duplication, status: frozen, pending_reason: null}
  - {id: cb_umi_collisions, status: frozen, pending_reason: null}
  - {id: ambient_index_hopping, status: frozen, pending_reason: null}
```

Verification (2026-10-13, interpreter `viralscan_test_full` Python 3.11, scratch copies only, nothing
committed). `apply.py` (appendix) extracts the `pending-edit`, `g3-freeze` and `sidecar` blocks
from this file, merges them into a scratch `protocol.yaml`, recomputes digests with the validator's
own functions, appends synthetic ledger records, then calls `validate_protocol_file`:

| Check | Result |
|---|---|
| `validate_v3_protocol.py` on `origin/main`, before any edit | valid |
| mode `pending`: all `pending-edit` blocks (defaults, workflow, failure, independent_review) | VALID, draft phase |
| mode `g3`: all blocks plus blocker removal, 4 section records and 1 `frozen_inputs` record | VALID, draft phase |
| negative controls on the g3 scratch: change `calibration.uncertainty.replicates`; unfreeze `viral_abundance` | both fail (digest mismatch; frozen partitions/calibration refuse an unfrozen factor), so the pass is not vacuous |
| sidecar skeleton | parses; required keys present; digest keys equal the frozen contract sections |

What G3 alone does not close (`--phase training` on the g3 scratch, 43 errors, all of these kinds):
unverified dataset and reference assets (the public FASTQ checksums, the synthetic and background
manifests, the curated index, the D-list and every native panel); 14 workflow rows without a
tool version and container digest, and `environment_pinning` not frozen (REL-03); the top-level
`status`, the freeze record and `training_allowed`; the remaining training blockers; and
`seeds.generation` (VAL-01). No error names a factor, a planned section, a ledger chain, a
partition or a calibration check, which is the point of the stanza.

The pending-phase pastes into `workflow_matrix` and `failure_and_deviation_reporting` each change
the digest of a section that already has a ledger chain. Each needs a re-baseline record at paste
time; otherwise `--phase training` reports pending-section drift (`pending_section_ledger_drift`,
`validate_v3_protocol.py:1221-1248`) until the G3 records supersede them. The draft phase does not
raise it, which is why the pending scratch run above is clean.

Recompute at apply time (never paste a digest by hand):

```bash
PYTHONPATH=src python - <<'PY'
import sys; sys.path.insert(0, "scripts")
import validate_v3_protocol as v
d = v.load_yaml(v.DEFAULT_PROTOCOL)
for n, f in (("partitions", v.partitions_sha256), ("calibration", v.calibration_sha256),
             ("workflow_matrix", v.workflow_matrix_sha256),
             ("failure_and_deviation_reporting", v.failure_reporting_sha256)):
    print(n, f(d[n]))
print("frozen_inputs", v.frozen_inputs_sha256(d))
PY
```

After pasting: append one ledger record per changed digest scope (the four sections and
`frozen_inputs`), declare their ids in `frozen_inputs.recorded_deviations`, and run
`scripts/validate_v3_protocol.py`, `scripts/check_ledger_append_only.py` and the protocol tests.
Number the records after the DEV-027 owner's records, not before.

## 10. Consolidated USER CHOICE list

Each line is the question and my recommendation. None has a recorded decision.

| ID | Question | Recommend |
|---|---|---|
| UC-P1 | Keep holdout fraction 0.3 and one holdout evaluation (agent-only, D:1078)? | keep |
| UC-P2 | Keep 2,000 cells and 1 technical replicate per sample (review round 8)? | keep |
| UC-P3 | Keep the three sibling pairs (review rounds 7-8)? | keep |
| UC-P4 | Is `pbmc_10x_healthy_v3` a planting background (grill 2, P:27) or external-only (protocol, D:1639)? | external-only |
| UC-C1 | The 45-point tier grid is unconsumed and no code emits probable/strong: freeze it as a contract on future tier code, drop it, or re-derive it? | freeze as contract |
| UC-C2 | Keep 2,000 bootstrap replicates, BCa, 95 % (agent-only)? | keep |
| UC-C3 | LOD grain: per virus x chemistry (draft) or per genus (grill 2)? | per virus, genus as descriptive roll-up |
| UC-C4 | Does grill 2's "no numeric specificity pass/fail" cover synthetic exact negatives (H2/E2 zero-call rule)? | no, real presumed negatives only |
| UC-C5 | Where do the T2/T3/T4 claim-gate minimums and scored T2 rows live: claim registry, new E11/D26, or dropped? | claim registry and Track B lock |
| UC-C6 | No eligible grid point: "no tier" (draft) or R3.2's "most conservative, screening only"? | "no tier" for tiers |
| UC-D1 | Defaults tuned jointly (6 candidates), coordinate-wise, or fully crossed? | joint, 6 candidates |
| UC-D2 | Accept the grid cost (about 19.9 core-h/sample, up to about 52; about 3.9k to 10.2k of 11k) against Q10d? | accept, measure multimap rerun cost first |
| UC-D3 | Enter the R2.5 within-sibling variant before G3? | defer |
| UC-D4 | Which defaults are chemistry-specific at G3? | none; confirm per chemistry |
| UC-W1 | Freeze the matrix before REL-03 pins, filling pins by one amendment? | yes |
| UC-W2 | Comparator version policy: latest tag at pin date, paper version, or per tool? | latest tag at pin date |
| UC-W3 | Extend W2 (ViralScan two-step) to every dataset per the 2026-10-06 same-arm-set rule? | yes, n/a reasons recorded |
| UC-F1 | Define "outcome" and hash the explore ledger into the sidecar (round-12 M6)? | yes, definition in section 6 |
| UC-F2 | Retries for scored rows: DSR-13 rule, none, or unlimited? | DSR-13 rule (OOM/TIMEOUT once at 2x) |
| UC-I1 | Pass bar: zero Blockers and Majors, plus Minors, or Majors accepted with reason? | zero Blockers and Majors |
| UC-I2 | Reviewer disagreement: either blocks, agent gates, or you adjudicate? | either blocks |
| UC-I3 | Sidecar location and hash convention? | `analysis/v3_validation/freeze_sidecar.yaml`, file bytes |
| UC-I4 | One freeze or two (top-level stays `draft` at G3)? | two |
| UC-I5 | G3 date and sign-off wording | yours |
| UC-X1 | Sign off the five agent-proposed factor level sets (B10-B14)? | accept |

## Appendix: apply.py (scratch verifier; run as `python apply.py REPO THIS_FILE pending|g3 WORKDIR`)

```python
"""Apply the fenced blocks of the pending-sections draft to a scratch copy and validate.

usage: apply.py REPO DOC MODE WORKDIR      (MODE = pending | g3)
Reads blocks tagged ```yaml pending-edit / ```yaml g3-freeze / ```yaml sidecar from DOC.
Never writes inside REPO.
"""
import copy
import pathlib
import re
import shutil
import sys

import yaml

repo, doc_md, mode, work = (pathlib.Path(a) for a in sys.argv[1:5])
mode = mode.name
sys.path.insert(0, str(repo / "scripts"))
import validate_v3_protocol as v  # noqa: E402

work.mkdir(parents=True, exist_ok=True)
text = doc_md.read_text()
blocks = {}
for tag, body in re.findall(r"^```yaml (pending-edit|g3-freeze|sidecar)\n(.*?)^```", text, re.S | re.M):
    blocks.setdefault(tag, []).append(yaml.safe_load(body))

protocol = v.load_yaml(repo / "analysis/v3_validation/protocol.yaml")
ledger = v.load_yaml(repo / "analysis/v3_validation/deviations.yaml")
schema = v.load_yaml(repo / "schemas/v3/validation_protocol.schema.json")
orig = copy.deepcopy(protocol)


def merge(dst, src):
    for key, val in src.items():
        if isinstance(dst.get(key), dict) and isinstance(val, dict):
            merge(dst[key], val)
        elif isinstance(dst.get(key), list) and isinstance(val, list) and val and all(
            isinstance(x, dict) and "id" in x for x in val
        ):
            by_id = {x["id"]: x for x in dst[key]}
            for item in val:
                merge(by_id[item["id"]], {k: x for k, x in item.items() if k != "id"})
        else:
            dst[key] = val


SECTION_FNS = {
    "partitions": v.partitions_sha256,
    "calibration": v.calibration_sha256,
    "workflow_matrix": v.workflow_matrix_sha256,
    "failure_and_deviation_reporting": v.failure_reporting_sha256,
}
counter = [0]


def add_record(scope, before, after, section):
    counter[0] += 1
    rid = f"DEV-SCRATCH-{counter[0]}"
    prev = ledger["deviations"][-1]["record_sha256"]
    rec = {
        "deviation_id": rid, "date": "2026-10-13", "affected_workflow_row_ids": [],
        "protocol_section": section, "description": "scratch", "reason": "scratch",
        "outcome_triggered": False, "outcome_triggered_justification": "scratch",
        "protocol_sha256_before": before, "protocol_sha256_after": after,
        "digest_scope": scope, "git_sha": "88f02a1", "review_reference": "scratch",
    }
    rec["record_sha256"] = v.ledger_record_sha256(rec, prev)
    ledger["deviations"].append(rec)
    return rid


for blk in blocks.get("pending-edit", []):
    merge(protocol, blk)

if mode == "g3":
    for blk in blocks["g3-freeze"]:
        merge(protocol, blk)
    # placeholder date -> scratch date
    for name in list(SECTION_FNS) + ["defaults_selection"]:
        protocol[name]["frozen_at"] = "2026-10-13"
    closed = {"partitions_metrics", "workflow_rows", "independent_review_round_2", "independent_review"}
    rd = protocol["execution_readiness"]
    rd["training_blockers"] = [b for b in rd["training_blockers"] if b["id"] not in closed]

new_ids = []
for name, fn in SECTION_FNS.items():
    section = protocol[name]
    old = section["contract_sha256"]
    section["contract_sha256"] = fn(section)
    if mode == "g3":
        chain = v._ledger_chain_for_section(ledger, name)
        before = chain[-1][1] if chain else section["contract_sha256"]
        new_ids.append(add_record(f"{name}.contract_sha256", before, section["contract_sha256"], name))

fi = protocol["frozen_inputs"]
if mode == "g3":
    fi["recorded_deviations"] += new_ids
changed_inputs = mode == "g3" or any(
    protocol[k] != orig[k] for k in ("defaults_selection", "execution_readiness", "factors")
)
if changed_inputs:
    before = fi["contract_sha256"]
    rid = f"DEV-SCRATCH-{counter[0] + 1}"
    fi["recorded_deviations"].append(rid)
    fi["contract_sha256"] = v.frozen_inputs_sha256(protocol)
    assert add_record("frozen_inputs.contract_sha256", before, fi["contract_sha256"], "frozen_inputs") == rid

(work / "protocol.yaml").write_text(yaml.safe_dump(protocol, sort_keys=False, width=100000, allow_unicode=True))
(work / "deviations.yaml").write_text(yaml.safe_dump(ledger, sort_keys=False, width=100000, allow_unicode=True))
shutil.copy(repo / "schemas/v3/validation_protocol.schema.json", work / "schema.json")

errs = v.validate_protocol_file(work / "protocol.yaml", work / "schema.json", ledger_path=work / "deviations.yaml")
print(f"[{mode}] draft phase:", "VALID" if not errs else errs)
if mode == "g3":
    errs = v.validate_protocol_file(
        work / "protocol.yaml", work / "schema.json", phase="training", ledger_path=work / "deviations.yaml"
    )
    print("[g3] training phase residual errors (expected: top-level status, blockers, assets, seeds.generation):")
    for e in errs:
        print("   -", e[:150])
    side = blocks["sidecar"][0]
    need = {"schema_version", "protocol_id", "gate", "frozen_at", "protocol_file_sha256", "ledger_sha256",
            "git_sha", "section_contract_sha256", "defaults_selection_sha256", "explore_ledger", "review",
            "user_signoff"}
    print("[g3] sidecar parses; missing keys:", sorted(need - set(side)) or "none",
          "| digest keys == frozen contract sections:",
          set(side["section_contract_sha256"]) == set(SECTION_FNS) | {"frozen_inputs"})
```
