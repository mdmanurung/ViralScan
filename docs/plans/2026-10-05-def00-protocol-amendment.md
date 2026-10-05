# DEF-00 — protocol amendment list (draft)

**Status:** written into `analysis/v3_validation/protocol.yaml` (v0.3.0) on
2026-10-05 and recorded as ledger record `DEV-020`. Every factor and section
stays `pending`. The route to G3 is now: independent review → user sign-off →
freeze.
**Target:** `analysis/v3_validation/protocol.yaml`
**Sources:** WP1E decisions (`PLAN.md:1682`), the VAL-01 design decisions D1–D7
(`docs/plans/2026-10-05-val01-generator-design.md`), and the REF-06 measurement
of 2026-10-05.

Each row names the section it edits, what changes, and why. A row marked
**blocked** cannot be written yet and says what it waits on.

---

## A. Pre-existing DEF-00 rows (PLAN.md:1725)

| # | Section | Amendment | Source |
|---|---|---|---|
| A1 | new `defaults_selection` | grid, per-default objective, hard constraint, 1-SE conservative tie-break, sample bootstrap | R2.1, R3.1, R3.2 |
| A2 | `supported_scope.chemistry` | add 10x 5′ as `10xv2`/`10xv3` geometry plus strand handling | Q5 |
| A3 | `counting_contract` | within-sibling-group allocation variant in D3; the conservation claim is amended **before** the variant ships | R2.5 |
| A4 | `defaults_selection` | per-chemistry defaults are allowed and are preregistered per chemistry | R2.7 |
| A5 | `datasets` | generator backgrounds: real PBMC primary plus a synthetic GRCh38 control | R2.9, superseded in scope by D1 below |
| A6 | `partitions.outcome_blinding` | the holdout blinding layout | R3.4, extended by D6 below |
| A7 | `sensitivity_only_callers` | the `knee` estimator is defective (SW-23): fix or replace it | SW-23 |

## B. VAL-01 decisions (user, 2026-10-05)

| # | Section | Amendment | Decision |
|---|---|---|---|
| B1 | `datasets[synthetic_factorial]`, `[synthetic_mixed_host_virus]` | background is a real healthy-donor PBMC library per chemistry; each planted sample has an **unplanted twin**; planting goes into that library's outcome-independent called cells; `pbmc_10x_healthy_v3` is excluded as a background because it is the VAL-04 evaluation negative | D1 |
| B2 | `partitions.unit`, `uncertainty.resampling_unit` | both read **background donor library**; ≥ 10 donor libraries per chemistry; generated samples nest within donors; the bootstrap is a donor-level cluster bootstrap | D1.1 |
| B3 | `factors[sibling_virus_similarity]` | levels `disjoint_cells`, `co_infected`, `one_member_only`; `one_member_only` falls on one sample per stratum and the dropped member is a seeded draw | D2 |
| B4 | `endpoints`, new denominator | **new metric** for absent-sibling false calls: calls for the unplanted pair member in a `one_member_only` sample, denominated on that sample's anchor cells. Without it `one_member_only` has no consumer, which is the frozen-value-with-no-call-site defect the 2026-07-27 review named | D2 |
| B5 | `factors[viral_abundance]` | redefine as **total planted molecules per virus per sample**; levels 1, 3, 10, 30, 100, 300, 1000. The current text says "per infected cell", which smears the E8 probit because MECH-B thresholds the summed count | D3 |
| B6 | `partitions.leakage_prohibitions` | "template" and "locus" mean **host and challenge** sources. Viral genomes are shared across partitions because they are the estimand; molecules, UMIs and barcodes stay disjoint | D4 |
| B7 | `generated_sample_structure` | `reads_per_cell` stays 25,000 (no amendment). `audit_artifacts` truth manifests change granularity: one row per read for every non-host read, one row per molecule for host molecules, with the host molecule id in the read name for training samples only | D5 |
| B8 | `partitions.outcome_blinding`, asset locators | holdout truth goes to an owner-only path; the protocol records only its hash. Holdout FASTQs carry opaque read names and a seeded shuffled read order, and their run manifests are owner-only | D6 |
| B9 | `seeds.generation` | holdout generation takes an owner-held secret salt; the protocol publishes `sha256(salt)` only, and the salt is revealed at unblinding. Training stays reproducible from the public seed alone | D6 |
| B10 | `factors[infected_cell_fraction]` | per virus: 0.01, 0.04, 0.15. Six disjoint sets at 0.15 use 1,800 of 2,000 cells; `co_infected` pairs share theirs. `is_infected_cell` means "received ≥ 1 planted molecule" | D3 consequence |
| B11 | `factors[low_complexity]` | 0, 1 %, 5 % of library reads, from the four measured artefact classes with tract ≥ 40 nt. Viral 3′ poly-A tails are always on and are not this factor | VAL-01 §4 |
| B12 | `factors[pcr_duplication]` | mean reads per molecule 2, 5, 10 (1 + geometric) | VAL-01 §4 |
| B13 | `factors[cb_umi_collisions]` | 0, 0.5 %, 2 % of viral molecules, half *same* key and half 1-Hamming *disjoint* | VAL-01 §4 |
| B14 | `factors[ambient_index_hopping]` | 0, 1 %, 5 %; hopped reads come only from a sample in the **same partition and chemistry** | VAL-01 §4, D4 |
| B15 | `cell_universe` | empty-droplet count and depth are generator parameters and must be frozen with the levels (the protocol requires empty droplets for D11 but declares no count) | VAL-01 §4 |

## C. From the REF-06 measurement (2026-10-05)

| # | Section | Amendment | Why |
|---|---|---|---|
| C1 | `cell_universe.source_precedence` | a kb-derived cell anchor must be frozen **after** the reference freeze (G4), or come from a reference-independent source | Swapping cat42b → cat42d moved emptyDrops calls by −22 % on the HPV16 rafts at the same seed, because intronic and intergenic host reads are now D-listed out. An anchor frozen before G4 would be invalidated by the reference freeze itself |
| C2 | `references[curated_human_virus]` | the candidate production reference is cat42d: cat42b plus a checksum-pinned GRCh38 genome D-list, 3,094,720 D-list k-mers | REF-06 gates passed |

## D. Blocked rows

| # | Section | Waits on |
|---|---|---|
| ~~D1~~ | `factors[host_virus_homology]` | **No longer blocked; the measurement changes the row.** REF-07 ran 2026-10-05: 98.8 % of the viral bases that align to GRCh38 are low-complexity against a 4.5 % panel background, 1 of 89 clusters carries a 31-mer outside low complexity, and 0 of those are exonic (finding F-023). A synthetic control recovers non-repeat decoys 100 % at 90-95 % identity and 80 % at 85 %, so this is a measurement rather than a sensitivity failure. The identity/aligned-length bands the factor describes have no material. See **F** below |
| D2 | `calibration.*` thresholds | REF-08, which calibrates on the training panel's own controls and therefore follows VAL-01, not precedes it |
| D3 | `tool_environments` | REL-03 |

## E. Consistency checks the amendment must pass

1. Every factor marked `frozen` has non-empty `levels`; the validator refuses a
   freeze otherwise (`scripts/validate_v3_protocol.py:194`, the frozen-factor check).
2. `viral_abundance` has ≥ 7 levels (`MINIMUM_ABUNDANCE_LEVELS`,
   `validate_v3_protocol.py:26`). B5 gives exactly 7.
3. The sample floor depends on the allocation unit (B2). For a
   synthetic-background dataset it is 4 × the stratum count: 7 × 2 × 3 = 42
   strata, 168 samples. For a real-background dataset it is ≥ 10 donor
   libraries per chemistry, each carrying the full 7 × 2 grid: 140 planted
   samples and 140 unplanted twins per chemistry, with 3 libraries in holdout,
   so every stratum has 3 holdout samples by construction.
4. Every new level has a consumer: B4 exists because B3 otherwise would not.
5. `partitions.unit` and `uncertainty.resampling_unit` agree (B2). They are
   separate keys today and could drift.
6. No section flips to `frozen` in this amendment. The freeze is G3, after
   review and sign-off.

## F. `host_virus_homology` after the REF-07 measurement (new, 2026-10-05)

F-023 found that viral/host nucleotide homology in this panel is 99.7 %
low-complexity. `host_virus_homology` and `low_complexity` are therefore not
separable axes, and the factor cannot carry identity bands. The protocol's
three-factor stratification (abundance x homology x chemistry) rests on this
factor, so the choice changes the panel size. It needs the user's decision.

**Decided: (a).** Decided by Claude at the user's instruction (2026-10-05).
(c) makes the holdout homology claim vacuous by the protocol's own
unsplittable-stratum rule, and (b) drops a stratification axis and forces an
endpoint renumbering that nothing requires. (a) keeps the protocol's structure
and records the confound per read.

- **(a) Chosen: two levels, `none` and `repeat_homology`.** Challenge reads
  are planted from the 89 measured clusters (70 training, 19 holdout, assigned
  whole), and each read carries its own `low_complexity_fraction` and
  `homology_cluster_id` in `host_homology_manifest.tsv`, so the confound is
  *recorded per read* rather than designed away. Strata become
  7 x 2 x 3 = 42. The floor is 168 for a synthetic-background dataset; a
  real-background dataset uses the donor-library floor in E3.
  `low_complexity` stays a separate covariate for the *artefact* classes
  (poly-G, poly-A, CAG, TSO joins), which are a different population from
  host-derived repeat reads.
- **(b) Merge the two factors into one `host_derived_challenge` factor** with
  levels none / low-complexity-artefact / repeat-homology. Honest, but it drops
  a stratification axis and forces a wider renumbering of the endpoints.
- **(c) Keep the factor as written and report it as unpopulated.** The strata
  would exist with one realizable level, which is the unsplittable-stratum case
  the protocol already routes to training, so the holdout homology claim would
  be vacuous.

Either way `REF-09` ("planted human-homology reads cannot reach probable or
strong") now tests mostly low-complexity reads, which `DEF-01`'s read-artefact
filter already targets. Say so in the protocol rather than letting the row read
as independent evidence.

## G. Where each row landed (2026-10-05, DEV-020)

| row | protocol location |
|---|---|
| A1, A3, A4 | new top-level `defaults_selection`. Every `grid` is empty and must be filled before G3. The section is covered by the `frozen_inputs` digest. |
| A2 | `supported_scope.five_prime` (geometry 10xv2/10xv3 and strand handling). The frozen `chemistry` factor is unchanged. |
| A5, B1 | `datasets[synthetic_factorial].planting_rule` and a new `pbmc_background_manifest` asset; `synthetic_host_only.source` |
| A6, B8 | `partitions.outcome_blinding`; the truth-manifest asset `pending_reason` |
| A7 | `cell_universe.sensitivity_only_callers: [none]`. The knee is **removed**, not fixed. A corrected knee re-enters only by amendment. |
| B2 | `partitions.unit`, `algorithm`, `stratification_rationale`, `minimum_samples_rationale`; `calibration.uncertainty.resampling_unit`, `prohibited_units`, `below_minimum_policy`; `tie_breaker_standard_error`; E9 reporting |
| B3, B4 | factor levels; new `D25_absent_sibling_false_calls`, added to E5; E5 acceptance; `sibling_pair_rationale` |
| B5, B10–B14, F | factor levels and descriptions, all `pending` |
| B6 | `leakage_prohibitions[1]` |
| B7 | `generated_sample_structure.file_layout`; new truth-manifest columns `background_library_id`, `twin_sample_id`, `record_granularity`, `origin` |
| B9 | `seeds.generation.pending_reason` |
| B15 | `cell_universe.synthetic_rule` |
| C1 | `cell_universe.public_data_rule` |
| C2 | `references[curated_human_virus].purpose`; the D-list asset is pinned to `e368c0e5…` as `recorded-reverify` |
| F caveat | `factors[host_virus_homology].description` and `principal_risks` (REF-09 tests mostly low-complexity reads) |

Calls made while writing, for review to check:
- **A7: remove the knee, don't fix it.** Fixing it is a code change with its own test, outside DEF-00.
- **B2: the allocation unit is per dataset.** Real-background datasets split by donor library and stratify by chemistry only, with the full abundance × homology grid nested in every library. Synthetic-background datasets keep the generated sample as the unit.
- **B15: a real library's empty droplets are its uncalled barcodes** with a host-only UMI total in [10, `lower`).
- **D25 excludes cells the unplanted twin also calls.** Those are background signal.

Open for review:
- D18 still says "exact truth-labelled" anchors. On a PBMC background, specificity is the planted-minus-twin difference, and D18's population text does not say so yet.
