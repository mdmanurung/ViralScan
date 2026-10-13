# Real-data sensitivity and specificity: minimal design (PLAN AUD-12, finding M-04)

> **Superseded in part (DEV-029, 2026-10-13, Q120).** Real-data validation is descriptive only: real-data labels are
> presumed, never known, so no arm below licenses a sensitivity, specificity, PPV or accuracy claim, and the T2/T3/T4
> library minimums gate descriptive wording only (the T4 minimum of 10 libraries no longer licenses "specificity
> assessed"). Real-data results are reported under protocol endpoint `E11_real_data_descriptive` with the number of
> units; accuracy claims come only from synthetic planted truth. Read the arms and the bounds below with that rule.

Status: design only, written 2026-10-12. No code, no tests, nothing run: `src/` and `tests/` are
untouched, and no outcome run is authorised by this note (the protocol still has
`freeze.training_execution_allowed: false`, `analysis/v3_validation/protocol.yaml`). Claims about
the repository cite the file or PLAN row they come from, read not run. Numbers I could not verify
are marked "unverified". Literature is cited only as the repository already cites it ("as cited in
PLAN"); no external reference was checked for this note.

Source of the finding: Biomni review, `analysis/audit_2026-10/findings.tsv` row `M-04` ("audits
used synthetic edge cases; no real-data specificity/sensitivity shown"; proposed fix "negative
libraries, known positives, contamination challenges, orthogonal assay; report LOD separately for
presence vs programme"). Raised by Biomni, not independently audited. Companion note for the
contamination arm: `docs/plans/ambient-and-barcode-swapping-controls.md` (M-03).

## 1. Problem

What exists, and what it supports:

| Evidence | Where | Truth | What it can say | What it cannot |
|---|---|---|---|---|
| Synthetic truth panel (tiny) | `scripts/generate_truth_panel.py`, `scripts/v3_score.py`, `VAL-01/02/03/05/06/07`, `VAL-08` (partly) | Exact, per read or molecule | Mechanism correctness; count conservation; recovery of planted molecules; false calls in synthetic host-only and planted host-homology negatives | Behaviour on real soup, real repeats, real chemistry artefacts. Its background is synthetic (generator docstring), although `D1` chose real PBMC (PLAN `VAL-01`, 2026-10-05) |
| Analytic LOD and capture | `sensitivity.py`, `results/sensitivity.tsv`, `SENS-01..05` | Model | A Poisson depth floor; `informative_negative` | Capture is unmeasured without a control, so `informative_negative` is false almost always by construction (`SENS-03`) |
| Downsampling of one real EBV run | `SENS-05` | Design-known presence | The observed detection floor equals the Poisson floor on that run | Anything beyond one virus, one library, one chemistry |
| Real public positives | protocol datasets `hhv6b_srr20710641`, `ebv_srr12682296`, `hsv1_srr8315713`, `kshv_gse190558`, `kshv_ebv_gse154900` | Contextual (`context-only`) | That the tool runs and finds the expected virus (`RUN-01..05`) | A sensitivity or specificity estimate: the protocol says "published recovery is contextual" |
| Real presumed negatives | `pbmc_10x_healthy_v3` (`VAL-04`, blocked until G3/G4), tonsil atlas subset (`DSR-07`), `reagent_control` and `empty_droplet_control` (unselected) | Presumed negative | Specificity evidence, once run | Absolute truth: healthy PBMC can hold genuine EBV or anellovirus (protocol, `synthetic_factorial` planting rule) |
| Dossier ladder | `docs/plans/2026-10-07-dsr-dossier-spec.md` | n/a | Round 1 tops out at claim C1 and no negative is informative (N1) | Needs a positive control and planted recovery (G01, G02) that do not yet exist |

So the audits' synthetic edge cases are real tests of mechanisms and tell nothing about operating
characteristics on data with real background. The anellovirus claim carries an explicit rule:
"without this sample, publish screening support only and make no real-sensitivity claim"
(`anellovirus_orthogonal_positive`).

The minimal design below answers one question: what is the smallest set of real-data experiments
whose results would let the release state a defensible, bounded sensitivity and specificity, and
what each result may and may not be used to claim.

## 2. Truth definition (the central decision)

Real data have no cell-level truth, and ViralScan cannot supply its own labels without circularity.
Truth therefore comes in tiers; every result carries its tier, and tiers are never pooled.

| Tier | Source | Unit of truth | Registered examples (repo) | Supports |
|---|---|---|---|---|
| T1 exact | Planted molecules, recorded per read or molecule | Molecule, cell | Generator (synthetic background); `scripts/plant_anello_10x.py` plants held-out anellovirus reads into a real 10x library with `planted/`, `unplanted/` (same background, nothing added) and `negative/` arms; DEV-025 twin exclusion (`.living/decisions.md`, 2026-10-06) | Recall, precision, LOD by abundance and depth; cell-level metrics |
| T2 design-known | The experiment: infected vs mock or isogenic parental; a cell line that carries the virus | Library (sample presence); cell-level unknown | GSE189670 rafts (SRR19537341 vs SRR19537339, `EXPL-HPV16`); iSLK.219 latent vs induced (GSE190558); BC-1 PEL (GSE154900); the public EBV, HHV-6B and HSV-1 runs | Sample-level sensitivity at native depth; programme direction |
| T3 orthogonal assay | Clinical typing, qPCR, in situ, or an independent method on the same sample | Library | HNSCC GSE164690, HN18 vs HN01 (`DSR-06`; the run plan gives 6 HPV+ and 12 HPV- tumours, whether each has an assay in the deposit is unverified); `anellovirus_orthogonal_positive` has none | The only tier that supports a real-sensitivity claim for a virus without a cell-line or design truth |
| T4 presumed negative | A source that should not hold the virus | Library | `pbmc_10x_healthy_v3`; tonsil atlas negative subset (the run plan expects about 0 to 2 EBV+ cells per sample, unverified); the parental raft | Specificity bound; never "absolute ground truth" (`VAL-04`) |
| T5 contextual | Published results | n/a | `context-only` protocol rows | Context, not scoring |

A computational orthogonal method (a STARsolo viral-genome arm, `scripts/compare_starsolo_viralscan.py`;
read-level replay, `viralscan evidence`, `DSR-02`) is useful and independent in *method* but shares
the reads, so it can contradict a call and cannot create a truth label.

Unit of analysis. The protocol's resampling unit is the library, and it prohibits treating cells,
molecules, reads, lanes or twins of one donor library as independent replicates
(`calibration.uncertainty`). Every denominator in this note is a count of libraries. Depth
strata made by downsampling one library are n = 1 biological unit, however many strata there are.

**Hopped and ambient molecules in a negative: a scoring decision.** The isogenic parental raft
carried 9 clean HPV16 CB-UMI pairs, 4 of which also occur in the infected raft (PLAN
`EXPL-HPV16`, by hand, not re-verified here). At the default gate (one molecule in one called cell)
that library would be a sample-level false call, although the molecules are real HPV16 sequence
that reached the library by contamination, not misread host reads. Three readings are
possible and the choice changes the specificity number:

1. Count it as a false call: the tool reports what it sees; users face the same thing.
2. Exclude libraries with documented swapping from the specificity denominator and report them
   separately, with the swap evidence from the M-03 tests.
3. Count it, but also score the tool's contamination flag (M-03): the claim becomes "false call
   with a flag" versus "false call without".

Recommended: report all three numbers, with (1) the headline, since a user cannot know a library's
history. It must be declared before any run (decision 3).

## 3. Minimal design

Four arms, each with its own claim; no arm is a substitute for another.

### Arm P: planted into a real background (T1)

Take real, checksum-pinned background libraries, plant target molecules at known abundance into
their called cells, keep the unplanted twin, and score planted molecules only. This is the existing
pattern of `plant_anello_10x.py` (planted, unplanted, negative) combined with the protocol's
`synthetic_factorial` planting rule (at least ten donor libraries per chemistry, twin per planted
sample, DEV-025 per-cell twin exclusion for specificity and PPV). Differences from today:

- The generator's background is synthetic. A real-background mode with real called cells is the
  `D1` decision and is not implemented. The cell anchor must be outcome-independent (the protocol's
  `cell_universe`), taken after the reference freeze.
- Real genomes are planted via the AUD-18 hook (`--viral-fasta`, `--gene-spans`, `--onlist`); the
  panel still plants uniformly across gene models, which fits presence but not programme (section
  4).
- Planted molecules are idealised poly-A 3-prime fragments with known truth class (unique or
  ambiguous). Their realism limits any claim about reads that are fragmentary, internally primed
  or antisense (`PROG`/strand rows).

What it gives: recall and LOD95 per virus and chemistry at known abundance and depth in the
presence of real soup, real host homology and real artefact reads; cell-level precision and recall
against exact truth; the false-call rate *in the unplanted twin*, reported as background signal
(DEV-025), not subtracted.

### Arm D: design-known positives at native depth and downsampled strata (T2, T3)

Run the registered real positives on the frozen reference. Depth strata come from molecule-level
downsampling in the style of `SENS-05`, on the quantified-molecule axis (`sensitivity.py` defines
depth as the sum of the count matrix, not raw reads). The proposal is to stop at strata that bracket
the observed floor; the protocol's "at least seven levels bracketing the knee" binds the planted
arm (E8), not this one, and no number of downsampled strata of one library is more than one unit.
Output: sample-level presence sensitivity per (virus, chemistry, library), reported as a count of
libraries with the virus detected at each stratum. Where a library is also T3, the orthogonal
status is stated.

Real positives are high-load libraries; they say little about behaviour at the detection knee
(that is Arm P). Their value is realism of what the true virus looks like (gene breadth, splice
structure, programme state), and a check that Arm P's recall is not an artefact of the planter.

### Arm N: matched negatives (T4, plus T2 negatives)

Matched on the variables that move the false-call rate: chemistry and geometry; reference and
index; host cell-type composition; quantified molecules (downsample the negative to the positive's
depth stratum, so the exposure to a false call is the same); called-cell count; read-filter and
host-filter settings. Pairs available from the registered datasets: HN18 vs HN01 (same study, same
chemistry), HPV16 raft vs parental raft, the tonsil atlas subset and `pbmc_10x_healthy_v3` against
10x v3 positives. Output: sample-level false-call rate (M2) per tier of evidence, plus the number
and nature of called cells (descriptive, since presumed negatives have no cell truth).

Specificity is defined over the whole panel: a sample-level false call is any probable-or-stronger
call for any virus outside the truth set (M2 reads "at least one probable-or-strong call").

### Arm X: contamination challenge (links to M-03)

Take real positive molecules and add them to negatives at known fractions as ambient-origin or
as swapped (CB-preserving) molecules, or use a same-pool negative library if one exists. This is
the arm that measures the contamination false-call rate and tests whether the M-03 flags separate
it from infection. It needs the generator's CB-preserving hop mode (M-03 section 7), which does not
exist, and it is the only arm with *exact* truth for a contamination event.

### Orthogonal assay

Any wet assay (qPCR on an aliquot, in situ) is a user-owned experiment. The design does not depend
on it, but it is the only route to a T3 truth for anellovirus, HHV-6A, HSV-2 and any other virus
without a design-known library, and the protocol records its absence as a blocker for a
real-sensitivity claim.

## 4. Presence and programme are separate questions

The finding asks for separate LODs. They have different truths and different instruments:

| | Presence | Programme |
|---|---|---|
| Output | `detected`, evidence tier | `gene_programs` states: latent, productive, `reactivating`, `indeterminate` (AUD-19: only EBV can reach `reactivating`; unresolved productive markers degrade absence-based states) |
| Instrument | Probit of detection on planted abundance (E8) | None yet. A programme call needs marker molecules in cells, and the markers sit at a different (higher) molecule count than a single detection |
| Real truth | T2/T3 library presence | T2 direction only: induced vs uninduced (GSE190558), latent in a latency-only line. No cell-level programme truth |
| Planted truth | Uniform over gene models | Needs expression-weighted planting (latent vs lytic gene sets); the VAL-01 design note defers this as "add only if review asks" |

Consequence: a programme LOD cannot be estimated from the current generator, and real data support
only a **direction** (the induced library scores more lytic than the uninduced one), not a limit.
The honest deliverable for the first release is a presence LOD with a stated capture, and for
programmes a per-library directional check and a statement that no programme LOD was estimated.

## 5. Metrics and acceptance

All metrics reuse protocol and scorer definitions; none is new. Acceptance thresholds are **not**
proposed here, because the protocol forbids post-outcome targets (`E3`, `E4`: "no post-outcome pass
target") and picking one now would be inventing it. They must be declared by the user before any
run (decision 5).

| Question | Metric | Unit and denominator | Existing definition |
|---|---|---|---|
| Sample-level specificity | False-call rate | Complete, comparable negative libraries; failed rows are non-estimable, not correct negatives | M2 / E2, `v3_score.negative_false_call_rate` |
| Tier behaviour | Specificity and PPV per tier, monotone candidate to probable to strong | Calls in that tier, success-conditioned and labelled so | M8 / E9, `v3_score.tier_specificity_ppv` |
| Molecule and cell recovery | Precision, recall, F1, AUPRC | Planted molecules or truth-infected cells (Arm P only) | M3, M4, `v3_score.precision_recall_f1`, `auprc` |
| Presence LOD | LOD95 per virus and chemistry, bootstrap interval over libraries, no extrapolation: if the curve does not reach 0.95 it is reported as not reached at the highest level | Planted abundance levels (at least seven, bracketing the knee) | E8, `limit_of_detection` |
| Sample-level presence sensitivity | Libraries detected / libraries truly positive, per stratum | Libraries of that tier and stratum | New use of the existing gate output |
| Contamination | False-call rate with and without an M-03 flag | Arm X libraries | M-03 flags (proposed) |
| Cross-method | Paired concordance | Matched samples | M9 / E10, descriptive |

Uncertainty follows the protocol: library as the unit, at least three units for an interval, below
that a point estimate and "not estimable". With few real libraries the interval is the message.
Rough bounds for the all-correct case (one-sided exact 95 %, `1 - 0.05^(1/n)` for 0 false calls in
`n` negative libraries and `0.05^(1/n)` for `n` of `n` detected; computed here, a one-sided bound
is not the two-sided 95 % interval, which would be wider):

| n libraries | Specificity: upper bound on false-call rate after 0 of n | Sensitivity: lower bound after n of n detected |
|---|---|---|
| 3 | 0.632 | 0.368 |
| 5 | 0.451 | 0.549 |
| 10 | 0.259 | 0.741 |
| 20 | 0.139 | 0.861 |
| 30 | 0.095 | 0.905 |
| 60 | 0.049 | 0.951 |
| 100 | 0.030 | 0.970 |

Reading: a real-data claim of the form "specificity of at least 95 %" needs about 60 clean negative
libraries of one matched kind; with the registered negatives (a handful) the supportable statement
is a wide bound, which is what should be reported. These are arithmetic consequences of the
sample size, not targets.

## 6. What each result may and may not support

| Result | May support | May not support |
|---|---|---|
| Arm P recall and LOD95 | Detection limit of the pipeline for planted molecules of that virus, chemistry, depth, in that real background | Limits for divergent natural strains (the planted sequences are the reference or held-out genomes of known identity; `SENS-CAL-01` is the capture calibration); a real infection's cell-level behaviour |
| Arm P twin false calls | Background signal that exists in the real library | A false-positive rate: the twin may hold a real virus (it is background signal, DEV-025) |
| Arm D, all detected | That the pipeline finds the virus in these libraries at native depth | Sensitivity at low load; sensitivity for other viruses, chemistries or tissues; any per-cell statement |
| Arm D, downsampled strata | How detection falls with depth in this one library | Replication: n = 1 per stratum |
| Arm N, 0 false calls in n libraries | Specificity at least the bound in section 5, for libraries like those, on this reference and settings | Absolute specificity; specificity in a sample type not represented; "no virus" (N2 is never claimable: the panel is the limit, dossier spec) |
| Arm N, a call | A candidate for read-level validation (`DSR-02`); not yet a false positive, since the library may hold the virus | A confirmed false positive without evidence of absence |
| Arm X | Whether contamination alone makes a call, and whether the flags see it | Rates in a real run, whose contamination level is unknown |
| T3 positive (orthogonal) | A real-sensitivity claim for that virus | Any other virus |
| No T3 anellovirus | Screening support only (protocol rule) | Any real-sensitivity claim for anellovirus |

This maps onto the dossier ladder (`2026-10-07-dsr-dossier-spec.md`): Arm P plus a measured capture
(`--positive-control-*`, `SENS-04`, `SENS-07`) is what makes N1 ("no panel virus above the LOD")
claimable for a library; without it every negative stays uninformative. Arm D and N feed C1/C2.
C3 (cell-type infection) needs labels these arms do not provide.

## 7. How it connects to the plan

- Truth panel generator (`scripts/generate_truth_panel.py`): rows `VAL-01` (generator), `VAL-02`
  (truth files), `VAL-03` (host-only and host-homology; reagent and empty-droplet controls missing),
  `VAL-05` (positives; real genomes via the AUD-18 hook), `VAL-06/07` (scorer and fixtures),
  `VAL-08` (golden run; the `kb count` arm is open), `VAL-09/10` (freeze and full panels, behind
  G3). Arm P is `VAL-01` with a real background; Arm X needs the M-03 generator change.
- `VAL-04` (`pbmc_10x_healthy_v3`): Arm N's first dataset; the run plan holds it until G3 and G4.
- `SENS-04/06/07`, `SENS-CAL-01`: the positive-control recipe and per-dataset planted-read
  sensitivity that give `capture_measured`; Arm P generalises `SENS-07`.
- `DSR-06` (HNSCC pilot), `DSR-07` (tonsil negatives), `DSR-10` (standard arm set on every
  dataset): the T3 pair and the T4 negatives. `DSR-16`: one reference cell set per sample for every
  arm, which the scoring needs.
- `DEF-00`: the protocol amendment (section 8).
- `AUD-21`, `AUD-07`: `detected` is the called-cell gate, so sensitivity and specificity as defined
  here are for that gate; a change to the gate changes the numbers.

## 8. Gating and order

1. Nothing in this note may run before the protocol freeze, except runs labelled exploratory as
   the `DSR` rows already are (they tune nothing and score nothing).
2. Upgrading any `context-only` positive to scoring (`external-evaluation` with a `design-known`
   truth status) is a protocol amendment with a deviation record, disclosed as outcome-triggered if
   made after any result is seen. The first useful, outcome-independent step is declaring the truth
   tier of each registered dataset *now*, from metadata only.
3. Order by information per cost: (a) declare tiers and the scoring rule for hopped or ambient
   molecules in negatives; (b) Arm P on the generator with the AUD-18 real-genome hook, then on
   real backgrounds once pinned; (c) Arm N on `pbmc_10x_healthy_v3` and the tonsil subset after G3
   and G4; (d) Arm D on the registered positives; (e) Arm X once the CB-preserving generator mode
   exists.

## 9. Failure modes

1. **Circular truth.** Labels derived from ViralScan output (or from the same reads) score the tool
   against itself. Only T1 to T3 are labels.
2. **Pseudo-replication.** Cells, lanes, downsampling strata or twins counted as replicates. The
   protocol forbids it; the bounds in section 5 count libraries.
3. **Presumed negatives that are positive.** A call in a healthy library may be real. The design
   reports "call needing read-level validation", not "false positive", until evidence says so.
4. **Contamination in the negatives** (the parental raft). Addressed by the scoring decision in
   section 2, with M-03 flags.
5. **Planter realism.** Idealised reads overstate recall for fragmentary or antisense signal.
6. **Matching failure.** A negative with a different depth, composition or chemistry is not a
   control for the positive. Matching is declared before running.
7. **Whole-panel specificity.** A sample-level false call over a panel of thousands of genomes
   (the shipped panel is 2,343 genomes per `EXPL-HPV16`, unverified here) is much easier to
   trigger than a per-virus one; report both the any-virus and the per-virus rates.
8. **Gate dependence.** Numbers belong to `--detection-threshold`, `--detection-min-cells`,
   `--infected-min-molecules`, the cell caller and the host-filter strategy; they are not
   comparable across settings (the `pct_infected_comparable` precedent).
9. **Reference dependence.** The frozen reference (`G4`) is part of the measurement; a
   reference change re-opens all arms (the cat42b to cat42d swap moved emptyDrops calls, DEF-00 C1).
10. **Single-library positives.** A real positive says the virus is found there, not how often.

## Documentation rows implied

None edited here. `docs/output_reference.md` (pinned) would need no new rows for Arm P, D or N,
which score existing outputs; any M-03 flag adds rows as in that note. A future `docs/` page or
vignette for the spike-in recipe closes `SENS-06`.

## Decisions needed from the user

1. Approve the truth-tier vocabulary (T1 to T5) and the rule that tiers are reported separately and
   never pooled, and which registered datasets are assigned to which tier.
2. Is a protocol amendment to promote selected `context-only` positives to scored `design-known`
   rows acceptable, and which ones (GSE189670 pair, GSE190558, GSE154900, the EBV, HHV-6B, HSV-1
   runs)?
3. Scoring rule for hopped or ambient molecules in a negative library: count as false calls,
   exclude with evidence, or report all three numbers (recommended: three numbers, with the
   unconditional one as the headline).
4. Real-background generator: implement the `D1` real-PBMC mode (planted plus unplanted twin) before
   any real-negative run, and the CB-preserving hop mode for Arm X.
5. Acceptance thresholds: you fix them before any run (sensitivity floor, specificity bound, minimum
   libraries per claim). This note proposes none.
6. Programme LOD: accept "direction only, no programme LOD" for the first release, or fund
   expression-weighted planting.
7. Wet inputs you own: an orthogonally confirmed positive (anellovirus, HHV-6A, HSV-2), a same-pool
   negative library, a reagent control. Without them those claims stay "screening support".
8. Specificity scope: any-virus sample-level false call as the headline, with per-virus rates
   alongside.

## Scorer implementation (PLAN M04-SCORE, 2026-10-12)

Decisions 1, 3 and 5 of the user's 2026-10-12 round, and the tier vocabulary of decision 1, are in
`scripts/v3_score.py`; tests in `tests/test_v3_scorer_m04.py`. No acceptance threshold is embedded: the
report carries counts and ratios and no verdict field.

- `TIERS` (T1 to T5) and `check_tier`. The scorer takes the tier as input and only validates it; no
  registered dataset is assigned to a tier (that waits for the G3 round). T5 is contextual and refused
  for scoring. `score_negative_libraries` refuses a list that mixes tiers; `score_negative_libraries_by_tier`
  returns one report per tier and has no pooled entry.
- `score_negative_libraries` reports a negative's hopped or ambient molecules three ways, always
  together: (i) the unconditional any-virus sample-level false call (the HEADLINE, per-virus rates
  alongside), (ii) the rate after leaving out libraries with documented contamination, listing their
  provenance (`hop_donor_sample_id`, `hop_cb_preserved`, `hop_landing`, or `origin` alone for a pre-M04
  manifest), (iii) the contamination-attributable false calls, over the same denominator as (i), beside
  the residual (`unconditional = attributable + residual`). `format_negative_report` prints all three.
- Conservative readings. A library leaves the denominator in (ii) only if its truth manifest holds a
  hopped or ambient viral molecule and its call, if any, is fully explained: every called virus has
  supporting `(cb, umi)` molecules and all are hopped or ambient. A call with one other supporting
  molecule, or none, stays a false call. The scorer does not re-run the detection gate, so the caller
  passes the supporting molecules of each call. A negative library holding cell-origin planted molecules
  is refused.
- `origin = unknown` (D1 background) is accepted, and `read_truth_manifest` refuses it if labelled
  planted; `truth_positive_molecules` never returns it. `hop_*` columns are validated when present
  (only on `hopped` rows; `hop_landing` is `cell`, `empty` or `other`).
- `realised_shortfalls(detail)` compares requested and realised constructions per sample, and a D1
  sample (`background_library_id` present) is exempt: it has `realised.hopped == 0` beside
  `requested_hopped > 0`.
