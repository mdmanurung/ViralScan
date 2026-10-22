# DISC-12 — discovery benchmark protocol (DRAFT, not frozen)

**Status:** decisions recorded 2026-10-22 (grilling Q1–Q6, Q12); nothing has been run. The protocol is frozen when the
user approves this PR after the viral-only DB is downloaded and the spike-in truth check (§2) is done, so §8 pins every
hash at once. Until then no discovery performance number may be quoted (gate G13 of
`docs/plans/2026-10-21-kraken2-discovery-plan.md`).

**Parent plan:** `docs/plans/2026-10-21-kraken2-discovery-plan.md` §12.2 (scientific performance), §9.2 (states) and
§9.3 (evidence adapter). PLAN rows: `DISC-12` (this), `DISC-13` (cut-offs that need this), `HOST-05`, `SW-07`.

Items marked **DECIDED** record the 2026-10-22 grilling answers (`.living/decisions.md`). Rows still marked
**PROPOSED** were accepted as written; they can change only before the freeze.

---

## 1. Questions the benchmark answers

1. **Host-subtraction loss.** What fraction of truth viral fragments does the STAR host filter remove before discovery
   sees them? This is `HOST-05` (open), measured with `SW-07`'s `lost_truth_counts()`
   (`src/viralscan/scripts/host_filter.py:191`, D15/D16). Discovery reuses the number; it does not re-measure it.
2. **Database choice.** Viral-only RefSeq versus a broad DB with human (PlusPF-class) on the same inputs.
3. **Kraken2 thresholds.** Whether the shipped defaults (confidence 0.05, minimum hit groups 4; user decision R7b,
   documented as not calibrated) should change.
4. **Candidate-state cut-offs.** The `k` of the `prioritized` state (§9.2) and the `DISC-13` route constants
   (3 units / 80 % / 1 shared cell, currently uncalibrated).
5. **Resources.** Wall time, peak RSS and I/O at 1e6, 1e7 and 1e8 fragments.

Out of scope: any change to the primary ViralScan counts, `detected` or tiers; absolute abundance; batch effects.

## 2. Inputs (each pinned by sha256 in the freeze manifest)

| Role | Library | Truth | Status |
|---|---|---|---|
| Positive, panel virus | EBV SRR12682296 (used in DISC-R10) | EBV; cells with primary unique EBV molecules | PROPOSED |
| Positive, panel virus | HHV-6B: the Track A library with the most primary unique HHV-6B molecules (rule fixed before any discovery run) | HHV-6B; primary unique molecules | DECIDED (Q2) |
| Positive, panel virus | HSV-1 Drop-seq (no on-list) | HSV-1 | EXCLUDED from v1 (Q3): no Track A run has `called_cells.tsv`; reported `not_assessed` |
| Positive, panel virus | covid PBMC | SARS-CoV-2 | PROPOSED; primary calls SARS-CoV-2 = 0 here (SENS-03), so it tests capture, not recovery |
| Negative | hsv1 37 °C uninfected | none expected | PROPOSED |
| Negative | L-94 | none expected | PROPOSED |
| Spike-in | host background + known viral reads incl. EVE-like and host-homologous segments | per-read truth label | DECIDED (Q4): reuse the spike-in detection-limit grid (array 25864985) inputs if its truth has per-read labels; a separate set only if it lacks them or lacks host-homologous segments |

**Truth for panel positives** is the primary run's per-cell unique molecules, which is independent of Kraken2 but not
of the kallisto reference. Recovery is therefore "agreement with ViralScan", not accuracy against ground truth. Only the
spike-in set has read-level truth.

## 3. Fixed configuration

- **Code:** one ViralScan commit for all arms (SHA in §8).
- **Databases:** PlusPF-16 2026-06-26 (as in DISC-R10) and RefSeq viral `k2_viral_20260626` (same release; Q1);
  `hash.k2d` sha256 recorded for both.
- **Host subtraction:** computed once per library (STAR `pinned`); `star-default` as the second arm for question 1
  only.
- **Grid:** confidence {0.0, 0.05, 0.1, 0.2} × minimum hit groups {2, 3, 4} = 12 configurations per library per DB.
  Only Kraken2 is rerun. PolyA/T trim on (default), plus one `--no-polya-trim` arm at the shipped defaults.
- **Compute:** Slurm partition `all`, 8 CPUs (as DISC-R10); memory-mapped DB. Estimated 12 × ~1 h CPU per library per
  DB (UNRESOLVED; this run measures it).

## 4. Metrics (computed by a script committed before any run)

Per library × DB × configuration:

- **Taxon sensitivity:** panel virus reached as `observation` at species or genus; fraction of truth cells (primary
  unique molecules > 0) with ≥ 1 resolved discovery group; per-cell Spearman of discovery groups vs primary molecules.
- **False-positive taxa:** taxa at `observation` and at `prioritized` in negatives; non-truth taxa in positives
  outside the host lineage.
- **Rank correctness:** resolved rank vs truth lineage (exact / ancestor / wrong branch).
- **Host-lineage leakage:** fraction of resolved groups in the 9606 subtree (broad DB only).
- **Off-cell enrichment:** called-cell vs off-cell group fraction per taxon.
- **Host-subtraction loss (spike-in only):** D15/D16 from `lost_truth_counts()`.
- **Resources:** wall time, peak RSS, bytes read.

## 5. Decision rules (fixed before results)

1. **Kraken2 defaults (Q5).** Keep 0.05 / 4 unless a configuration has no more `prioritized` false-positive taxa than
   0.05 / 4 in **both** negatives **and** a higher median per-cell recovery of the panel virus in **every** positive. If
   several qualify: fewest false-positive taxa, then higher confidence.
2. **`k` for `prioritized` (Q6).** The smallest `k` at which neither negative has any taxon at `prioritized` while the
   panel virus still reaches `prioritized` in **every** benchmark positive. If no `k` meets both, report it and leave
   `prioritized` unbuilt.
3. **DISC-13 route constants.** Re-derived from the spike-in set only; the panel positives can confirm but not set them.
4. **DB.** Reported side by side; no automatic choice. The viral-only DB stays capped at `observation` (G15) whatever
   it scores.
5. A rule that cannot be evaluated (e.g. a library fails) is reported as `not_assessed`; it never defaults to "keep".

## 6. Outputs and receipts

- `results/discovery_benchmark/<date>/`: per-arm `qc.json`, `candidates.tsv`, `candidate_states.tsv`, metrics TSV,
  `freeze_manifest.json` (input, DB and code hashes), Slurm job IDs and `sacct` snapshot.
- A results note under `.living/findings/` citing the receipts; PLAN `DISC-12` flipped to `[x]` only then.

## 7. Known limits

- Panel-positive truth is ViralScan itself (§2), so agreement does not establish accuracy.
- Two negative libraries cannot estimate a false-positive rate; they only show whether one occurs.
- Kitome and reagent contamination have no built-in control (§9.1); the negatives are the only check.

## 8. Freeze record

| Field | Value |
|---|---|
| Approved by | — |
| Approval date | — |
| Protocol commit | — |
| Code commit for runs | — |
| Input/DB manifest | — |
