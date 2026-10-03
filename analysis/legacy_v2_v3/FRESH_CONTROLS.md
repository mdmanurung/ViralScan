# Fresh matched-FASTQ controls: ViralScan 2.2.0 versus frozen v3 (LVC-11 to LVC-13)

Scope: diagnostic only, outcome-ineligible. Historical outcomes were already
viewed (`confirmatory_holdout_eligible: false`). Nothing here closes a truth-panel,
comparator, calibration, release, or publication gate, and nothing here claims
viral absence, BCC association, calibrated performance, or package superiority.

## Arms and disclosures

- v2 arm: ViralScan 2.2.0, original index and t2g, its packaged 195 GTFs, no HIV GTF.
- v3 arm: frozen v3 wheel (3.0.0.dev0), same index, t2g, whitelist and FASTQs, with
  explicit `-gtf` = the v2-arm GTFs plus one HIV GTF (protocol 1.2.0, outcome-triggered).
  **The arms therefore differ on HIV**: the v2 install cannot call HIV, so its HIV
  result is structurally negative and is not a like-for-like sensitivity comparison.
- The residual GTF/t2g parity gap (89 GTF-only genes, 190 non-Ensembl t2g-only genes)
  was accepted fail-closed (digest `cad05a71...`).
- The frozen v3 wheel predates SW-13 (no barcode correction without `-w`; `-w` was
  supplied) and carries SW-24 (absolute tolerance in the multimap conservation check).
  Neither was patched. No row failed because of either.
- Tier note: v2 `SRR6825024` ran on the highmem tier after the attempt-2 out-of-memory
  failure (retained, protocol 1.1.0).

## Run status (all ten fresh rows)

All ten rows have status `success`, exit 0 and no artifact validation errors
(v3 rows also `validate-run` exit 0). No row is excluded. Earlier failed attempts
(attempt 1 spool-root, attempt 2 validator defect, missing cache, v2 OOM) are retained
in the tracker. Evidence attempt 1 (five records) failed because `kallisto` and
`bustools` were not on `PATH` in the job; those records are retained as
`*.evidence.attempt1_failed.json` and the retry is the table below.

## 1. Did fresh v2 reproduce its archive?

Yes. All five v2 rows match the archived legacy tables on every gene, virus and
total-viral-load record (tolerance 1e-6): 99/99, 97/97 and 95/95 records for the
three EBV controls, and 12/12 and 9/9 for the two HIV controls
(`fresh_vs_archive.json`). Archived v2 total viral load equals fresh v2 for every
control (1,792,862; 262,394; 31,016.3; 38; 12).

## 2. Fresh v3 versus archived v3 (identical-BUS arm)

The fresh v3 host-conservative per-virus totals match the archived v3 host-conservative
values for 104-105 of 106-116 virus records per control. All 26 mismatching records
(tolerance 0, so any difference counts), classified from the report without further
metrics:

| Sample | Record | Fresh v3 | Archived v3 | Reading |
|---|---|---:|---:|---|
| SRR12682296, 7 | Epstein-Barr virus | 874,039 / 121,749 | 874,039.0000000 (-6e-9) / 121,749 (-7e-11) | floating-point noise |
| SRR12682298 | Epstein-Barr virus | 18,343 | 18,345 | -2 molecules (0.011 percent); unresolved |
| SRR12682296, 6825024, 6825025 | Cercopithecine herpesvirus | 0 | 1 | one molecule not recovered; unresolved |
| SRR6825024, SRR6825025 | IMMUNO_HIV1gp1-10 (20 records) | 744,569 / 90,464 total | 0 | reference difference: archived v3 lacked an HIV feature, fresh v3 carries the HIV GTF |

No mismatch was attributed to SW-13 or SW-24; the 1-2 molecule differences are
left `unresolved` rather than guessed.

## 3. Fresh v2 versus fresh v3 on identical FASTQs (full stacks)

Descriptive paired counts (`fresh_control_stack_comparison.tsv`). v2 is the archived
estimator; v3 is host-conservative molecules. They are different quantities
(see `REPORT.md`) and are not rescaled.

| Control | Target | v2 | v3 host-conservative | v3 / v2 |
|---|---|---:|---:|---:|
| SRR12682296 | EBV | 1,792,860 | 874,039 | 0.49 |
| SRR12682297 | EBV | 262,394 | 121,749 | 0.46 |
| SRR12682298 | EBV | 31,011 | 18,343 | 0.59 |
| SRR6825024 | HIV (IMMUNO_HIV1gp1-10) | not callable | 744,569 summed over ten genes (gp1 338,344) | not defined |
| SRR6825025 | HIV (IMMUNO_HIV1gp1-10) | not callable | 90,464 summed over ten genes (gp1 28,225) | not defined |

Non-target calls in the two HIV controls (v2 versus v3): HHV-6b 6 vs 2 and 6 vs 3;
vesicular stomatitis 20 vs 6 (SRR6825024); monkeypox 4 vs 1; molluscum
contagiosum 2 vs 1 (SRR6825025); HHV-1 (2 and 2), HHV-8 (6, SRR6825024), and human
parainfluenza (2, SRR6825025) are present in v2 and absent (0) in v3. Cercopithecine
herpesvirus is 5 in both stacks in SRR12682298 only.

## 4. Expected-target recovery (qualitative)

- EBV (three controls): recovered by both stacks; abundant. This is qualitative
  positive-control evidence only.
- HIV (two controls): recovered by v3 (exact-read evidence below) and not called by
  v2. Because the v2 install carries no HIV annotation, this is a reference
  difference between arms, not evidence that v3 is more sensitive. Per the
  interpretation rule, failure of the v2 arm to recover HIV is recorded as a
  sensitivity failure of that arm as configured.

## 5. Exact-read evidence (read extraction only, no new reference)

See `fresh_control_evidence.tsv`. All expected-target steps succeeded (read-lineage
records: 1,555,661; 229,228; 26,875; 1,158,876; 336,555). Largest non-target step:
skipped for two EBV controls (no positive non-target call); Cercopithecine
herpesvirus in SRR12682298 (5 records, a persistent low-count call: specificity
concern). For the HIV controls the script's "largest non-target" was another HIV
gene (IMMUNO_HIV1gp9; 301,796 and 9,713 records), so no non-HIV candidate was
traced there; this is an open limitation of the helper, not a result.
Read-lineage records are diagnostic and do not by themselves confirm infection.

## 6. Persistent unexpected calls

Low-molecule, non-target calls (HHV-6b, vesicular stomatitis, monkeypox, molluscum
contagiosum in the HIV controls; Cercopithecine herpesvirus in SRR12682298) persist
in both stacks, at reduced counts in v3 except Cercopithecine herpesvirus. They are recorded as specificity concerns. Calls present in v2 only
(HHV-1, HHV-8, parainfluenza) did not persist in v3. A lost call supports an artifact
explanation but does not prove viral absence.

## Not claimed

Viral absence, BCC association, calibrated diagnostic performance, general package
superiority, and any gate closure. Technical repeats (SRR6825025 and similar) are
not independent samples. v3-unique and v3-equal endpoints were not compared in the
fresh arm (exported only as AnnData layers).
