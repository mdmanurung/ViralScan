# DSR round dossier: what must exist to support a finding (DOSSIER-01)

Status: implemented. Machine spec: `analysis/dsr_round1/dossier_spec.tsv` (46 items); roles and acceptance per
sample: `analysis/dsr_round1/dossier_roles.tsv`; checker: `scripts/dsr05_dossier.py` (read-only).

```
python scripts/dsr05_dossier.py <round_dir>      # writes <round_dir>/dossier/
```

## Why

A call, or the absence of one, is only as good as the files behind it. The dossier lists those files for every
sample, checks each against the round on disk and states which claims the files can carry. It never regenerates
anything. Interim rule (user, 2026-10-07): an unexpected call is **unverified** until its evidence items are `ok`.

## Sections

| § | Section | Unit | Required items (ids in the TSV) |
|---|---|---|---|
| A | Identity and provenance | sample, arm | chemistry JSON (A01); FASTQ sha256 R1/R2 (A02-3); deposit md5 check (A04); code SHA (A05); reference fingerprint (A06); STAR host index fingerprint (A07); tool versions (A08); command line (A09); `run_complete.json` (A10) |
| B | Input QC and counting | arm | `run_info.json` (B01), `count_audit.tsv` (B02), host-filter audit (B03), read-filter audit (B04) |
| C | Cell calling | arm, sample | `emptydrops_cells.tsv` (C01), `called_cells.tsv` (C02), `host_called_cells.tsv` (C03), parameters (C04), `cell_calling_summary.tsv` (C05), vendor barcodes (C06, only if a role names them), frozen reference cells (C07), cell-level table (C08) |
| D | Detection per arm | arm | `viral_summary`, `per_cell_viral`, `multimap_evidence`, `virus_identity`, `sensitivity` (D01-5); alignment columns or an explicit n/a (D06); `common_cells_summary.tsv` (D07) |
| E | Read-level evidence | call | evidence directory files (E01-6), `molecule_verdict_summary.tsv` (E07), evidence code SHA (E08), a row in `verdicts.tsv` (E09) |
| F | Specificity and artefacts | round, sample | `recurrence.tsv` (F01), `arm_concordance.tsv` (F02) |
| G | Sensitivity and negative claims | sample | certified positive control (G01), planted-read recovery (G02) |
| H | Biology | sample | cell-type or donor labels (H01) |
| I, J | Role and reporting | dataset | role row (I01), findings recorded (J01), claim-to-artifact map (J02) |

A **call** is a virus with at least 3 estimated molecules, or any anellovirus call, in an operative arm
(`scripts/dsr02_enumerate_calls.py`). Operative arms: `combined_off`, `combined_artefact`, `twostep_v2`.

## Statuses

`ok` present and non-empty; `missing`; `failed` (the arm directory exists but has no `run_complete.json`);
`stale` (older than the file it derives from); `n_a` (recorded in `status/n_a.tsv`, or an explicit reason such as
"--anello-align not run"); `not_started` (AIFI, Terekhova, GEO aging PBMC: access blockers).

## Claim ladder

An item gates a rung only on the sample's reference arm, `combined_off`.

| Rung | Claim | Needs |
|---|---|---|
| C0 | N cells were called | A, B, C items marked C0 |
| C1 | virus X has at least 1 UMI in n of N cells (kallisto level, not verified) | C0 and D |
| C2 | virus X is supported at read level | C1 and, for every call, E01-3, E06, E07 and E09; at least one call must exist |
| C3 | virus X infects cells of type T | C2 and H01 |
| N1 | no panel virus above the LOD | C1 and G01, G02 |
| N2 | no virus | never claimable: the panel is the limit |

## First result on `dsr_round1` (2026-10-07)

2,593 item rows: 1,480 ok, 887 missing, 194 n_a, 17 stale, 12 failed, 3 not started. Every sample tops out at
**C1** and no negative is informative (N1). The dossier shows why:

- `molecule_verdict_summary.tsv`, the evidence code SHA and the `verdicts.tsv` row are missing for all 140 calls;
  53 of the 140 calls have no evidence directory at all (calls are enumerated from the runs, `calls.tsv` is stale).
- `common_cells_summary.tsv` is stale for all 17 samples (older than the `twostep_v2` redetects).
- No deposit md5 check (A04), no STAR host index fingerprint (A07), and no `time.txt` for any `twostep_v2` run (A09).
- `positive_control.json` is `not-configured` and there are no planted-read recoveries (G01, G02) for any sample.
- No cell-type or donor labels (H01), and no claim-to-artifact map (J02) for any sample.
- 9 of 17 samples have no recorded finding (J01); x223's Cell Ranger barcodes are not located (C06).
- Derived `cell_calling_summary.tsv`: the twostep host-matrix cell sets overlap `combined_off` poorly (Jaccard
  hhv6b 0.03, covid x213 0.40) and hhv6b's knee equals its inflection (253.6).
- Derived `recurrence.tsv`: Cercopithecine herpesvirus (9 samples, 5 datasets), HPV77 (7, 4), HHV-8 (7, 4, in
  datasets with no KSHV), HSV-1 (5, 4), HPV118 (4, 4), HPV29 (7, 3).
