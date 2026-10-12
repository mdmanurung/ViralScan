# CAT-30 retraction notice: DRAFT, NOT PUBLISHED

Status: awaiting the user's sign-off on the wording (PLAN `DEF-07`, R3.5). Nothing here is applied: `claims/registry.json` is unchanged and no file has been placed next to the published results.
Drafted 2026-10-12. Every number was re-read from its source; line references are in this worktree.

## 0. Corrections to the request that produced this draft (read first)

1. **The 27,098 figure is not a covid number.** "0 of 27,098 KP343824.1-hitting reads carry a non-homopolymer viral 31-mer" comes from the HPV16 positive-control dataset (GSE189670, SRR19537339, first 20M R2 reads): `.living/findings/hpv16-positive-control-and-polya-anellovirus-sink.md:27` and `:57`, PLAN.md:1573-1575. That is F-021, a poly-A sink in a different library. It supports "the same class of artefact" and must not be cited as read evidence about the covid samples. The covid evidence is F-019 (section 3).
2. **`1,167,103.0` is x213, not x216-g.jul2.** x213-g and x213-g.jul2 publish 1,167,103.0; both x216 directories (x216-g and x216-g.jul2) publish 1,605,631.0. PLAN.md:4365-4366 groups them loosely ("x213-g, x216-g.jul2"). Checked against the four `viral_summary.tsv` files (table in section 1).
3. **"17.63-19.05 % of cells" is `pct_infected` over `total_cells`** (all barcodes in the matrix, 578,938 and 586,305), not over called cells. Over called cells the same columns read 56.64 % (x213-g) and 64.77 % (x216-g). The notice should name the denominator.
4. F-005's registry line says "no genuine viral infection in these COVID PBMC samples" (FINDINGS_REGISTRY.md:11 and PLAN.md:4361-4362 quote it). F-019 says absence of true anellovirus is not established (`.living/findings/covid-anellovirus-signal-is-polyg-reads.md:31`), and PLAN `ANELLO-PRIOR` (PLAN.md:1565-1568) says anelloviruses are expected commensals. The draft therefore retracts the published *magnitude*, not the possible presence of a virus.

## 1. What was published

`covid_viralscan/results/<sample>/results/viral_summary.tsv` (four sample directories; verified by reading the files):

| Directory | `Alphatorquevirus` total_umi | infected_cells / total_cells | pct_infected | Also in |
|---|---:|---:|---:|---|
| `LUM-SJ-x213-g` | 1,167,103.0 | 102,065 / 578,938 | 17.6297 | `viral_summary.knee.tsv`, `per_cell_viral*.tsv` |
| `LUM-SJ-x213-g.jul2` | 1,167,103.0 | 102,065 / 578,938 | 17.6297 | `per_cell_viral.tsv` |
| `LUM-SJ-x216-g` | 1,605,631.0 | 111,672 / 586,305 | 19.0467 | `viral_summary.knee.tsv`, `per_cell_viral*.tsv` |
| `LUM-SJ-x216-g.jul2` | 1,605,631.0 | 111,672 / 586,305 | 19.0467 | `per_cell_viral.tsv` |

`covid_viralscan/results/SURVEY_SUMMARY.md` and `SURVEY_SUMMARY.jul2.md` repeat the figure as a combined 2,772,734.0 UMI (= 1,167,103 + 1,605,631) in 213,737 infected cells. The earlier "~90 % of cells at >= 5 UMI" statement is already closed (F-005, "Do NOT cite", `covid-viralscan-no-sars2-anellovirus-dominant.md:23`); this notice covers the numbers still sitting in the tables.

## 2. Notice text (proposed, for sign-off)

> **Retraction notice (CAT-30), 2026-10-DD.**
> The `Alphatorquevirus` rows in the `viral_summary.tsv`, `viral_summary.knee.tsv` and `per_cell_viral*.tsv` files of this directory (1,167,103.0 UMI in x213; 1,605,631.0 UMI in x216; 17.63 % and 19.05 % of all barcodes in the matrix) are **retracted**. They are not evidence of anellovirus infection or of a viral load in these samples.
> Re-extracting every read the index placed in an anellovirus equivalence class (1,269,176 reads, x213) showed that 99.999 % carry a homopolymer of 15 nt or more. Overall, 90.1 % of the reads were best-aligned to a viral contig (almost all to the poly(C) tract of encephalomyocarditis virus) and 9.8 % to the human genome; in a 200,000-read sample, 91 % were at least 80 % poly-G, the signature of empty clusters on two-colour Illumina chemistry. Only 9 reads were free of homopolymer and template-switch oligo, and all 9 aligned best to the human genome. No clean read had a viral best hit (F-005, F-019).
> What this does not show: low-level, divergent anellovirus is not excluded. The reference captures only a small part of anellovirus sequence space, and anelloviruses are common in human blood. SARS-CoV-2 was not detected in either sample and that result stands (F-005).
> The original files are unchanged and kept for the record. Do not cite their `Alphatorquevirus` rows.

Placement: a new file `RETRACTION_NOTICE_CAT-30.md` (the text above) is added **beside** each `viral_summary.tsv` in the four `results/` directories, plus one at `covid_viralscan/results/` next to `SURVEY_SUMMARY*.md`. No published file is edited, renamed or regenerated ("beside, never over"). Whether `covid_viralscan/` (outside the repo, untracked) can hold these files is the user's call; this draft does not write there.

## 3. Why it is retracted, with sources

| Statement | Source |
|---|---|
| 1,269,176 R2 reads in anellovirus ECs (x213); 1,143,775 (90.1 %) viral-best with a homopolymer >= 15 nt; 124,918 (9.8 %) host-best; 474 unmapped; 9 clean, all host-best | `.living/findings/covid-anellovirus-signal-is-polyg-reads.md:12-19` |
| 1,143,452 viral-best reads hit NC_001479.1 (EMCV) poly(C); 91 % of a 200,000-read sample >= 80 % poly-G | same, `:21-22` |
| Absence of true anellovirus not established | same, `:31` |
| F-005: 0 / 4,500,299 primary alignments on a viral contig (first 5M reads, primary only; mechanism revised by F-019) | `covid-viralscan-no-sars2-anellovirus-dominant.md:13`; `covid-anellovirus-signal-is-polyg-reads.md:30` |
| Current code reproduces 1,083,687 and 1,499,051 against the published 1,167,103 and 1,605,631 | `covid-anellovirus-signal-is-polyg-reads.md:44` |
| Forward-strand specific: 9-40x lower under `--strand reverse` | `:45`, FINDINGS_REGISTRY.md:20 (F-020) |
| The 19,785 aligned reads: 56.5 % contain the 10x TSO verbatim | `:51`, `:61` |
| Same artefact class in another library (poly-A, not poly-G): 27,098 / 20M reads, 0 with a non-homopolymer viral 31-mer (HPV16 control, not covid) | `hpv16-positive-control-and-polya-anellovirus-sink.md:27`, `:57`; PLAN.md:1573-1575 |
| Published wording agreed in grill Q11 | PLAN.md:1927; `covid-anellovirus-signal-is-polyg-reads.md:35` |

## 4. What replaces it

- Reference side, done: homopolymer mask >= 20 nt and the anellovirus low-complexity 31-mer mask (CAT-42/CAT-17). On the HPV16 control this took Gamma/Betatorquevirus from 10,037 to 0 molecules (`hpv16-...md` update 2026-10-03).
- Covid side, **not yet done**: no replacement `Alphatorquevirus` number exists. `DSR-04` (covid x213/x216 rerun on cat42d, `--strand reverse`, artefact-aware read filter, emptydrops) is open (PLAN.md:5158-5159). The notice says "retracted", not "replaced by X".

## 5. Proposed `claims/registry.json` change (NOT applied)

The registry already holds the claim as `rejected` (`claims/registry.json:94-110`, id `ttv-biological-signal`, `public: false`). It needs no new status; the enum (`schemas/v3/claim_registry.schema.json:92`) has no `retracted`, and the schema forbids extra fields (`additionalProperties: false`). What is stale is `reason`: it names the mechanism F-005 proposed ("cDNA-only host-homology false-positive risk"), which F-019 revised. Proposed one-field edit:

```json
{
  "id": "ttv-biological-signal",
  "reason": "Retracted (CAT-30): the published COVID Alphatorquevirus counts (1,167,103.0 UMI x213; 1,605,631.0 UMI x216) are 99.999 % homopolymer reads (about 90 % poly-G or poly(C), about 10 % host-best) with no clean viral read (F-005, F-019); low-level divergent anellovirus is not excluded."
}
```

Exact diff against the current line 110:

```diff
-      "reason": "Existing evidence demonstrates a cDNA-only host-homology false-positive risk."
+      "reason": "Retracted (CAT-30): the published COVID Alphatorquevirus counts (1,167,103.0 UMI x213; 1,605,631.0 UMI x216) are 99.999 % homopolymer reads (about 90 % poly-G or poly(C), about 10 % host-best) with no clean viral read (F-005, F-019); low-level divergent anellovirus is not excluded."
```

Applying it changes the sha256 of `claims/registry.json`, so `analysis/v3_artifact_inventory.tsv:5` (pin `claim-registry`) and any governance pin of that file must be re-pinned in the same commit, as earlier commits did. Not done here.

## 6. UNVERIFIED

- Whether the date and byline of the notice should be the retraction decision date (placeholder `2026-10-DD`).
- Read counts and percentages are taken from the F-019 finding text, not recomputed from the raw reads (the job directory `viralscan_work/f005_readorigin/` is outside the repo). The table at `covid-anellovirus-signal-is-polyg-reads.md:15-19` (rows `:16-19`) is internally consistent: 9 clean of 1,269,176 is 0.0007 %, so 99.999 % carrying a homopolymer follows from it. The EMCV count (1,143,452, `:21`) is a full count; only the 91 % poly-G figure (`:22`) comes from a 200,000-read sample.
- That `covid_viralscan/` is the location the user means by "the published results" and that adding files there is acceptable (directory is untracked, outside the repo's governance pins).
- Whether `per_cell_viral*.tsv` / `viral_summary.knee.tsv` carry the same figures: they exist and the knee summary lists `Alphatorquevirus`, but only the four `viral_summary.tsv` files were read value by value.
- F-019 reproduced x213 only (job 25670680); x216 is covered by the same index and the reproduction run (job 25672275) but its read-level extraction was not done.
