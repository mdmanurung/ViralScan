# F-017 — Adding EBV type 2 moves 16 % of EBV molecules onto it through the equal split of shared k-mers

**Status:** confirmed (2026-09-29) · **Tags:** reference-panel, EBV, sibling, multimapping, dilution, max-panel

## Observation

Run: SRR12682296, EBV LCL, 10xv2, combined host+virus index.

The run was repeated on the 4,127-genome max panel (`viral_panel_max_2026-09-28/`, kb-ref job 25662863). The panel adds EBV type 2 `NC_009334.1`, the 35 F-015 accessions, 174 HPV types, influenza, SARS-CoV-2 lineages and 1,396 anelloviruses. The comparison baseline is the v1 2,215-genome panel. Both runs used the same code and the same `-w None` barcode handling, so the comparison is like for like.

| | v1 panel | max panel |
|---|---:|---:|
| EBV type 1 molecules | 906,202 | **759,972 (−16.1 %)** |
| EBV type 2 molecules | n/a | **146,634**, reported as 80 per-gene raw-ID rows |
| EBV 1 + 2 total | 906,202 | 906,606 (+0.04 %) |
| HSV-1, SRR8315713 | 32,404 | 32,419 (+0.05 %) |
| HHV-2 `possible_em_bleed` | fires | fires (591:1) |
| anellovirus false positives | 0 | 0 |

Layer split for the max-panel EBV run (`adata_multimap.h5ad`):

| layer | EBV-1 | EBV-2 |
|---|---:|---:|
| `counts_unique` | 183,964 | 1,357 |
| `counts_ambiguous_allocated` | 576,008 | 145,277 |

## Interpretation

- **The EBV-2 signal is allocation, not evidence.** 99 % of EBV-2 mass comes from EBV-1/EBV-2 shared equivalence classes. Unique support is 136:1 in favour of type 1.
- **Cause.** Under the default `host-conservative` method, virus–virus ambiguity is split *equally* (`docs/v3_counting_contract.md`). Two near-identical sibling genomes therefore each receive half of the shared mass, whatever their unique support.
- **Why no flag fired.**
  - EBV-2 resolves to raw gene IDs, not a virus name, so the output shows 80 rows.
  - `SIBLING_VIRUS_PAIRS` is keyed on legacy names, so `possible_em_bleed` cannot fire for EBV-2.

## Consequence

- **Do not ship the max panel as the default** until two things are in place:
  - sibling genomes group through the Virus Identity table (`sibling_group`, consolidation candidate A);
  - detection aggregates at virus level (candidate B).
- **Reconsider the equal split for virus–virus siblings.** A unique-weighted or EM allocation *within* a sibling group would give EBV-2 about 0.7 % instead of 16 %. That is a multimap-policy question for the D3 experiments.
- **Measure near-identical additions before relying on them.** The max panel also carries two copies each of HPV45/58/39/51 (96–98 % containment), OC43 ×2 and 18 SARS-CoV-2 genomes. See `near_duplicates.tsv`.

## Provenance

- v1 run: `ebv_latest_ref_2026-09-27/runs/combined/SRR12682296/`
- max-panel runs: `viral_panel_max_2026-09-28/runs/combined_max/{SRR12682296,SRR8315713}/`
- max-panel runs launched from `scripts/slurm_quant_max_regression.sh` (job 25662864)
- HEAD at run time: `759ce04`, with a dirty tree

## Update 2026-09-30 — size confirmed with barcode correction on

Job 25672273 reran the max panel with barcode correction on (SW-13/SW-20 fixed):

- EBV-1: 743,325 molecules. EBV-2: 142,963, which is **16.1 %** of 886,288, against 16.2 % uncorrected.
- The EBV total is +0.04 % against the corrected v1 run (885,933).
- HSV-1: 23,181 molecules (+0.05 % against the corrected v1 run); the HHV-2 bleed is 34.8.

The F-017 sizes therefore do not depend on barcode correction.
