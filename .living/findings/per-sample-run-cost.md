# F-024: Per-sample run cost is about 1.3 core-h (combined) and 5.6 core-h (STAR two-step) per 50 M pairs

**Status:** preliminary. Based on 7 completed full-depth jobs on 3 real libraries, measured 2026-10-05.
**Task:** `COST-01` (VAL-01 D5: measure cost before sizing the `defaults_selection` grids)
**Script:** `scripts/cost01_per_sample_cost.py`
**Output:** `/exports/para-lipg-hpc/mdmanurung/ViralScan/cost01/{per_sample_cost.tsv,sacct_snapshot.txt}`

## Claim
Normalised to the VAL-01 sample (50 M read pairs = 2,000 cells × 25 k reads), using Slurm allocated core-hours (CPUTimeRAW):

| arm | jobs | alloc core-h / 50 M | used CPU-h / 50 M | peak RAM |
|---|---|---|---|---|
| combined (incl. corrected) | 5 (8 cores) | 0.94–1.79, mean 1.31 | 0.41–0.71 | 1.7–9.0 GB |
| STAR two-step | 2 (16 cores) | 5.58–5.62 | 1.33–1.42 | 15.1 GB |

These are below the val01 §5 estimates (2 and 8 core-h).

## Projection for grid sizing
- `rerun-multimap` skips `kb count`, so the 4 multimap modes reuse one quantification per sample and host strategy (`menu.py:311`).
- Per training sample, both host strategies cost 1.8 + 5.6 = **7.4 core-h** (upper combined value).
- The upper bound is **18.2 core-h**, assuming each of the 6 extra multimap reruns costs as much as a full combined run. Rerun cost is unmeasured; EM modes reprocess the BUS file.
- The shared-default share of the budget is 55 % of 20 k = 11 k core-h. That covers about **600–1,490 training samples**, counting each unplanted twin as a sample (the twins double the run count).

## Caveats
- Real 10xv2/dropseq libraries, not synthetic panels. Planted viral reads are a tiny share of reads, so the per-read cost should carry over.
- No 10xv3 two-step row: job 25652882_1 failed.
- Node variance is about 1.5×: the same hsv1 library took 2.89 and 1.88 alloc core-h.
- CPU efficiency is about 25–40 %. Requesting 4 cores instead of 8–16 would cut allocated cost roughly in half; this was not tested.
