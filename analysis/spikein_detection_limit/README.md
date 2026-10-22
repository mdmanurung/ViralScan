# Spike-in detection limit (PLAN SPIKE-01 / SPIKE-02)

Purpose: measure the empirical detection limit of three arms on real background libraries with real viral reads
spiked in at known levels: the old kb index, the flank index v2, and STARsolo host+virus (the `star_confirm` rule).
Output: a limit per background x arm x layout, a null (false-positive) floor, the flank keep/retire verdict, and the
STAR gate evidence (STAR-GATE-01 is decided only after this).

Nothing here has been run for real. Every command below is for a human to submit.

## Design (user decisions 2026-10-21/22)

- Levels 3, 10, 30, 100, 300 spiked UMIs x 3 seeds x 2 layouts (5 cells | 100 cells) x 2 backgrounds = 60 conditions,
  plus 3 unspiked replicates per background = 66 (`grid.py`, `conditions.tsv`).
- Backgrounds: D03 SRR29320037 (5' 10xv2, reverse) with HHV-6B source GSM6416022; R01_SOC4 (3' 10xv3, forward) with EBV
  source GSM4796271. Each ~45 M pairs by seeded random subsample (`spikein.py bg --frac --seed`).
- Call = >= 3 kb molecules AND >= 1 STAR read on the same truth contig. kb arm = the product's multimap step
  (unique + allocated ambiguous molecules), not the raw gene matrix.
- Limit = lowest level met in >= 2 of 3 seeds AND with every higher level also met (monotone rule).
- Flank index kept only if, in every background x layout, its limit <= old limit AND it adds <= 5 unsupported UMIs per
  null library (per-replicate difference vs old, maximum). Retire if any background fails; undetermined if data are
  missing. See the `aggregate.py` docstring for the exact definitions.

## Run order

All from `analysis/spikein_detection_limit/` with `PATH=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_test_full/bin:$PATH`.

1. Tests (foreground, no cluster): `PYTHONPATH=../../src python -m pytest ../../tests/test_spikein_*.py -q`
2. Smoke: `sbatch smoke.sbatch` (D03, first 5 M pairs, levels 10/100, seed 1; done 2026-10-21, 20 min, see PLAN SPIKE-01).
3. Write the conditions table and log directory: `python grid.py` (writes `$SPIKE_GRID_ROOT/conditions.tsv`; default root
   `viralscan_work/spikein_grid_2026-10`).
4. Dry run one condition (prints every command, touches nothing): `python run_condition.py --cond-id D03-L100-s1-c005 --dry-run`
5. Pilot condition for real, ONE array task, before the grid: find the row (`grep -n D03-L100-s1-c005 conditions.tsv`, row =
   line - 1) and `sbatch --array=<row> grid.sbatch`. Check its `cond/<id>/bundle/scores.tsv`, the STAR spike names surviving
   the host filter (the runner fails loudly at level >= 100 with 0 STAR reads), and the flags passed to `viralscan run`.
   These were UNVERIFIED when written.
6. Grid: `sbatch --array=1-66%6 grid.sbatch`.
7. Resume after any failure: submit the same array again. Finished conditions (`DONE`) and arms (`stage_<arm>.done`) are
   skipped; half-finished arms are wiped and redone. To redo only the failed ones: `sbatch --array=<comma list>%6 grid.sbatch`.
8. Completeness: `python aggregate.py check` (exit 0 = all 66 DONE with valid scores; lists every problem otherwise).
9. Aggregate: `python aggregate.py report [--out DIR]` (default `ROOT/aggregate`). It runs the check first and refuses to
   write anything, exit 2, if any condition is missing or failed. Outputs: `report.md`, `recovery_by_condition.tsv`,
   `lod.tsv`, `fp_floor.tsv`, `verdict.tsv`, `star_gate.tsv`.

## Budget

From `grid.sbatch` (the explorer report itself was not available to this task, so these are the sbatch settings, not
measurements): per task 8 cores, 80 GB RAM, 12 h limit, `highmem` partition; peak scratch ~13 GB per task (STAR tmp +
host-filtered reads) + ~6 GB inputs, so at `%6` about 115 GB concurrent plus cached backgrounds (~3 GB each, 6 in all,
~18 GB). Final bundles are small. Runtime at 45 M pairs is not measured: the only data point is the 5 M-pair smoke (20 min for
two levels and two kb arms), and the host filter repeats per arm and dominates. Treat 12 h as a ceiling, not an estimate.
Total = 66 tasks; wall time depends on `%6` and queue.

## What the report does NOT claim

Single source virus per background; NH policy for EBV (`--max-nh 2`) unvalidated; no error-free simulation arm in the
grid; spikes come from one real library per virus (seeds resample reads, not libraries); kb arms have no per-class
recovery; STAR "absent" lumps unmapped-too-short with host-filter loss; limits are grid levels. The report's final section
states this again.
