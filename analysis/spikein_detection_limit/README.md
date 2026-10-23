# Spike-in computational detection threshold (PLAN SPIKE-01 / SPIKE-02 / SPIKE-03)

Purpose: measure the computational detection threshold (source, background and layout-specific) of three arms on real background libraries with real viral reads
spiked in at known levels: the old kb index, the flank index v2, and STARsolo host+virus (the `star_confirm` rule).
Output: a threshold per background x arm x layout, a null (false-positive) floor, the flank keep/retire verdict, and the
STAR gate evidence (STAR-GATE-01 is decided only after this).

**Status (2026-10-23): the grid has been run** (Slurm array 25864985, 66/66 conditions) and rescored with the viral-only
scorer; results and verdict (flank index retired, F-028) are in `results_2026-10-23/` with `SHA256SUMS`. Do not rerun the
grid to reproduce them; the commands below document how they were produced.

**Known limit of the STAR fates:** `star_fates` reads `star_hostvirus/viral.bam`, which `star_confirm` prunes to `VIRUS|`
records, so a spike whose primary placement is on the host appears as `absent`, not `host`. `fate_host` in `star_gate.tsv`
therefore does not measure host competition; STAR detection (`on_truth`, any record on a truth contig) is unaffected.
Measuring host competition needs the fates extracted before pruning (a rerun).

## Design (user decisions 2026-10-21/22)

- Levels 3, 10, 30, 100, 300 spiked UMIs x 3 seeds x 2 layouts (5 cells | 100 cells) x 2 backgrounds = 60 conditions,
  plus 3 unspiked replicates per background = 66 (`grid.py`, `conditions.tsv`).
- Backgrounds: D03 SRR29320037 (5' 10xv2, reverse) with HHV-6B source GSM6416022; R01_SOC4 (3' 10xv3, forward) with EBV
  source GSM4796271. Each ~45 M pairs by seeded random subsample (`spikein.py bg --frac --seed`).
- Call = >= 3 kb molecules AND >= 1 STAR read on the same truth contig. kb arm = the product's multimap step
  (unique + allocated ambiguous molecules), not the raw gene matrix.
- Threshold = lowest level met in >= 2 of 3 seeds AND with every higher level also met (monotone rule). Labels: a level;
  `<=3` (left-censored: the lowest level is met, the true threshold may be lower); `>300 (not established)` (top level not
  met); `undetermined` (missing seeds or conditions leave a level undecided). A threshold that is not established is NOT a
  negative result. `non_monotone` marks a met level below an unmet one; the report warns.
- Flank index kept only if, in every background x layout, its threshold <= old threshold AND it adds <= 5 unsupported UMIs per
  null library (per-replicate difference vs old, maximum). Retire if any background fails; undetermined if data are
  missing. See the `aggregate.py` docstring for the exact definitions.

## Run order

All from `analysis/spikein_detection_limit/` with `PATH=$SITE_ROOT/conda/envs/viralscan_test_full/bin:$PATH`.

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
9. Rescore (SPIKE-03): `python rescore.py` re-scores every DONE condition with the current scorer into
   `bundle/scores_v2.tsv` (never overwrites `scores.tsv`; refuses conditions without DONE). Needed because bundles scored
   before the viral-only fix counted ~41k host genes as "viral" molecules in the floor columns
   (`kb_unsupported_molecules`, `kb_total_viral_molecules`, `n_false_calls`); truth-side columns are unchanged. Then
   `python aggregate.py check --scores-name scores_v2.tsv`.
10. Report: `python aggregate.py report --scores-name scores_v2.tsv [--out DIR]` (default `ROOT/aggregate`; `--scores-name`
   defaults to `scores.tsv`). It runs the check first and refuses to write anything, exit 2, if any condition is missing
   or failed. Outputs: `report.md`, `recovery_by_condition.tsv`, `lod.tsv` (the thresholds), `fp_floor.tsv`,
   `verdict.tsv`, `star_gate.tsv`.
11. Recovery ratio: `python recovery_ratio.py --scores-name scores_v2.tsv` writes `recovery_ratio.tsv` / `.md` (flank / old
   net kb recovery per background x layout x level; RATIO_MIN 0.80; flags `needs_user_decision` when the Q15 verdict is keep
   and the ratio is below it; never changes the verdict). Same refusal rule as the report.

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
recovery; STAR "absent" lumps unmapped-too-short with host-filter loss; thresholds are grid levels. The report's final section
states this again.

## Host-homologous arm and host-subtraction loss (PLAN DISC-12-HH; status: tooling only, NOT RUN)

Protocol `docs/plans/2026-10-22-discovery-benchmark.md` §1 Q1 / §4 and decision Q16. `spikein.py spike --pool outside_truth`
spikes ONLY the source reads the default mode drops as `alignment_outside_truth_set` (real HHV-6B / EBV reads that also align to
`HOST|` or another virus; use `--max-nh 20`, NH=1 sources cannot have such reads). Truth tables gain a trailing `pool` column
in that mode only; default output is byte-identical. An empty pool, or one smaller than the top level, exits loudly.
`hosthomolog.py` runs ONLY the product STAR host filter (`host_filter._starsolo_filter`, `pinned` then `star-default`) on
the grid's cached background plus the spikes and writes `cond/<id>/host_loss.tsv` from `lost_truth_counts()` (D15/D16).
30 conditions per pool (5 levels x 3 seeds x c100 x D03/R01); the ordinary pool is a second, separate root:

    PATH=<env>/bin:$PATH VIRALSCAN_SITE_ROOT=... python hosthomolog.py grid --pool outside_truth   # HOSTLOSS_ROOT default
    cd $HOSTLOSS_ROOT/slurm && sbatch --array=1-30%6 .../hosthomolog.sbatch
    HOSTLOSS_ROOT=<other> python hosthomolog.py grid --pool default   # baseline number, own 30 conditions

Caveat: the source BAMs are STARsolo alignments of host-UNMAPPED reads, so reads the host filter already removes are absent
from the pool by construction; the measured loss is a lower bound for the real library. Whether the EBV pool is non-empty is
checked by the run itself (it fails loudly).
