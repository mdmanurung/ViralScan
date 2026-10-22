#!/usr/bin/env bash
# Spike-in detection-limit grid (PLAN SPIKE-02): Slurm ARRAY, task i = row i of conditions.tsv, run end to end by
# run_condition.py (resumable: finished conditions and arms are skipped; fails loudly; DONE marker last).
# Before submitting (a human does this; no commit submits anything):
#   PATH=<env>/bin:$PATH python grid.py              # writes $SPIKE_GRID_ROOT/conditions.tsv and the slurm/ log dir
#   sbatch --array=1-66%6 grid.sbatch                # array size = number of rows; %6 caps concurrent disk use
# Peak scratch per task ~13 GB (STAR tmp + host-filtered reads) + ~6 GB inputs; cached backgrounds ~3 GB each (6 in all).
# Re-submit the same array to resume after a failure.
#SBATCH -J spikegrid
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -p highmem
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/spikein_grid_2026-10/slurm/%x_%A_%a.log
set -euo pipefail
M=/exports/archive/hg-funcgenom-research/mdmanurung
A=${SPIKE_GRID_CODE:-$M/vs_work/spikein/analysis/spikein_detection_limit}
exec "$M/conda/envs/viralscan_test_full/bin/python" "$A/run_condition.py" --row "${SLURM_ARRAY_TASK_ID:?run as a Slurm array}"
