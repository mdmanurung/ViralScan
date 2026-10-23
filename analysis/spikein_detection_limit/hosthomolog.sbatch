#!/usr/bin/env bash
# DISC-12-HH host-subtraction loss: Slurm ARRAY, task i = row i of $HOSTLOSS_ROOT/conditions.tsv (hosthomolog.py run).
# Before submitting (a human does this; no commit submits anything):
#   PATH=<env>/bin:$PATH python hosthomolog.py grid --pool outside_truth   # 30 rows + slurm/ log dir
#   cd $HOSTLOSS_ROOT/slurm && sbatch --array=1-30%6 <this>/hosthomolog.sbatch
# Default (ordinary) pool = separate 30 conditions: HOSTLOSS_ROOT=<other root> ... grid --pool default.
# Needs the SPIKE-02 grid's cached backgrounds (read-only). Two STAR host filters per task (pinned, star-default).
# Scratch per task ~13 GB per filter run (deleted after each) + ~6 GB inputs. Re-submit the same array to resume.
#SBATCH -J hosthomolog
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -p highmem
#SBATCH -o %x_%A_%a.log  # relative to the submit directory
set -euo pipefail
M=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
A=${SPIKE_GRID_CODE:-$M/vs_work/disc12-hosthomolog/analysis/spikein_detection_limit}
exec "$M/conda/envs/viralscan_test_full/bin/python" "$A/hosthomolog.py" run --row "${SLURM_ARRAY_TASK_ID:?run as a Slurm array}"
