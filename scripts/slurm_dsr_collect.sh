#!/usr/bin/env bash
# Status collector for DSR round-1 arrays. Submit with --dependency=afterany:<array jobid>:
#   JOBS="123 124" TAG=S_off sbatch --dependency=afterany:123:124 scripts/slurm_dsr_collect.sh
# Writes $ROUND_DIR/status/<TAG>.tsv (ROUND_DIR defaults to round 1; for round 2 also pass
# `-o $ROUND_DIR/logs/%x_%j.log`): one row per array task (state, elapsed, MaxRSS, exit).
#SBATCH -J dsr_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH -o slurm_%x_%j.log
set -uo pipefail
: "${JOBS:?}" "${TAG:?}"
# >>> dsr-env (tests/test_dsr_round_env.py evaluates this block alone)
M=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
R=${ROUND_DIR:-$M/viralscan_work/dsr_round1}
# <<< dsr-env
mkdir -p "$R/status"
O=$R/status/$TAG.tsv
printf 'JobID\tJobName\tState\tElapsed\tMaxRSS\tReqMem\tExitCode\n' > "$O"
for j in $JOBS; do
  sacct -j "$j" -P -n --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,ExitCode | awk -F'|' '$1 !~ /\.(extern|batch)$/ || $1 ~ /\.batch$/' | tr '|' '\t' >> "$O"
done
