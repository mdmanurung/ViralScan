#!/usr/bin/env bash
# Status collector for DSR round-1 arrays. Submit with --dependency=afterany:<array jobid>:
#   JOBS="123 124" TAG=S_off sbatch --dependency=afterany:123:124 scripts/slurm_dsr_collect.sh
# Writes status/<TAG>.tsv: one row per array task (state, elapsed, MaxRSS, exit).
#SBATCH -J dsr_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/logs/%x_%j.log
set -uo pipefail
: "${JOBS:?}" "${TAG:?}"
O=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/status/$TAG.tsv
printf 'JobID\tJobName\tState\tElapsed\tMaxRSS\tReqMem\tExitCode\n' > "$O"
for j in $JOBS; do
  sacct -j "$j" -P -n --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,ExitCode | awk -F'|' '$1 !~ /\.(extern|batch)$/ || $1 ~ /\.batch$/' | tr '|' '\t' >> "$O"
done
