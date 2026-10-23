#!/usr/bin/env bash
# Real exit-139 fallback check on a tiny synthetic Drop-seq library (no on-list, random reads, zero mapped).
#SBATCH -J star139
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH -t 01:00:00
#SBATCH -p highmem
#SBATCH -o slurm_%x_%j.log
set -euo pipefail
M=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
export PATH=$M/conda/envs/viralscan_test_full/bin:$PATH PYTHONPATH=$M/vs_work/star-confirm/src
date; STAR --version
python $M/vs_work/star-confirm/analysis/star_confirm_acceptance/fallback139.py $M/viralscan_work/star_confirm_accept_2026-10/fallback139 \
  --index $M/viralscan_work/star_confirm_accept_2026-10/index
date; echo STAR139_DONE
