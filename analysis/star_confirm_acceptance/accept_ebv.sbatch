#!/usr/bin/env bash
# usage: sbatch accept_ebv.sbatch <worktree-src> ; output dir is created next to the other work dirs
#SBATCH -J star_accept
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 04:00:00
#SBATCH -p highmem
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/star_confirm_accept_2026-10/slurm_%x_%j.log
set -euo pipefail
M=/exports/archive/hg-funcgenom-research/mdmanurung
export PATH=$M/conda/envs/viralscan_test_full/bin:$PATH PYTHONPATH=${1:-$M/vs_work/star-confirm/src}
ACC=$M/viralscan_work/star_confirm_accept_2026-10
mkdir -p $ACC
date; STAR --version
python -I -c 'import sys' && python $M/vs_work/star-confirm/analysis/star_confirm_acceptance/accept_ebv.py $ACC
date; echo STAR_ACCEPT_DONE
